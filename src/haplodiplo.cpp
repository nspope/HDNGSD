// TODO: missing data handler should rescale genotype likelihoods for haploids in log space. Otherwise can get underflow if heterozygote GL is very high. -> DONE, but not all methods are using it

// TODO: I think I have to explicitly initialize "do not copy" members in parallel workers

// TODO: How about... various methods for getting genotype probs. And covariances. For example, SAF is one. "Admixture" ala PCANGSD another.
//  -basically right now "Covar" estimates genotype probs based on given admixture coefficients.
//  -instead, put this functionality in a new method for calculating genotype probs
//  -and have one method for Covar that always uses genotype probs

// TODO: ShfWindow takes way to long to calculate variant site pairs. E.g. with 400K sites, it would take a very very very long time to iterate through all pairs. Alternative strategies:
//    -if we have a hard distance threshold, we can ignore pairings beyond this threshold. This requires initial sorting based on pos (good idea in any case)
//    -could maintain a list of "current distances" for all sites, that would be decremented uniformly as we move variant to variant. Then, for each variant, we can calculate # invariant by a simple query. The complexity is then nvariant * total rather than total^2.
//    -combine these two ideas: track site, move forward and back until threshold reached
//    -what about tracking boundary positions (in terms of indices). Basically:
//       -move to next variant; set ideal boundaries as ibound[j] = cM[i] + bins[j]
//       -take current boundary points, move right until cM[k] > ibound[j], set bound[j] = k-1
//       -within slices defined by bound[j], get variants, count total number of sites
// (this is partially complete??)

#include <RcppArmadillo.h> 
#include <RcppParallel.h> 
#include <string>
#include <sys/stat.h>
#include <vector>
#include <iterator>

// [[Rcpp::depends(RcppArmadillo,RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
arma::cube cycle_cube (const arma::cube& inp)
{
  arma::cube out (inp.n_cols, inp.n_slices, inp.n_rows);
  for (unsigned i=0; i<inp.n_rows; ++i)
    for (unsigned j=0; j<inp.n_cols; ++j)
      for (unsigned k=0; k<inp.n_slices; ++k)
        out.at(j,k,i) = inp.at(i,j,k);
  return out;
}

bool multi_loop (arma::uvec& bin, arma::uvec& nchr, int i)
{
  if(++bin[i] <= nchr[i]) return true;
  if(i==0) return false;
  if(multi_loop(bin, nchr, i-1)) 
  {
    bin[i] = 0;
    return true;
  }
  return false;
}

// [[Rcpp::export]]
arma::mat hypergeometric_basis (unsigned N, unsigned n)
{
  if (n > N)
    Rcpp::stop("[hypergeometric_basis] Can only down-project");

  // 0 <= K <= N are bins in the original SFS
  // 0 <= k <= n are bins in the projected SFS
  arma::mat out (n + 1, N + 1, arma::fill::zeros);
  for (unsigned K=0; K<=N; ++K)
    for (unsigned k=0; k<=n; ++k)
      out.at(k,K) += 
        exp(R::lchoose(K, k) + R::lchoose(N-K,n-k) - R::lchoose(N,n));

  return out;
}

namespace utils 
{

  arma::umat all_haplotype_configurations (const unsigned n)
  {
    // generate all haplotype configurations summing to n chromosomes,
    // with the reference state omitted. Ordered decreasing in n.

    arma::umat out (3, 0);
    arma::uvec c (3);
    for (c[0]=0; c[0]<=n; ++c[0])
      for (c[1]=0; c[1]<=n; ++c[1])
        for (c[2]=0; c[2]<=n; ++c[2])
          if (arma::accu(c)<=n)
            out.insert_cols(out.n_cols, c);
    out = out.cols(arma::stable_sort_index(arma::trans(arma::sum(out, 0)), "descend"));
    return out;
  }

  double multichoose (const arma::uvec vals)
  {
    unsigned n = arma::accu(vals);
    double out = std::lgamma(n + 1);
    for (auto v : vals)
      out -= std::lgamma(v + 1);
    return exp(out);
  }

  std::tuple<arma::sp_mat, arma::umat> haplotype_projection (const arma::umat cfg, const unsigned prj)
  {
    // generate a hypergeometric projection of a set of haplotype configurations onto
    // a set of configurations of smaller cardinality

    if (cfg.n_rows != 3)
      Rcpp::stop("[haplotype_projection] Each haplotype configuration (column) must have 3 entries");
    if (prj < 1 || prj > cfg.max()) 
      Rcpp::stop("[haplotype_projection] Projection size must be positive and less than the number of chromosomes");

    const double errtol = std::pow(10., -9);

    // generate new configuration
    auto config = utils::all_haplotype_configurations(prj);

    // loop over old configuration; for each entry in new configuration, get hypergeometric probability
    unsigned n = cfg.max();
    arma::mat proj (config.n_cols, cfg.n_cols);

    for (unsigned c1 = 0; c1 < cfg.n_cols; ++c1)
      for (unsigned c2 = 0; c2 < config.n_cols; ++c2)
        proj.at(c2, c1) = exp(
            R::lchoose(n-arma::accu(cfg.col(c1)), prj-arma::accu(config.col(c2))) + 
            R::lchoose(cfg.at(0,c1), config.at(0,c2)) + 
            R::lchoose(cfg.at(1,c1), config.at(1,c2)) + 
            R::lchoose(cfg.at(2,c1), config.at(2,c2)) - 
            R::lchoose(n, prj)); // multivariate hypergeometric pmf

    // sparsify projection matrix
    proj.clean(errtol);
    proj.each_col([](arma::vec& x){ return x/arma::accu(x); }); //normalize

    return std::make_tuple(arma::sp_mat(proj), config);
  }

  std::tuple<arma::vec, std::vector<arma::umat>>
    hfs1d_downproject (arma::vec hfs, arma::umat configuration, const unsigned projection)
  {
    if (hfs.n_elem != configuration.n_cols || configuration.n_rows != 3)
      Rcpp::stop("[downproject_hfs1d] Dimension mismatch");
    if (projection > configuration.max())
      Rcpp::stop("[downproject_hfs1d] Size of projection must be less than current HFS");

    auto projector0 = utils::haplotype_projection(configuration, projection);
    arma::sp_mat proj0 = std::get<0>(projector0);
    arma::umat config0 = std::get<1>(projector0);

    arma::vec hfsp (config0.n_cols, arma::fill::zeros);

    for (unsigned i=0; i<hfs.n_elem; ++i)
      for (arma::sp_mat::const_col_iterator val0 = proj0.begin_col(i); val0 != proj0.end_col(i); ++val0)
        hfsp.at(val0.row()) += (*val0) * hfs.at(i);

    std::vector<arma::umat> config = {config0};

    return std::make_tuple(hfsp, config);
  }

  std::tuple<arma::mat, std::vector<arma::umat>>
    hfs2d_downproject (arma::mat hfs, 
                       arma::umat configuration0, 
                       const unsigned projection0,
                       arma::umat configuration1,
                       const unsigned projection1)
  {
    if (hfs.n_rows != configuration0.n_cols || hfs.n_cols != configuration1.n_cols ||
        configuration0.n_rows != 3          || configuration1.n_rows != 3          )
      Rcpp::stop("[downproject_hfs2d] Dimension mismatch");
    if (projection0 > configuration0.max() || projection1 > configuration1.max())
      Rcpp::stop("[downproject_hfs2d] Size of projection must be less than current HFS");

    auto projector0 = utils::haplotype_projection(configuration0, projection0);
    arma::sp_mat proj0 = std::get<0>(projector0);
    arma::umat config0 = std::get<1>(projector0);

    auto projector1 = utils::haplotype_projection(configuration1, projection1);
    arma::sp_mat proj1 = std::get<0>(projector1);
    arma::umat config1 = std::get<1>(projector1);

    arma::mat hfsp (config0.n_cols, config1.n_cols, arma::fill::zeros);

    for (unsigned i=0; i<hfs.n_rows; ++i)
      for (unsigned j=0; j<hfs.n_cols; ++j)
        for (arma::sp_mat::const_col_iterator val0 = proj0.begin_col(i); val0 != proj0.end_col(i); ++val0)
          for (arma::sp_mat::const_col_iterator val1 = proj1.begin_col(j); val1 != proj1.end_col(j); ++val1)
            hfsp.at(val0.row(), val1.row()) += (*val0) * (*val1) * hfs.at(i,j);

    std::vector<arma::umat> config = {config0, config1};

    return std::make_tuple(hfsp, config);
  }

  std::tuple<arma::cube, std::vector<arma::umat>>
    hfs3d_downproject (arma::cube hfs, 
                       arma::umat configuration0, 
                       const unsigned projection0,
                       arma::umat configuration1,
                       const unsigned projection1,
                       arma::umat configuration2,
                       const unsigned projection2)
  {
    if (hfs.n_rows != configuration0.n_cols || hfs.n_cols != configuration1.n_cols || hfs.n_slices != configuration2.n_cols ||
        configuration0.n_rows != 3          || configuration1.n_rows != 3          || configuration2.n_rows != 3            )
      Rcpp::stop("[downproject_hfs3d] Dimension mismatch");
    if (projection0 > configuration0.max() || projection1 > configuration1.max() || projection2 > configuration2.max())
      Rcpp::stop("[downproject_hfs3d] Size of projection must be less than current HFS");

    auto projector0 = utils::haplotype_projection(configuration0, projection0);
    arma::sp_mat proj0 = std::get<0>(projector0);
    arma::umat config0 = std::get<1>(projector0);

    auto projector1 = utils::haplotype_projection(configuration1, projection1);
    arma::sp_mat proj1 = std::get<0>(projector1);
    arma::umat config1 = std::get<1>(projector1);

    auto projector2 = utils::haplotype_projection(configuration2, projection2);
    arma::sp_mat proj2 = std::get<0>(projector2);
    arma::umat config2 = std::get<1>(projector2);

    arma::cube hfsp (config0.n_cols, config1.n_cols, config2.n_cols, arma::fill::zeros);

    for (unsigned i=0; i<hfs.n_rows; ++i)
      for (unsigned j=0; j<hfs.n_cols; ++j)
        for (unsigned k=0; k<hfs.n_slices; ++k)
          for (arma::sp_mat::const_col_iterator val0 = proj0.begin_col(i); val0 != proj0.end_col(i); ++val0)
            for (arma::sp_mat::const_col_iterator val1 = proj1.begin_col(j); val1 != proj1.end_col(j); ++val1)
              for (arma::sp_mat::const_col_iterator val2 = proj2.begin_col(k); val2 != proj2.end_col(k); ++val2)
                hfsp.at(val0.row(), val1.row(), val2.row()) += (*val0) * (*val1) * (*val2) * hfs.at(i,j,k);

    std::vector<arma::umat> config = {config0, config1, config2};

    return std::make_tuple(hfsp, config);
  }
}

struct MultiAllelicGenotypeLikelihood
{
  const arma::uvec ploidy;
  const unsigned samples;
  unsigned loci;

  std::vector<arma::cube> like;
  arma::umat pos;

  MultiAllelicGenotypeLikelihood (const arma::uvec ploidy)
    : ploidy (ploidy)
    , samples (ploidy.n_elem)
    , like (0)
    , pos (4, 0)
    , loci (0)
  {}

  void add_locus (const unsigned chrom, const unsigned start, const unsigned end, arma::imat genotypes, const double dropout, const double err)
  {
    if (genotypes.n_cols != samples || 
        genotypes.n_rows != 2       )
      Rcpp::stop("[MultiAllelicGenotypeLikelihood::add_locus] Dimension mismatch");

    if (err < 0 || err > 1)
      Rcpp::stop("[MultiAllelicGenotypeLikelihood::add_locus] Genotyping error rate must be in [0,1]");
    if (dropout < 0 || dropout > 1)
      Rcpp::stop("[MultiAllelicGenotypeLikelihood::add_locus] Dropout error rate must be in [0,1]");
    
    // convert genotypes to unique codes
    std::vector<unsigned> nonmissing;
    for (unsigned i=0; i<samples; ++i)
    {
      if (genotypes.at(0,i) <= 0 && genotypes.at(1,i) <= 0)
      {
        genotypes.at(0,i) = -1;
        genotypes.at(1,i) = -1;
      } 
      else if (genotypes.at(0,i) > 0 && genotypes.at(1,i) <= 0)
      {
        if (ploidy.at(i) != 1)
          Rcpp::stop("[MultiAllelicGenotypeLikelihood::add_locus] One genotype must be missing (coded as 0) for haploids");
        genotypes.at(1,i) = -1;
        nonmissing.push_back(genotypes.at(0,i));
      }
      else if (genotypes.at(0,i) <= 0 && genotypes.at(1,i) > 0)
      {
        if (ploidy.at(i) != 1)
          Rcpp::stop("[MultiAllelicGenotypeLikelihood::add_locus] One genotype must be missing (coded as 0) for haploids");
        genotypes.at(0,i) = genotypes.at(1,i);
        genotypes.at(1,i) = -1;
        nonmissing.push_back(genotypes.at(0,i));
      }
      else
      {
        if (genotypes.at(1,i) < genotypes.at(0,i))
        {
          int tmp = genotypes.at(1,i);
          genotypes.at(1,i) = genotypes.at(0,i);
          genotypes.at(0,i) = tmp;
        }
        nonmissing.push_back(genotypes.at(0,i));
        nonmissing.push_back(genotypes.at(1,i));
      }
    }

    arma::uvec unique_alleles = arma::unique(arma::uvec(nonmissing));
    unsigned alleles = unique_alleles.n_elem;
    if (alleles == 0)
      Rcpp::stop("[MultiAllelicGenotypeLikelihood::add_locus] No observed alleles");
    for (unsigned i = 0; i < alleles; ++i)
      genotypes.replace(unique_alleles[i], i);

    // calculate likelihoods
    arma::cube likelihood (alleles, alleles, samples, arma::fill::zeros);
    for (unsigned i=0; i<samples; ++i)
      if (ploidy.at(i) == 1)
        likelihood.slice(i) = wang_error_model_haploid(alleles, genotypes.at(0,i), err);
      else if (ploidy.at(i) == 2)
        likelihood.slice(i) = wang_error_model_diploid(alleles, genotypes.at(0,i), genotypes.at(1,i), dropout, err);
    
    // store
    like.push_back(likelihood);
    pos = arma::join_horiz(pos, arma::uvec({chrom, start, end, alleles}));
    loci++;
    fprintf(stderr, "[MultiAllelicGenotypeLikelihood::add_locus] Inserted locus at %u:%u-%u with %u alleles\n", chrom, start, end, alleles);
  }

  arma::mat wang_error_model_haploid (const unsigned k, const int u, const double E2)
  {
    double e2 = E2/(double(k) - 1.);

    arma::mat L (k, k);
    L.fill(arma::datum::nan);

    if (u < 0)
    {
      for (unsigned w=0; w<k; ++w)
        L.at(w,w) = 1./double(k);
    }
    else
    {
      for (unsigned w=0; w<k; ++w)
        L.at(w,w) = w == u ? (1 - E2) : e2;
    }

    return L;
  }

  arma::mat wang_error_model_diploid (const unsigned k, const int u, const int v, const double E1, const double E2)
  {
    double e1 = E1/(1. + E1),
           e2 = E2/(double(k) - 1.);

    arma::mat L (k, k);
    L.fill(arma::datum::nan);

    if (u < 0 || v < 0)
    {
      for (unsigned w=0; w<k; ++w)
        for (unsigned x=w; x<k; ++x)
          L.at(w,x) = 1./std::pow(double(k),2);
    }
    else
    {
      double l;
      for (unsigned w=0; w<k; ++w)
        for (unsigned x=w; x<k; ++x)
        {
          if (w == x)
          {
            if (u == v && v == w)
              l = std::pow(1 - E2, 2);
            else if ((u == w && v != w) || (v == w && u != w))
              l = 2. * e2 * (1 - E2); 
            else
              l = (2 - int(u == v)) * std::pow(e2, 2);
          }
          else
          {
            if ((u == w && v == x) || (u == x && v == w)) // correct genotype, both match
              l = std::pow(1 - E2, 2) + std::pow(e2, 2) - 2 * e1 * std::pow(1 - E2 - e2, 2);
            else if ((u == v && v == w) || (u == v && v == x)) // genotype appears homozygous with one match
              l = e2 * (1 - E2) + e1 * std::pow(1 - E2 - e2, 2);
            else if (u != w && u != x && v != w && v != x) // incorrect genotype, zero match
              l = (2 - int(u == v)) * std::pow(e2, 2);
            else // genotype appears heterozygous with one match T[1,2] O[1,3]
              l = e2 * (1 - E2 + e2);
          }
          L.at(w,x) = l;
        }
      // step 1 ...
      // aa -> [aa; 1]
      // ab -> [aa; e1] [bb; e1] [ab; 1-2e1]
      // ba -> [aa; e1] [bb; e1] [ba; 1-2e1]
      // step 2 ... [ok I get this. (1-E2+e2) is prob that a stays a OR mutates to b
      // aa -> [aa: (1-E2)(1-E2)] [[ab: e2(1-E2)] [ba: (1-E2)e2]] [bb: e2 e2] [[bc: e2 e2] [cb: e2 e2]]
      // ab -> [aa: (1-E2)e2] [ab: (1-E2)(1-E2)] [bb: e2(1-E2)] 
      //       [ac: (1-E2)e2] [bc: e2 e2] [cb: e2 (1-E2)] [ca: e2 e2]
      //       [cc: e2 e2] [cd: e2 e2] [dc: e2 e2]
    }

    return L;
  }
};

struct GenotypeLikelihood
{
  // Read in beagle binary format >like
  // parse ploidy, find index for haploids / diploids >haploid, diploid, ploidy
  // replace missing (0.3333 0.3333 0.3333) with 0. and count number of missing chroms per site >missing
  // count number of chroms, samples, sites >chromosomes, samples, sites

  const std::string filepath;
  const arma::uvec ploidy;

  arma::uvec haploids, diploids; 

  unsigned chromosomes, samples; 
  size_t sites;

  arma::Mat<short>  missing; // missing data pattern: samples x sites
  arma::Cube<float> like; // genotype likelihoods: 3 x samples x sites
  arma::Cube<short> cnts; // major/minor allele counts: 2 x samples x sites
  arma::Cube<float> post; // genotype probabilities: 3 x samples x sites

  arma::umat pos; // chromosome, position : 2 x sites

  //FILE* handle;
  //unsigned buffer_size;
  //arma::cube like;
  //size_t offset;

//  GenotypeLikelihood (const std::string filepath, const arma::uvec& ploidy) 
//    : ploidy (ploidy)
//    , filepath (filepath)
//  {
//    if (arma::max(ploidy) > 2 || arma::min(ploidy) < 1)
//      Rcpp::stop ("GLReader: allowed ploidy are (1, 2)\n");
//    haploids = arma::find(ploidy == 1);
//    diploids = arma::find(ploidy == 2);
//    chromosomes = haploids.n_elem + diploids.n_elem * 2;
//    samples = ploidy.n_elem;
//
//    //get size of file and calculate number of columns
//    //n_chr[i] = filesize(saf_file[i]) / (n_site * sizeof(double)) - 1;
//    FILE *handle = std::fopen(filepath.c_str(), "r");
//    if (handle == NULL)
//      Rcpp::stop ("GLReader: could not open file " + filepath + "\n");
//
//    //check if gzipped
//    unsigned char magic[2];
//    std::fread(magic, sizeof(magic), 1, handle);
//    if (magic[0] == 0x1f && magic[1] == 0x8b)
//      Rcpp::stop ("GLReader: cannot read gzipped files\n");
//
//    //read in data
//    sites = filesize(filepath) / (sizeof(double) * samples * 3.);
//    like = arma::ones<arma::cube>(3, samples, sites);
//    std::cout << std::fread(like.memptr(), sizeof(double), sites * samples * 3, handle) << std::endl;
//    std::fclose(handle);
//
//    // normalize by site/sample (just in case), set to zero if GLs are equal, and calculate missingness
//    // TODO why not row/colsum missing
//    missing_samples = arma::zeros<arma::uvec>(sites);
//    missing_sites = arma::zeros<arma::uvec>(samples);
//    missing = arma::zeros<arma::umat>(samples, sites);
//    for (size_t i=0; i<sites; ++i)
//    {
//      for (unsigned j=0; j<samples; ++j)
//      {
//        like.slice(i).col(j) /= arma::accu(like.slice(i).col(j)); //pretty sure this is pointless, better to check if any are negative
//        if (fabs(like.at(0,i,j) - like.at(2,i,j)) < 1e-8)
//        {
//          like.slice(i).col(j) *= 0.;
//          missing_samples[i] += ploidy[j];
//          missing_sites[j] += 1;
//          missing.at(j,i) = 1;
//        }
//      }
//    }
//    std::fprintf(stderr, "GLReader: read %u sites across %u chromosomes in %u samples\n", sites, chromosomes, samples);
//  }
//
//  size_t filesize(const std::string& filename) {
//    struct stat st;
//    if(stat(filename.c_str(), &st) != 0) {
//      return 0;
//    }
//    return st.st_size;   
//  }
  
  GenotypeLikelihood (const arma::cube& lmat, const arma::ucube& cmat, const arma::umat& pmat, const arma::uvec& pldy) 
    : ploidy (pldy)
    , like (arma::conv_to<arma::Cube<float>>::from(lmat))
    , cnts (2, cmat.n_cols, cmat.n_slices, arma::fill::zeros)
    , post (arma::size(like))
    , pos (pmat)
  {
    if (arma::max(ploidy) > 2 || arma::min(ploidy) < 1)
      Rcpp::stop ("GLReader: allowed ploidy are (1, 2)\n");
    if (lmat.n_cols   != ploidy.n_elem ||
        lmat.n_rows   != 3             || 
        cmat.n_cols   != ploidy.n_elem ||
        cmat.n_slices != lmat.n_slices || 
        cmat.n_rows   != 4             || 
        pos.n_cols    != lmat.n_slices || 
        pos.n_rows    != 4             || 
        pos.row(2).max() > 3           ||
        pos.row(3).max() > 3           )
      Rcpp::stop ("GLReader: dimension mismatch\n");
    if (arma::any(arma::vectorise(like) > 1e-12))
      Rcpp::stop ("GLReader: genotype likelihoods should be negative (e.g. log-scaled)\n");

    haploids = arma::find(ploidy == 1);
    diploids = arma::find(ploidy == 2);
    chromosomes = haploids.n_elem + diploids.n_elem * 2;
    samples = ploidy.n_elem;
    sites = like.n_slices;

    // normalize by site/sample (just in case), set to zero if GLs are equal, and calculate missingness
    missing = arma::zeros<arma::Mat<short>>(samples, sites);
    for (size_t i=0; i<sites; ++i)
    {
      // figure out major and minor; 0123 ~ ACGT
      const unsigned major = pos.at(2,i), minor = pos.at(3,i);

      for (unsigned j=0; j<samples; ++j)
      {
        // rescale
        like.slice(i).col(j) -= like.slice(i).col(j).max();

        // recognize missing data, get counts
        if((fabs(like.at(0,j,i) - like.at(2,j,i)) < 1e-8) && 
           (fabs(like.at(0,j,i) - like.at(1,j,i)) < 1e-8)) // missing 
        {
          missing.at(j,i) = 1;
        } else {
          cnts.at(0,j,i) = cmat.at(major,j,i);
          cnts.at(1,j,i) = cmat.at(minor,j,i);
        }

        // transform to [0,1]
        like.slice(i).col(j)  =  arma::exp(like.slice(i).col(j));
        like.slice(i).col(j) /= arma::accu(like.slice(i).col(j)); 

        // posterior with uniform prior
        post.slice(i).col(j)  = arma::conv_to<arma::Col<float>>::from(handler(j,i));
      }
    }

    //// set genotype probabilities to NaN
    //post.fill(arma::datum::nan);

    // sort by chromosome then position
    arma::uvec index = arma::regspace<arma::uvec>(0, sites-1);
    std::sort (index.begin(), index.end(), 
      [&](int i, int j) { 
        return pos.at(0,i)==pos.at(0,j) ? pos.at(1,i)<pos.at(1,j) : pos.at(0,i)<pos.at(0,j); 
      });
    missing = missing.cols(index);
    like    = like.slices(index);
    cnts    = cnts.slices(index);
    post    = post.slices(index);
    pos     = pos.cols(index);

    std::fprintf(stderr, "[GLReader] Read %lu sites across %u chromosomes in %u samples\n", sites, chromosomes, samples);
    std::fprintf(stderr, "[GLReader] NB: sites sorted according to linkage group then position\n");
  }

  size_t remove_sites (const arma::uvec& drop)
  {
    if (drop.n_elem == 0)
      return 0;

    if (drop.max() >= sites)
      Rcpp::stop("[GenotypeLikelihood::remove_sites] invalid indices");

    arma::uvec udrop = arma::unique(drop);
    size_t filtered = udrop.n_elem;

    like.shed_slices(udrop);
    cnts.shed_slices(udrop);
    post.shed_slices(udrop);
    missing.shed_cols(udrop);
    pos.shed_cols(udrop);
    sites -= udrop.n_elem;

    return filtered;
  }

  size_t retain_sites (const arma::uvec& keep)
  {
    if (keep.n_elem == 0)
      return 0;

    if (keep.max() >= sites)
      Rcpp::stop("[GenotypeLikelihood::retain_sites] invalid indices");

    arma::uvec ukeep = arma::unique(keep);
    size_t filtered = sites - ukeep.n_elem;

    like = like.slices(ukeep);
    cnts = cnts.slices(ukeep);
    post = post.slices(ukeep);
    missing = missing.cols(ukeep);
    pos = pos.cols(ukeep);
    sites = ukeep.n_elem;

    return filtered;
  }

  size_t polarize (const arma::uvec& ref)
  {
    if (ref.n_elem != sites)
      Rcpp::stop ("[GenotypeLikelihood::polarize] Reference sequence is wrong length");

    size_t polarized = 0;
    for (size_t i=0; i<sites; ++i)
      if (pos.at(3,i) == ref[i]) // minor == ref
      {
        polarized++;

        // flip counts, glf, pos
        pos.at(3,i) = pos.at(2,i);
        pos.at(2,i) = ref[i];
        like.slice(i) = arma::reverse(like.slice(i));
        cnts.slice(i) = arma::reverse(cnts.slice(i));
        post.slice(i) = arma::reverse(post.slice(i));
      }

    return polarized;
  }

  arma::vec::fixed<3> handler (const unsigned j, const size_t i) const
  {
    // return genotype likelihood for sample/site after processing for ploidy/missing data
    arma::vec::fixed<3> pj;

    if (missing(j, i))
      pj.fill(1./3.);
    else
    {
      pj[0] = like.at(0,j,i);
      pj[1] = like.at(1,j,i);
      pj[2] = like.at(2,j,i);
    }
    pj[1] *= double(ploidy[j] - 1);
    pj    /= arma::accu(pj);

    if (!arma::is_finite(pj)) // check for underflow that can happen with haploid "heterozygotes"
    {
      pj.fill(1./3.);
      pj[1] *= double(ploidy[j] - 1);
      pj    /= arma::accu(pj);
    }

    return pj;
  }
};

struct Genotypes : public RcppParallel::Worker
{ //DONE

  // estimate genotype frequencies

  const GenotypeLikelihood &GL;
  const arma::uvec site_index, sample_index;
  arma::mat freq;

  private:
  const arma::uvec &rsite_index, &rsample_index;
  arma::mat &rfreq;

  public:
  Genotypes (const GenotypeLikelihood& GL, const arma::uvec site_index, const arma::uvec sample_index)
    : GL (GL)
    , site_index (site_index)
    , sample_index (sample_index)
    , freq (site_index.n_elem, 5, arma::fill::zeros)
    // references accessed by slaves
    , rsite_index (site_index)
    , rsample_index (sample_index)
    , rfreq (freq)
  {
    if (site_index.max() > GL.sites)
      Rcpp::stop("[Genotypes] Invalid entries in site index");
    if (sample_index.max() > GL.samples)
      Rcpp::stop("[Genotypes] Invalid entries in sample index");

    RcppParallel::parallelFor(0, site_index.n_elem, *this);
    //(*this)(0, site_index.n_elem);
  }

  Genotypes (const Genotypes& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    // references accessed by slaves
    , rsite_index (rhs.rsite_index)
    , rsample_index (rhs.rsample_index)
    , rfreq (rhs.rfreq)
  {}

  double EM (arma::vec& eta, const size_t i)
  {
    const double tol = 1e-8;

    const size_t i2 = rsite_index[i];

    double loglik = 0;

    arma::vec upd (5, arma::fill::zeros),
                p (3, arma::fill::zeros),
              num (5, arma::fill::zeros);

    for (unsigned j = 0; j < rsample_index.n_elem; ++j)
    {
      const unsigned j2 = rsample_index[j];

      p = GL.handler (j2, i2);
      num.zeros();

      if (GL.ploidy[j2] == 1)
      {
        num[3] = p[0] * eta[3];
        num[4] = p[2] * eta[4];
      } 
      else if (GL.ploidy[j2] == 2)
      {
        num[0] = p[0] * eta[0];
        num[1] = p[1] * eta[1];
        num[2] = p[2] * eta[2];
      }
      double den = arma::accu(num);
      upd += num / den;
      loglik += log(den);
    }

    arma::vec sr2 = eta;

    eta  = upd / double(rsample_index.n_elem);
    sr2 -= eta;

    if (!arma::is_finite(sr2) || sqrt(arma::accu(arma::pow(sr2,2))) < tol)
      return true;

    return false;
  }

  void operator () (const size_t start, const size_t end)
  {
    const unsigned maxiter = 100;
    arma::vec eta (5);

    for (size_t i = start; i != end; ++i)
    {
      eta.fill(1./5.);
      for (unsigned iter = 0; iter < maxiter; ++iter)
        if (EM (eta, i)) break;
      rfreq.row(i) = eta.t() * double(rsample_index.n_elem);
    }
  }
};

struct Paralogs : public RcppParallel::Worker
{ // DONE
  public:
  const GenotypeLikelihood &GL;
  const arma::uvec site_index, sample_index;
  arma::mat diphet, haphet;

  private:
  const arma::uvec &rsite_index, &rsample_index;
  arma::mat &rdiphet, &rhaphet;

  public:
  Paralogs (const GenotypeLikelihood& GL, const arma::uvec site_index, const arma::uvec sample_index)
    : GL (GL)
    , site_index (site_index)
    , sample_index (sample_index)
    , diphet (0, 0)
    , haphet (0, 0)
    // references accessed by slaves
    , rsite_index (site_index)
    , rsample_index (sample_index)
    , rdiphet (diphet)
    , rhaphet (haphet)
  {
    if (site_index.max() > GL.sites)
      Rcpp::stop("Paralogs: invalid entries in site index");
    if (sample_index.max() > GL.samples)
      Rcpp::stop("Paralogs: invalid entries in sample index");

    // only allocate if haploids/diploids present
    if (arma::any(GL.ploidy.elem(sample_index) == 2))
    {
      diphet.set_size (site_index.n_elem, 3);
      diphet.fill(arma::datum::nan);
    }
    if (arma::any(GL.ploidy.elem(sample_index) == 1))
    {
      haphet.set_size (site_index.n_elem, 3);
      haphet.fill(arma::datum::nan);
    }

    RcppParallel::parallelFor(0, site_index.n_elem, *this);
    //(*this)(0, site_index.n_elem);
  }

  Paralogs (const Paralogs& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    // references accessed by slaves
    , rsite_index (rhs.rsite_index)
    , rsample_index (rhs.rsample_index)
    , rdiphet (rhs.rdiphet)
    , rhaphet (rhs.rhaphet)
  {}

  void operator() (const size_t start, const size_t end)
  {
    arma::uvec keep (GL.samples);
    arma::fvec tmp  (GL.samples);
    for (size_t i=start; i!=end; ++i)
    {
      unsigned i2 = rsite_index[i];
      keep.zeros();
      for (unsigned j=0; j<rsample_index.n_elem; ++j) 
      {
        unsigned j2 = rsample_index[j];
        if (!GL.missing.at(j2,i2))
          keep[j2] = GL.ploidy[j2];
      }

      tmp = arma::trans(GL.like.slice(i2).row(1));
      if (arma::any(keep == 1))
      { 
        PoissonBinomial poibin (tmp.elem(arma::find(keep == 1)));
        rhaphet.at(i,0) = double(poibin.dim);
        rhaphet.at(i,1) = poibin.expectation();
        rhaphet.at(i,2) = poibin.cdf(0.);
      }
      if (arma::any(keep == 2))
      {
        PoissonBinomial poibin (tmp.elem(arma::find(keep == 2)));
        rdiphet.at(i,0) = double(poibin.dim);
        rdiphet.at(i,1) = poibin.expectation();
        rdiphet.at(i,2) = poibin.cdf(0.5);
      }
    }
  }

  struct PoissonBinomial
  {
    const unsigned dim;
    const arma::fvec PDF;

    PoissonBinomial (const arma::fvec& P)
      : dim (P.n_elem)
      , PDF (pdf(P))
    {}

    arma::fvec pdf (const arma::fvec& P)
    {
      if (arma::any(P < 0. || P > 1.))
        Rcpp::stop ("PoissonBinomial: probabilities must be in [0,1]\n");

      // convolution algorithm from Biscarri
      arma::fvec out (P.n_elem+1, arma::fill::zeros);
      arma::fvec::fixed<2> p;
      if (out.n_elem == 1)
        out(0) = 1.;
      else
      {
        out(0) = 1.-P(0); out(1) = P(0);
        for (unsigned i=1; i<P.n_elem; ++i)
        {
          p(0) = 1.-P(i); p(1) = P(i);
          for (int j=i; j>=0; --j)
            out(j+1) = out(j) * p(1) + out(j+1) * p(0);
          out(0) *= p(0);
        }
      }
      return out;
    }

    double cdf (const double& thresh)
    {
      // always round threshold down to nearest integer:
      // [0,1,2,3,4,5,6] -> heterozygote majority would be 3 or more
      // [0,1,2,3,4,5] -> heterozygote majority would be 3 or more
      return arma::accu(PDF.subvec(floor(thresh * double(dim)), dim));
    }

    double expectation (void)
    {
      double out = 0.;
      for (unsigned i=0; i<PDF.n_elem; ++i)
        out += double(i) * PDF[i];
      return out;
    }
  };
};

struct Frequencies : public RcppParallel::Worker
{ //DONE
  
  // This struct estimates allele frequencies from genotype likelihoods using the EM algorithm described in Li, H. 2011. Bioinformatics
  // It is essentially the same as that implemented in ANGSD, with the exception that haplodiploid samples are allowed.

  public:
  const GenotypeLikelihood &GL;
  const arma::uvec site_index, sample_index;
  const unsigned chromosomes;
  arma::vec freq, lrt, pval;

  private:
  const arma::uvec &rsite_index, &rsample_index;
  arma::vec &rfreq, &rlrt, &rpval;
  
  public:
  Frequencies (const GenotypeLikelihood& GL, const arma::uvec site_index, const arma::uvec sample_index)
    : site_index (site_index)
    , sample_index (sample_index)
    , GL (GL)
    , chromosomes (arma::accu(GL.ploidy.elem(sample_index)))
    , freq (site_index.n_elem)
    , lrt (site_index.n_elem)
    , pval (site_index.n_elem)
    // references that slave interacts with
    , rsite_index (site_index)
    , rsample_index (sample_index)
    , rfreq (freq)
    , rlrt (lrt)
    , rpval (pval)
  {
    if (arma::max(site_index) >= GL.sites)
      Rcpp::stop ("Frequencies: invalid entries in site index");
    if (arma::max(sample_index) >= GL.samples)
      Rcpp::stop ("Frequencies: invalid entries in sample index");

    freq.fill(arma::datum::nan);
    lrt.fill(arma::datum::nan);
    pval.fill(arma::datum::nan);

    RcppParallel::parallelFor(0, site_index.n_elem, *this);
    //(*this)(0, site_index.n_elem);
  }

  Frequencies (const Frequencies& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , chromosomes (rhs.chromosomes)
    // references that slave interacts with
    , rsite_index (rhs.site_index)
    , rsample_index (rhs.sample_index)
    , rfreq (rhs.rfreq)
    , rlrt (rhs.rlrt)
    , rpval (rhs.rpval)
  {}

  void operator() (const size_t start, const size_t end)
  {
    // should only interact with references
    double loglik;
    for (size_t i = start; i != end; ++i)
    {
      // check that there's data
      unsigned missing = 0;
      for (auto j2 : rsample_index)
        missing += GL.missing.at(j2, rsite_index[i]) * GL.ploidy[j2];

      // estimate allele frequencies via EM
      if (chromosomes - missing > 0) 
      {
        loglik = EM (rfreq[i], i);
        // test against null of monomorphic ancestral or monomorphic fixed
        rlrt[i]  = -2.*null(i, rfreq[i]>0.5) + 2.*loglik;
        rpval[i] = R::pchisq(rlrt[i], 1., false, false); 
      }
    }
  }

  double null (const size_t i, bool fixed)
  {
    // likelihood under null model (freq = 0 || 1)
    double loglik = 0;
    for (unsigned j=0; j<rsample_index.n_elem; ++j)
    {
      unsigned j2 = rsample_index[j],
               i2 = rsite_index[i];
      if (!GL.missing(j2,i2))
        loglik += log(GL.like.at(int(fixed)*2,j2,i2));
    }
    return loglik;
  }

  bool freaky (double& f, double& loglik, const size_t i)
  {
    // missing samples have likelihood set to zero, so do not increment num/den
    // just remember to divide by the number of nonmissing samples!
    const double tol  = 0.00001, //same as in angsd
                 rtol = 0.;

    const unsigned i2 = rsite_index[i];

    double num = 0.,
           den = 0.,
           of  = f;
    unsigned n = 0;

    f      = 0.; 
    loglik = 0.;

    for (unsigned j=0; j<rsample_index.n_elem; ++j)
    {
      const unsigned j2 = rsample_index[j];
      if (!GL.missing(j2,i2))
      {
        if (GL.ploidy[j2] == 1)
        {
          num = GL.like.at(2,j2,i2) * of;      
          den = GL.like.at(2,j2,i2) * of + 
                GL.like.at(0,j2,i2) * (1.-of);
        } 
        else if (GL.ploidy[j2] == 2)
        {
          // sum_{samples} sum_{0 \leq g \leq ploidy} g lik(g) dbinom(g, ploidy, freq) / (nsamples sum_{samples} sum_{0 \leq g \leq ploidy} lik(g) dbinom(g, ploidy, freq))
          num = 2. * GL.like.at(2,j2,i2) * 1. * of * of +     // hom
                1. * GL.like.at(1,j2,i2) * 2. * of * (1.-of); // het
          den = GL.like.at(2,j2,i2) * 1. * of * of +          // hom
                GL.like.at(1,j2,i2) * 2. * of * (1.-of) +     // het
                GL.like.at(0,j2,i2) * 1. * (1.-of) * (1.-of); // hom
        }
        n      += GL.ploidy[j2];
        f      += num/den;
        loglik += log(den);
      }
    }
    f /= double(n);

    return (f - of < tol && of - f < tol) || (f/of < 1.+rtol && f/of > 1.-rtol) || !arma::is_finite(f);
  }

  double EM (double& freak, const size_t i)
  {
    // starting conditions, iters, stoping conditions are same as angsd
    const unsigned maxiter = 100; 
    double loglik = 0;
    freak = 0.001; // could use better guess
    for (unsigned iter=0; iter<maxiter; ++iter)
      if (freaky(freak, loglik, i))
        break;
    return loglik;
  }
};

struct Admixture : public RcppParallel::Worker
{ 
  public:
  const GenotypeLikelihood &GL;				// contains likelihoods, ploidy, dimensions
  const arma::uvec sample_index,      // only use these samples
                   site_index; 			  // only use these sites
  const unsigned K;   								// number of clusters, number of chromosomes to hold out
  arma::mat Qmat, Fmat, loglik;       // admixture, frequencies, per site loglik and cross-validation
  arma::cube Acoef, Bcoef;						// per-site coefficients (could move to private eventually)
	double loglikelihood = 0.;					// total loglikelihood
  unsigned iter = 0;

  private:
  bool   acceleration = false;       // use EM acceleration?
  double errtol = 1e-8;              // not using dynamic bounds
  double stepmax = 1., stepmin = 1.; // adaptive steplength

	// references for slaves
  const arma::uvec &rsite_index, &rsample_index;
  const arma::mat &rQmat, &rFmat;
  arma::mat &rloglik;
  arma::cube &rAcoef, &rBcoef;

	//temporaries
	arma::mat Q0, Q1, Q2, F0, F1, F2, Qd1, Qd2, Qd3, Fd1, Fd2, Fd3;

  public:
  Admixture (const GenotypeLikelihood& GL, const arma::uvec site_index, const arma::uvec sample_index, const arma::mat& Qstart, const bool fixQ)
    : GL (GL)
    , site_index (site_index)
    , sample_index (sample_index)
    , K (Qstart.n_rows)
    , Qmat (K, sample_index.n_elem)
    , Fmat (K, site_index.n_elem, arma::fill::randu)
    , loglik (sample_index.n_elem, site_index.n_elem)
    , Acoef (K, sample_index.n_elem, site_index.n_elem)
    , Bcoef (K, sample_index.n_elem, site_index.n_elem)
    // read-write references for workers 
    , rsite_index (site_index)
    , rsample_index (sample_index)
    , rQmat (Qmat)
    , rFmat (Fmat)
    , rloglik (loglik)
    , rAcoef (Acoef)
    , rBcoef (Bcoef)
		// temporaries
		, Q0 (arma::size(Qmat), arma::fill::zeros)
		, Q1 (arma::size(Qmat), arma::fill::zeros)
		, Q2 (arma::size(Qmat), arma::fill::zeros)
		, Qd1 (arma::size(Qmat), arma::fill::zeros)
		, Qd2 (arma::size(Qmat), arma::fill::zeros)
		, Qd3 (arma::size(Qmat), arma::fill::zeros)
		, F0 (arma::size(Fmat), arma::fill::zeros)
		, F1 (arma::size(Fmat), arma::fill::zeros)
		, F2 (arma::size(Fmat), arma::fill::zeros)
		, Fd1 (arma::size(Fmat), arma::fill::zeros)
		, Fd2 (arma::size(Fmat), arma::fill::zeros)
		, Fd3 (arma::size(Fmat), arma::fill::zeros)
  {
		// everything here should interface with originals not references
    if (K < 1)
      Rcpp::stop ("[Admixture] must have at least one cluster");
    if (K > 999)
      Rcpp::stop ("[Admixture] number of clusters capped at 999");
    if (arma::max(site_index) >= GL.sites)
      Rcpp::stop ("[Admixture] Invalid entries in site index");
    if (arma::max(sample_index) >= GL.samples)
      Rcpp::stop ("[Admixture] Invalid entries in sample index");
    if (Qstart.n_rows != Qmat.n_rows || Qstart.n_cols != Qmat.n_cols)
      Rcpp::stop ("[Admixture] Starting values for Q are of wrong dimension");

    const unsigned maxiter = 5000;

    Qmat = Qstart;
    Qmat.each_col([](arma::vec& x) { x /= arma::accu(x); });

    Fmat.randu();
    Fmat = arma::clamp(Fmat, errtol, 1.-errtol);
		
		for (iter=0; iter<maxiter; ++iter)
    {
      bool converged = EMaccel(fixQ ? Q0 : Qmat, Fmat); // inputs get updated, but only Qmat is used in likelihood calc (Q0 is like /dev/null)
      //TODO make mutables clearer; should always be passed as arguments (e.g. loglikelihood, cvscore!
      if((iter+1) % 10 == 0) 
        fprintf(stderr, "[Admixture] Iteration %u, loglik = %f\n", iter+1, loglikelihood);
			if(converged)
				break;
    }
  
		if (iter == maxiter)
			Rcpp::warning("[Admixture] Did not converge in maximum number of iterations");
  }

  Admixture (const Admixture& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , K (rhs.K)
    // references that interface with slaves
    , rsite_index (rhs.rsite_index)
    , rsample_index (rhs.rsample_index)
    , rQmat (rhs.rQmat)
    , rFmat (rhs.rFmat)
    , rloglik (rhs.rloglik)
    , rAcoef (rhs.rAcoef)
    , rBcoef (rhs.rBcoef)
  {}

  // TODO: I would prefer to simplify so that only arguments are mutables or indices (e.g. get rid of Q/F in below)
  void Diploid (arma::cube& A, arma::cube& B, arma::mat& ll, const size_t i, const unsigned j)
  {
    // everything in here should interact with references
    // j = sample, i = site, k = cluster
    unsigned j2 = rsample_index[j],
             i2 = rsite_index[i];

    // deal with missing data or holdout set
    arma::vec::fixed<3> p;
    if (GL.missing.at(j2, i2))
      p.fill(1./3.);
    else
    {
      p[0] = GL.like(0,j2,i2);
      p[1] = GL.like(1,j2,i2);
      p[2] = GL.like(2,j2,i2);
    }
    p /= arma::accu(p);

    // check
    const double h  = arma::accu(Qmat.col(j) % Fmat.col(i)),
                 hn = arma::accu(Qmat.col(j) % (1.-Fmat.col(i)));
    const arma::vec::fixed<3> hwe = {p[0] * (1.-h) * (1.-h),
     																 p[1] * 2. * h * (1.-h),
     																 p[2] * h * h};
    double num = 0., den = 0.;
    for (unsigned g=0; g<3; ++g)
    {
      num += double(g) * hwe[g];
      den += hwe[g];
    }
			
    ll.at(j,i) = log(den);
    for (unsigned k=0; k<K; ++k)
    {
      A.at(k,j,i) = 0.5 * num/den * Qmat.at(k,j)*Fmat.at(k,i)/h;
      B.at(k,j,i) = 0.5 * (2.-num/den) * Qmat.at(k,j)*(1.-Fmat.at(k,i))/hn;
    }
  }

  void Haploid (arma::cube& A, arma::cube& B, arma::mat& ll, const size_t i, const unsigned j)
  {
    // everything in here should interact with references
    // j = sample, i = site, k = cluster
    unsigned j2 = rsample_index[j],
             i2 = rsite_index[i];

    // TODO: if converting to GL.handler, need to change p[1] to p[2]
    // deal with missing data and holdouts
    arma::vec::fixed<2> p;
    if (GL.missing.at(j2, i2))
      p.fill(1./2.);
    else
    {
      p[0] = GL.like.at(0,j2,i2);
      p[1] = GL.like.at(2,j2,i2);
    }
    p /= arma::accu(p);

    //check
    const double h  = arma::accu(Qmat.col(j) % Fmat.col(i)),
                 hn = arma::accu(Qmat.col(j) % (1.-Fmat.col(i)));
    double den = p[0] * hn + p[1] * h;

		ll.at(j,i) = log(den);
		for (unsigned k=0; k<K; ++k)
		{
			A.at(k,j,i) = 0.5 * p[1] * Qmat.at(k,j)*Fmat.at(k,i)/den;
			B.at(k,j,i) = 0.5 * p[0] * Qmat.at(k,j)*(1.-Fmat.at(k,i))/den;
		}
  }

  void operator () (const size_t start, const size_t end)
  {
    // every large container passed should be a reference
    for (size_t i=start; i!=end; ++i)
    {
      //TODO handle totally missing data? Or totally fixed sites? Necessary?
      for (unsigned j=0; j<rsample_index.n_elem; ++j)
      {
        unsigned j2 = rsample_index[j];
        if (GL.ploidy[j2] == 1)
          Haploid(rAcoef, rBcoef, rloglik, i, j);
        else if (GL.ploidy[j2] == 2)
          Diploid(rAcoef, rBcoef, rloglik, i, j);
      }
    }
  }

  double EM (arma::mat& Q, arma::mat& F)
  {
    // this should interact with the original containers
    // this uses Qmat/Fmat in to calculate likelihood updates, and stores the output in Q/F

    loglik.zeros();
    Acoef.zeros();
    Bcoef.zeros();

    RcppParallel::parallelFor(0, site_index.n_elem, *this); 
    //(*this)(0, site_index.n_elem);
    
    if (!(arma::is_finite(Acoef) && arma::is_finite(Bcoef)))
      Rcpp::stop("[Admixture] Numeric underflow, remove nonsegregating sites");

    F = arma::sum(Acoef, 1) / (arma::sum(Acoef, 1) + arma::sum(Bcoef, 1));
    Q = (arma::sum(Acoef, 2) + arma::sum(Bcoef, 2));
    Q.each_col([](arma::vec& x) { x /= arma::accu(x); }); // make sum_k q_{jk} = 1

    return arma::accu(loglik); 
  }

  void project (arma::mat& Q, arma::mat& F)
  {
    F = arma::clamp(F, errtol, 1.-errtol);
    Q = arma::clamp(Q, errtol, arma::datum::inf);
    Q.each_col([](arma::vec& x) { x /= arma::accu(x); }); // make sum_k q_{jk} = 1
  }

  bool EMaccel (arma::mat& Q, arma::mat& F)
  {
    // this should interact with the original containers
    // this uses Qmat/Fmat in to calculate likelihood updates, and stores the output in Q/F

		const double tol = 1e-5;
		const double mstep = 4;

    double ll0 = loglikelihood;

    Q0 = Q;	F0 = F;

		// first secant condition
    loglikelihood = EM (Q, F); project (Q, F);
	  Q1  = Q;		  	F1  = F;
		Qd1 = Q - Q0;   Fd1 = F - F0;
    double ll1 = loglikelihood;
	  double sr2 = arma::accu(arma::pow(Qd1,2)) + arma::accu(arma::pow(Fd1,2));
	  if (!arma::is_finite(sr2) || fabs(ll1 - ll0) < tol || sqrt(sr2) < tol)
			return true;
    
    if (!acceleration) // vanilla EM
      return false;

		// second secant condition
    loglikelihood = EM (Q, F); project (Q, F);
	  Q2  = Q;		  	F2  = F;
		Qd2 = Q - Q1;   Fd2 = F - F1;
    double em  = loglikelihood;
	  double sq2 = arma::accu(arma::pow(Qd2,2)) + arma::accu(arma::pow(Fd2,2));
	  if (!arma::is_finite(sq2) || fabs(em - ll1) < tol || sqrt(sq2) < tol)
			return true;

		// the magic happens
		Qd3 = Qd2 - Qd1; Fd3 = Fd2 - Fd1;
	  double sv2 = arma::accu(arma::pow(Qd3,2)) + arma::accu(arma::pow(Fd3,2));
	  double alpha = sqrt(sr2/sv2);
		alpha = std::max(stepmin, std::min(stepmax, alpha));

    F = F0 + 2.*alpha*Fd1 + alpha*alpha*Fd3;
    Q = Q0 + 2.*alpha*Qd1 + alpha*alpha*Qd3;
    project (Q, F);

		// stabilize
		loglikelihood = EM (Q, F); project (Q, F);

    // revert to simple EM iteration if loglikelihood has decreased after acceleration
    if (!arma::is_finite(loglikelihood) || loglikelihood < em)
    {
      loglikelihood = em;
      Q = Q2;
      F = F2;
    }

		if (alpha == stepmax)
			stepmax *= mstep;
		
		return false;
  }
};

struct AdmixtureHybrid : public RcppParallel::Worker
{ 
  // This supercedes Admixture, remove that
  
  struct Microsat
  {
    const MultiAllelicGenotypeLikelihood &MAGL;
    const arma::uvec sample_index, eval_msat;
    const arma::cube L;
    const arma::uword num_alleles, num_samples, num_clusters;

    arma::mat Fmat, Emat, Smat,          // frequencies, expected frequencies, expected counts
              F0, F1, F2, Fd1, Fd2, Fd3; // temporaries for acceleration
    arma::vec loglikelihood;

    Microsat (const MultiAllelicGenotypeLikelihood& MAGL, const arma::uvec sample_index, const arma::uvec eval_msat, arma::mat Q, const arma::cube like)
      : MAGL (MAGL)
      , sample_index (sample_index)
      , eval_msat (eval_msat)
      , L (like)
      , num_alleles (like.n_rows)
      , num_samples (sample_index.n_elem)
      , num_clusters (Q.n_rows)
      , loglikelihood (sample_index.n_elem)
    {
      if (Q.n_cols      != num_samples             ||
          num_samples   != eval_msat.n_elem        ||
          MAGL.samples  != like.n_slices           ||
          like.n_cols   != num_alleles             ||
          arma::max(sample_index) >= like.n_slices )
        Rcpp::stop ("[Microsat] Dimension mismatch");

      Q.each_col([](arma::vec& x) { x /= arma::accu(x); });
      //Fmat = arma::randu<arma::mat>(num_alleles, num_clusters);
      Fmat = arma::ones<arma::mat>(num_alleles, num_clusters);
      Fmat.each_col([](arma::vec& x) { x /= arma::accu(x); });
      Emat = expected_frequencies(Q, Fmat);
      Smat = expected_counts(Q, Fmat, Emat);
    }

    double EM (const arma::mat& Q, arma::mat& Qup)
    {
      if (Q.n_cols   != num_samples  ||
          Q.n_rows   != num_clusters ||
          Qup.n_cols != num_samples  ||
          Qup.n_rows != num_clusters )
        Rcpp::stop ("[Microsat] Dimension mismatch");

      arma::mat F = Fmat;
      Smat = expected_frequencies(Q, F);
      Emat = expected_counts(Q, F, Smat);
      Fmat = update_frequencies(Q, F, Smat, Emat);
      Qup += update_admixture(Q, F, Smat, Emat);

      loglikelihood = loglik(Q, F, Smat);

      return arma::accu(loglikelihood);
    }

    arma::mat expected_frequencies (const arma::mat& Q, const arma::mat& F)
    {
      // compute HWE frequencies per sample
      arma::mat out (num_alleles, num_samples, arma::fill::zeros);
      for (unsigned s=0; s<num_samples; ++s) //samples
      {
        if (eval_msat[s])
        {
          for (unsigned k=0; k<num_clusters; ++k) //clusters
            for (unsigned i=0; i<num_alleles; ++i) //alleles
              out.at(i,s) += Q.at(k,s) * F.at(i,k);
        }
      }

      return out;
    }

    arma::mat expected_counts (const arma::mat& Q, const arma::mat& F, const arma::mat& S)
    {
      // posterior expectation of allele counts for each sample
      arma::mat out (num_alleles, num_samples, arma::fill::zeros);
      for (unsigned s=0; s<num_samples; ++s) //samples
      {
        if (eval_msat[s])
        {
          double den = 0., p = 0.;
          arma::uword s2 = sample_index[s];
          if (MAGL.ploidy.at(s2)==2)
          {
            for (unsigned i=0; i<num_alleles; ++i) //alleles
              for (unsigned j=i; j<num_alleles; ++j) //alleles again
              {
                if (i==j)
                {
                  p = L.at(i,i,s2) * std::pow(S.at(i,s), 2);
                  den += p;
                  out.at(i,s) += 2. * p;
                }
                else
                {
                  p = L.at(i,j,s2) * S.at(j,s) * S.at(i,s);
                  den += 2. * p;
                  // 8/30/20 adding the factor of 2 below seems to have fixed the EM updates
                  out.at(i,s) += 2. * p;
                  out.at(j,s) += 2. * p;
                }
              }
          } 
          else 
          {
            for (unsigned i=0; i<num_alleles; ++i) //alleles
            {
              p = L.at(i,i,s2) * S.at(i,s);
              den += p;
              out.at(i,s) += p;
            }
          }
          out.col(s) /= 2.*den; //scales haploid expectation to 0.5 as for snps
        }
      }
      return out;
    }

    arma::mat update_admixture (const arma::mat& Q, const arma::mat& F, const arma::mat& S, const arma::mat& E)
    {
      // contribution to admixture coeff update 
      arma::mat Qup = arma::zeros<arma::mat>(arma::size(Q));
      for (unsigned s=0; s<num_samples; ++s) // samples
      {
        if (eval_msat[s])
        {
          for (unsigned k=0; k<num_clusters; ++k) // clusters
            for (unsigned i=0; i<num_alleles; ++i) // alleles
              Qup.at(k,s) += E.at(i,s) * Q.at(k,s) * F.at(i,k) / S.at(i,s);
        }
      }
      return Qup;
    }

    arma::mat update_frequencies (const arma::mat& Q, const arma::mat& F, const arma::mat& S, const arma::mat& E)
    {
      // compute new frequencies 
      arma::mat Fup = arma::zeros<arma::mat>(arma::size(F));
      for (unsigned k=0; k<num_clusters; ++k) // clusters
      {
        arma::vec A (num_alleles, arma::fill::zeros);
        for (unsigned i=0; i<num_alleles; ++i)
        {
          for (unsigned s=0; s<num_samples; ++s)
          {
            if (eval_msat[s])
            {
              A.at(i) += E.at(i,s) * Q.at(k,s) * F.at(i,k) / S.at(i,s);
            }
          }
        }
        Fup.col(k) = A / arma::accu(A);
      }

      return Fup;
    }

    arma::vec loglik (const arma::mat& Q, const arma::mat& F, const arma::mat& S)
    {
      // calculate per-sample loglikelihood
      arma::vec ll (num_samples, arma::fill::zeros);

      for (unsigned s=0; s<num_samples; ++s)
      {
        if (eval_msat[s])
        {
          double den = 0., p = 0.;
          arma::uword s2 = sample_index[s];
          if (MAGL.ploidy.at(s2)==2)
          {
            for (unsigned i=0; i<num_alleles; ++i)
              for (unsigned j=i; j<num_alleles; ++j)
              {
                if (i == j)
                {
                  p = L.at(i,i,s2) * std::pow(S.at(i,s), 2);
                  den += p;
                }
                else
                {
                  p = L.at(i,j,s2) * S.at(j,s) * S.at(i,s);
                  den += 2. * p;
                }
              }
          } 
          else 
          {
            for (unsigned i=0; i<num_alleles; ++i)
            {
              p = L.at(i,i,s2) * S.at(i,s);
              den += p;
            }
          }
          ll[s] = log(den);
        }
      }
      return ll;
    }
  };
  
  const GenotypeLikelihood &GL;				        // contains SNP likelihoods, ploidy, dimensions
  const MultiAllelicGenotypeLikelihood &MAGL; // contains microsat likelihoods
  const arma::uvec sample_index,              // only use these samples
                   site_index;                // only use these sites
  const unsigned K;   								        // number of clusters, number of chromosomes to hold out
  arma::mat D, Qdeme, Qmat, Fmat, loglik;     // admixture, frequencies, per site loglik
  arma::cube Acoef, Bcoef;						        // per-site coefficients (could move to private eventually)
  std::vector<Microsat> Msat;                 // microsatellite genotype likelihoods and frequencies
	double loglikelihood = 0.;					        // total loglikelihood
  unsigned iter = 0;

  private:
  bool   fix_Fmat = false,
         fix_Qmat = false,
         verbose  = true;
  bool   acceleration = true;        // use EM acceleration?
  double errtol = 1e-16;             // not using dynamic bounds
  double stepmax = 1., stepmin = 1.; // adaptive steplength

  arma::uvec eval_msat,
             eval_snp;

	// references for workers
  const arma::uvec &rsample_index, &rsite_index; 
  arma::uvec &reval_snp;
  const arma::mat &rQmat, &rFmat;
  arma::mat &rloglik;
  arma::cube &rAcoef, &rBcoef;

	//temporaries
	arma::mat Q0, Q1, Q2, F0, F1, F2, Qd1, Qd2, Qd3, Fd1, Fd2, Fd3;

  public:
  AdmixtureHybrid (const GenotypeLikelihood& GL, const MultiAllelicGenotypeLikelihood& MAGL, 
                   const arma::uvec site_index, const arma::uvec sample_index, 
                   const arma::uvec deme_index, const arma::uvec marker_type, const arma::mat Qstart, const arma::mat Fstart, const unsigned maxiter)
    : GL (GL)
    , MAGL (MAGL)
    , sample_index (sample_index) 
    , site_index (site_index)
    , K (Qstart.n_rows)
    , D (arma::max(deme_index) + 1, sample_index.n_elem)
    , Qdeme (K, arma::max(deme_index) + 1)
    , Qmat (K, sample_index.n_elem)
    , Fmat (K, site_index.n_elem)
    , loglik (sample_index.n_elem, site_index.n_elem)
    , Acoef (K, sample_index.n_elem, site_index.n_elem)
    , Bcoef (K, sample_index.n_elem, site_index.n_elem)
    , eval_msat (sample_index.n_elem, arma::fill::zeros)
    , eval_snp (sample_index.n_elem, arma::fill::zeros)
    // read-write references for workers 
    , rsample_index (sample_index)
    , rsite_index (site_index)
    , reval_snp (eval_snp)
    , rQmat (Qmat)
    , rFmat (Fmat)
    , rloglik (loglik)
    , rAcoef (Acoef)
    , rBcoef (Bcoef)
		// temporaries
		, Q0 (arma::size(Qdeme), arma::fill::zeros)
		, Q1 (arma::size(Qdeme), arma::fill::zeros)
		, Q2 (arma::size(Qdeme), arma::fill::zeros)
		, Qd1 (arma::size(Qdeme), arma::fill::zeros)
		, Qd2 (arma::size(Qdeme), arma::fill::zeros)
		, Qd3 (arma::size(Qdeme), arma::fill::zeros)
		, F0 (arma::size(Fmat), arma::fill::zeros)
		, F1 (arma::size(Fmat), arma::fill::zeros)
		, F2 (arma::size(Fmat), arma::fill::zeros)
		, Fd1 (arma::size(Fmat), arma::fill::zeros)
		, Fd2 (arma::size(Fmat), arma::fill::zeros)
		, Fd3 (arma::size(Fmat), arma::fill::zeros)
  {
		// everything here should interface with originals not references
    if (K < 1)
      Rcpp::stop ("[AdmixtureHybrid] Must have at least one cluster");
    if (K > 999)
      Rcpp::stop ("[AdmixtureHybrid] Number of clusters capped at 999");
    if (arma::max(site_index) >= GL.sites)
      Rcpp::stop ("[AdmixtureHybrid] Invalid entries in site index");
    if (arma::max(sample_index) >= GL.samples)
      Rcpp::stop ("[AdmixtureHybrid] Invalid entries in sample index");
    if (Qstart.n_rows != Qdeme.n_rows || Qstart.n_cols != Qdeme.n_cols)
      Rcpp::stop ("[AdmixtureHybrid] Starting values for Q are of wrong dimension");
    if (Fstart.n_rows != Qdeme.n_rows || Fstart.n_cols != site_index.n_elem)
      Rcpp::stop ("[AdmixtureHybrid] Starting values for F are of wrong dimension");
    if (deme_index.n_elem != sample_index.n_elem)
      Rcpp::stop ("[AdmixtureHybrid] Length of deme index does not match length of sample index");
    if (marker_type.n_elem != sample_index.n_elem)
      Rcpp::stop ("[AdmixtureHybrid] Length of marker type does not match length of sample index");
    if (arma::max(deme_index) >= deme_index.n_elem)
      Rcpp::stop ("[AdmixtureHybrid] Deme indices must be 0-based contiguous integers");
    if (arma::max(marker_type) >= 3)
      Rcpp::stop ("[AdmixtureHybrid] Marker type must be in [0, 1, 2]");
    if (Fstart.max() > 1.0 || Fstart.min() < 0.0)
      Rcpp::stop ("[AdmixtureHybrid] Starting allele frequencies must be in [0,1]");
    if (Qstart.min() < 0.0)
      Rcpp::stop ("[AdmixtureHybrid] Starting admixture proportions must be positive");

    fprintf(stderr,
        "[AdmixtureHybrid] [verbose=%d fix_Qmat=%d fix_Fmat=%d acceleration=%d errtol=%f stepmax=%f stepmin=%f maxiter=%d]\n",
        verbose, fix_Qmat, fix_Fmat, acceleration, errtol, stepmax, stepmin, maxiter);

    eval_snp.elem(arma::find(marker_type <= 1)).ones();
    eval_msat.elem(arma::find(marker_type >= 1)).ones();

    // indicator matrix for samples in demes
    D.zeros();
    for (unsigned i = 0; i < sample_index.n_elem; ++i)
      D.at(deme_index[i], i) = 1.;

    // initialize
    Qdeme = Qstart;
    Qdeme.each_col([](arma::vec& x) { x /= arma::accu(x); });

    //Fmat.randu();
    //Fmat = arma::clamp(Fmat, errtol, 1.-errtol);
    Fmat = Fstart;

    for (auto& ms : MAGL.like)
      Msat.emplace_back(MAGL, sample_index, eval_msat, Qmat, ms);
		
    // EM
		for (iter=0; iter<maxiter; ++iter)
    {
      bool converged = EMaccel(Qdeme, Fmat, Msat);
      if((iter+1) % 100 == 0) 
      {
        fprintf(stderr, "[AdmixtureHybrid] [%u] log-likelihood = %f\n", iter+1, loglikelihood);
        if (verbose)
          Qdeme.raw_print("Dumping admixture coefficients:");
      }
			if(converged)
				break;
    }
  
		if (iter == maxiter)
			Rcpp::warning("[AdmixtureHybrid] Did not converge in maximum number of iterations");
    else
      fprintf(stderr, "[AdmixtureHybrid] Converged in %u iterations\n", iter);
  }

  AdmixtureHybrid (const AdmixtureHybrid& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , MAGL (rhs.MAGL)
    , K (rhs.K)
    // references that interface with slaves
    , rsite_index (rhs.rsite_index)
    , rsample_index (rhs.rsample_index)
    , reval_snp (rhs.reval_snp)
    , rQmat (rhs.rQmat)
    , rFmat (rhs.rFmat)
    , rloglik (rhs.rloglik)
    , rAcoef (rhs.rAcoef)
    , rBcoef (rhs.rBcoef)
  {}

  void Diploid (arma::cube& A, arma::cube& B, arma::mat& ll, const size_t i, const unsigned j)
  {
    // everything in here should interact with references
    // j = sample, i = site, k = cluster
    unsigned j2 = rsample_index[j],
             i2 = rsite_index[i];

    arma::vec::fixed<3> p = GL.handler(j2,i2);

    // check
    const double h  = arma::accu(Qmat.col(j) % Fmat.col(i)),
                 hn = arma::accu(Qmat.col(j) % (1.-Fmat.col(i)));
    const arma::vec::fixed<3> hwe = {p[0] * (1.-h) * (1.-h),
     																 p[1] * 2. * h * (1.-h),
     																 p[2] * h * h};
    double num = 0., den = 0.;
    for (unsigned g=0; g<3; ++g)
    {
      num += double(g) * hwe[g];
      den += hwe[g];
    }
			
    ll.at(j,i) = log(den);
    for (unsigned k=0; k<K; ++k)
    {
      A.at(k,j,i) = 0.5 * num/den * Qmat.at(k,j)*Fmat.at(k,i)/h;
      B.at(k,j,i) = 0.5 * (2.-num/den) * Qmat.at(k,j)*(1.-Fmat.at(k,i))/hn;
    }
  }

  void Haploid (arma::cube& A, arma::cube& B, arma::mat& ll, const size_t i, const unsigned j)
  {
    // everything in here should interact with references
    // j = sample, i = site, k = cluster
    unsigned j2 = rsample_index[j],
             i2 = rsite_index[i];

    arma::vec::fixed<3> p = GL.handler(j2,i2);

    //check
    const double h  = arma::accu(Qmat.col(j) % Fmat.col(i)),
                 hn = arma::accu(Qmat.col(j) % (1.-Fmat.col(i)));
    double den = p[0] * hn + p[2] * h;

		ll.at(j,i) = log(den);
		for (unsigned k=0; k<K; ++k)
		{
			A.at(k,j,i) = 0.5 * p[2] * Qmat.at(k,j)*Fmat.at(k,i)/den;
			B.at(k,j,i) = 0.5 * p[0] * Qmat.at(k,j)*(1.-Fmat.at(k,i))/den;
      // this is equivalent to:
      //   0.5 * (p[2]*h)/den * Qmat.at(k,j)*Fmat.at(k,i)/h
      // which is analogous to the diploid case
		}
  }

  void operator () (const size_t start, const size_t end)
  {
    // every large container passed should be a reference
    for (size_t i=start; i!=end; ++i)
    {
      for (unsigned j=0; j<rsample_index.n_elem; ++j)
      {
        if (eval_snp[j])
        {
          unsigned j2 = rsample_index[j];
          if (GL.ploidy[j2] == 1)
            Haploid(rAcoef, rBcoef, rloglik, i, j);
          else if (GL.ploidy[j2] == 2)
            Diploid(rAcoef, rBcoef, rloglik, i, j);
        }
      }
    }
  }

  double EM (arma::mat& Q, arma::mat& F, std::vector<Microsat>& Ms)
  {
    loglik.zeros();
    Acoef.zeros();
    Bcoef.zeros();

    Qmat = Q * D;
    Fmat = F;

    arma::mat Qmat_update = arma::zeros<arma::mat>(arma::size(Qmat));
    double    loglik_msat = 0.;
    for (auto& ms : Ms)
      loglik_msat += ms.EM(Qmat, Qmat_update);

    RcppParallel::parallelFor(0, site_index.n_elem, *this); //uses Qmat,Fmat
    //(*this)(0, site_index.n_elem);

    if (!(arma::is_finite(Acoef) && arma::is_finite(Bcoef) && arma::is_finite(Qmat_update)))
    {
      //instead: warning and return NaN?
      Rcpp::stop("[AdmixtureHybrid] Numeric underflow, remove nonsegregating sites");
    }
    Qmat_update += arma::sum(Acoef, 2) + arma::sum(Bcoef, 2);

    if (!fix_Fmat)
      F = arma::sum(Acoef, 1) / (arma::sum(Acoef, 1) + arma::sum(Bcoef, 1));

    if (!fix_Qmat)
    {
      Q = Qmat_update * D.t(); 
      Q.each_col([](arma::vec& x) { x /= arma::accu(x); }); // make sum_k q_{jk} = 1
    }
    
    return arma::accu(loglik) + loglik_msat;
  }

  void project (arma::mat& Q, arma::mat& F, std::vector<Microsat>& Ms)
  {
    for (auto& ms : Ms)
    {
      ms.Fmat = arma::clamp(ms.Fmat, errtol, arma::datum::inf);
      ms.Fmat.each_col([](arma::vec& x) { x /= arma::accu(x); });
    }
    F = arma::clamp(F, errtol, 1.-errtol);
    Q = arma::clamp(Q, errtol, arma::datum::inf);
    Q.each_col([](arma::vec& x) { x /= arma::accu(x); }); // make sum_k q_{jk} = 1
  }

  bool EMaccel (arma::mat& Q, arma::mat& F, std::vector<Microsat>& Ms)
  {
		const double tol = 1e-8;
		const double mstep = 4;

    double ll0 = loglikelihood;

    Q0 = Q;	F0 = F;
    for (auto& ms : Ms)
      ms.F0 = ms.Fmat;

		// first secant condition
    loglikelihood = EM(Q, F, Ms); //project(Q, F, Ms);
	  Q1  = Q;		  	F1  = F;
		Qd1 = Q - Q0;   Fd1 = F - F0;
    double ll1 = loglikelihood;
	  double sr2 = arma::accu(arma::pow(Qd1,2)) + arma::accu(arma::pow(Fd1,2));
    for (auto& ms : Ms)
    {
      ms.F1 = ms.Fmat; ms.Fd1 = ms.F1 - ms.F0;
      sr2  += arma::accu(arma::pow(ms.Fd1,2));
    }
	  if (!arma::is_finite(sr2) || fabs(ll1 - ll0) < tol || sqrt(sr2) < tol)
			return true;
    
    if (!acceleration) // vanilla EM
      return false;

		// second secant condition
    loglikelihood = EM(Q, F, Ms); //project(Q, F, Ms);
	  Q2  = Q;		  	F2  = F;
		Qd2 = Q - Q1;   Fd2 = F - F1;
    double em  = loglikelihood;
	  double sq2 = arma::accu(arma::pow(Qd2,2)) + arma::accu(arma::pow(Fd2,2));
    for (auto& ms : Ms)
    {
      ms.F2 = ms.Fmat; ms.Fd2 = ms.F2 - ms.F1;
      sq2  += arma::accu(arma::pow(ms.Fd2,2));
    }
	  if (!arma::is_finite(sq2) || fabs(em - ll1) < tol || sqrt(sq2) < tol)
			return true;

		// the magic happens
		Qd3 = Qd2 - Qd1; Fd3 = Fd2 - Fd1;
	  double sv2 = arma::accu(arma::pow(Qd3,2)) + arma::accu(arma::pow(Fd3,2));
    for (auto& ms : Ms)
    {
      ms.Fd3  = ms.Fd2 - ms.Fd1;
      sv2    += arma::accu(arma::pow(ms.Fd3,2));
    }
	  double alpha = sqrt(sr2/sv2);
		alpha = std::max(stepmin, std::min(stepmax, alpha));

    F = F0 + 2.*alpha*Fd1 + alpha*alpha*Fd3;
    Q = Q0 + 2.*alpha*Qd1 + alpha*alpha*Qd3;
    for (auto& ms : Ms)
      ms.Fmat = ms.F0 + 2.*alpha*ms.Fd1 + alpha*alpha*ms.Fd3;
    project (Q, F, Ms);

		// stabilize
		loglikelihood = EM(Q, F, Ms); project(Q, F, Ms);

    // revert to simple EM iteration if loglikelihood has decreased after acceleration
    if (!arma::is_finite(loglikelihood) || loglikelihood < em)
    {
      loglikelihood = em; Q = Q2; F = F2;
      for (auto& ms : Ms) ms.Fmat = ms.F2;
    }

		if (alpha == stepmax)
			stepmax *= mstep;
		
		return false;
  }

  // getters 

  std::vector<arma::mat> microsat_frequencies (void) const
  {
    std::vector<arma::mat> out;
    for (auto& ms : Msat)
      out.push_back(ms.Fmat);
    return out;
  }

  arma::mat microsat_loglikelihood (void) const
  {
    arma::mat out (sample_index.n_elem, Msat.size());
    for (arma::uword locus = 0; locus < Msat.size(); ++locus)
      out.col(locus) = Msat[locus].loglikelihood;
    return out;
  }

  arma::mat snp_frequencies (void) const
  {
    return Fmat;
  }

  arma::mat snp_loglikelihood (void) const
  {
    return loglik;
  }

  arma::mat admixture_proportions (void) const
  {
    return Qdeme;
  }
};

struct Saf : public RcppParallel::Worker
{
// On rescaling: this seems like a good idea in general and is necessary if bounds are dynamically updated. For diploids, we have
//   h[j,x] = 1/choose(chrom[j], x) * 
//        (choose(chrom[j-1], x) * L[0,j] * h[j-1,x] + 2 * choose(chrom[j-1], x-1) * L[1,j] * h[j-1,x-1] + choose(chrom[j-1],x-2) * L[2,j] * h[j-1,x-2])
//          = 1/(chrom[j]*(chrom[j]-1)) *
//            ( (chrom[j] - x) * (chrom[j] - x - 1) * L[0,j] * h[j-1,x] +
//              2 * x * (chrom[j] - x) * L[1,j] * h[j-1,x-1] +
//              x * (x - 1) * L[2,j] * h[j-1,x-2] )
// For haploids: ... 
//   h[j,x] = 1/choose(chrom[j], x) * 
//        (choose(chrom[j-1], x) * L[0,j] * h[j-1,x] + choose(chrom[j-1], x-1) * L[1,j] * h[j-1,x-1])
//          = 1/chrom[j] *
//            ( (chrom[j] - x) * L[0,j] * h[j-1,x] +
//              x * L[1,j] * h[j-1,x-1])

  const GenotypeLikelihood &GL;
  const arma::uvec site_index, sample_index;
  const unsigned chromosomes;
  arma::sp_mat SAF;

  private:
  const arma::uvec &rsite_index, &rsample_index;
  std::vector<unsigned> rowind, colind;
  std::vector<double>  values;

  public:
  Saf (const GenotypeLikelihood& GL, const arma::uvec site_index, const arma::uvec sample_index)
    : GL (GL)
    , sample_index (sample_index)
    , site_index (site_index)
    , chromosomes (arma::accu(GL.ploidy.elem(sample_index)))
    // references used by slaves
    , rsite_index (site_index)
    , rsample_index (sample_index)
  {
    if (site_index.max() > GL.sites)
      Rcpp::stop("Saf: invalid entries in site index");
    if (sample_index.max() > GL.samples)
      Rcpp::stop("Saf: invalid entries in sample index");

    // TODO: have to check if join works in same order as split (so that order of sites is preserved on merger)
    // Does it matter? It's all going into a matrix anyway

    RcppParallel::parallelReduce(0, site_index.n_elem, *this);
    //(*this)(0, site_index.n_elem);

    // make sparse matrix
    SAF = arma::sp_mat(arma::join_vert(arma::urowvec(rowind), arma::urowvec(colind)), arma::vec(values), chromosomes+1, site_index.n_elem);
  }

  Saf (const Saf& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , chromosomes (rhs.chromosomes)
    // references used by slaves
    , rsite_index (rhs.rsite_index)
    , rsample_index (rhs.rsample_index)
  {}

  void operator() (const size_t start, const size_t end)
  {
    // heuristic for memory allocation
    rowind.reserve((end-start+1)*5);
    colind.reserve((end-start+1)*5);
    values.reserve((end-start+1)*5);
    for (size_t i=start; i!=end; ++i)
    {
      // TODO totally missing data? is that an issue?
      score_limited(i, rowind, colind, values);
    }
  }

  void join (const Saf& rhs)
  {
    values.insert(values.end(), 
                  std::make_move_iterator(rhs.values.begin()), 
                  std::make_move_iterator(rhs.values.end()));
    rowind.insert(rowind.end(), 
                  std::make_move_iterator(rhs.rowind.begin()), 
                  std::make_move_iterator(rhs.rowind.end()));
    colind.insert(colind.end(), 
                  std::make_move_iterator(rhs.colind.begin()), 
                  std::make_move_iterator(rhs.colind.end()));
  }

  double likelihood0 (const arma::vec::fixed<3>& p, const arma::vec& h, const int pos, const unsigned nchr, const unsigned ploidy)
  { // not rescaled, as would be done in ANGSD
    double out = 0.;
    if (ploidy == 2)
    {
      if (pos >= 0 && pos <= nchr)
      {
        out += p[0] * h[pos];
        if (pos-1 >= 0)
        {
          out += 2 * p[1] * h[pos-1];
          if (pos-2 >= 0)
            out += p[2] * h[pos-2];
        }
      }
    } else if (ploidy == 1)
    {
      if (pos >= 0 && pos <= nchr)
      {
        out += p[0] * h[pos];
        if (pos-1 >= 0)
          out += p[2] * h[pos-1];
      }
    }
    return out;
  }

  double likelihood (const arma::vec::fixed<3>& p, const arma::vec& h, const int pos, const unsigned nchr, const unsigned ploidy)
  {
    double out = 0.;
    if (ploidy == 2)
    {
      if (pos >= 0 && pos <= nchr)
      {
        out += double(nchr-pos) * double(nchr-pos-1) * p[0] * h[pos];
        if (pos-1 >= 0)
        {
          out += 2. * double(pos) * double(nchr-pos) * p[1] * h[pos-1];
          if (pos-2 >= 0)
            out += double(pos) * double(pos-1) * p[2] * h[pos-2];
        }
      }
    } else if (ploidy == 1)
    {
      if (pos >= 0 && pos <= nchr)
      {
        out += double(nchr-pos) * p[0] * h[pos];
        if (pos-1 >= 0)
          out += double(pos) * p[2] * h[pos-1];
      }
    }
    return out;
  }

  void score_limited (const unsigned i, std::vector<unsigned>& rind, std::vector<unsigned>& cind, std::vector<double>& vals) 
  {
    const double epsilon = std::pow(10., -9);
    const unsigned i2 = rsite_index[i];

    arma::vec saf (chromosomes + 1, arma::fill::zeros);
    arma::vec::fixed<3> p;
    unsigned nchr, bestguess;
    int left, right;
    double check;

    // initialize with first sample
    unsigned j2 = rsample_index[0];
    p = GL.handler(j2,i2);
    //if (GL.missing.at(j2,i2))
    //  p.fill(1./3.);
    //else {
    //  p[0]  = GL.like.at(0,j2,i2);
    //  p[1]  = GL.like.at(1,j2,i2);
    //  p[2]  = GL.like.at(2,j2,i2);
    //}
    //p[1] *= double(GL.ploidy[j2]-1); // drop hets for haploids
    //p    /= arma::accu(p); // the amazing transformation of likelihood to posterior

    nchr  = GL.ploidy[j2];
    left  = 0;
    right = GL.ploidy[j2];
    if (GL.ploidy[j2] == 2)
    {
      saf[0]  = p[0];
      saf[1]  = p[1];
      saf[2]  = p[2];
    } 
    else if (GL.ploidy[j2] == 1)
    {
      saf[0]  = p[0];
      saf[1]  = p[2];
    }

    // main loop
    for (unsigned j=1; j<rsample_index.n_elem; ++j)
    {
      j2    = rsample_index[j];
      nchr += GL.ploidy[j2];

      p = GL.handler(j2,i2);
      //// assign equal likelihoods to missing data
      //// drop "haploid hets" and normalize
      //if (GL.missing.at(j2,i2))
      //  p.fill(1./3.);
      //else {
      //  p[0]  = GL.like.at(0,j2,i2);
      //  p[1]  = GL.like.at(1,j2,i2);
      //  p[2]  = GL.like.at(2,j2,i2);
      //}
      //p[1] *= double(GL.ploidy[j2]-1); 
      //p    /= arma::accu(p); // the amazing transformation of likelihood to posterior

      // MLE genotype
      if (p[0] > p[1] && p[0] > p[2])
        bestguess = 0;
      else if (p[1] > p[2] && GL.ploidy[j2] == 2)
        bestguess = 1;
      else // if missing this will always be the case
        bestguess = GL.ploidy[j2];

      // update left bound
      left  += bestguess;
      while (true)
      {
        check = likelihood (p, saf, left, nchr, GL.ploidy[j2]);
        if (check < epsilon || left == 0)
          break;
        left -= 1;
      }

      // update right bound
      right += bestguess;
      while (true)
      {
        check  = likelihood (p, saf, right, nchr, GL.ploidy[j2]);
        if (check < epsilon || right == nchr)
          break;
        right += 1;
      }

      // normalize
      for (int k=nchr; k>=0; --k)
        saf.at(k) = (k >= left && k <= right) ? likelihood(p, saf, k, nchr, GL.ploidy[j2]) : 0.;
      saf /= saf.max();
    }

    // store as a sparse matrix
    // (can't have concurrent write access to sparse matrix)
    for (unsigned k=left; k<=right; ++k)
    {
      rind.push_back(k);
      cind.push_back(i);
      vals.push_back(saf[k]);
    }
  }
};

struct Shf : public RcppParallel::Worker
{
  const GenotypeLikelihood &GL;
  const arma::uvec sample_index;
  const arma::umat site_pair;
  const unsigned chromosomes;

  const arma::umat configurations;

  arma::sp_mat SHF;

  private:
  const arma::uvec &rsample_index;
  const arma::umat &rsite_pair;
  std::vector<unsigned> rowind, colind;
  std::vector<double> values;
  arma::uvec ptr;

  public:
  Shf (const GenotypeLikelihood& GL, const arma::umat pairs, const arma::uvec samples)
    : GL (GL)
    , sample_index (samples)
    , site_pair (pairs)
    , chromosomes (arma::accu(GL.ploidy.elem(sample_index)))
    , configurations (utils::all_haplotype_configurations(chromosomes))
    // references used by slaves
    , rsample_index (sample_index)
    , rsite_pair (site_pair)
  {
    if (site_pair.max() > GL.sites)
      Rcpp::stop("Haf: invalid entries in site index");
    if (sample_index.max() > GL.samples)
      Rcpp::stop("Haf: invalid entries in sample index");

    RcppParallel::parallelReduce(0, site_pair.n_cols, *this);
    //(*this)(0, site_pair.n_cols);

    // make sparse matrix
    SHF = arma::sp_mat(arma::join_vert(arma::urowvec(rowind), arma::urowvec(colind)), arma::vec(values), configurations.n_cols, site_pair.n_cols);
  }

  Shf (const Shf& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , chromosomes (rhs.chromosomes)
    , configurations (rhs.configurations)
    // don't allocate me
    , sample_index (0, arma::fill::zeros)
    , site_pair (0, 0, arma::fill::zeros)
    // references used by slaves
    , rsample_index (rhs.rsample_index)
    , rsite_pair (rhs.rsite_pair)
  {}

  void operator() (const size_t start, const size_t end)
  {
    // starting points for each iteration
    arma::uvec sum = arma::trans(arma::sum(configurations,0));
    ptr = arma::zeros<arma::uvec>(chromosomes+1); 
    for (unsigned i=0; i<=chromosomes; ++i)
    {
      arma::uvec tmp = arma::find(sum==i);
      ptr[i] = tmp.at(0);
    }
    // heuristic for memory allocation, assuming average of 5 nonzero values per pair
    rowind.reserve((end-start+1)*5);
    colind.reserve((end-start+1)*5);
    values.reserve((end-start+1)*5);
    for (size_t i=start; i!=end; ++i)
    {
      dyn_algo (i, rowind, colind, values);
    }
  }

  void join (const Shf& rhs)
  {
    values.insert(values.end(), 
                  std::make_move_iterator(rhs.values.begin()), 
                  std::make_move_iterator(rhs.values.end()));
    rowind.insert(rowind.end(), 
                  std::make_move_iterator(rhs.rowind.begin()), 
                  std::make_move_iterator(rhs.rowind.end()));
    colind.insert(colind.end(), 
                  std::make_move_iterator(rhs.colind.begin()), 
                  std::make_move_iterator(rhs.colind.end()));
  }

  void dyn_algo (const unsigned index, std::vector<unsigned>& rind, std::vector<unsigned>& cind, std::vector<double>& vals)
  {
    // calculate coniguration likelihoods using dynamic programming algorithm
    // "returns" a compressed sparse column format
    
    const double epsilon = std::pow(10., -9);

    const unsigned site0 = rsite_pair.at(0, index),
                   site1 = rsite_pair.at(1, index);

    unsigned nchr = 0;
    int i, j, k;

    // initialize
    arma::cube buffer (chromosomes+1, chromosomes+1, chromosomes+1, arma::fill::zeros);
    buffer.at(0,0,0) = 1.;

    // calculate
    for (unsigned s=0; s<rsample_index.n_elem; ++s) {

      const unsigned sample = rsample_index[s];

      nchr += GL.ploidy[sample];

      arma::vec p0 = GL.handler(sample, site0),
                p1 = GL.handler(sample, site1);

      arma::mat lik = likelihood (p0, p1, GL.ploidy[sample]); // haplotype likelihoods

      for (unsigned l=ptr[nchr]; l<configurations.n_cols; ++l)
      {
        i = configurations.at(0,l);
        j = configurations.at(1,l);
        k = configurations.at(2,l);
        if (GL.ploidy[sample] == 1)
        {
          buffer.at(i,j,k) *= lik.at(0,0);
          buffer.at(i,j,k) += i-1 < 0 ? 0. : buffer.at(i-1,j,k) * lik.at(1,1);
          buffer.at(i,j,k) += j-1 < 0 ? 0. : buffer.at(i,j-1,k) * lik.at(2,2);
          buffer.at(i,j,k) += k-1 < 0 ? 0. : buffer.at(i,j,k-1) * lik.at(3,3);
        } 
        else if (GL.ploidy[sample] == 2)
        {
          buffer.at(i,j,k) *= lik.at(0,0);
          buffer.at(i,j,k) += i-2 < 0 ? 0. : buffer.at(i-2,j,k) * lik.at(1,1);
          buffer.at(i,j,k) += j-2 < 0 ? 0. : buffer.at(i,j-2,k) * lik.at(2,2);
          buffer.at(i,j,k) += k-2 < 0 ? 0. : buffer.at(i,j,k-2) * lik.at(3,3);
          buffer.at(i,j,k) += i-1 < 0 ? 0. : buffer.at(i-1,j,k) * (lik.at(0,1) + lik.at(1,0));
          buffer.at(i,j,k) += j-1 < 0 ? 0. : buffer.at(i,j-1,k) * (lik.at(0,2) + lik.at(2,0));
          buffer.at(i,j,k) += k-1 < 0 ? 0. : buffer.at(i,j,k-1) * (lik.at(0,3) + lik.at(3,0));
          buffer.at(i,j,k) += i-1 < 0 || j-1 < 0 ? 0. : buffer.at(i-1,j-1,k) * (lik.at(1,2) + lik.at(2,1));
          buffer.at(i,j,k) += i-1 < 0 || k-1 < 0 ? 0. : buffer.at(i-1,j,k-1) * (lik.at(1,3) + lik.at(3,1));
          buffer.at(i,j,k) += j-1 < 0 || k-1 < 0 ? 0. : buffer.at(i,j-1,k-1) * (lik.at(2,3) + lik.at(3,2));
        }
      }
    }

    // transfer, sparsify, and store
    arma::vec out (configurations.n_cols, arma::fill::zeros);
    for (unsigned l=0; l<configurations.n_cols; ++l)
    {
      i = configurations.at(0,l);
      j = configurations.at(1,l);
      k = configurations.at(2,l);
      out.at(l) = buffer.at(i,j,k) / utils::multichoose(arma::uvec({chromosomes-unsigned(i)-unsigned(j)-unsigned(k), unsigned(i), unsigned(j), unsigned(k)}));
    }
    out.clean(epsilon);
    out /= arma::accu(out);

    arma::uvec nz = arma::find(out);
    for (unsigned l=0; l<nz.n_elem; ++l)
    {
      rind.push_back(nz[l]);
      cind.push_back(index);
      vals.push_back(out[nz[l]]);
    }

    return;
  }

  arma::mat likelihood (const arma::vec::fixed<3>& p0, const arma::vec::fixed<3>& p1, const unsigned ploidy)
  {
    // haplotype likelihoods from genotype likelihoods
    
    arma::mat::fixed<4,4> out;
    if (ploidy == 1)
    {
      out.at(0,0) = p0.at(0) * p1.at(0); // ab
      out.at(1,1) = p0.at(2) * p1.at(0); // Ab
      out.at(2,2) = p0.at(0) * p1.at(2); // aB
      out.at(3,3) = p0.at(2) * p1.at(2); // AB
    }
    else if (ploidy == 2)
    {
      out.at(0,0) = p0.at(0) * p1.at(0); // ab.ab
      out.at(1,0) = p0.at(1) * p1.at(0); // Ab.ab
      out.at(2,0) = p0.at(0) * p1.at(1); // aB.ab
      out.at(3,0) = p0.at(1) * p1.at(1); // AB.ab
      out.at(0,1) = p0.at(1) * p1.at(0); // ab.Ab
      out.at(1,1) = p0.at(2) * p1.at(0); // Ab.Ab
      out.at(2,1) = p0.at(1) * p1.at(1); // aB.Ab
      out.at(3,1) = p0.at(2) * p1.at(1); // AB.Ab
      out.at(0,2) = p0.at(0) * p1.at(1); // ab.aB
      out.at(1,2) = p0.at(1) * p1.at(1); // Ab.aB
      out.at(2,2) = p0.at(0) * p1.at(2); // aB.aB
      out.at(3,2) = p0.at(1) * p1.at(2); // AB.aB
      out.at(0,3) = p0.at(1) * p1.at(1); // ab.AB
      out.at(1,3) = p0.at(2) * p1.at(1); // Ab.AB
      out.at(2,3) = p0.at(1) * p1.at(2); // aB.AB
      out.at(3,3) = p0.at(2) * p1.at(2); // AB.AB
    }
    return out;
  }
};

struct ShfInvariant : public RcppParallel::Worker
{
  // for a set of sites, calculate the SHF for a left- and right- invariant pairing
  
  const GenotypeLikelihood &GL;
  const arma::uvec sample_index, site_index;
  const unsigned chromosomes;

  const arma::umat configurations;

  arma::sp_mat lSHF, rSHF;

  private:
  const arma::uvec &rsample_index, &rsite_index;
  std::vector<unsigned> rowindl, colindl, rowindr, colindr;
  std::vector<double> valuesl, valuesr;
  arma::uvec ptr;

  public:
  ShfInvariant (const GenotypeLikelihood& GL, const arma::uvec sites, const arma::uvec samples)
    : GL (GL)
    , sample_index (samples)
    , site_index (sites)
    , chromosomes (arma::accu(GL.ploidy.elem(sample_index)))
    , configurations (utils::all_haplotype_configurations(chromosomes))
    // references used by slaves
    , rsample_index (sample_index)
    , rsite_index (site_index)
  {
    if (site_index.max() > GL.sites)
      Rcpp::stop("[ShfInvariant] invalid entries in site index");
    if (sample_index.max() > GL.samples)
      Rcpp::stop("[ShfInvariant] invalid entries in sample index");

    RcppParallel::parallelReduce(0, site_index.n_elem, *this);
    //(*this)(0, site_index.n_elem);

    // make sparse matrices
    lSHF = arma::sp_mat(arma::join_vert(arma::urowvec(rowindl), arma::urowvec(colindl)), arma::vec(valuesl), configurations.n_cols, site_index.n_elem);
    rSHF = arma::sp_mat(arma::join_vert(arma::urowvec(rowindr), arma::urowvec(colindr)), arma::vec(valuesr), configurations.n_cols, site_index.n_elem);
  }

  ShfInvariant (const ShfInvariant& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , chromosomes (rhs.chromosomes)
    , configurations (rhs.configurations)
    // don't allocate me
    , sample_index (0, arma::fill::zeros)
    , site_index (0, arma::fill::zeros)
    // references used by slaves
    , rsample_index (rhs.rsample_index)
    , rsite_index (rhs.rsite_index)
  {}

  void operator() (const size_t start, const size_t end)
  {
    // starting points for each iteration
    arma::uvec sum = arma::trans(arma::sum(configurations,0));
    ptr = arma::zeros<arma::uvec>(chromosomes+1); 
    for (unsigned i=0; i<=chromosomes; ++i)
    {
      arma::uvec tmp = arma::find(sum==i);
      ptr[i] = tmp.at(0);
    }
    // heuristic for memory allocation, assuming average of 5 nonzero values per pair
    rowindl.reserve((end-start+1)*5);
    colindl.reserve((end-start+1)*5);
    valuesl.reserve((end-start+1)*5);
    rowindr.reserve((end-start+1)*5);
    colindr.reserve((end-start+1)*5);
    valuesr.reserve((end-start+1)*5);
    for (size_t i=start; i!=end; ++i)
    {
      dyn_algo (i, rowindl, colindl, valuesl, true);
      dyn_algo (i, rowindr, colindr, valuesr, false);
    }
  }

  void join (const ShfInvariant& rhs)
  {
    valuesl.insert(valuesl.end(), 
                  std::make_move_iterator(rhs.valuesl.begin()), 
                  std::make_move_iterator(rhs.valuesl.end()));
    rowindl.insert(rowindl.end(), 
                  std::make_move_iterator(rhs.rowindl.begin()), 
                  std::make_move_iterator(rhs.rowindl.end()));
    colindl.insert(colindl.end(), 
                  std::make_move_iterator(rhs.colindl.begin()), 
                  std::make_move_iterator(rhs.colindl.end()));
    valuesr.insert(valuesr.end(), 
                  std::make_move_iterator(rhs.valuesr.begin()), 
                  std::make_move_iterator(rhs.valuesr.end()));
    rowindr.insert(rowindr.end(), 
                  std::make_move_iterator(rhs.rowindr.begin()), 
                  std::make_move_iterator(rhs.rowindr.end()));
    colindr.insert(colindr.end(), 
                  std::make_move_iterator(rhs.colindr.begin()), 
                  std::make_move_iterator(rhs.colindr.end()));
  }

  void dyn_algo (const unsigned index, std::vector<unsigned>& rind, std::vector<unsigned>& cind, std::vector<double>& vals, const bool left)
  {
    // calculate configuration likelihoods using dynamic programming algorithm
    // "returns" a compressed sparse column format
    
    const double epsilon = std::pow(10., -9);

    const unsigned site = rsite_index.at(index);

    unsigned nchr = 0;
    int i, j, k;

    // initialize
    arma::cube buffer (chromosomes+1, chromosomes+1, chromosomes+1, arma::fill::zeros);
    buffer.at(0,0,0) = 1.;

    // calculate
    for (unsigned s=0; s<rsample_index.n_elem; ++s) {

      const unsigned sample = rsample_index[s];

      nchr += GL.ploidy[sample];

      arma::vec p0, p1;
      if (left) // hypothetical site to the left is invariant
      {
        p0 = arma::vec({1., 0., 0.}); 
        p1 = GL.handler(sample, site);
      } 
      else // hypothetical site to the right is invariant
      {
        p0 = GL.handler(sample, site);
        p1 = arma::vec({1., 0., 0.}); 
      }

      arma::mat lik = likelihood (p0, p1, GL.ploidy[sample]); // haplotype likelihoods

      for (unsigned l=ptr[nchr]; l<configurations.n_cols; ++l)
      {
        i = configurations.at(0,l);
        j = configurations.at(1,l);
        k = configurations.at(2,l);
        if (GL.ploidy[sample] == 1)
        {
          buffer.at(i,j,k) *= lik.at(0,0);
          buffer.at(i,j,k) += i-1 < 0 ? 0. : buffer.at(i-1,j,k) * lik.at(1,1);
          buffer.at(i,j,k) += j-1 < 0 ? 0. : buffer.at(i,j-1,k) * lik.at(2,2);
          buffer.at(i,j,k) += k-1 < 0 ? 0. : buffer.at(i,j,k-1) * lik.at(3,3);
        } 
        else if (GL.ploidy[sample] == 2)
        {
          buffer.at(i,j,k) *= lik.at(0,0);
          buffer.at(i,j,k) += i-2 < 0 ? 0. : buffer.at(i-2,j,k) * lik.at(1,1);
          buffer.at(i,j,k) += j-2 < 0 ? 0. : buffer.at(i,j-2,k) * lik.at(2,2);
          buffer.at(i,j,k) += k-2 < 0 ? 0. : buffer.at(i,j,k-2) * lik.at(3,3);
          buffer.at(i,j,k) += i-1 < 0 ? 0. : buffer.at(i-1,j,k) * (lik.at(0,1) + lik.at(1,0));
          buffer.at(i,j,k) += j-1 < 0 ? 0. : buffer.at(i,j-1,k) * (lik.at(0,2) + lik.at(2,0));
          buffer.at(i,j,k) += k-1 < 0 ? 0. : buffer.at(i,j,k-1) * (lik.at(0,3) + lik.at(3,0));
          buffer.at(i,j,k) += i-1 < 0 || j-1 < 0 ? 0. : buffer.at(i-1,j-1,k) * (lik.at(1,2) + lik.at(2,1));
          buffer.at(i,j,k) += i-1 < 0 || k-1 < 0 ? 0. : buffer.at(i-1,j,k-1) * (lik.at(1,3) + lik.at(3,1));
          buffer.at(i,j,k) += j-1 < 0 || k-1 < 0 ? 0. : buffer.at(i,j-1,k-1) * (lik.at(2,3) + lik.at(3,2));
        }
      }
    }

    // transfer, sparsify, and store
    arma::vec out (configurations.n_cols, arma::fill::zeros);
    for (unsigned l=0; l<configurations.n_cols; ++l)
    {
      i = configurations.at(0,l);
      j = configurations.at(1,l);
      k = configurations.at(2,l);
      out.at(l) = buffer.at(i,j,k) / utils::multichoose(arma::uvec({chromosomes-unsigned(i)-unsigned(j)-unsigned(k), unsigned(i), unsigned(j), unsigned(k)}));
    }
    out.clean(epsilon);
    out /= arma::accu(out);

    arma::uvec nz = arma::find(out);
    for (unsigned l=0; l<nz.n_elem; ++l)
    {
      rind.push_back(nz[l]);
      cind.push_back(index);
      vals.push_back(out[nz[l]]);
    }

    return;
  }

  arma::mat likelihood (const arma::vec::fixed<3>& p0, const arma::vec::fixed<3>& p1, const unsigned ploidy)
  {
    // haplotype likelihoods from genotype likelihoods
    
    arma::mat::fixed<4,4> out;
    if (ploidy == 1)
    {
      out.at(0,0) = p0.at(0) * p1.at(0); // ab
      out.at(1,1) = p0.at(2) * p1.at(0); // Ab
      out.at(2,2) = p0.at(0) * p1.at(2); // aB
      out.at(3,3) = p0.at(2) * p1.at(2); // AB
    }
    else if (ploidy == 2)
    {
      out.at(0,0) = p0.at(0) * p1.at(0); // ab.ab
      out.at(1,0) = p0.at(1) * p1.at(0); // Ab.ab
      out.at(2,0) = p0.at(0) * p1.at(1); // aB.ab
      out.at(3,0) = p0.at(1) * p1.at(1); // AB.ab
      out.at(0,1) = p0.at(1) * p1.at(0); // ab.Ab
      out.at(1,1) = p0.at(2) * p1.at(0); // Ab.Ab
      out.at(2,1) = p0.at(1) * p1.at(1); // aB.Ab
      out.at(3,1) = p0.at(2) * p1.at(1); // AB.Ab
      out.at(0,2) = p0.at(0) * p1.at(1); // ab.aB
      out.at(1,2) = p0.at(1) * p1.at(1); // Ab.aB
      out.at(2,2) = p0.at(0) * p1.at(2); // aB.aB
      out.at(3,2) = p0.at(1) * p1.at(2); // AB.aB
      out.at(0,3) = p0.at(1) * p1.at(1); // ab.AB
      out.at(1,3) = p0.at(2) * p1.at(1); // Ab.AB
      out.at(2,3) = p0.at(1) * p1.at(2); // aB.AB
      out.at(3,3) = p0.at(2) * p1.at(2); // AB.AB
    }
    return out;
  }
};

// problem: if done separately per stratum, sites won't line up
// solution: (a) do all strata in a single run
//           (b) call snps using all samples in GL (easy)
struct ShfWindow : public RcppParallel::Worker
{
  struct RecombinationMap
  {
    // recombination map from Hapmap format: 
    //   chromosome
    //   position (bp)
    //   rate (cM/Mb)
    //   position (cM)

    arma::uvec chrom, begin, end;
    arma::vec cM, rate;

    RecombinationMap (arma::mat hapmap)
      : chrom (hapmap.n_rows-1, arma::fill::zeros)
      , begin (hapmap.n_rows-1, arma::fill::zeros)
      , end (hapmap.n_rows-1, arma::fill::zeros)
      , cM (hapmap.n_rows-1, arma::fill::zeros)
      , rate (hapmap.n_rows-1, arma::fill::zeros)
    {
      if (hapmap.n_cols != 4)
        Rcpp::stop("[RecombinationMap] Hapmap format has four columns");
      if (hapmap.min() < 0.)
        Rcpp::stop("[RecombinationMap] Invalid map");

      arma::uvec chr = arma::conv_to<arma::uvec>::from(hapmap.col(0)),
        pos = arma::conv_to<arma::uvec>::from(hapmap.col(1));
      arma::vec  rat = hapmap.col(2),
        cen = hapmap.col(3);

      // sort by chromosome and then position
      arma::uvec index = arma::regspace<arma::uvec>(0, hapmap.n_rows-1);
      std::sort (index.begin(), index.end(), 
          [&](int i, int j) { 
          return chr[i]==chr[j] ? pos[i]<pos[j] : chr[i]<chr[j]; 
          });
      chr = chr.elem(index);
      pos = pos.elem(index);
      rat = rat.elem(index);
      cen = cen.elem(index);

      unsigned j = 0;
      for (unsigned i=0; i<pos.n_elem-1; ++i)
      {
        if (chr.at(i) == chr.at(i+1))
        {
          chrom.at(j) = chr.at(i);
          begin.at(j) = pos.at(i);
          end.at(j)   = pos.at(i+1);
          cM.at(j)    = cen.at(i);
          rate.at(j)  = rat.at(i);
          j++;
        }
      }
      chrom = chrom.head(j);
      begin = begin.head(j);
      end = end.head(j);
      cM = cM.head(j);
      rate = rate.head(j);
    }

    arma::vec interpolate (const arma::umat& position) const
    {
      // given a set of positions and a piecewise constant recombination map (rows 0 & 1 of input)
      // interpolate location from start of chromosome in centimorgans

      auto chr = position.row(0);
      auto pos = position.row(1);

      arma::vec out (pos.n_elem);
      out.fill(arma::datum::nan);

      for (unsigned i=0; i<pos.n_elem; ++i)
      {
        arma::uvec j = arma::find(chrom == chr[i] && begin <= pos[i] && end > pos[i]);
        if (j.n_elem == 1)
          out[i] = cM[j[0]] + double(pos[i] - begin[j[0]]) * rate[j[0]] * 1.e-6;
      }
      return out;
    }
  };

  const GenotypeLikelihood& GL;

  const arma::uvec sample_index;
  const arma::mat bins; // "distance bins" (cM) in which to calculate SHF 
  const RecombinationMap rmap;
  const arma::vec  cM; // positions of sites in cM
  const arma::uvec variant_index; // 1-based SNP id, 0 if invariant

  // outputs
  std::vector<arma::sp_mat> SHF;
  std::vector<arma::uvec> weights;
  std::vector<arma::umat> configurations;

  private:
  const arma::vec &rcM;
  const arma::uvec &rvariant_index;
  std::vector<arma::umat> variant; // indices of pairs of variant sites
  std::vector<arma::uvec> left_invariant, right_invariant; // for each site, number of invariant sites to right/left
  std::vector<unsigned>   invariant; // number of invariant pairs

  public:
  ShfWindow (const GenotypeLikelihood& GL,
             const arma::uvec samples,
             const arma::mat  hapmap,   // each row is a breakpoint, cols are: chr, pos(bp), rate, pos(cM)
             const arma::mat  bin,      // start/end points of bins (cM)
             const double     thresh)   // P-value to call SNP
    : GL (GL)
    , sample_index (samples)
    , rmap (hapmap)
    , bins (buffer_bins(bin))
    , cM (rmap.interpolate(GL.pos))
    , variant_index (call_variants(thresh))
    // references for workers
    , rvariant_index (variant_index)
    , rcM (cM)
  {
    // determine indices of variant pairs, and number of left/right/dual- invariant pairs.
    std::fprintf(stderr, "[ShfWindow] Counting variant pairs per window\n");
    RcppParallel::parallelReduce(0, GL.sites-1, *this);
    //(*this)(0, GL.sites-1);

    // calculate the left- and right- invariant SHF for each variant site
    std::fprintf(stderr, "[ShfWindow] Calculating SHF for variant-invariant pairs\n");
    ShfInvariant shfinv (GL, arma::find(variant_index), sample_index); 

    // invariant-invariant SHF
    arma::sp_mat iSHF (shfinv.lSHF.n_rows, 1);
    arma::uvec invariant_bin = arma::find(
        shfinv.configurations.row(0) == 0 && 
        shfinv.configurations.row(1) == 0 && 
        shfinv.configurations.row(2) == 0 );
    iSHF.at(invariant_bin[0]) = 1.;

    // assemble SHF, weights per bin
    std::fprintf(stderr, "[ShfWindow] Calculating SHF for variant-variant pairs\n");
    for (unsigned i=0; i<bins.n_rows; ++i)
    {
      if (variant[i].n_cols)
      {
        // ensure that variant index pairs are sorted
        arma::uvec index = arma::regspace<arma::uvec>(0, variant[i].n_cols-1);
        std::sort (index.begin(), index.end(), 
            [&](int j, int k){ return
            variant[i].at(0,j)==variant[i].at(0,k) ? 
            variant[i].at(1,j)<variant[i].at(1,k)  : 
            variant[i].at(0,j)<variant[i].at(0,k); });
        variant[i] = variant[i].cols(index);

        // calculate the variant-variant SHF for a bin
        Shf shf (GL, variant[i], sample_index);

        SHF.push_back(arma::join_horiz(arma::join_horiz(arma::join_horiz(shf.SHF, shfinv.lSHF), shfinv.rSHF), iSHF));
        weights.push_back(arma::join_vert(arma::join_vert(arma::join_vert(arma::ones<arma::uvec>(variant[i].n_cols), left_invariant[i]), right_invariant[i]), arma::uvec({invariant[i]})));
        configurations.push_back(shf.configurations);
      }
      else
      {
        SHF.push_back(arma::join_horiz(arma::join_horiz(shfinv.lSHF, shfinv.rSHF), iSHF));
        weights.push_back(arma::join_vert(arma::join_vert(left_invariant[i], right_invariant[i]), arma::uvec({invariant[i]})));
        configurations.push_back(shfinv.configurations);
      }
    }
  }

  ShfWindow (const ShfWindow& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , sample_index (rhs.sample_index)
    , rmap (rhs.rmap)
    , bins (rhs.bins)
    // don't allocate
    , cM (0, arma::fill::zeros)
    , variant_index (0, arma::fill::zeros)
    // references for workers
    , rvariant_index (rhs.rvariant_index)
    , rcM (rhs.rcM)
  {}

  void operator() (const size_t begin, const size_t end)
  {
    unsigned number_variants = rvariant_index.max();

    // allocate storage
    for(unsigned i=0; i<bins.n_rows; ++i)
    {
      variant.push_back(arma::zeros<arma::umat>(2, 0));
      left_invariant.push_back(arma::zeros<arma::uvec>(number_variants));
      right_invariant.push_back(arma::zeros<arma::uvec>(number_variants));
      invariant.push_back(0);
    }
        
    // in each window, assembly lists of sites to calculate SHF for and tallies of invariant sites
    for (unsigned i = begin; i != end; ++i)
    {
      // moving right
      unsigned bin = 0;
      for (unsigned j=i+1; j<GL.sites; ++j)
      {
        if (GL.pos.at(0,j) != GL.pos.at(0,i)) // reached end of chromosome; exit
          break;
        double dis = std::fabs(rcM[j] - rcM[i]);
        if (dis >= bins.at(bin,0)) // distance is greater that start of initial bin
        {
          if (dis > bins.at(bin,1)) // distance exceeds end point of current bin
            bin++;
          if (bin == bins.n_rows) // past final bin
            break;
          if (rvariant_index[i] && rvariant_index[j]) 
            // append to list of variant pairs in window
            variant[bin].insert_cols(variant[bin].n_cols, arma::uvec({i, j}));
          else if (!rvariant_index[i] && rvariant_index[j]) 
            // increment # left-invariants of site j
            left_invariant[bin].at(rvariant_index[j]-1)++;
          else if (rvariant_index[i] && !rvariant_index[j]) 
            // increment # right-invariants of site i
            right_invariant[bin].at(rvariant_index[i]-1)++;
          else 
            // increment # invariant pairs
            invariant[bin]++;
        }
      }
    }
  }

  void join (const ShfWindow& rhs)
  {
    for(unsigned i=0; i<bins.n_rows; ++i)
    {
      variant[i].insert_cols(variant[i].n_cols, rhs.variant[i]);
      left_invariant[i] += rhs.left_invariant[i];
      right_invariant[i] += rhs.right_invariant[i];
      invariant[i] += rhs.invariant[i];
    }
  }
  
  arma::mat buffer_bins (arma::mat bin)
  {
    if (bin.min() < 0.)
      Rcpp::stop("[ShfWindow::buffer_bins] Bin positions cannot be negative");
    if (bin.n_cols != 2)
      Rcpp::stop("[ShfWindow::buffer_bins] Bin matrix must have two columns (start and end in cM)");
    for (unsigned i=0; i<bin.n_rows; ++i)
      if (bin.at(i,0) >= bin.at(i,1))
        Rcpp::stop("[ShfWindow::buffer_bins] Bin start points must be strictly less than end points");

    // not a good idea to buffer with 0. and Inf.

    return bin;
  }

  arma::uvec call_variants (const double thresh)
  {
    // call and store indices of variants
    // NB: currently use all samples available. This is less storage efficient but ensures that
    //     all SHF that generated from the same GL have equivalent positions
    Frequencies freq (GL, arma::regspace<arma::uvec>(0, GL.sites-1), arma::regspace<arma::uvec>(0, GL.samples-1));
    arma::uvec out (GL.sites, arma::fill::zeros); // note this will be 1-based
    unsigned k = 1;
    for (unsigned i=0; i<GL.sites; ++i)
      if (freq.pval[i] < thresh)
        out[i] = k++;
    return out;
  }
};

struct CovarGP : public RcppParallel::Worker
{
  // expected number of variants per site/sample, using precomputed genotype probabilities

  const GenotypeLikelihood &GL;
  arma::vec global; 
  arma::mat covar, zscore;

  private:
  arma::vec &rglobal;
  arma::mat &rzscore;

  public:
  CovarGP (const GenotypeLikelihood& GL)
    : GL (GL)
    , global (GL.sites, arma::fill::zeros)
    , zscore (GL.sites, GL.samples, arma::fill::zeros)
    , covar (GL.samples, GL.samples, arma::fill::zeros)
    // references used by workers 
    , rglobal (global)
    , rzscore (zscore)
  {
    RcppParallel::parallelReduce (0, GL.sites, *this);
    //(*this)(0, GL.sites);
  }

  CovarGP (const CovarGP& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , global (0, arma::fill::zeros)
    , zscore (0, 0, arma::fill::zeros)
    , covar (arma::size(rhs.covar), arma::fill::zeros)
    // references used by workers
    , rglobal (rhs.rglobal)
    , rzscore (rhs.rzscore)
  {}

  double Z (const unsigned i, const unsigned j)
  {
    double pi = rglobal.at(i),
           e = double(GL.ploidy[j]) * pi,
           g = GL.ploidy[j] == 1 ? GL.post.at(2,j,i) : GL.post.at(1,j,i) + 2.*GL.post.at(2,j,i),
           v = sqrt(e * (1.-pi));

    return (g-e)/v;
  }

  double variance (const unsigned i, const unsigned j)
  {
    double pi = rglobal.at(i),
           e = double(GL.ploidy[j]) * pi,
           v = e * (1.-pi),
           g = GL.ploidy[j] == 1 ?
                 std::pow(0.-e,2) * GL.post.at(0,j,i) + std::pow(1.-e,2) * GL.post.at(2,j,i) :
                 std::pow(0.-e,2) * GL.post.at(0,j,i) + std::pow(1.-e,2) * GL.post.at(1,j,i) +
                 std::pow(2.-e,2) * GL.post.at(2,j,i);
    return g / v;
  }

  void operator() (const size_t start, const size_t end)
  {
    arma::vec var (GL.samples);
    arma::mat outer (GL.samples, GL.samples);

    for (size_t i=start; i!=end; ++i)
    {
      // global allele frequencies
      rglobal.at(i) = 0.;
      for (unsigned j=0; j<GL.samples; ++j)
        rglobal.at(i) += GL.ploidy[j] == 1 ? GL.post.at(2,j,i) : GL.post.at(1,j,i) + 2.*GL.post.at(2,j,i);
      rglobal.at(i) /= double(arma::accu(GL.ploidy));

      // Z scores and covariance
      for (unsigned j=0; j<GL.samples; ++j)
      {
        var.at(j)       = variance(i, j);
        rzscore.at(i,j) = Z(i, j);
      }
      outer        = rzscore.row(i).t() * rzscore.row(i);
      outer.diag() = var;
      covar       += outer / double(GL.sites - 1);
    }
  }

  void join (const CovarGP& rhs)
  {
    covar += rhs.covar;
  }
};

struct Covar : public RcppParallel::Worker
{// DONE except include estimation of global allele frequencies

  // I am basically using Anders' method but providing the individual allele freqs rather than iteratively estimating them
  // so, I would get something like the ngsTools method if I first use EM to get global freqs and use these
  // and I would get something like pcangsd if I use Admixture to get structure-aware freqs

  // for diploid genotype we have var(g) = 2p(1-p), e[g] = 2p, g = 2 (1-p_i) p_i l(1) 1 + p_i^2 l(2) 2
  // for haploid genotype we have var(g) = p(1-p), e[g] = p, g = l(2) p_i
  // so correlation for haploid-haploid:
  //    (l_i(2) * p - p) (l_j(2) * p - p) / (p * (1-p))
  // so correlation for haploid-diploid:
  //    (l_i(2) * p - p) (2 * l_j(1) * (1-p) * p + 2 * l_j(2) * p * p - 2 * p) / (sqrt(2) * p * (1-p)) 
  // and diploid diploid:
  //    (2 * l_i(1) * (1-p) * p + 2 * l_i(2) * p * p - 2 * p) (2 * l_j(1) * (1-p) * p + 2 * l_j(2) * p * p - 2 * p) / (2 * p * (1-p)) 

  const GenotypeLikelihood &GL;
  const arma::mat  freq;
  const arma::vec  global; //global allele frequency estimates: I could include estimation as part of this routine
  const arma::uvec site_index, sample_index;
  arma::mat        covar;
  arma::umat       sites;

  private:
  const arma::mat  &rfreq;
  const arma::vec  &rglobal;
  const arma::uvec &rsite_index, &rsample_index;

  public:
  Covar (const GenotypeLikelihood& GL, const arma::uvec site_index, const arma::uvec sample_index, const arma::mat freq, const arma::vec global)
    : GL (GL)
    , freq (freq)
    , global (global)
    , site_index (site_index)
    , sample_index (sample_index)
    , covar (sample_index.n_elem, sample_index.n_elem, arma::fill::zeros)
    , sites (sample_index.n_elem, sample_index.n_elem, arma::fill::zeros)
    // references used by slaves
    , rfreq (freq)
    , rglobal (global)
    , rsite_index (site_index)
    , rsample_index (sample_index)
  {
    if (site_index.max() > GL.sites)
      Rcpp::stop("Covar: invalid entries in site index");
    if (sample_index.max() > GL.samples)
      Rcpp::stop("Covar: invalid entries in sample index");
    if (freq.n_rows != sample_index.n_elem) 
      Rcpp::stop("Covar: individual site frequency matrix must have as many rows as entries of the sample index");
    if (freq.n_cols != site_index.n_elem)
      Rcpp::stop("Covar: individual site frequency matrix must have as many cols as entries of the site index");
    if (global.n_elem != site_index.n_elem)
      Rcpp::stop("Covar: global site frequency vector must have as many elements as entries of the site index");

    RcppParallel::parallelReduce (0, site_index.n_elem, *this);
    //(*this)(0, site_index.n_elem);
    
    covar /= arma::conv_to<arma::mat>::from(sites);
  }

  Covar (const Covar& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , covar (arma::size(rhs.covar), arma::fill::zeros)
    , sites (arma::size(rhs.sites), arma::fill::zeros)
    // references used by slaves
    , rfreq (rhs.rfreq)
    , rglobal (rhs.rglobal)
    , rsite_index (rhs.rsite_index)
    , rsample_index (rhs.rsample_index)
  {}

  double expectation (const double f, const arma::vec::fixed<3>& p, const unsigned ploidy)
  {
    // expected genotype given HWE and likelihoods
    double num = 0., den = 0.;
    if (ploidy == 1)
    {
      num = p[2] * f;
      den = p[0] * (1.-f) + p[2] * f;
    }
    else if (ploidy == 2)
    {
      num = p[1] * 2. * f * (1.-f) + 
            2. * p[2] * f * f;
      den = p[0] * (1.-f) * (1.-f) +
            p[1] * (1.-f) * f * 2. +
            p[2] * f * f;
    }
    return num/den;
  }

  double covariance (const unsigned i, const unsigned j, const unsigned k)
  {
    const unsigned i2 = rsite_index[i],
                   j2 = rsample_index[j],
                   k2 = rsample_index[k];
    arma::vec::fixed<3> pj, pk;
    double fj, ej, vj, gj, fk, ek, vk, gk, pi;

    // deal with missing data
    if (GL.missing(j2, i2))
      pj.fill(1./3.);
    else
    {
      pj[0] = GL.like(0,j2,i2);
      pj[1] = GL.like(1,j2,i2);
      pj[2] = GL.like(2,j2,i2);
    }
    pj[1] *= double(GL.ploidy[j2] - 1);
    pj    /= arma::accu(pj);

    if (GL.missing(k2, i2))
      pk.fill(1./3.);
    else
    {
      pk[0] = GL.like(0,k2,i2);
      pk[1] = GL.like(1,k2,i2);
      pk[2] = GL.like(2,k2,i2);
    }
    pk[1] *= double(GL.ploidy[k2] - 1);
    pk    /= arma::accu(pk);

    // covariance calculation
    fj = rfreq.at(j,i);
    fk = rfreq.at(k,i);
    pi = rglobal.at(i);
    ej = double(GL.ploidy[j2]) * pi;
    ek = double(GL.ploidy[k2]) * pi;
    gj = expectation(fj, pj, GL.ploidy[j2]);
    gk = expectation(fk, pk, GL.ploidy[k2]);
    vj = sqrt(double(GL.ploidy[j2]) * pi * (1.-pi));
    vk = sqrt(double(GL.ploidy[k2]) * pi * (1.-pi));

    return (gj-ej)/vj * (gk-ek)/vk;
  }

  double variance (const unsigned i, const unsigned j)
  {
    const unsigned i2 = rsite_index[i],
                   j2 = rsample_index[j];
    arma::vec::fixed<3> pj;
    double fj, ej, vj, gj, pi, num = 0., den = 0.;
    
    // deal with missing data
    if (GL.missing(j2, i2))
      pj.fill(1./3.);
    else
    {
      pj[0] = GL.like(0,j2,i2);
      pj[1] = GL.like(1,j2,i2);
      pj[2] = GL.like(2,j2,i2);
    }
    pj[1] *= double(GL.ploidy[j2] - 1);
    pj    /= arma::accu(pj);

    // variance calculation
    fj = rfreq.at(j,i);
    pi = rglobal.at(i);
    ej = double(GL.ploidy[j2]) * pi;
    vj = double(GL.ploidy[j2]) * pi * (1.-pi);

    if (GL.ploidy[j2] == 1)
    {
      num = std::pow(0.-ej,2) * pj[0] * (1.-fj) + 
            std::pow(1.-ej,2) * pj[2] * fj;
      den = pj[0] * (1.-fj) + pj[2] * fj;
    }
    else if (GL.ploidy[j2] == 2)
    {
      num = std::pow(0.-ej,2) * pj[0] * (1.-fj) * (1.-fj) +
            std::pow(1.-ej,2) * pj[1] * fj * (1.-fj) * 2. +
            std::pow(2.-ej,2) * pj[2] * fj * fj; 
      den = pj[0] * (1.-fj) * (1.-fj) + 
            pj[1] * fj * (1.-fj) * 2. +
            pj[2] * fj * fj;
    }
    gj = num/den;
    
    return gj / vj;
  }

  void operator() (const size_t start, const size_t end)
  {
    double out;
    for (size_t i=start; i!=end; ++i)
      for (unsigned j=0; j<rsample_index.n_elem; ++j)
        for (unsigned k=j; k<rsample_index.n_elem; ++k)
        {
          if (k == j)
          {
            out = variance(i, j);
            if (arma::is_finite(out))
            {
              sites.at(j,j) += 1;
              covar.at(j,j) += out;
            }
          }
          else
          {
            out = covariance(i, j, k);
            if (arma::is_finite(out))
            {
              sites.at(j,k) += 1;
              sites.at(k,j) += 1;
              covar.at(j,k) += out;
              covar.at(k,j) += out;
            }
          }
        }
  }

  void join (const Covar& rhs)
  {
    covar += rhs.covar;
    sites += rhs.sites;
  }
};

struct Linkage : public RcppParallel::Worker
{
  // this struct implements pairwise haplotype frequency estimation, much like in ngsLD, except
  // allowing for mixed haploid-diploid samples.

  const GenotypeLikelihood &GL;
  const arma::umat sitepair_index; 
  const arma::uvec sample_index;
  const arma::vec freq;
  arma::mat hap, ld; 

  private:
  const arma::umat &rsitepair_index; 
  const arma::uvec &rsample_index;
  const arma::vec  &rfreq;
  arma::mat &rhap, &rld;

  public:
  Linkage (const GenotypeLikelihood& GL, const arma::umat sitepair_index, const arma::uvec sample_index)
    : GL (GL)
    , sitepair_index (sitepair_index)
    , sample_index (sample_index)
    , freq (frequencies(sample_index))
    , hap (sitepair_index.n_cols, 4)
    , ld (sitepair_index.n_cols, 4)
    // references shared across slaves
    , rsitepair_index (sitepair_index)
    , rsample_index (sample_index)
    , rfreq (freq)
    , rhap (hap)
    , rld (ld)
  {
    if (sitepair_index.n_rows != 2)
      Rcpp::stop ("[Linkage::Linkage] Must provide a pair of site indices in each row");
    if (sitepair_index.max() >= GL.sites)
      Rcpp::stop ("[Linkage::Linkage] Invalid entries in site indices");
    if (sample_index.max() >= GL.samples)
      Rcpp::stop ("[Linkage::Linkage] Invalid entries in sample indices");

    hap.fill(arma::datum::nan);
    ld.fill(arma::datum::nan);

    RcppParallel::parallelFor(0, sitepair_index.n_cols, *this);
    //(*this)(0, sitepair_index.n_cols);
  }

  Linkage (const Linkage& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , rsitepair_index (rhs.rsitepair_index)
    , rsample_index (rhs.rsample_index)
    , rfreq (rhs.rfreq)
    , rhap (rhs.rhap)
    , rld (rhs.rld)
  {}

  arma::vec frequencies (const arma::uvec& samp_ind)
  {
    const double errtol = 1e-6;
    Frequencies freq (GL, arma::regspace<arma::uvec>(0, GL.sites-1), samp_ind);
    return arma::clamp(freq.freq, errtol, 1.-errtol);
  }

  bool EM (arma::mat::fixed<2,2>& eta, double& loglik, const size_t i2, const size_t k2)
  {
    // EM step updating haplotype frequencies, from Li 2011 Bioinformatics.

    const double tol = 1e-5; // as per ngsLD
    const bool verbose = true; // for debugging, don't use in parallel

    arma::mat::fixed<2,2> Num; Num.zeros();
    double Den = 0.;
    loglik = 0.;

    for (auto j2 : rsample_index)
    {
      arma::vec::fixed<3> pj, Pj;
      arma::mat::fixed<2,2> num; num.zeros();
      double den = 0.;

      pj = GL.handler(j2, i2);
      Pj = GL.handler(j2, k2);

      // calculate conditional expectation of haplotype count in sample
      if (GL.ploidy[j2] == 2)
      {
        for (unsigned a=0; a<2; ++a) //strand 1, locus 1
          for (unsigned A=0; A<2; ++A) //strand 1, locus 2
            for (unsigned b=0; b<2; ++b) //strand 2, locus 1
              for (unsigned B=0; B<2; ++B) //strand 2, locus 2
              {
                num.at(a, A) += eta.at(b, B) * pj.at(a+b) * Pj.at(A+B);
                den += eta.at(a, A) * eta.at(b, B) * pj.at(a+b) * Pj.at(A+B);
              }
        num %= eta;
        loglik += log(den);
      }
      else if (GL.ploidy[j2] == 1)
      {
        for (unsigned a=0; a<2; ++a)
          for (unsigned A=0; A<2; ++A)
            num.at(a, A) = eta.at(a, A) * pj.at(2*a) * Pj.at(2*A);
        den = arma::accu(num);
        loglik += log(den);
      }

      Num += double(GL.ploidy[j2]) * num/den;
      Den += double(GL.ploidy[j2]);
    }

    double eps = arma::abs(eta - Num/Den).max();

    if (verbose)
      std::fprintf(stderr, "[Linkage::EM]\t(%lu,%lu) loglik: %f, |eta'-eta|: %f\n", i2, k2, loglik, eps);

    eta  = Num/Den;
    eta /= arma::accu(eta);

    return (eps < tol) || !arma::is_finite(eps);
  }

  void operator() (const size_t start, const size_t end)
  {
    const unsigned maxiter = 100;

    for (size_t p = start; p != end; ++p)
    {
      size_t i2 = rsitepair_index.at(0,p),
             k2 = rsitepair_index.at(1,p);

      arma::mat::fixed<2,2> eta, eta0;
      double loglik, loglik0;

      // initialize with HWE frequencies
      eta.at(0,0) = (1. - rfreq.at(i2)) * (1. - rfreq.at(k2));
      eta.at(1,0) = rfreq.at(i2) * (1. - rfreq.at(k2));
      eta.at(0,1) = (1. - rfreq.at(i2)) * rfreq.at(k2);
      eta.at(1,1) = rfreq.at(i2) * rfreq.at(k2);

      // optimize
      for (unsigned iter=0; iter<maxiter; ++iter)
        if (EM(eta, loglik, i2, k2))
          break;

      // store estimated haplotype frequencies
      rhap.at(p,0) = eta.at(0,0);
      rhap.at(p,1) = eta.at(1,0);
      rhap.at(p,2) = eta.at(0,1);
      rhap.at(p,3) = eta.at(1,1);

      // null likelihood
      double af1 = eta.at(1,0) + eta.at(1,1),
             af2 = eta.at(0,1) + eta.at(1,1);
      eta0.at(0,0) = (1. - af1) * (1. - af2);
      eta0.at(1,0) = af1 * (1. - af2);
      eta0.at(0,1) = (1. - af1) * af2;
      eta0.at(1,1) = af1 * af2;
      EM(eta0, loglik0, i2, k2);

      // calculate LD and related stats
      // D = P_BA * P_ba - P_Ba * P_bA = P_BA - p_A p_B
      rld.at(p,0) = eta.at(1,1) * eta.at(0,0) - eta.at(1,0) * eta.at(0,1); 
      rld.at(p,1) = rld.at(p,0) / (rld.at(p,0) < 0 ? 
            std::max(-(1.-af1)*(1.-af2), -af1*af2) : 
              std::min((1.-af1)*af2, af1*(1.-af2))); // D' TODO this doesn't seem to work?
      rld.at(p,2) = rld.at(p,0) / sqrt(af1 * (1.-af1) * af2 * (1.-af2)); // r
      rld.at(p,3) = 2.*loglik - 2.*loglik0;
    }
  }
};

struct SNPCaller : public RcppParallel::Worker
{
  // Calls SNPs using an Empirical Bayes method similar to Nielsen et al 2012 PLoS One
  // The fundamental idea is to calculate HWE probabilities from SAF + prior SFS, and use
  // these as priors when calculating genotype posteriors.
  
  // Filters sites down to SNPs with a posterior probability of being monomorphic < threshold,
  // or expected frequency > threshold. Works with samples stratified into multiple demes.

  GenotypeLikelihood& GL; // genotype posteriors are saved in GL.post

  const arma::uvec sample_index;
  const arma::uvec strata;

  arma::uvec site_index;

  const std::vector<arma::vec>    sfs_list; // important to initialize prior to saf_list
  const std::vector<arma::sp_mat> saf_list;

  size_t filtered = 0;

  private:
  const std::vector<arma::vec>    &rsfs_list;
  const std::vector<arma::sp_mat> &rsaf_list;
  const arma::uvec                &rsite_index;

  public:
  SNPCaller (GenotypeLikelihood& GL, 
             const arma::uvec& samples,
             const arma::uvec& samples_strata, 
             const std::vector<arma::vec>& sfs, 
             const double minmaf,
             const double snp_pval)
    : GL (GL)
    , sample_index (check_samples(samples))
    , strata (check_strata(samples_strata, sample_index))
    , sfs_list (fold_sfs(sfs, strata, sample_index))
    , saf_list (call_snps(site_index, sfs_list, minmaf, snp_pval)) // sets site_index
    // access-references for workers
    , rsfs_list (sfs_list)
    , rsaf_list (saf_list)
    , rsite_index (site_index)
  {
    // calculate genotype posteriors
    fprintf(stderr, "[SNPCaller] Calculating genotype posteriors\n");

    RcppParallel::parallelFor(0, site_index.n_elem, *this);
    //(*this)(0, site_index.n_elem);

    // filter
    filtered = GL.retain_sites(site_index); 
    fprintf(stderr, 
            "[SNPCaller] Filtered %lu sites that were monomorphic in %u samples across %u strata\n", 
            filtered, sample_index.n_elem, strata.max()+1);
  }

  SNPCaller (const SNPCaller& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , sample_index (rhs.sample_index)
    , strata (rhs.strata)
    // access references for workers
    , rsfs_list (rhs.rsfs_list)
    , rsaf_list (rhs.rsaf_list)
    , rsite_index (rhs.rsite_index)
  {}

  void operator () (const size_t start, const size_t end)
  {
    for (unsigned i = start; i != end; ++i)
    {
      const unsigned i2 = rsite_index.at(i);
      for (unsigned j = 0; j < sample_index.n_elem; ++j)
      {
        const unsigned j2 = sample_index.at(j);

        arma::vec gl   = GL.handler(j2,i2), 
                  post = arma::zeros<arma::vec>(arma::size(rsfs_list[strata[j]]));

        for (arma::sp_mat::const_col_iterator k = rsaf_list[strata[j]].begin_col(i2); 
             k != rsaf_list[strata[j]].end_col(i2); ++k) 
          post.at(k.row()) = (*k) * rsfs_list[strata[j]].at(k.row());
        post /= arma::accu(post);

        GL.post.slice(i2).col(j2) = GL.ploidy[j2] == 2 ? 
          probabilities (gl[0], gl[1], gl[2], post) :
          probabilities (gl[0], gl[2], post);
      }
    }
  }

  arma::uvec check_samples (const arma::uvec& sample)
  {
    fprintf(stderr, "[SNPCaller::check_samples]\n");//DEBUG

    if (sample.max() >= GL.samples)
      Rcpp::stop("[SNPCaller] Sample indices out of bounds");

    return sample;
  }

  arma::uvec check_strata (const arma::uvec& strat, const arma::uvec& sample)
  {
    fprintf(stderr, "[SNPCaller::check_strata]\n");//DEBUG

    arma::uvec uniq_strata = arma::unique(strat);

    if (strat.n_elem != sample.n_elem)
      Rcpp::stop("[SNPCaller] Must provide stratum for each sample");
    if (strat.max() != uniq_strata.n_elem-1)
      Rcpp::stop("[SNPCaller] Stratum ID must be 0-based and contiguous");

    return strat;
  }

  std::vector<arma::vec> fold_sfs (std::vector<arma::vec> sfs, const arma::uvec& strat, const arma::uvec& sample)
  {
    fprintf(stderr, "[SNPCaller::fold_sfs]\n");//DEBUG

    if (sfs.size() != strat.max() + 1)
      Rcpp::stop("[SNPCaller] Must provide a folded or unfolded SFS for each stratum");

    // check dimensions of SFS
    for (unsigned j=0; j<sfs.size(); ++j)
    {
      unsigned n = arma::accu(GL.ploidy.elem(sample(arma::find(strat == j))));

      if (sfs[j].n_elem != n+1)
      {
        fprintf(stderr, "[SNPCaller] SFS for stratum %u has %u chromosomes but data has %u chromosomes\n", j, sfs[j].n_elem-1, n);
        Rcpp::stop("[SNPCaller] Exiting");
      }

      if (sfs[j].min() < 0)
        Rcpp::stop("[SNPCaller] Cannot have negative entries of the SFS");

      sfs[j].replace(arma::datum::nan, 0.);
      sfs[j]  = 0.5 * (sfs[j] + arma::reverse(sfs[j]));
    }

    return sfs;
  }

  std::vector<arma::sp_mat> call_snps (arma::uvec& sites, const std::vector<arma::vec>& sfs, const double minmaf, const double pval)
  {
    // calculate SAF and call SNPs; indices of retained SNPs stored in "sites"

    fprintf(stderr, "[SNPCaller] Calculating SAF and calling SNPs across strata\n");

    std::vector<arma::sp_mat> saf;

    arma::vec post = arma::ones<arma::vec>(GL.sites);
    arma::vec freq = arma::zeros<arma::vec>(GL.sites);

    unsigned ns = 0;
    for (unsigned j=0; j<sfs.size(); ++j)
    {
      Saf SAF (GL, arma::regspace<arma::uvec>(0, GL.sites-1), sample_index.elem(arma::find(strata==j)));
      saf.push_back (SAF.SAF);
      SNPCaller::Henchman igor (saf[j], sfs[j]); // "what hump?"

      // total frequency
      freq += igor.stats.col(0) * double(sfs[j].n_elem - 1);
      ns += sfs[j].n_elem - 1;

      // assuming independence, posterior probability that site is globally monomorphic is product across strata
      post %= igor.stats.col(2);
    }

    freq /= double(ns);
    freq.transform([](double val){ return std::min(val, 1.-val); });

    // keep sites that are SNPs in at least one stratum and meet MAF filter
    sites = arma::find(post < pval && freq >= minmaf);

    return saf;
  }

  arma::Col<float>::fixed<3> probabilities (const double gl0, 
                                            const double gl1, 
                                            arma::vec post) 
  {
    // genotype probabilities for a haploid

    const unsigned n = post.n_elem - 1;

    arma::Col<float>::fixed<3> out;
    out.zeros();

    post /= arma::accu(post);

    out.at(2) = arma::dot(arma::regspace(0.,n), post)/double(n);
    out.at(0) = 1.-out.at(2);

    out.at(2) *= gl1;
    out.at(0) *= gl0;

    out /= arma::accu(out);

    return out;
  }

  arma::Col<float>::fixed<3> probabilities (const double gl0, 
                                            const double gl1, 
                                            const double gl2, 
                                            arma::vec post) 
  {
    // genotype probabilities for a diploid

    const unsigned n = post.n_elem - 1;

    arma::Col<float>::fixed<3> out;
    out.zeros();

    post /= arma::accu(post);

    ////I believe these expectations create false differentiation b/w haps and dips
    //arma::vec f = arma::regspace(0., n)/double(n);
    //out.at(2) = gl2 * arma::dot(arma::pow(f,2), post);
    //out.at(1) = gl1 * arma::dot(2.* f % (1.-f), post);
    //out.at(0) = gl0 * arma::dot(arma::pow(1.-f,2), post);
    
    //instead, use same expectation as for haploids
    double f = arma::dot(arma::regspace(0.,n), post)/double(n);    
    out.at(2) = gl2 * f * f;
    out.at(1) = gl1 * 2. * f * (1.-f);
    out.at(0) = gl0 * (1.-f) * (1.-f);

    out /= arma::accu(out);

    return out;
  }

  struct Henchman : public RcppParallel::Worker
  {
    // does SNP calling using Emp Bayes with SFS

    const arma::sp_mat& saf;
    arma::mat stats;
    
    private:
    arma::vec fsfs, y;
    arma::mat &rstats;

    public:
    Henchman (const arma::sp_mat& saf, 
              const arma::vec& sfs)
      : saf (saf)
      , fsfs (sfs)
      , y (sfs.n_elem)
      // don't copy these
      , stats (saf.n_cols, 3, arma::fill::zeros)
      // access references used by workers
      , rstats (stats)
    {
      if (sfs.n_elem != saf.n_rows)
        Rcpp::stop("[SNPCaller::Henchman] SFS dimension does not match SAF dimension");

      fsfs.replace(arma::datum::nan, 0.);
      fsfs = 0.5 * (fsfs + arma::reverse(fsfs));
      for (unsigned k=0; k<fsfs.n_elem; ++k)
        y.at(k) = k;
      //  y.at(k) = std::min(k, fsfs.n_elem-k-1); // this would estimate MAF


      RcppParallel::parallelFor(0, stats.n_rows, *this);
      //(*this)(0, stats.n_rows);
    }

    Henchman (const Henchman& rhs, RcppParallel::Split)
      : saf (rhs.saf)
      , fsfs (rhs.fsfs)
      , y (rhs.y)
      // don't copy these
      , stats (0, 0)
      // references used by workers
      , rstats (rhs.rstats)
    {}

    void operator () (const size_t start, const size_t end)
    {
      arma::vec buffer (fsfs.n_elem);

      for (size_t i = start; i != end; ++i)
      {
        buffer.zeros();
        for (arma::sp_mat::const_col_iterator k = saf.begin_col(i); 
             k != saf.end_col(i); ++k) 
          buffer.at(k.row()) = (*k) * fsfs.at(k.row());
        buffer /= arma::accu(buffer);

        rstats.at(i,2) = buffer.at(0) + buffer.at(fsfs.n_elem-1);        // posterior probability of monomorphic
        rstats.at(i,1) = log(rstats.at(i,2)/(1-rstats.at(i,2)));         // log Bayes' factor
        rstats.at(i,0) = arma::accu(y % buffer) / double(fsfs.n_elem-1); // posterior mean
      }
    }
  };
};

struct Filter
{
  // various in-place filtering operations

  struct MaxHet
  {
    // removes sites with a heterozygote majority greater than 0.5, given cutoff of estimated frequency

    size_t filtered = 0;

    MaxHet (GenotypeLikelihood& GL, const arma::uvec sample_index, const double freq)
    {
      if (sample_index.max() >= GL.samples)
        Rcpp::stop("[Filter::MaxHet] invalid sample indices");

      arma::uvec dips = arma::find(GL.ploidy.elem(sample_index) == 2);

      if (dips.n_elem > 0)
      {
        Genotypes geno (GL, arma::regspace<arma::uvec>(0, GL.sites-1), sample_index);
        arma::uvec drop = arma::find(geno.freq.col(1) > double(dips.n_elem) * freq);
        filtered = GL.remove_sites(drop);
      }
      fprintf(stderr, "[Filter::MaxHet] Removed %lu sites using %u diploids\n", filtered, dips.n_elem);
    }
  };

  struct MajHet
  {
    // removes sites with a heterozygote majority, given cutoff of probability

    size_t filtered = 0;
    
    MajHet (GenotypeLikelihood& GL, const arma::uvec sample_index, const double prob)
    {
      if (sample_index.max() >= GL.samples)
        Rcpp::stop("[Filter::MajHet] invalid sample indices");
      if (arma::any(GL.ploidy.elem(sample_index) == 2))
      {
        Paralogs para (GL, arma::regspace<arma::uvec>(0, GL.sites-1), sample_index);
        arma::uvec drop = arma::find(para.diphet.col(2) > prob);
        filtered = GL.remove_sites(drop);
      }
      fprintf(stderr, "[Filter::MajHet] Removed %lu sites\n", filtered);
    }
  };

  struct MinMAF
  {
    // removes sites with a minimum allele frequency lower than cutoff, and/or with SNP pval lower than cutoff

    size_t filtered;
    
    MinMAF (GenotypeLikelihood& GL, const arma::uvec sample_index, const double minmaf, const double pval)
    {
      if (sample_index.max() >= GL.samples)
        Rcpp::stop("[Filter::MinMAF] invalid sample indices");

      Frequencies freq (GL, arma::regspace<arma::uvec>(0, GL.sites-1), sample_index);
      freq.freq.transform ( [](double x) { return x > 0.5 ? 1.-x : x; } ); // to MAF
      arma::uvec keep = arma::find(freq.pval < pval && freq.freq >= minmaf);
      filtered = GL.retain_sites(keep);
      fprintf(stderr, "[Filter::MinMAF] Removed %lu sites\n", filtered);
    }

  };

  struct Minchrom
  {
    // removes sites with less than a minimum number of chromosomes per stratum
    size_t filtered;

    Minchrom (GenotypeLikelihood& GL, const arma::uvec sample_index, const arma::uvec strata, const arma::uvec mchr)
    {
      arma::uvec uniq_strata = arma::unique(strata);

      //TODO count chroms per strata, don't allow mchr to exceed this

      if (sample_index.max() >= GL.samples)
        Rcpp::stop("[Filter::Minchrom] Invalid sample indices");
      if (strata.n_elem != sample_index.n_elem)
        Rcpp::stop("[Filter::Minchrom] Strata dimensions don't match indices");
      if (mchr.n_elem != uniq_strata.n_elem)
        Rcpp::stop("[Filter::Minchrom] Must provide minimum chromosomes per stratum");
      if (uniq_strata.n_elem != uniq_strata.max()+1)
        Rcpp::stop("[Filter::Minchrom] Strata must be 0-based and contiguous");

      arma::Col<short> nonmiss (mchr.n_elem, arma::fill::zeros),
                       pass (GL.sites, arma::fill::ones);

      for (unsigned i=0; i<GL.sites; ++i)
      {
        nonmiss.zeros();
        for (unsigned j=0; j<sample_index.n_elem; ++j)
          nonmiss[strata[j]] += GL.missing.at(sample_index[j],i) ? 0 : GL.ploidy[sample_index[j]];
        for (unsigned k=0; k<mchr.n_elem; ++k)
          pass[i] *= int(nonmiss[k] >= mchr[k]);
      }

      filtered = GL.retain_sites (arma::find(pass));
      fprintf(stderr, "[Filter::Minchrom] Removed %lu sites\n", filtered);
    }
  };
  
  struct Minsample
  {
    // removes sites with less than a minimum number of samples per stratum
    size_t filtered;

    Minsample (GenotypeLikelihood& GL, const arma::uvec sample_index, const arma::uvec strata, const arma::uvec mind)
    {
      arma::uvec uniq_strata = arma::unique(strata);

      //TODO count samples per strata, don't allow mind to exceed this

      if (sample_index.max() >= GL.samples)
        Rcpp::stop("[Filter::Minsample] Invalid sample indices");
      if (strata.n_elem != sample_index.n_elem)
        Rcpp::stop("[Filter::Minsample] Strata dimensions don't match indices");
      if (mind.n_elem != uniq_strata.n_elem)
        Rcpp::stop("[Filter::Minsample] Must provide minimum samples per stratum");
      if (uniq_strata.n_elem != uniq_strata.max()+1)
        Rcpp::stop("[Filter::Minsample] Strata must be 0-based and contiguous");

      arma::Col<short> nonmiss (mind.n_elem, arma::fill::zeros),
                       pass (GL.sites, arma::fill::ones);

      for (unsigned i=0; i<GL.sites; ++i)
      {
        nonmiss.zeros();
        for (unsigned j=0; j<sample_index.n_elem; ++j)
          nonmiss[strata[j]] += GL.missing.at(sample_index[j],i) ? 0 : 1;
        for (unsigned k=0; k<mind.n_elem; ++k)
          pass[i] *= int(nonmiss[k] >= mind[k]);
      }

      filtered = GL.retain_sites (arma::find(pass));
      fprintf(stderr, "[Filter::Minsample] Removed %lu sites\n", filtered);
    }
  };

  struct Hetbias : public RcppParallel::Worker
  {
    // attempts to remove "heterozygous" sites that
    // diverge greatly from expectations
    
    GenotypeLikelihood& GL;
    const arma::uvec sample_index;
    const arma::vec::fixed<2> pval;
    const arma::vec err;
    arma::Row<short> missing_samples;
    arma::Col<size_t> missing_sites;

    size_t filtered;

    Hetbias (GenotypeLikelihood& GL, const arma::uvec sample_index, const arma::vec::fixed<2> pval, const arma::vec err)
      : GL (GL)
      , sample_index (sample_index)
      , pval (pval)
      , err  (err)
      , filtered (0)
    {
      if (err.n_elem != sample_index.n_elem)
        Rcpp::stop ("[Filter::Hetbias] Dimension mismatch");

      RcppParallel::parallelReduce(0, GL.sites, *this);
      //(*this)(0, GL.sites);

      fprintf(stderr, "[Filter::Hetbias] Removed %lu individual genotypes\n", filtered);
    }

    Hetbias (Hetbias& rhs, RcppParallel::Split)
      : GL (rhs.GL)
      , sample_index (rhs.sample_index)
      , pval (rhs.pval)
      , err  (rhs.err)
      , filtered (0)
    {}

    void operator () (const size_t start, const size_t end)
    {
      for (size_t i=start; i!=end; ++i)
        for (unsigned j=0; j<sample_index.n_elem; ++j)
          filter (i, j);
    }

    void filter (const unsigned i2, const unsigned j)
    {
      const unsigned j2 = sample_index[j];

      // if not heterozygote call, skip
      if (GL.missing.at(j2,i2) || (GL.like.at(1,j2,i2) < std::max(GL.like.at(0,j2,i2), GL.like.at(2,j2,i2))))
        return;

      // major and minor count
      const short major = GL.cnts.slice(i2).col(j2).max(),
                  minor = GL.cnts.slice(i2).col(j2).min();

      // filter _diploid heterozygote calls_ based on departure from 50/50 expectation
      // filter _haploid heterozygote calls_ based on departure from 100/0 expectation
      bool set_missing = pval[GL.ploidy[j2]-1] > (GL.ploidy[j2] == 1 ?
        binomial_cdf (minor, minor+major, err[j]) :
        binomial_cdf (major, minor+major, 0.5)) ;

      if (set_missing)
      {
        filtered += 1;
        GL.missing.at(j2,i2) = 1;
        GL.like.slice(i2).col(j2).fill(1./3.);
        GL.cnts.slice(i2).col(j2).zeros();
      }
    }

    void join (const Hetbias& rhs)
    {
      filtered += rhs.filtered;
    }

    double binomial_cdf (const unsigned y, const unsigned n, const double x)
    {
      double p = 0;
      for (unsigned k=y; k<=n; ++k)
        p += R::dbinom(k, n, x, false);
      return p;
    }
  };
};

struct Haplodiplo
{
  GenotypeLikelihood GL;
  MultiAllelicGenotypeLikelihood MAGL;

  Haplodiplo (const arma::cube& like, const arma::ucube& cnts, const arma::umat& pos, arma::uvec ploidy)
    : GL (like, cnts, arma::trans(pos), ploidy)
    , MAGL (ploidy)
  {}

  void multiallelic (const unsigned chrom, const unsigned start, const unsigned end, const arma::imat& genotypes, const double dropout, const double error)
  {
    MAGL.add_locus(chrom, start, end, genotypes, dropout, error);
  }

  Rcpp::List multiallelic_likelihoods (void) const
  {
    Rcpp::List like (MAGL.loci);
    for (int i=0; i<MAGL.loci; ++i)
    {
      Rcpp::NumericVector tmp = Rcpp::wrap(MAGL.like[i]);
      tmp.attr("dim") = Rcpp::Dimension(MAGL.like[i].n_rows, MAGL.like[i].n_cols, MAGL.like[i].n_slices);
      like[i] = tmp; 
    }
    return like;
  }

  arma::umat multiallelic_positions (void) const
  {
    return MAGL.pos;
  }

  Rcpp::List paralogs (arma::uvec site_index, arma::uvec sample_index)
  {
    Paralogs para (GL, site_index, sample_index);
    return Rcpp::List::create(
        Rcpp::_["haphet"] = para.haphet,
        Rcpp::_["diphet"] = para.diphet
        );
  }

  Rcpp::List poibin (arma::fvec prob)
  {
    Paralogs::PoissonBinomial pb (prob);
    return Rcpp::List::create(
        Rcpp::_["pdf"] = pb.PDF,
        Rcpp::_["cdf"] = pb.cdf(0.5),
        Rcpp::_["mean"] = pb.expectation()
        );
  }

  Rcpp::List frequencies (arma::uvec site_index, arma::uvec sample_index)
  {
    arma::arma_rng::set_seed(1);

    Frequencies freq (GL, site_index, sample_index);
    return Rcpp::List::create(
        Rcpp::_["freq"] = freq.freq,
        Rcpp::_["lrt"] = freq.lrt,
        Rcpp::_["pval"] = freq.pval
        );
  }

  Rcpp::List postsnp (arma::uvec site_index, arma::uvec sample_index, arma::vec sfs)
  {
      Saf saf (GL, site_index, sample_index);
      SNPCaller::Henchman igor (saf.SAF, sfs); // "Not the third switch!"
      return Rcpp::List::create(
          Rcpp::_["mean"]  = igor.stats.col(0),
          Rcpp::_["logbf"] = igor.stats.col(1),
          Rcpp::_["post"]  = igor.stats.col(2)
          );
  }

  arma::mat genofreq (arma::uvec site_index, arma::uvec sample_index)
  {
    arma::arma_rng::set_seed(1);

    Genotypes geno (GL, site_index, sample_index);
    return geno.freq;
  }

  Rcpp::List admixture (arma::uvec site_index, arma::uvec sample_index, arma::mat Qstart)
  {
    arma::arma_rng::set_seed(1);

    //arma::mat Qstart = arma::randu<arma::mat>(K, sample_index.n_elem);

    fprintf(stderr, "[Haplodiplo::admixture] Using all data\n");    
    Admixture admix (GL, site_index, sample_index, Qstart, false);

    return Rcpp::List::create(
        Rcpp::_["K"] = Qstart.n_rows,
        Rcpp::_["Loglik"] = admix.loglikelihood,    
        Rcpp::_["F"] = admix.Fmat,
        Rcpp::_["Q"] = admix.Qmat,
        Rcpp::_["AIC"] = -2. * admix.loglikelihood + 2 * double(admix.Fmat.n_elem + admix.Qmat.n_elem)
        );
  }

  Rcpp::List admixture_hybrid (arma::uvec site_index, arma::uvec sample_index, arma::uvec deme_index, arma::uvec marker_type, arma::mat Q, arma::mat F, unsigned maxiter = 1000)
  {
    fprintf(stderr, "[Haplodiplo::admixture_hybrid] Using %lu biallelic and %u multiallelic loci\n", GL.sites, MAGL.loci);    

    AdmixtureHybrid admix (GL, MAGL, site_index, sample_index, deme_index, marker_type, Q, F, maxiter);

    // correctly setting dimensions
    std::vector<arma::mat> Fms = admix.microsat_frequencies();
    Rcpp::List Fms_out (Fms.size());
    unsigned pars_ms = 0;
    for (int i=0; i<Fms.size(); ++i)
    {
      Rcpp::NumericVector tmp = Rcpp::wrap(Fms[i]);
      tmp.attr("dim") = Rcpp::Dimension(Fms[i].n_rows, Fms[i].n_cols);
      Fms_out[i] = tmp; 
      pars_ms += Q.n_rows * (Fms[i].n_rows - 1);
    }

    return Rcpp::List::create(
        Rcpp::_["K"] = Q.n_rows,
        Rcpp::_["Loglik"] = admix.loglikelihood,    
        Rcpp::_["F_snp"] = admix.snp_frequencies(),
        Rcpp::_["F_msat"] = Fms_out, 
        Rcpp::_["Q"] = admix.admixture_proportions(),
        Rcpp::_["AIC"] = -2. * admix.loglikelihood + 
          2 * double(admix.snp_frequencies().n_elem + (Q.n_rows - 1) * admix.admixture_proportions().n_cols + pars_ms)
        );
  }

  arma::sp_mat saf (arma::uvec site_index, arma::uvec sample_index)
  {
    Saf SAF (GL, site_index, sample_index);
    return SAF.SAF;
  }

  Rcpp::List shf (arma::umat sitepairs_index, arma::uvec sample_index)
  {
    Shf SHF (GL, arma::trans(sitepairs_index), sample_index);
    return Rcpp::List::create(
        Rcpp::_["SHF"] = SHF.SHF,
        Rcpp::_["config"] = SHF.configurations
        );
  }

  Rcpp::List shfwindow (arma::uvec sample_index, arma::mat hapmap, arma::mat bins, double pval)
  {
    ShfWindow SHF (GL, sample_index, hapmap, bins, pval);

    Rcpp::List SHFbin (SHF.SHF.size()),
               weight (SHF.SHF.size()),
               config (SHF.SHF.size());

    for (unsigned i=0; i<SHF.SHF.size(); ++i)
    {
      SHFbin[i] = SHF.SHF[i];
      weight[i] = SHF.weights[i];
      config[i] = SHF.configurations[i];
    }

    return Rcpp::List::create(
        Rcpp::_["SHF"] = SHFbin,
        Rcpp::_["weights"] = weight,
        Rcpp::_["config"] = config,
        Rcpp::_["bins"] = SHF.bins);
  }

  arma::umat configurations (arma::uvec sample_index)
  {
    unsigned nchr = arma::accu(GL.ploidy.elem(sample_index));
    return utils::all_haplotype_configurations(nchr);
  }

  Rcpp::List covar (arma::uvec site_index, arma::uvec sample_index, arma::mat freq)
  {
    //TODO: currently this will blow up if any site has global frequency of 0 or 1. I suppose it also might blow up if genotype likelihoods are 0/1 and do not align with HWE predictions from 'freq' matrix.
    //one way to avoid this is to filter out on basis of "global"
    arma::arma_rng::set_seed(1);

    Frequencies global (GL, site_index, sample_index);
    Covar cov (GL, site_index, sample_index, freq, global.freq);
    return Rcpp::List::create(
        Rcpp::_["covar"] = cov.covar,
        Rcpp::_["sites"] = cov.sites
        );
  }

  Rcpp::List covarGP (void)
  {
    CovarGP cov (GL);
    return Rcpp::List::create(
        Rcpp::_["freq"] = cov.global,
        Rcpp::_["covar"] = cov.covar,
        Rcpp::_["zscore"] = cov.zscore
        );
  }

  Rcpp::List linkage (arma::umat sitepairs_index, arma::uvec sample_index)
  {
    arma::arma_rng::set_seed(1);

    Linkage link (GL, arma::trans(sitepairs_index), sample_index);
    return Rcpp::List::create(
        Rcpp::_["hap"] = link.hap,
        Rcpp::_["ld"] = link.ld
        );
  }

  // filtering operations
  
  size_t retain (arma::uvec keep)
  {
    return GL.retain_sites (keep);
  }

  size_t remove (arma::uvec drop)
  {
    return GL.remove_sites (drop);
  }

  size_t minsample (arma::uvec sample_index, arma::uvec strata, arma::uvec mind)
  {
    Filter::Minsample msample (GL, sample_index, strata, mind);

    return msample.filtered;
  }

  size_t minchrom (arma::uvec sample_index, arma::uvec strata, arma::uvec mchr)
  {
    Filter::Minchrom mchrom (GL, sample_index, strata, mchr);

    return mchrom.filtered;
  }

  size_t minmaf (arma::uvec sample_index, double minmaf, double snp_pval)
  {
    Filter::MinMAF mmaf (GL, sample_index, minmaf, snp_pval);

    return mmaf.filtered;
  }

  size_t maxhet (arma::uvec sample_index, double freq)
  {
    Filter::MaxHet mhet (GL, sample_index, freq);

    return mhet.filtered;
  }

  size_t majhet (arma::uvec sample_index, double prob)
  {
    Filter::MajHet mhet (GL, sample_index, prob);

    return mhet.filtered;
  }

  size_t hetbias (arma::uvec sample_index, arma::vec err, arma::vec pval)
  {
    Filter::Hetbias hbias (GL, sample_index, arma::clamp(err, 0., 1.), pval);

    return hbias.filtered;
  }

  size_t callsnps (arma::uvec sample_index, arma::uvec strata, std::vector<arma::vec> sfs, double minmaf, double snp_pval)
  {
    SNPCaller caller (GL, sample_index, strata, sfs, minmaf, snp_pval);

    return caller.filtered;
  }

  size_t polarize (arma::uvec ref)
  {
    return GL.polarize(ref);
  }

  // getters
  
  arma::Cube<float> likelihoods (void) const
  {
    return GL.like;
  }

  arma::Cube<float> posteriors (void) const
  {
    return GL.post;
  }

  arma::Mat<short> genotypes (bool impute) const
  {
    // TODO: this uses posterior. Won't do shit if posterior hasn't been calculated by callsnps
    // Could instead initially set posterior to likelihoods, corrected for ploidy (e.g. by GL.handler)
    // I'll need to think about whether this is desirable behavior or not.
    // TODO: I think I did this^^^ ??

    arma::Mat<short> geno (GL.samples, GL.sites);
    geno.fill(9);

    for (unsigned i=0; i<GL.sites; ++i)
    {
      for (unsigned j=0; j<GL.samples; ++j)
      {
        if (GL.post.slice(i).col(j).is_finite() && (impute || !GL.missing(j,i)))
          geno.at(j,i) = arma::index_max(GL.post.slice(i).col(j)) / 
            (GL.ploidy[j] == 2 ? 1 : 2);
      }
    }

    return geno;
  }

  arma::Cube<short> majorminor (void) const
  {
    return GL.cnts;
  }

  arma::Cube<short> counts (void) const
  {
    arma::Cube<short> cnts (4, GL.samples, GL.sites, arma::fill::zeros);

    for (unsigned i=0; i<GL.sites; ++i)
    {
      unsigned major = GL.pos(2,i),
               minor = GL.pos(3,i);
      for (unsigned j=0; j<GL.samples; ++j)
      {
        cnts.at(major,j,i) = GL.cnts.at(0,j,i);
        cnts.at(minor,j,i) = GL.cnts.at(1,j,i);
      }
    }

    return cnts;
  }

  arma::umat position (void) const
  {
    return arma::trans(GL.pos);
  }

  arma::uvec ploidy (void) const
  {
    return GL.ploidy;
  }

  arma::uvec haploids (void) const
  {
    return GL.haploids;
  }

  arma::uvec diploids (void) const
  {
    return GL.diploids;
  }

  arma::Col<size_t> missing_sites (void) const
  {
    arma::Col<size_t> miss (GL.samples, arma::fill::zeros);
    for (unsigned i=0; i<GL.samples; ++i)
      miss[i] = arma::accu(arma::conv_to<arma::urowvec>::from(GL.missing.row(i)));
    return miss;
  }

  arma::urowvec missing_samples (void) const
  {
    return arma::conv_to<arma::urowvec>::from(arma::sum(GL.missing, 0));
  }

  arma::sp_umat missing (void) const
  {
    return arma::conv_to<arma::sp_umat>::from(arma::SpMat<short>(GL.missing));
  }

  size_t sites (void) const
  {
    return GL.sites;
  }

  unsigned samples (void) const
  {
    return GL.samples;
  }

  unsigned chromosomes (void) const
  {
    return GL.chromosomes;
  }

  void vcfbody (const std::string outfile, const std::vector<std::string> chrnames, const bool impute = false, const bool append = false) const
  {
    // dump vcf body to file

    if (chrnames.size() != arma::max(GL.pos.row(0)))
      Rcpp::stop("[[Haplodiplo::vcfbody]] Wrong length for chromosome names");

    std::ofstream out;
    out.open(outfile, append ? std::ios::app : std::ios::out);

    std::string acgtn = "ACGTN";

    if (append)
    {
    //std::string out = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    //for (unsigned i=0; i<GL.samples; ++i)
    //  out += "\t" + smpnames[i];
    //out += "\n";
    }

    // call SNPs
    auto gno = genotypes (impute);

    std::string gt, ad, gl;

    for (unsigned i=0; i<GL.sites; ++i)
    {
      unsigned chr = GL.pos.at(0,i),
               pos = GL.pos.at(1,i),
               ref = GL.pos.at(2,i),
               alt = GL.pos.at(3,i);

      out << chrnames[chr] + "\t" + 
       std::to_string(pos) + "\t" + 
                       "." + "\t" + 
                acgtn[ref] + "\t" + 
                acgtn[alt] + "\t" +
                       "." + "\t" +
                       "." + "\t" +
                       "." + "\t" +
                       "GT:AD:GL";
      for (unsigned j=0; j<GL.samples; ++j)
      {
        if (GL.ploidy[j] == 1)
        {
          if (gno.at(j,i) == 9)
          {
            gt = ".";
            ad = ".,.";
            gl = ".,.";
          }
          else
          {
            gt = std::to_string(gno.at(j,i));
            ad = std::to_string(GL.cnts.at(ref,j,i)) + "," + std::to_string(GL.cnts.at(alt,j,i));
            gl = std::to_string(GL.like.at(0,j,i)) + "," + std::to_string(GL.like.at(2,j,i));
          }
        } 
        else 
        {
          if (gno.at(j,i) == 9)
          {
            gt = "./.";
            ad = ".,.";
            gl = ".,.,.";
          }
          ad = std::to_string(GL.cnts.at(ref,j,i)) + "," + std::to_string(GL.cnts.at(alt,j,i));
          gl = std::to_string(GL.like.at(0,j,i)) + "," + std::to_string(GL.like.at(1,j,i)) + "," + std::to_string(GL.like.at(2,j,i));
        }
        out << "\t" + gt + ":" + ad + ":" + gl;
      }
      out << "\n";
    }
  }
};

RCPP_EXPOSED_CLASS_NODECL(Haplodiplo)

RCPP_MODULE(Haplodiplo) {
  using namespace Rcpp;
  class_<Haplodiplo>("Haplodiplo")
//    .constructor<std::string, arma::uvec>()
    .constructor<arma::cube, arma::ucube, arma::umat, arma::uvec>()
    .method("multiallelic", &Haplodiplo::multiallelic)
    .method("multiallelic_likelihoods", &Haplodiplo::multiallelic_likelihoods)
    .method("multiallelic_positions", &Haplodiplo::multiallelic_positions)
    .method("admixture", &Haplodiplo::admixture)
    .method("admixture_hybrid", &Haplodiplo::admixture_hybrid)
    .method("frequencies", &Haplodiplo::frequencies)
    .method("postsnp", &Haplodiplo::postsnp)
    .method("genofreq", &Haplodiplo::genofreq)
    .method("saf", &Haplodiplo::saf)
    .method("shf", &Haplodiplo::shf)
    .method("shfwindow", &Haplodiplo::shfwindow)
    .method("configurations", &Haplodiplo::configurations)
    .method("covar", &Haplodiplo::covar)
    .method("covarGP", &Haplodiplo::covarGP)
    .method("paralogs", &Haplodiplo::paralogs)
    .method("poibin", &Haplodiplo::poibin)
    .method("linkage", &Haplodiplo::linkage)
    .method("remove", &Haplodiplo::remove)
    .method("retain", &Haplodiplo::retain)
    .method("minsample", &Haplodiplo::minsample)
    .method("minchrom", &Haplodiplo::minchrom)
    .method("minmaf", &Haplodiplo::minmaf)
    .method("maxhet", &Haplodiplo::maxhet)
    .method("majhet", &Haplodiplo::majhet)
    .method("hetbias", &Haplodiplo::hetbias)
    .method("polarize", &Haplodiplo::polarize)
    .method("callsnps", &Haplodiplo::callsnps)
    .method("likelihoods", &Haplodiplo::likelihoods)
    .method("posteriors", &Haplodiplo::posteriors)
    .method("genotypes", &Haplodiplo::genotypes)
    .method("counts", &Haplodiplo::counts)
    .method("majorminor", &Haplodiplo::majorminor)
    .method("position", &Haplodiplo::position)
    .method("ploidy", &Haplodiplo::ploidy)
    .method("diploids", &Haplodiplo::diploids)
    .method("haploids", &Haplodiplo::haploids)
    .method("missing_sites", &Haplodiplo::missing_sites)
    .method("missing_samples", &Haplodiplo::missing_samples)
    .method("missing", &Haplodiplo::missing)
    .method("sites", &Haplodiplo::sites)
    .method("chromosomes", &Haplodiplo::chromosomes)
    .method("samples", &Haplodiplo::samples)
    .method("vcfbody", &Haplodiplo::vcfbody)
    ;
}

// 1D SFS ... TODO make into a module 

struct SFS1d : public RcppParallel::Worker
{
  const bool fold = false;
  const arma::uvec &multiplier, &block;
  const arma::sp_mat &saf; // these must be exponentiated!
  arma::vec sfs;
  double loglikelihood;

  const std::array<arma::uvec,3> folder;

  private:
  arma::vec upd, buffer;
  double loglik;
  size_t sites;

  // settings
  bool   acceleration = true;        // use EM acceleration?
  double errtol = 1e-16;             // not using dynamic bounds
  double stepmax = 1., stepmin = 1.; // adaptive steplength

  public:

  SFS1d (const arma::sp_mat& saf, const arma::uvec& block, const arma::uvec& multiplier, const bool fold)
    : fold (fold)
    , saf (saf)
    , multiplier (multiplier)
    , block (block)
    , sfs (saf.n_rows, arma::fill::ones)
    , folder (fold_mapper(sfs))
    // accumulated quantities
    , loglik (0.)
    , upd (arma::size(sfs), arma::fill::zeros)
    , sites (0)
    // temporaries
    , buffer (arma::size(sfs))
  {
    const unsigned maxiter = 1000;

    // TODO when module, move checks outside
    if (saf.n_cols != block.n_elem || block.max() >= multiplier.n_elem)
      Rcpp::stop("[Saf::SFS1d] Dimension mismatch");
    if (saf.min() < 0.)
      Rcpp::stop("[Saf::SFS1d] SAF must be non-negative");

    bool converged;
    for (unsigned iter=0; iter<maxiter; ++iter)
    {
      converged = EMaccel();
      if (converged) 
        break;
      fprintf(stderr, "[Saf::SFS1d]\titer %u, loglik = %f\n", iter, loglikelihood);
    }

    if (!converged)
      Rcpp::warning("[Saf::SFS1d] EM did not converge in maximum number of iterations");

    sfs.elem(arma::find(folder[2] == 0)).fill(arma::datum::nan);
    sfs *= double(sites);
  }

  SFS1d (const SFS1d& rhs, RcppParallel::Split)
    : fold (rhs.fold)
    , saf (rhs.saf)
    , multiplier (rhs.multiplier)
    , block (rhs.block)
    , sfs (rhs.sfs)
    , folder (rhs.folder)
    // accumulated quantities
    , loglik (0.)
    , upd (arma::size(sfs), arma::fill::zeros)
    , sites (0)
    // temporaries
    , buffer (arma::size(sfs))
  {}

  std::array<arma::uvec,3> fold_mapper (const arma::vec& sf)
  {
    // generate indices mapping SFS entries onto (un)folded spectrum
    // and set "out-of-bounds" entries to zero

    std::array<arma::uvec,3> f = 
    {
      arma::zeros<arma::uvec> (arma::size(sf)), // remapped indices
       arma::ones<arma::uvec> (arma::size(sf)), // weights
       arma::ones<arma::uvec> (arma::size(sf)), // in-bounds
    };

    for (unsigned s = 0; s < sf.n_elem; ++s)
    {
      if (fold && 2*s > sf.n_elem-1)
      {
        f[0].at(s) = sf.n_elem-s-1;
        f[2].at(s) = 0;
      } 
      else
      {
        f[0].at(s) = s;
        if (fold && 2*s == sf.n_elem-1)
          f[1].at(s) = 2;
      }
    }

    return f;
  }

  double likelihood (const arma::vec& inp, arma::vec& out, size_t& ns, const unsigned i2)
  {
    const unsigned mult = multiplier[block[i2]];

    if (mult == 0)
      return 0.;

    //arma::vec buffer (arma::size(inp), arma::fill::zeros);

    double den = 0.;

    // taking care not to iterate over entire SFS
    for (arma::sp_mat::const_col_iterator val0 = saf.begin_col(i2); 
         val0 != saf.end_col(i2); ++val0)
    {
      buffer.at(val0.row()) = folder[1].at(val0.row()) * 
              inp.at(folder[0].at(val0.row())) * (*val0);
      den += buffer.at(val0.row());
    }

    if (fabs(den) < arma::datum::eps) // SAF empty for this site
      return 0.;

    for (arma::sp_mat::const_col_iterator val0 = saf.begin_col(i2); 
         val0 != saf.end_col(i2); ++val0)
    {
      out.at(folder[0].at(val0.row())) += mult * buffer.at(val0.row())/den;
    }

    ns += mult; // track number of sites

    return mult * log(den);
  }

  void operator() (size_t start, size_t end)
  {
    for (size_t i = start; i != end; ++i)
      loglik += likelihood (sfs, upd, sites, i);
  }

  void join (const SFS1d& rhs)
  {
    loglik += rhs.loglik;
    upd += rhs.upd;
    sites += rhs.sites;
  }

  double EM (void)
  {
    upd.zeros();
    loglik = 0.;
    sites = 0;

    sfs %= arma::conv_to<arma::vec>::from(folder[2]);
    sfs /= arma::accu(sfs);

    RcppParallel::parallelReduce(0, block.n_elem, *this);
    //(*this)(0, multiplier.n_elem);

    sfs = upd / double(sites);

    return loglik;
  }

  bool EMaccel (void)
  {
    const double tol = 1e-5;
    const double mstep = 4.;

    arma::vec sfs0 = sfs;
    double    ll0  = loglikelihood;

    loglikelihood = EM ();
    arma::vec sfs1 = sfs;
    double    ll1  = loglikelihood,
              sr2  = arma::accu(arma::pow(sfs1 - sfs0,2));
    if (!arma::is_finite(sr2) || fabs(ll1 - ll0) < tol || sqrt(sr2) < tol)
      return true;

    if (!acceleration)
      return false;

    loglikelihood  = EM ();
    arma::vec sfs2 = sfs;
    double    ll2  = loglikelihood,
              sq2  = arma::accu(arma::pow(sfs2 - sfs1,2));
    if (!arma::is_finite(sq2) || fabs(ll2 - ll1) < tol || sqrt(sq2) < tol)
      return true;

    // accelerate
    double sv2   = arma::accu(arma::pow(sfs2 - 2*sfs1 + sfs0,2)),
           alpha = sqrt(sr2/sv2);
    alpha = std::max(stepmin, std::min(stepmax, alpha));
    sfs   = sfs0 + 2.*alpha*(sfs1 - sfs0) + alpha*alpha*(sfs2 - 2*sfs1 + sfs0);

    // project
    sfs  = arma::clamp(sfs, errtol, 1.-errtol); 
    sfs %= arma::conv_to<arma::vec>::from(folder[2]);
    sfs /= arma::accu(sfs);

    // stabilize and reset if needed
    loglikelihood = EM ();
    if (!arma::is_finite(loglikelihood) || loglikelihood < ll2)
    {
      loglikelihood = ll2;
      sfs = sfs2;
    }

    if (alpha == stepmax)
      stepmax *= mstep;

    return false;
  }
};

// [[Rcpp::export()]]
arma::mat sfs1d (const arma::sp_mat saf, const arma::uvec block, const unsigned num_boot, const bool fold)
{
  arma::uvec uniq_block = arma::unique(block);

  if (block.max() >= uniq_block.n_elem)
    Rcpp::stop("Block indices must be 0-based and contiguous");
  if (block.n_elem != saf.n_cols)
    Rcpp::stop("Dimension mismatch");

  arma::mat sfs (saf.n_rows, num_boot+1);
  arma::uvec multiplier (uniq_block.n_elem, arma::fill::ones);

  for (unsigned boot = 0; boot <= num_boot; ++boot)
  {
    if (boot)
    {
      multiplier.zeros();
      arma::uvec resample = arma::randi<arma::uvec>(uniq_block.n_elem, arma::distr_param(0, uniq_block.n_elem-1));
      for (auto b : resample)
        multiplier[b] += 1;
    }
    SFS1d estimate (saf, block, multiplier, fold);
    sfs.col(boot) = estimate.sfs;
  }

  return sfs;
}

// 2d SFS ... TODO make this a module

struct SFS2d : public RcppParallel::Worker
{
  const bool fold = false;
  const arma::uvec &multiplier, &block;
  const arma::sp_mat &saf0, &saf1; // these must be exponentiated!
  arma::mat sfs;
  double loglikelihood;

  const std::array<arma::umat,3> folder;

  private:
  arma::mat upd, buffer;
  double loglik;
  size_t sites;

  // settings
  bool   acceleration = true;        // use EM acceleration?
  double errtol = 1e-16;             // not using dynamic bounds
  double stepmax = 1., stepmin = 1.; // adaptive steplength

  public:

  SFS2d (const arma::sp_mat& saf0, const arma::sp_mat& saf1, const arma::uvec& block, const arma::uvec& multiplier, const bool fold)
    : fold (fold)
    , saf0 (saf0)
    , saf1 (saf1)
    , block (block)
    , multiplier (multiplier)
    , sfs (saf0.n_rows, saf1.n_rows, arma::fill::ones)
    , folder (fold_mapper(sfs))
    // accumulated quantities
    , loglik (0.)
    , upd (arma::size(sfs), arma::fill::zeros)
    , sites (0)
    // temporaries
    , buffer (arma::size(sfs))
  {
    const unsigned maxiter = 1000;

    if (saf0.n_cols != saf1.n_cols || block.n_elem != saf0.n_cols || block.max() >= multiplier.n_elem)
      Rcpp::stop("[Saf::SFS2d] Dimension mismatch"); 
    if (saf0.min() < 0. || saf1.min() < 0.)
      Rcpp::stop("[Saf::SFS2d] SAF must be non-negative");

    bool converged;
    for (unsigned iter=0; iter<maxiter; ++iter)
    {
      converged = EMaccel();
      if (converged) 
        break;
      fprintf(stderr, "[Saf::SFS2d]\titer %u, loglik = %f\n", iter, loglikelihood);
    }

    if (!converged)
      Rcpp::warning("[Saf::SFS2d] EM did not converge in maximum number of iterations");

    sfs.elem(arma::find(folder[2] == 0)).fill(arma::datum::nan);
    sfs *= double(sites);
  }

  SFS2d (const SFS2d& rhs, RcppParallel::Split)
    : fold (rhs.fold)
    , saf0 (rhs.saf0)
    , saf1 (rhs.saf1)
    , block (rhs.block)
    , multiplier (rhs.multiplier)
    , sfs (rhs.sfs)
    , folder (rhs.folder)
    // accumulated quantities
    , loglik (0.)
    , upd (arma::size(sfs), arma::fill::zeros)
    , sites (0)
    // temporaries
    , buffer (arma::size(sfs))
  {}

  std::array<arma::umat,3> fold_mapper (const arma::mat& sf)
  {
    // generate indices mapping SFS entries onto (un)folded spectrum
    // and set "out-of-bounds" entries to zero

    std::array<arma::umat,3> f = 
    {
      arma::zeros<arma::umat> (arma::size(sf)), // remapped indices
       arma::ones<arma::umat> (arma::size(sf)), // weights
       arma::ones<arma::umat> (arma::size(sf)), // is zero or not
    };

    for (unsigned s0 = 0; s0 < sf.n_rows; ++s0)
      for (unsigned s1 = 0; s1 < sf.n_cols; ++s1)
      {
        if (fold && 2*(s0+s1) > sf.n_rows+sf.n_cols-2)
        {
          f[0].at(s0,s1) = arma::sub2ind(arma::size(sf), sf.n_rows - s0 - 1, sf.n_cols - s1 - 1);
          f[2].at(s0,s1) = 0;
        } 
      else
      {
        f[0].at(s0,s1) = arma::sub2ind(arma::size(sf), s0, s1);
        if (fold && 2*s0 == sf.n_rows-1 && 2*s1 == sf.n_cols-1)
          f[1].at(s0,s1) = 2;
      }
    }

    return f;
  }

  double likelihood (const arma::mat& inp, arma::mat& out, size_t& ns, const unsigned i2)
  {
    const unsigned mult = multiplier[block[i2]];

    if (mult == 0)
      return 0.;

    //arma::mat buffer (arma::size(inp), arma::fill::zeros);

    double den = 0.;

    for (arma::sp_mat::const_col_iterator val0 = saf0.begin_col(i2); 
         val0 != saf0.end_col(i2); ++val0)
      for (arma::sp_mat::const_col_iterator val1 = saf1.begin_col(i2); 
          val1 != saf1.end_col(i2); ++val1)
    {
      buffer.at(val0.row(), val1.row()) = 
        folder[1].at(val0.row(),val1.row()) *
        inp.at(folder[0].at(val0.row(), val1.row())) * 
        (*val0) * (*val1);
      den += buffer.at(val0.row(), val1.row());
    }

    if (fabs(den) < arma::datum::eps) // SAF empty for this site
      return 0.;

    // taking care not to iterate over entire SFS
    for (arma::sp_mat::const_col_iterator val0 = saf0.begin_col(i2); 
         val0 != saf0.end_col(i2); ++val0)
      for (arma::sp_mat::const_col_iterator val1 = saf1.begin_col(i2); 
          val1 != saf1.end_col(i2); ++val1)
    {
      out.at(folder[0].at(val0.row(), val1.row())) += 
        mult * buffer.at(val0.row(), val1.row())/den;
    }

    ns += mult;

    return mult * log(den);
  }

  void operator() (size_t start, size_t end)
  {
    for (size_t i = start; i != end; ++i)
      loglik += likelihood (sfs, upd, sites, i);
  }

  void join (const SFS2d& rhs)
  {
    loglik += rhs.loglik;
    upd += rhs.upd;
    sites += rhs.sites;
  }

  double EM (void)
  {
    upd.zeros();
    loglik = 0.;
    sites = 0;

    sfs %= arma::conv_to<arma::mat>::from(folder[2]);
    sfs /= arma::accu(sfs);

    RcppParallel::parallelReduce(0, block.n_elem, *this);
    //(*this)(0, multiplier.n_elem);

    sfs = upd / double(sites);

    return loglik;
  }

  bool EMaccel (void)
  {
    const double tol = 1e-5;
    const double mstep = 4.;

    arma::mat sfs0 = sfs;
    double    ll0  = loglikelihood;

    loglikelihood  = EM ();
    arma::mat sfs1 = sfs;
    double    ll1  = loglikelihood,
              sr2  = arma::accu(arma::pow(sfs1 - sfs0,2));
    if (!arma::is_finite(sr2) || fabs(ll1 - ll0) < tol || sqrt(sr2) < tol)
      return true;

    if (!acceleration)
      return false;

    loglikelihood  = EM ();
    arma::mat sfs2 = sfs;
    double    ll2  = loglikelihood,
              sq2  = arma::accu(arma::pow(sfs2 - sfs1,2));
    if (!arma::is_finite(sq2) || fabs(ll2 - ll1) < tol || sqrt(sq2) < tol)
      return true;

    // accelerate
    double sv2   = arma::accu(arma::pow(sfs2 - 2*sfs1 + sfs0,2)),
           alpha = sqrt(sr2/sv2);
    alpha = std::max(stepmin, std::min(stepmax, alpha));
    sfs   = sfs0 + 2.*alpha*(sfs1 - sfs0) + alpha*alpha*(sfs2 - 2*sfs1 + sfs0);

    // project
    sfs  = arma::clamp(sfs, errtol, 1.-errtol);
    sfs %= arma::conv_to<arma::mat>::from(folder[2]);
    sfs /= arma::accu(sfs);

    // stabilize and reset if needed
    loglikelihood = EM ();
    if (!arma::is_finite(loglikelihood) || loglikelihood < ll2)
    {
      loglikelihood = ll2;
      sfs = sfs2;
    }

    if (alpha == stepmax)
      stepmax *= mstep;

    return false;
  }
};

// [[Rcpp::export()]]
arma::cube sfs2d (const arma::sp_mat saf0, const arma::sp_mat saf1, const arma::uvec block, const unsigned num_boot, const bool fold)
{
  arma::uvec uniq_block = arma::unique(block);

  if (block.max() >= uniq_block.n_elem)
    Rcpp::stop("Block indices must be 0-based and contiguous");
  if (block.n_elem != saf0.n_cols || saf0.n_cols != saf1.n_cols)
    Rcpp::stop("Dimension mismatch");

  arma::cube sfs (saf0.n_rows, saf1.n_rows, num_boot+1);
  arma::uvec multiplier (uniq_block.n_elem, arma::fill::ones);

  for (unsigned boot = 0; boot <= num_boot; ++boot)
  {
    if (boot)
    {
      multiplier.zeros();
      arma::uvec resample = arma::randi<arma::uvec>(uniq_block.n_elem, arma::distr_param(0, uniq_block.n_elem-1));
      for (auto b : resample)
        multiplier[b] += 1;
    }
    SFS2d estimate (saf0, saf1, block, multiplier, fold);
    sfs.slice(boot) = estimate.sfs;
  }

  return sfs;
}

// 3d SFS ...

struct SFS3d : public RcppParallel::Worker
{
  const bool fold = false;
  const arma::uvec &multiplier, &block;
  const arma::sp_mat &saf0, &saf1, &saf2; // these must be exponentiated!
  arma::cube sfs;
  double loglikelihood;

  const std::array<arma::ucube,3> folder;

  private:
  arma::cube upd, buffer;
  double loglik;
  size_t sites;

  // settings
  bool   acceleration = true;        // use EM acceleration?
  double errtol = 1e-16;             // not using dynamic bounds
  double stepmax = 1., stepmin = 1.; // adaptive steplength

  public:

  SFS3d (const arma::sp_mat& saf0, const arma::sp_mat& saf1, const arma::sp_mat& saf2, const arma::uvec& block, const arma::uvec& multiplier, const bool fold)
    : fold (fold)
    , saf0 (saf0)
    , saf1 (saf1)
    , saf2 (saf2)
    , block (block)
    , multiplier (multiplier)
    , sfs (saf0.n_rows, saf1.n_rows, saf2.n_rows, arma::fill::ones)
    , folder (fold_mapper(sfs))
    // accumulated quantities
    , loglik (0.)
    , upd (arma::size(sfs), arma::fill::zeros)
    , sites (0)
    // temporaries
    , buffer (arma::size(sfs))
  {
    const unsigned maxiter = 1000;

    if (saf0.n_cols != saf1.n_cols || saf0.n_cols != saf2.n_cols || block.n_elem != saf0.n_cols || block.max() >= multiplier.n_elem)
      Rcpp::stop("[Saf::SFS3d] Dimension mismatch"); 
    if (saf0.min() < 0. || saf1.min() < 0. || saf2.min() < 0.)
      Rcpp::stop("[Saf::SFS3d] SAF must be non-negative");

    bool converged;
    for (unsigned iter=0; iter<maxiter; ++iter)
    {
      converged = EMaccel();
      if (converged) 
        break;
      fprintf(stderr, "[Saf::SFS3d]\titer %u, loglik = %f\n", iter, loglikelihood);
    }

    if (!converged)
      Rcpp::warning("[Saf::SFS3d] EM did not converge in maximum number of iterations");

    sfs.elem(arma::find(folder[2] == 0)).fill(arma::datum::nan);
    sfs *= double(sites);
  }

  SFS3d (const SFS3d& rhs, RcppParallel::Split)
    : fold (rhs.fold)
    , saf0 (rhs.saf0)
    , saf1 (rhs.saf1)
    , saf2 (rhs.saf2)
    , block (rhs.block)
    , multiplier (rhs.multiplier)
    , sfs (rhs.sfs)
    , folder (rhs.folder)
    // accumulated quantities
    , loglik (0.)
    , upd (arma::size(sfs), arma::fill::zeros)
    , sites (0)
    // temporaries
    , buffer (arma::size(sfs))
  {}

  std::array<arma::ucube,3> fold_mapper (const arma::cube& sf)
  {
    // generate indices mapping SFS entries onto (un)folded spectrum
    // and set "out-of-bounds" entries to zero

    std::array<arma::ucube,3> f = 
    {
      arma::zeros<arma::ucube> (arma::size(sf)), // remapped indices
       arma::ones<arma::ucube> (arma::size(sf)), // weights
       arma::ones<arma::ucube> (arma::size(sf)), // is zero or not
    };

    for (unsigned s0 = 0; s0 < sf.n_rows; ++s0)
      for (unsigned s1 = 0; s1 < sf.n_cols; ++s1)
        for (unsigned s2 = 0; s2 < sf.n_slices; ++s2)
      {
        if (fold && 2*(s0+s1+s2) > sf.n_rows+sf.n_cols+sf.n_slices-3)
        {
          f[0].at(s0,s1,s2) = arma::sub2ind(arma::size(sf), sf.n_rows - s0 - 1, sf.n_cols - s1 - 1, sf.n_slices - s2 - 1);
          f[2].at(s0,s1,s2) = 0;
        } 
      else
      {
        f[0].at(s0,s1,s2) = arma::sub2ind(arma::size(sf), s0, s1, s2);
        if (fold && 2*s0 == sf.n_rows-1 && 2*s1 == sf.n_cols-1 && 2*s2 == sf.n_slices-1)
          f[1].at(s0,s1,s2) = 2;
      }
    }

    return f;
  }

  double likelihood (const arma::cube& inp, arma::cube& out, size_t& ns, const unsigned i2)
  {
    const unsigned mult = multiplier[block[i2]];

    if (mult == 0)
      return 0.;

    //arma::cube buffer (arma::size(inp), arma::fill::zeros);

    double den = 0.;

    for (arma::sp_mat::const_col_iterator val0 = saf0.begin_col(i2); 
         val0 != saf0.end_col(i2); ++val0)
      for (arma::sp_mat::const_col_iterator val1 = saf1.begin_col(i2); 
          val1 != saf1.end_col(i2); ++val1)
        for (arma::sp_mat::const_col_iterator val2 = saf2.begin_col(i2); 
            val2 != saf2.end_col(i2); ++val2)
    {
      buffer.at(val0.row(), val1.row(), val2.row()) = 
        folder[1].at(val0.row(), val1.row(), val2.row()) *
        inp.at(folder[0].at(val0.row(), val1.row(), val2.row())) * 
        (*val0) * (*val1) * (*val2);
      den += buffer.at(val0.row(), val1.row(), val2.row());
    }

    if (fabs(den) < arma::datum::eps) // SAF empty for this site
      return 0.;

    // taking care not to iterate over entire SFS
    for (arma::sp_mat::const_col_iterator val0 = saf0.begin_col(i2); 
         val0 != saf0.end_col(i2); ++val0)
      for (arma::sp_mat::const_col_iterator val1 = saf1.begin_col(i2); 
          val1 != saf1.end_col(i2); ++val1)
        for (arma::sp_mat::const_col_iterator val2 = saf2.begin_col(i2); 
            val2 != saf2.end_col(i2); ++val2)
    {
      out.at(folder[0].at(val0.row(), val1.row(), val2.row())) += 
        mult * buffer.at(val0.row(), val1.row(), val2.row())/den;
    }

    ns += mult;

    return mult * log(den);
  }

  void operator() (size_t start, size_t end)
  {
    for (size_t i = start; i != end; ++i)
      loglik += likelihood (sfs, upd, sites, i);
  }

  void join (const SFS3d& rhs)
  {
    loglik += rhs.loglik;
    upd += rhs.upd;
    sites += rhs.sites;
  }

  double EM (void)
  {
    upd.zeros();
    loglik = 0.;
    sites = 0;

    sfs %= arma::conv_to<arma::cube>::from(folder[2]);
    sfs /= arma::accu(sfs);

    RcppParallel::parallelReduce(0, block.n_elem, *this);
    //(*this)(0, multiplier.n_elem);

    sfs = upd / double(sites);

    return loglik;
  }

  bool EMaccel (void)
  {
    const double tol = 1e-5;
    const double mstep = 4.;

    arma::cube sfs0 = sfs;
    double    ll0  = loglikelihood;

    loglikelihood  = EM ();
    arma::cube sfs1 = sfs;
    double    ll1  = loglikelihood,
              sr2  = arma::accu(arma::pow(sfs1 - sfs0,2));
    if (!arma::is_finite(sr2) || fabs(ll1 - ll0) < tol || sqrt(sr2) < tol)
      return true;

    if (!acceleration)
      return false;

    loglikelihood  = EM ();
    arma::cube sfs2 = sfs;
    double    ll2  = loglikelihood,
              sq2  = arma::accu(arma::pow(sfs2 - sfs1,2));
    if (!arma::is_finite(sq2) || fabs(ll2 - ll1) < tol || sqrt(sq2) < tol)
      return true;

    // accelerate
    double sv2   = arma::accu(arma::pow(sfs2 - 2*sfs1 + sfs0,2)),
           alpha = sqrt(sr2/sv2);
    alpha = std::max(stepmin, std::min(stepmax, alpha));
    sfs   = sfs0 + 2.*alpha*(sfs1 - sfs0) + alpha*alpha*(sfs2 - 2*sfs1 + sfs0);

    // project
    sfs  = arma::clamp(sfs, errtol, 1.-errtol);
    sfs %= arma::conv_to<arma::cube>::from(folder[2]);
    sfs /= arma::accu(sfs);

    // stabilize and reset if needed
    loglikelihood = EM ();
    if (!arma::is_finite(loglikelihood) || loglikelihood < ll2)
    {
      loglikelihood = ll2;
      sfs = sfs2;
    }

    if (alpha == stepmax)
      stepmax *= mstep;

    return false;
  }
};

// [[Rcpp::export()]]
Rcpp::NumericVector sfs3d (const arma::sp_mat saf0, const arma::sp_mat saf1, const arma::sp_mat saf2, const arma::uvec block, const unsigned num_boot, const bool fold)
{
  arma::uvec uniq_block = arma::unique(block);

  if (block.max() >= uniq_block.n_elem)
    Rcpp::stop("Block indices must be 0-based and contiguous");
  if (block.n_elem != saf0.n_cols || saf0.n_cols != saf1.n_cols || saf0.n_cols != saf2.n_cols)
    Rcpp::stop("Dimension mismatch");

  std::vector<int> dims = { int(saf0.n_rows), int(saf1.n_rows), int(saf2.n_rows), int(num_boot+1) };
  Rcpp::NumericVector sfs = Rcpp::NumericVector(Rcpp::Dimension(Rcpp::wrap(dims))); // wrap needed to convert to SEXP

  arma::uvec multiplier (uniq_block.n_elem, arma::fill::ones);

  for (unsigned boot = 0; boot <= num_boot; ++boot)
  {
    if (boot)
    {
      multiplier.zeros();
      arma::uvec resample = arma::randi<arma::uvec>(uniq_block.n_elem, arma::distr_param(0, uniq_block.n_elem-1));
      for (auto b : resample)
        multiplier[b] += 1;
    }
    SFS3d estimate (saf0, saf1, saf2, block, multiplier, fold);

    // copying into R array ... linear index for sfs[i,j,k,l] is i + j*rows + k*rows*cols + l*rows*cols*slices
    for (unsigned i=0; i<saf0.n_rows; ++i)
      for (unsigned j=0; j<saf1.n_rows; ++j)
        for (unsigned k=0; k<saf2.n_rows; ++k)
          sfs[i + j*saf0.n_rows + k*saf0.n_rows*saf1.n_rows + 
                  boot*saf0.n_rows*saf1.n_rows*saf2.n_rows] = estimate.sfs.at(i,j,k);
  }

  return sfs;
}

// Fst

struct Fst : public RcppParallel::Worker
{
  const arma::sp_mat &saf0, &saf1; // these must be exponentiated!

  arma::mat fsfs, numerator, denominator, fst;

  const std::array<arma::umat,3> folder;

  private:
  arma::mat &rfst;
  arma::mat buffer;

  public:

  Fst (const arma::sp_mat& saf0, const arma::sp_mat& saf1, const arma::mat& sfs)
    : saf0 (saf0)
    , saf1 (saf1)
    , fsfs (sfs)
    , numerator (arma::size(fsfs), arma::fill::zeros)
    , denominator (arma::size(fsfs), arma::fill::zeros)
    , fst (saf0.n_cols, 2)
    , folder (fold_mapper(fsfs, numerator, denominator))
    // references for workers
    , rfst (fst)
    // temporaries
    , buffer (arma::size(fsfs))
  {
    const unsigned maxiter = 1000;

    if (saf0.n_cols != saf1.n_cols  || 
        saf0.n_rows != fsfs.n_rows  || 
        saf1.n_rows != fsfs.n_cols  )
      Rcpp::stop("[Fst] Dimension mismatch"); 
    if (saf0.min() < 0. || saf1.min() < 0. || fsfs.min() < 0.)
      Rcpp::stop("[Fst] SAF/SFS must be non-negative");

    RcppParallel::parallelFor(0, fst.n_rows, *this);
    //(*this)(0, fst.n_rows);
  }
  
  Fst (const Fst& rhs, RcppParallel::Split)
    : saf0 (rhs.saf0)
    , saf1 (rhs.saf1)
    , fsfs (rhs.fsfs)
    , numerator (rhs.numerator)
    , denominator (rhs.denominator)
    , fst (0, 0) //don't allocate
    , folder (rhs.folder)
    // references for workers
    , rfst (rhs.rfst)
    // temporaries
    , buffer (arma::size(fsfs))
  {}

  std::array<arma::umat,3> fold_mapper (arma::mat& sfs, arma::mat& num, arma::mat& den)
  {
    // -generate indices mapping SFS entries onto folded spectrum
    // -calculate numerator and denominator of Fst estimator
    // -fold input SFS

    const bool fold = true;

    std::array<arma::umat,3> f = 
    {
      arma::zeros<arma::umat> (arma::size(sfs)), // remapped indices
       arma::ones<arma::umat> (arma::size(sfs)), // weights
       arma::ones<arma::umat> (arma::size(sfs)), // is zero or not
    };

    for (unsigned s0 = 0; s0 < sfs.n_rows; ++s0)
      for (unsigned s1 = 0; s1 < sfs.n_cols; ++s1)
      {
        if (fold && 2*(s0+s1) > sfs.n_rows+sfs.n_cols-2)
        {
          f[0].at(s0,s1) = arma::sub2ind(arma::size(sfs), sfs.n_rows-s0-1, sfs.n_cols-s1-1);
          f[2].at(s0,s1) = 0;
        } 
        else
        {
          f[0].at(s0,s1) = arma::sub2ind(arma::size(sfs), s0, s1);
          if (fold && 2*s0 == sfs.n_rows-1 && 2*s1 == sfs.n_cols-1)
            f[1].at(s0,s1) = 2;

          // precompute numerator and denominator of Bhatia estimator
          double p0 = double(s0)/double(sfs.n_rows),
                 p1 = double(s1)/double(sfs.n_cols);
          num.at(s0,s1) = std::pow(p0-p1, 2) - 
             p0*(1.-p0)/double(sfs.n_rows-1) - 
             p1*(1.-p1)/double(sfs.n_cols-1);
          den.at(s0,s1) = p0*(1.-p1) + p1*(1.-p0);
        }
    }

    // fold sfs 
    sfs.replace(arma::datum::nan, 0.);
    arma::mat sfs_fold (arma::size(sfs), arma::fill::zeros);
    for (unsigned s0 = 0; s0 < sfs.n_rows; ++s0)
      for (unsigned s1 = 0; s1 < sfs.n_cols; ++s1)
        sfs_fold.at(f[0].at(s0,s1)) += sfs.at(s0,s1);
    sfs = sfs_fold;

    return f;
  }

  arma::rowvec::fixed<2> bhatia (const arma::mat& sfs, const unsigned i2)
  {
    // numerator and denominator of Bhatia et al 2015 Genome Research estimator
    
    arma::rowvec::fixed<2> out = arma::zeros<arma::rowvec>(2);

    double den = 0.;

    for (arma::sp_mat::const_col_iterator val0 = saf0.begin_col(i2); 
         val0 != saf0.end_col(i2); ++val0)
      for (arma::sp_mat::const_col_iterator val1 = saf1.begin_col(i2); 
          val1 != saf1.end_col(i2); ++val1)
    {
      buffer.at(val0.row(), val1.row()) = 
        folder[1].at(val0.row(),val1.row()) *
        sfs.at(folder[0].at(val0.row(), val1.row())) * 
        (*val0) * (*val1);
      den += buffer.at(val0.row(), val1.row());
    }

    // taking care not to iterate over entire SFS
    for (arma::sp_mat::const_col_iterator val0 = saf0.begin_col(i2); 
         val0 != saf0.end_col(i2); ++val0)
      for (arma::sp_mat::const_col_iterator val1 = saf1.begin_col(i2); 
          val1 != saf1.end_col(i2); ++val1)
    {
      out.at(0) += 
        numerator.at(folder[0].at(val0.row(), val1.row())) *
        buffer.at(val0.row(), val1.row())/den;
      out.at(1) += 
        denominator.at(folder[0].at(val0.row(), val1.row())) *
        buffer.at(val0.row(), val1.row())/den;
    }

    if (fabs(den) < arma::datum::eps) // SAF empty for this site
      out.fill(arma::datum::nan);

    return out;
  }

  void operator() (size_t start, size_t end)
  {
    for (size_t i = start; i != end; ++i)
      fst.row(i) = bhatia(fsfs, i);
  }
};

// [[Rcpp::export()]]
arma::mat FST (const arma::sp_mat saf0, const arma::sp_mat saf1, const arma::mat sfs)
{
  Fst estimate (saf0, saf1, sfs);

  return estimate.fst;
}

// Theta

struct Theta : public RcppParallel::Worker
{
  const arma::sp_mat &saf; // these must be exponentiated!

  arma::vec fsfs, segregating, watterson, pairwise;
  arma::mat theta;

  const std::array<arma::uvec,3> folder;

  private:
  arma::mat &rtheta;
  arma::vec buffer;

  public:

  Theta (const arma::sp_mat& saf, const arma::vec& sfs)
    : saf (saf)
    , fsfs (sfs)
    , segregating (arma::size(fsfs), arma::fill::zeros)
    , watterson (arma::size(fsfs), arma::fill::zeros)
    , pairwise (arma::size(fsfs), arma::fill::zeros)
    , theta (saf.n_cols, 3)
    , folder (fold_mapper(fsfs, segregating, watterson, pairwise))
    // references for workers
    , rtheta (theta)
    // temporaries
    , buffer (arma::size(fsfs))
  {
    if (saf.n_rows != fsfs.n_elem)
      Rcpp::stop("[Theta] Dimension mismatch"); 
    if (saf.min() < 0. || fsfs.min() < 0.)
      Rcpp::stop("[Theta] SAF/SFS must be non-negative");

    RcppParallel::parallelFor(0, theta.n_rows, *this);
    //(*this)(0, fst.n_rows);
  }
  
  Theta (const Theta& rhs, RcppParallel::Split)
    : saf (rhs.saf)
    , fsfs (rhs.fsfs)
    , segregating (rhs.segregating)
    , watterson (rhs.watterson)
    , pairwise (rhs.pairwise)
    , theta (0, 0) //don't allocate
    , folder (rhs.folder)
    // references for workers
    , rtheta (rhs.rtheta)
    // temporaries
    , buffer (arma::size(fsfs))
  {}

  std::array<arma::uvec,3> fold_mapper (arma::vec& sfs, arma::vec& seg, arma::vec& wat, arma::vec& pairs)
  {
    // -generate indices mapping SFS entries onto folded spectrum
    // -pre calculate thetas
    // -fold input SFS

    const bool fold = true;

    std::array<arma::uvec,3> f = 
    {
      arma::zeros<arma::uvec> (arma::size(sfs)), // remapped indices
       arma::ones<arma::uvec> (arma::size(sfs)), // weights
       arma::ones<arma::uvec> (arma::size(sfs)), // is zero or not
    };

    //double a1 = arma::accu(1./arma::regspace(1.,double(sfs.n_elem-1))),
    //       a2 = arma::accu(1./arma::pow(arma::regspace(1.,double(sfs.n_elem-1)),2)),
    //       b1 = double(sfs.n_elem + 1)/(3. * double(sfs.n_elem - 1)),
    //       b2 = 2.*double(sfs.n_elem*sfs.n_elem + sfs.n_elem + 3)/double(9 * sfs.n_elem * (sfs.n_elem - 1)),
    //       c1 = b1 - 1/a1,
    //       c2 = b2 - double(sfs.n_elem+2)/(a1*double(sfs.n_elem)) + a2/(a1*a1),
    //       e1 = c1/a1,
    //       e2 = c2/(a1*a1 + a2);
    //    seg == 0 ? 0. : pi - watterson / sqrt(e1*seg + e2*seg*(seg-1)); // is this additive? no it is not.
    
    double a1 = arma::accu(1./arma::regspace(1.,double(sfs.n_elem-1)));

    for (unsigned s = 0; s < sfs.n_elem; ++s)
    {
      if (fold && 2*s > sfs.n_elem-1)
      {
        f[0].at(s) = sfs.n_elem-s-1;
        f[2].at(s) = 0;
      } 
      else
      {
        f[0].at(s) = s;
        if (fold && 2*s == sfs.n_elem-1)
          f[1].at(s) = 2;

        seg.at(s) = s > 0 && s < sfs.n_elem ? 1. : 0.;
        wat.at(s) = s > 0 && s < sfs.n_elem ? 1./a1 : 0.;
        pairs.at(s) = s * (sfs.n_elem - s);
      }
    }

    // fold sfs 
    sfs.replace(arma::datum::nan, 0.);
    arma::vec sfs_fold (arma::size(sfs), arma::fill::zeros);
    for (unsigned s = 0; s < sfs.n_elem; ++s)
      sfs_fold.at(f[0].at(s)) += sfs.at(s);
    sfs = sfs_fold;

    return f;
  }

  arma::rowvec::fixed<3> theta_stats (const arma::mat& sfs, const unsigned i2)
  {
    arma::rowvec::fixed<3> out = arma::zeros<arma::rowvec>(3);

    double den = 0.;

    for (arma::sp_mat::const_col_iterator val0 = saf.begin_col(i2); 
         val0 != saf.end_col(i2); ++val0)
    {
      buffer.at(val0.row()) = 
        folder[1].at(val0.row()) *
        sfs.at(folder[0].at(val0.row())) * 
        (*val0);
      den += buffer.at(val0.row());
    }

    // taking care not to iterate over entire SFS
    for (arma::sp_mat::const_col_iterator val0 = saf.begin_col(i2); 
         val0 != saf.end_col(i2); ++val0)
    {
      out.at(0) += 
        segregating.at(folder[0].at(val0.row())) *
        buffer.at(val0.row())/den;
      out.at(1) += 
        watterson.at(folder[0].at(val0.row())) *
        buffer.at(val0.row())/den;
      out.at(2) += 
        pairwise.at(folder[0].at(val0.row())) *
        buffer.at(val0.row())/den;
    }

    if (fabs(den) < arma::datum::eps) // SAF empty for this site
      out.fill(arma::datum::nan);

    return out;
  }

  void operator() (size_t start, size_t end)
  {
    for (size_t i = start; i != end; ++i)
      theta.row(i) = theta_stats(fsfs, i);
  }
};

// [[Rcpp::export()]]
arma::mat theta (const arma::sp_mat saf, const arma::mat sfs)
{
  Theta estimate (saf, sfs);

  return estimate.theta;
}

// [[Rcpp::export()]]
arma::mat slider (arma::mat inp, arma::uvec coord, const unsigned window, const unsigned step)
{
  // sum over a sliding window. Pretty inefficient at the moment

  arma::uvec order = arma::stable_sort_index(coord);

  coord = coord.elem(order);
  inp   = inp.rows(order);

  const unsigned coord_max = coord.max(),
                 coord_min = coord.min();

  unsigned lower = coord_min,
           upper = coord_min + window;

  arma::mat out (3 + inp.n_cols, 0);

  while (lower <= coord_max)
  {
    arma::uvec sites = arma::find(coord >= lower && coord < upper); // has to be more efficient way, only need min/max

    arma::vec window (out.n_rows);

    window.at(0) = double(lower);
    window.at(1) = double(upper);
    window.at(2) = double(sites.n_elem);
    for (unsigned i=0; i<inp.n_cols; ++i)
      window.at(3+i) = sites.n_elem ? arma::accu(inp.submat(sites.min(), i, sites.max(), i)) : arma::datum::nan;

    out.insert_cols (out.n_cols, window); //this is very inefficient

    lower += step;
    upper += step;
  }

  return arma::trans(out);
}

namespace RobinsonHill
{
  // these are horrendously bloated but whatever
  
  double Di (arma::vec hfs, arma::umat config)
  {
    // Unbiased estimate of expected linkage disequillibrium 'E[D]'
    //
    // Projection: 1pop to 2chr
    // D_i = E[p_AB p_ab - p_Ab p_aB]
    
    if (config.max() < 2)
      Rcpp::stop("[RobinsonHill::Di] Insufficient number of chromosomes");

    auto dp = utils::hfs1d_downproject (hfs, config, 2);

    hfs    = std::get<0>(dp);
    config = std::get<1>(dp)[0];

    arma::uvec bin1001 = arma::find (config.row(0) == 0 && config.row(1) == 0 && config.row(2) == 1),
               bin0110 = arma::find (config.row(0) == 1 && config.row(1) == 1 && config.row(2) == 0);

    return 1./2. * hfs.at(bin1001[0]) - 
           1./2. * hfs.at(bin0110[0]);
  }

  double DiDi (arma::vec hfs, arma::umat config)
  {
    // Unbiased estimate of variance of linkage disequillibrium 'E[D^2]'
    //
    // Projection: 1pop to 4chr
    // D_i^2 = E[p_AB^2 p_ab^2 - 2 p_AB p_Ab p_aB p_ab + p_Ab^2 p_aB^2]
    
    if (config.max() < 4)
      Rcpp::stop("[RobinsonHill::DiDi] Insufficient number of chromosomes");

    auto dp = utils::hfs1d_downproject (hfs, config, 4);

    hfs    = std::get<0>(dp);
    config = std::get<1>(dp)[0];

    arma::uvec bin2002 = arma::find (config.row(0) == 0 && config.row(1) == 0 && config.row(2) == 2),
               bin1111 = arma::find (config.row(0) == 1 && config.row(1) == 1 && config.row(2) == 1),
               bin0220 = arma::find (config.row(0) == 2 && config.row(1) == 2 && config.row(2) == 0);

    return 1./6. * hfs.at(bin2002[0]) - 
           2./24. * hfs.at(bin1111[0]) + 
           1./6. * hfs.at(bin0220[0]);
  }

  double DiDj (arma::mat hfs, arma::umat config0, arma::umat config1)
  {
    // Unbiased estimate of covariance of linkage disequillibrium 'E[D_i D_j]'
    //
    // Projection: 2pop to 2x2chr
    // D_i D_j   = E[p_AB p_ab q_AB q_ab - p_AB p_ab q_Ab q_aB - p_Ab p_aB q_AB q_ab + p_Ab p_aB q_Ab q_aB]
    
    if (config0.max() < 2 || config1.max() < 2)
      Rcpp::stop("[RobinsonHill::DiDj] Insufficient number of chromosomes");

    auto dp = utils::hfs2d_downproject (hfs, config0, 2, config1, 2);

    hfs     = std::get<0>(dp);
    config0 = std::get<1>(dp)[0];
    config1 = std::get<1>(dp)[1];

    arma::uvec bin1001a = arma::find (config0.row(0) == 0 && config0.row(1) == 0 && config0.row(2) == 1),
               bin1001b = arma::find (config1.row(0) == 0 && config1.row(1) == 0 && config1.row(2) == 1),
               bin0110a = arma::find (config0.row(0) == 1 && config0.row(1) == 1 && config0.row(2) == 0),
               bin0110b = arma::find (config1.row(0) == 1 && config1.row(1) == 1 && config1.row(2) == 0);

    return 1./4. * hfs.at(bin1001a[0], bin1001b[0]) - 
           1./4. * hfs.at(bin1001a[0], bin0110b[0]) - 
           1./4. * hfs.at(bin0110a[0], bin1001b[0]) +
           1./4. * hfs.at(bin0110a[0], bin0110b[0]);
  }

  double DiZii (arma::mat hfs, arma::umat config)
  {
    // Unbiased estimate of single population, frequency-weighted linkage disequillibrium 'E[D_i z_ii]'
    //
    // Projection: 1pop to 2,3,4chr
    // D_i z_ii  = E[ (p_AB p_ab - p_Ab p_aB) (1 - 2 p_AB - 2 p_Ab) (1 - 2 p_AB - 2 p_aB) ] =
    //             E[ (p_AB p_ab - p_Ab p_aB) (2 p_AB + 2 p_ab - 4 (p_AB p_ab - p_Ab p_aB) - 1) ] =
    //             E[2 p_AB^2 p_ab + 2 p_AB p_ab^2 - 2 p_AB p_Ab p_aB - 2 p_ab p_Ab p_aB - 4 D_i^2 - D_i]
    
    if (config.max() < 4)
      Rcpp::stop("[RobinsonHill::DiZii] Insufficient number of chromosomes");

    auto dp = utils::hfs1d_downproject (hfs, config, 3);

    double di = RobinsonHill::Di (hfs, config),
           didi = RobinsonHill::DiDi (hfs, config);

    hfs     = std::get<0>(dp);
    config  = std::get<1>(dp)[0];

    arma::uvec bin1002 = arma::find (config.row(0) == 0 && config.row(1) == 0 && config.row(2) == 2),
               bin2001 = arma::find (config.row(0) == 0 && config.row(1) == 0 && config.row(2) == 1),
               bin0111 = arma::find (config.row(0) == 1 && config.row(1) == 1 && config.row(2) == 1),
               bin1110 = arma::find (config.row(0) == 1 && config.row(1) == 1 && config.row(2) == 0);

    return 2./3. * hfs.at(bin1002[0]) +
           2./3. * hfs.at(bin2001[0]) - 
           2./6. * hfs.at(bin0111[0]) - 
           2./6. * hfs.at(bin1110[0]) -
           di - 4. * didi;
  }

  double DiZij (arma::mat hfs2d, arma::umat config0, arma::umat config1)
  {
    // Unbiased estimate of two population, frequency-weighted linkage disequillibrium 'E[D_i z_ij]'
    //
    // Projection: 2pop to 2,3x1
    // D_i z_ij  = E[ (p_AB p_ab - p_Ab p_aB) (1 - 2 p_AB - 2 p_Ab) (1 - 2 q_AB - 2 q_aB) ] =
    //             E[ (p_AB p_ab - p_Ab p_aB) (1 - 2 p_AB - 2 p_Ab - 2 q_AB + 4 p_AB q_AB + 4 p_Ab q_AB - 2 q_aB + 4 p_AB q_aB + 4 p_Ab q_aB)] =
    //             D_i + 
    //             - 2 p_AB^2 p_ab + 2 p_AB p_Ab p_aB - 2 p_Ab p_AB p_ab + 2 p_Ab^2 p_aB 
    //             - 2 q_AB p_AB p_ab + 2 q_AB p_Ab p_aB - 2 q_aB p_AB p_ab + 2 q_aB p_Ab p_aB 
    //             + 4 p_AB^2 q_AB p_ab - 4 p_AB q_AB p_Ab p_aB + 4 p_Ab q_AB p_AB p_ab - 4 p_Ab^2 q_AB p_aB + 4 p_AB^2 q_aB p_ab - 4 p_AB q_aB p_Ab p_aB + 4 p_Ab q_aB p_AB p_ab - 4 p_Ab^2 q_aB p_aB

    if (config0.max() < 3 || config1.max() < 1)
      Rcpp::stop("[RobinsonHill::DiZij] Insufficient number of chromosomes");

    auto dp3 = utils::hfs2d_downproject (hfs2d, config0, 3, config1, 1),
         dp2 = utils::hfs2d_downproject (hfs2d, config0, 2, config1, 1);

    arma::mat hfs3 = std::get<0>(dp3),
              hfs2 = std::get<0>(dp2);

    arma::umat config0a = std::get<1>(dp3)[0],
               config1a = std::get<1>(dp3)[1],
               config0b = std::get<1>(dp2)[0],
               config1b = std::get<1>(dp2)[1];

    // marginalize over second population by summing over columns
    arma::vec mhfs3 = arma::sum(hfs3, 1),
              mhfs2 = arma::sum(hfs2, 1);

    arma::uvec bin1002 = arma::find (config0a.row(0) == 0 && config0a.row(1) == 0 && config0a.row(2) == 2),
               bin0111 = arma::find (config0a.row(0) == 1 && config0a.row(1) == 1 && config0a.row(2) == 1),
               bin1101 = arma::find (config0a.row(0) == 1 && config0a.row(1) == 0 && config0a.row(2) == 1),
               bin0210 = arma::find (config0a.row(0) == 2 && config0a.row(1) == 1 && config0a.row(2) == 0),
               bin1001 = arma::find (config0b.row(0) == 0 && config0b.row(1) == 0 && config0b.row(2) == 1),
               bin0110 = arma::find (config0b.row(0) == 1 && config0b.row(1) == 1 && config0b.row(2) == 0),
               bin0010 = arma::find (config1a.row(0) == 0 && config1a.row(1) == 1 && config1a.row(2) == 0),
               bin0001 = arma::find (config1a.row(0) == 0 && config1a.row(1) == 0 && config1a.row(2) == 1);

    double di = 0.5 * (mhfs2.at(bin1001[0]) - mhfs2.at(bin0110[0]));

    return di +
      // 3x0
      -2./3. * mhfs3.at(bin1002[0]) +
       2./6. * mhfs3.at(bin0111[0]) +
      -2./6. * mhfs3.at(bin1101[0]) +
       2./3. * mhfs3.at(bin0210[0]) +
     // 2x1
      -2./2. * hfs2.at(bin1001[0], bin0001[0]) +
       2./2. * hfs2.at(bin0110[0], bin0001[0]) +
      -2./2. * hfs2.at(bin1001[0], bin0010[0]) +
       2./2. * hfs2.at(bin0110[0], bin0010[0]) +
     // 3x1
       4./3. * hfs3.at(bin1002[0], bin0001[0]) +
      -4./6. * hfs3.at(bin0111[0], bin0001[0]) +
       4./6. * hfs3.at(bin1101[0], bin0001[0]) +
      -4./3. * hfs3.at(bin0210[0], bin0001[0]) +
       4./3. * hfs3.at(bin1002[0], bin0010[0]) +
      -4./6. * hfs3.at(bin0111[0], bin0010[0]) +
       4./6. * hfs3.at(bin1101[0], bin0010[0]) +
      -4./3. * hfs3.at(bin0210[0], bin0010[0]);
  }

  double DiZji (arma::mat hfs2d, arma::umat config0, arma::umat config1)
  {
    // Unbiased estimate of two population, frequency-weighted linkage disequillibrium 'E[D_i z_ji]'
    //
    // Projection: 2pop to 2,3x1
    // D_i z_ji  = E[ (p_AB p_ab - p_Ab p_aB) (1 - 2 q_AB - 2 q_Ab) (1 - 2 p_AB - 2 p_aB) ] =
    //             E[ (p_AB p_ab - p_Ab p_aB) (1 - 2 q_AB - 2 q_Ab - 2 p_AB + 4 q_AB p_AB + 4 q_Ab p_AB - 2 p_aB + 4 q_AB p_aB + 4 q_Ab p_aB)] =
    //             D_i + 
    //             - 2 p_AB^2 p_ab + 2 p_Ab p_aB p_AB - 2 p_aB p_AB p_ab + 2 p_aB^2 p_Ab 
    //             - 2 p_AB p_ab q_AB + 2 p_Ab p_aB q_AB - 2 p_AB p_ab q_Ab + 2 p_Ab p_aB q_Ab
    //             + 4 q_AB p_AB^2 p_ab - 4 q_AB p_AB p_Ab p_aB + 4 q_Ab p_AB^2 p_ab - 4 q_Ab p_AB p_Ab p_aB
    //             + 4 q_AB p_aB p_AB p_ab - 4 q_AB p_aB^2 p_Ab + 4 q_Ab p_aB p_AB p_ab - 4 q_Ab p_aB^2 p_Ab

    if (config0.max() < 3 || config1.max() < 1)
      Rcpp::stop("[RobinsonHill::DiZji] Insufficient number of chromosomes");

    auto dp3 = utils::hfs2d_downproject (hfs2d, config0, 3, config1, 1);
    auto dp2 = utils::hfs2d_downproject (hfs2d, config0, 2, config1, 1);

    arma::mat hfs3 = std::get<0>(dp3),
              hfs2 = std::get<0>(dp2);

    arma::umat config0a = std::get<1>(dp3)[0],
               config1a = std::get<1>(dp3)[1],
               config0b = std::get<1>(dp2)[0],
               config1b = std::get<1>(dp2)[1];

    // marginalize over second population by summing over columns
    arma::vec mhfs3 = arma::sum(hfs3, 1),
              mhfs2 = arma::sum(hfs2, 1);

    arma::uvec bin1002 = arma::find (config0a.row(0) == 0 && config0a.row(1) == 0 && config0a.row(2) == 2),
               bin0111 = arma::find (config0a.row(0) == 1 && config0a.row(1) == 1 && config0a.row(2) == 1),
               bin1011 = arma::find (config0a.row(0) == 0 && config0a.row(1) == 1 && config0a.row(2) == 1),
               bin0120 = arma::find (config0a.row(0) == 1 && config0a.row(1) == 2 && config0a.row(2) == 0),
               bin1001 = arma::find (config0b.row(0) == 0 && config0b.row(1) == 0 && config0b.row(2) == 1),
               bin0110 = arma::find (config0b.row(0) == 1 && config0b.row(1) == 1 && config0b.row(2) == 0),
               bin0100 = arma::find (config1a.row(0) == 1 && config1a.row(1) == 0 && config1a.row(2) == 0),
               bin0001 = arma::find (config1a.row(0) == 0 && config1a.row(1) == 0 && config1a.row(2) == 1);

    double di = 0.5 * (mhfs2.at(bin1001[0]) - mhfs2.at(bin0110[0]));

    return di +
      // 3x0
      -2./3. * mhfs3.at(bin1002[0]) +
       2./6. * mhfs3.at(bin0111[0]) +
      -2./6. * mhfs3.at(bin1011[0]) +
       2./3. * mhfs3.at(bin0120[0]) +
     // 2x1
      -2./2. * hfs2.at(bin1001[0], bin0001[0]) +
       2./2. * hfs2.at(bin0110[0], bin0001[0]) +
      -2./2. * hfs2.at(bin1001[0], bin0100[0]) +
       2./2. * hfs2.at(bin0110[0], bin0100[0]) +
     // 3x1
       4./3. * hfs3.at(bin1002[0], bin0001[0]) +
      -4./6. * hfs3.at(bin0111[0], bin0001[0]) +
       4./6. * hfs3.at(bin1011[0], bin0001[0]) +
      -4./3. * hfs3.at(bin0120[0], bin0001[0]) +
       4./3. * hfs3.at(bin1002[0], bin0100[0]) +
      -4./6. * hfs3.at(bin0111[0], bin0100[0]) +
       4./6. * hfs3.at(bin1011[0], bin0100[0]) +
      -4./3. * hfs3.at(bin0120[0], bin0100[0]);
  }

  double DiZjk (arma::cube hfs3d, arma::umat config0, arma::umat config1, arma::umat config2)
  {
    // Unbiased estimate of two population, frequency-weighted linkage disequillibrium 'E[D_i z_jk]'
    //
    // Projection: 3pop to 2x1x1
    // D_i z_jk  = E[ (p_AB p_ab - p_Ab p_aB) (1 - 2 q_AB - 2 q_Ab) (1 - 2 r_AB - 2 r_aB) ] =
    //             E[ (p_AB p_ab - p_Ab p_aB) (1 - 2 q_AB - 2 q_Ab - 2 r_AB - 2 r_aB + 4 q_AB r_AB + 4 q_Ab r_AB + 4 q_AB r_aB + 4 q_Ab r_aB) ] =    
    //             D_i +
    //             - 2 p_AB p_ab q_AB + 2 p_Ab p_aB q_AB - 2 p_AB p_ab q_Ab + 2 p_Ab p_aB q_Ab 
    //             - 2 p_AB p_ab r_AB + 2 p_Ab p_aB r_AB - 2 p_AB p_ab r_aB + 2 p_Ab p_aB r_aB
    //             + 4 q_AB r_AB p_AB p_ab - 4 q_AB r_AB p_Ab p_aB + 4 q_Ab r_AB p_AB p_ab - 4 q_Ab r_AB p_Ab p_aB 
    //             + 4 q_AB r_aB p_AB p_ab - 4 q_AB r_aB p_Ab p_aB + 4 q_Ab r_aB p_AB p_ab - 4 q_Ab r_aB p_Ab p_aB

    if (config0.max() < 3 || config1.max() < 1 || config2.max() < 1)
      Rcpp::stop("[RobinsonHill::DiZjk] Insufficient number of chromosomes");

    auto dp = utils::hfs3d_downproject (hfs3d, config0, 2, config1, 1, config2, 1);

    arma::cube hfs = std::get<0>(dp);

    config0 = std::get<1>(dp)[0];
    config1 = std::get<1>(dp)[1];
    config2 = std::get<1>(dp)[2];

    // marginalize over third populations
    arma::mat mhfs_ij (config0.n_cols, config1.n_cols, arma::fill::zeros),
              mhfs_ik (config0.n_cols, config2.n_cols, arma::fill::zeros);
    for (unsigned i=0; i<config0.n_cols; ++i)
      for (unsigned j=0; j<config1.n_cols; ++j)
        for (unsigned k=0; k<config2.n_cols; ++k)
        {
          mhfs_ij.at(i,j) += hfs.at(i,j,k);
          mhfs_ik.at(i,k) += hfs.at(i,j,k);
        }
    arma::vec mhfs = arma::sum(mhfs_ij, 1);

    arma::uvec bin1001 = arma::find (config0.row(0) == 0 && config0.row(1) == 0 && config0.row(2) == 1),
               bin0110 = arma::find (config0.row(0) == 1 && config0.row(1) == 1 && config0.row(2) == 0),
               bin0100a = arma::find (config1.row(0) == 1 && config1.row(1) == 0 && config1.row(2) == 0),
               bin0001a = arma::find (config1.row(0) == 0 && config1.row(1) == 0 && config1.row(2) == 1),
               bin0010b = arma::find (config2.row(0) == 0 && config2.row(1) == 1 && config2.row(2) == 0),
               bin0001b = arma::find (config2.row(0) == 0 && config2.row(1) == 0 && config2.row(2) == 1);

    double di = 0.5 * (mhfs.at(bin1001[0]) - mhfs.at(bin0110[0]));

    return di +
      // 2x1x0
      -2./2. * mhfs_ij.at(bin1001[0], bin0001a[0]) +
       2./2. * mhfs_ij.at(bin0110[0], bin0001a[0]) +
      -2./2. * mhfs_ij.at(bin1001[0], bin0100a[0]) +
       2./2. * mhfs_ij.at(bin0110[0], bin0100a[0]) +
     // 2x0x1
      -2./2. * mhfs_ik.at(bin1001[0], bin0001b[0]) +
       2./2. * mhfs_ik.at(bin0110[0], bin0001b[0]) +
      -2./2. * mhfs_ik.at(bin1001[0], bin0010b[0]) +
       2./2. * mhfs_ik.at(bin0110[0], bin0010b[0]) +
     // 2x1x1
       4./2. * hfs.at(bin1001[0], bin0001a[0], bin0001b[0]) +
      -4./2. * hfs.at(bin0110[0], bin0001a[0], bin0001b[0]) +
       4./2. * hfs.at(bin1001[0], bin0100a[0], bin0001b[0]) +
      -4./2. * hfs.at(bin0110[0], bin0100a[0], bin0001b[0]) +
       4./2. * hfs.at(bin1001[0], bin0001a[0], bin0010b[0]) +
      -4./2. * hfs.at(bin0110[0], bin0001a[0], bin0010b[0]) +
       4./2. * hfs.at(bin1001[0], bin0100a[0], bin0010b[0]) +
      -4./2. * hfs.at(bin0110[0], bin0100a[0], bin0010b[0]);
  }

  double DiZkj (arma::cube hfs3d, arma::umat config0, arma::umat config1, arma::umat config2)
  {
    // Unbiased estimate of two population, frequency-weighted linkage disequillibrium 'E[D_i z_kj]'
    //
    // Projection: 3pop to 2x1x1
    // D_i z_kj  = E[ (p_AB p_ab - p_Ab p_aB) (1 - 2 r_AB - 2 r_Ab) (1 - 2 q_AB - 2 q_aB) ] =
    //             E[ (p_AB p_ab - p_Ab p_aB) (1 - 2 r_AB - 2 r_Ab - 2 q_AB - 2 q_aB + 4 r_AB q_AB + 4 r_Ab q_AB + 4 r_AB q_aB + 4 r_Ab q_aB) ] =    
    //             D_i +
    //             - 2 p_AB p_ab r_AB + 2 p_Ab p_aB r_AB - 2 p_AB p_ab r_Ab + 2 p_Ab p_aB r_Ab 
    //             - 2 p_AB p_ab q_AB + 2 p_Ab p_aB q_AB - 2 p_AB p_ab q_aB + 2 p_Ab p_aB q_aB
    //             + 4 r_AB q_AB p_AB p_ab - 4 r_AB q_AB p_Ab p_aB + 4 r_Ab q_AB p_AB p_ab - 4 r_Ab q_AB p_Ab p_aB 
    //             + 4 r_AB q_aB p_AB p_ab - 4 r_AB q_aB p_Ab p_aB + 4 r_Ab q_aB p_AB p_ab - 4 r_Ab q_aB p_Ab p_aB

    if (config0.max() < 3 || config1.max() < 1 || config2.max() < 1)
      Rcpp::stop("[RobinsonHill::DiZkj] Insufficient number of chromosomes");

    auto dp = utils::hfs3d_downproject (hfs3d, config0, 2, config1, 1, config2, 1);

    arma::cube hfs = std::get<0>(dp);

    config0 = std::get<1>(dp)[0];
    config1 = std::get<1>(dp)[1];
    config2 = std::get<1>(dp)[2];

    // marginalize over third populations
    arma::mat mhfs_ij (config0.n_cols, config1.n_cols, arma::fill::zeros),
              mhfs_ik (config0.n_cols, config2.n_cols, arma::fill::zeros);
    for (unsigned i=0; i<config0.n_cols; ++i)
      for (unsigned j=0; j<config1.n_cols; ++j)
        for (unsigned k=0; k<config2.n_cols; ++k)
        {
          mhfs_ij.at(i,j) += hfs.at(i,j,k);
          mhfs_ik.at(i,k) += hfs.at(i,j,k);
        }
    arma::vec mhfs = arma::sum(mhfs_ij, 1);

    arma::uvec bin1001 = arma::find (config0.row(0) == 0 && config0.row(1) == 0 && config0.row(2) == 1),
               bin0110 = arma::find (config0.row(0) == 1 && config0.row(1) == 1 && config0.row(2) == 0),
               bin0010a = arma::find (config1.row(0) == 0 && config1.row(1) == 1 && config1.row(2) == 0),
               bin0001a = arma::find (config1.row(0) == 0 && config1.row(1) == 0 && config1.row(2) == 1),
               bin0100b = arma::find (config2.row(0) == 1 && config2.row(1) == 0 && config2.row(2) == 0),
               bin0001b = arma::find (config2.row(0) == 0 && config2.row(1) == 0 && config2.row(2) == 1);

    double di = 0.5 * (mhfs.at(bin1001[0]) - mhfs.at(bin0110[0]));

    return di +
      // 2x1x0
      -2./2. * mhfs_ij.at(bin1001[0], bin0001a[0]) +
       2./2. * mhfs_ij.at(bin0110[0], bin0001a[0]) +
      -2./2. * mhfs_ij.at(bin1001[0], bin0010a[0]) +
       2./2. * mhfs_ij.at(bin0110[0], bin0010a[0]) +
     // 2x0x1
      -2./2. * mhfs_ik.at(bin1001[0], bin0001b[0]) +
       2./2. * mhfs_ik.at(bin0110[0], bin0001b[0]) +
      -2./2. * mhfs_ik.at(bin1001[0], bin0100b[0]) +
       2./2. * mhfs_ik.at(bin0110[0], bin0100b[0]) +
     // 2x1x1
       4./2. * hfs.at(bin1001[0], bin0001a[0], bin0001b[0]) +
      -4./2. * hfs.at(bin0110[0], bin0001a[0], bin0001b[0]) +
       4./2. * hfs.at(bin1001[0], bin0010a[0], bin0001b[0]) +
      -4./2. * hfs.at(bin0110[0], bin0010a[0], bin0001b[0]) +
       4./2. * hfs.at(bin1001[0], bin0001a[0], bin0100b[0]) +
      -4./2. * hfs.at(bin0110[0], bin0001a[0], bin0100b[0]) +
       4./2. * hfs.at(bin1001[0], bin0010a[0], bin0100b[0]) +
      -4./2. * hfs.at(bin0110[0], bin0010a[0], bin0100b[0]);
  }

}

// [[Rcpp::export]]
arma::mat hfs1d (const arma::sp_mat shf, const arma::uvec weights)
{
  arma::uvec block = arma::regspace<arma::uvec>(0, shf.n_cols-1);
  SFS1d estimate (shf, block, weights, false);
  return estimate.sfs;
}

// [[Rcpp::export]]
Rcpp::NumericVector basis (arma::cube hfs3d, arma::umat config0, arma::umat config1, arma::umat config2)
{
  // calculate terms in the Robinson Hill basis involving E[D_i]

  // marginalize
  arma::mat hfs_ij (config0.n_cols, config1.n_cols, arma::fill::zeros),
            hfs_ik (config0.n_cols, config2.n_cols, arma::fill::zeros);
  arma::vec hfs_i (config0.n_cols, arma::fill::zeros);
  for (unsigned i=0; i<config0.n_cols; ++i)
    for (unsigned j=0; j<config1.n_cols; ++j)
      for (unsigned k=0; k<config2.n_cols; ++k)
      {
        hfs_ij.at(i,j) += hfs3d.at(i,j,k);
        hfs_ik.at(i,k) += hfs3d.at(i,j,k);
        hfs_i.at(i)    += hfs3d.at(i,j,k);
      }

  // various statistics; D_i, D_i^2, D_i D_j, D_i D_k, D_i z_ii, D_i z_ij, D_i z_ji,  D_i z_ik, D_i z_ki, D_i z_jk, D_i z_kj
  return Rcpp::NumericVector::create(
      Rcpp::_["Di"] = RobinsonHill::Di(hfs_i, config0),
      Rcpp::_["DiDi"] = RobinsonHill::DiDi(hfs_i, config0),
      Rcpp::_["DiDj"] = RobinsonHill::DiDj(hfs_ij, config0, config1),
      Rcpp::_["DiDk"] = RobinsonHill::DiDj(hfs_ik, config0, config2),
      Rcpp::_["DiZii"] = RobinsonHill::DiZii(hfs_i, config0),
      Rcpp::_["DiZij"] = RobinsonHill::DiZij(hfs_ij, config0, config1),
      Rcpp::_["DiZji"] = RobinsonHill::DiZji(hfs_ij, config0, config1),
      Rcpp::_["DiZik"] = RobinsonHill::DiZij(hfs_ik, config0, config2),
      Rcpp::_["DiZki"] = RobinsonHill::DiZji(hfs_ik, config0, config2),
      Rcpp::_["DiZjk"] = RobinsonHill::DiZjk(hfs3d, config0, config1, config2),
      Rcpp::_["DiZkj"] = RobinsonHill::DiZkj(hfs3d, config0, config1, config2)
      );
}
