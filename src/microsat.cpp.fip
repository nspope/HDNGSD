// to use this effectively I need to use only variable sites
// and "augment" likelihoods/counts for samples that are missing snps but have msats

struct Microsat
{
  // for use with AdmixtureMkII
  
  const GenotypeLikelihood& GL;
  const arma::uvec sample_index;
  const arma::cube L;
  arma::mat oldfreq, freq, hwe, counts;

  Microsat (const GenotypeLikelihood& GL, const arma::uvec samples, arma::mat Qstart, const arma::cube like)
    : GL (GL),
    , sample_index (samples)
    , L (like)
  {
    if (sample_index.n_elem != L.n_slices)
      Rcpp::stop ("[Microsat] Dimension mismatch");

    Qstart.each_col([](arma::vec& x) { x /= arma::accu(x); });
    freq = arma::randu<arma::mat>(L.n_rows, Qstart.n_rows);
    freq.each_col([](arma::vec& x) { x /= arma::accu(x); });
    oldfreq = freq;
    hwe = expected_frequencies(Qstart, freq);
    counts = expected_counts(Qstart, freq, hwe);
  }

  double EM (const arma::mat& Q, arma::mat& Qup)
  {
    oldfreq = freq;
    hwe = expected_frequencies(Q, oldfreq);
    counts = expected_counts(Q, oldfreq, hwe);
    freq = update_frequencies(Q, oldfreq, hwe, counts);
    Qup += update_admixture(Q, oldfreq, hwe, counts);

    return loglik (Q, oldfreq, hwe);
  }

  arma::mat expected_frequencies (const arma::mat& Q, const arma::mat& F)
  {
    // compute new frequencies
    if (Q.n_cols != sample_index.n_elem)
      Rcpp::stop ("[Microsat::expected_frequencies] Dimension mismatch");
  
    // compute HWE frequencies per sample (input const Q&)
    arma::mat out (F.n_rows, Q.n_cols, arma::fill::zeros);
    for (unsigned s=0; s<L.n_slices; ++s)
      for (unsigned k=0; k<Q.n_cols; ++k)
        for (unsigned i=0; i<F.n_rows; ++i)
          out.at(i,s) += Q.at(k,s) * F.at(i,k);
    return out;
  }
  
  arma::mat expected_counts (const arma::mat& Q, const arma::mat& F, const arma::mat& S)
  {
    // compute new frequencies

    // posterior expectation of allele counts for each sample
    arma::mat out (F.n_rows, L.n_slices, arma::fill::zeros);
    for (unsigned s=0; s<L.n_slices; ++s)
    {
      double den = 0., p = 0.;
      if (GL.ploidy(sample_index[s])==2)
      {
        for (unsigned i=0; i<F.n_rows; ++i)
        {
          p = L.at(i,i,s) * std::pow(S.at(i,s), 2);
          den += p;
          out.at(i,s) += 2. * p;
          for (unsigned j=i+1; j<F.n_rows; ++j)
          {
            p = L.at(i,j,s) * S.at(j,s) * S.at(i,s);
            den += 2. * p;
            out.at(i,s) += p;
            out.at(j,s) += p;
          }
        }
      } 
      else 
      {
        for (unsigned i=0; i<F.n_rows; ++i)
        {
          p = L.at(i,i,s) * S.at(i,s);
          den += p;
          out.at(i,s) += p;
        }
      }
      out.col(s) /= 2.*den; //scales haploid expectation to 0.5 as for snps
    }
    return out;
  }
  
  arma::mat update_admixture (const arma::mat& Q, const arma::mat& F, const arma::mat& S, const arma::mat& E)
  {
    // contribute to admixture coeff update (input Qup&, const Q&)
    if (Q.n_cols != F.n_cols || Q.n_rows != E.n_rows)
      Rcpp::stop ("[Microsat::update_admixture] Dimension mismatch");
  
    Qup = arma::zeros<arma::mat>(arma::size(Q));
    for (unsigned s=0; s<Q.n_cols; ++s)
      for (unsigned k=0; k<Q.n_rows; ++k)
        for (unsigned i=0; i<E.n_rows; ++i)
          Qup.at(k,s) += E.at(i,s) * Q.at(k,s) * F.at(i,k) / S.at(i,s);
    return Qup;
  }
  
  arma::mat update_frequencies (const arma::mat& Q, const arma::mat& F, const arma::mat& S, const arma::mat& E)
  {
    // compute new frequencies (input const Q&)
    if (Q.n_cols != F.n_cols || Q.n_rows != E.n_rows)
      Rcpp::stop ("[Microsat::update_frequencies] Dimension mismatch");
  
    arma::mat Fup = arma::zeros<arma::mat>(arma::size(F));
    for (unsigned k=0; k<F.n_cols; ++k)
    {
      arma::vec A (F.n_rows, arma::fill::zeros);
      for (unsigned i=0; i<F.n_rows; ++i)
      {
        for (unsigned s=0; s<Q.n_cols; ++s)
        {
          A.at(i) += E.at(i,s) * Q.at(k,s) * F.at(i,k) / S.at(i,s);
        }
      }
      Fup.col(k) = A / arma::accu(A);
    }
    return Fup;
  }

  double loglik (const arma::mat& Q, const arma::mat& F, const arma::mat& S)
  {
    // loglikelihood
    double ll = 0.;
    for (unsigned s=0; s<L.n_slices; ++s)
    {
      double den = 0., p = 0.;
      if (GL.ploidy(sample_index[s])==2)
      {
        for (unsigned i=0; i<F.n_rows; ++i)
        {
          p = L.at(i,i,s) * std::pow(S.at(i,s), 2);
          den += p;
          for (unsigned j=i+1; j<F.n_rows; ++j)
          {
            p = L.at(i,j,s) * S.at(j,s) * S.at(i,s);
            den += 2. * p;
          }
        }
      } 
      else 
      {
        for (unsigned i=0; i<F.n_rows; ++i)
        {
          p = L.at(i,i,s) * S.at(i,s);
          den += p;
        }
      }
    }
    return ll;
  }
};

struct AdmixtureMkII : public RcppParallel::Worker
{ 
  //TODO for now, no acceleration for msats
  
  public:
  const GenotypeLikelihood &GL;				// contains likelihoods, ploidy, dimensions
  const arma::uvec sample_index;      // only use these samples
  const unsigned K;   								// number of clusters, number of chromosomes to hold out
  arma::mat Qmat, Fmat, loglik;       // admixture, frequencies, per site loglik and cross-validation
  arma::cube Acoef, Bcoef;						// per-site coefficients (could move to private eventually)
  std::vector<Microsat> msat;         // microsatellite loci
	double loglikelihood = 0.;					// total loglikelihood
  unsigned iter = 0;

  private:
  bool   acceleration = true; // use EM acceleration?
  double errtol = 1e-8;// not using dynamic bounds
  double stepmax = 1., stepmin = 1.; // adaptive steplength

	// references for slaves
  const arma::uvec &rsample_index;
  const arma::mat &rQmat, &rFmat;
  arma::mat &rloglik;
  arma::cube &rAcoef, &rBcoef;

	//temporaries
	arma::mat Q0, Q1, Q2, F0, F1, F2, Qd1, Qd2, Qd3, Fd1, Fd2, Fd3;

  public:
  AdmixtureMkII (const GenotypeLikelihood& GL, const arma::uvec samples, const arma::mat Qstart, const std::vector<arma::cube> microsat)
    : GL (GL)
    , sample_index (samples)
    , K (Qstart.n_rows)
    , Qmat (K, sample_index.n_elem)
    , Fmat (K, GL.sites, arma::fill::randu)
    , loglik (sample_index.n_elem, GL.sites)
    , Acoef (K, sample_index.n_elem, GL.sites)
    , Bcoef (K, sample_index.n_elem, GL.sites)
    // read-write references for workers 
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
      Rcpp::stop ("Admixture: must have at least one cluster");
    if (K > 999)
      Rcpp::stop ("Admixture: number of clusters capped at 999");
    if (arma::max(sample_index) >= GL.samples)
      Rcpp::stop ("Admixture: invalid entries in sample index");
    if (Qstart.n_rows != Qmat.n_rows || Qstart.n_cols != Qmat.n_cols)
      Rcpp::stop ("Admixture: starting values for Q are of wrong dimension");

    const unsigned maxiter = 5000;

    Qmat = Qstart;
    Qmat.each_col([](arma::vec& x) { x /= arma::accu(x); });

    Fmat.randu();
    Fmat = arma::clamp(Fmat, errtol, 1.-errtol);

    for (auto ms : microsat)
      msat.emplace_back(GL, sample_index, Qmat, ms);
		
		for (iter=0; iter<maxiter; ++iter)
    {
      bool converged = EMaccel(Qmat, Fmat);
      //TODO make mutables clearer; should always be passed as arguments (e.g. loglikelihood, cvscore!
      if((iter+1) % 10 == 0) 
        std::cout << "Admixture: iteration " << iter+1 << ", loglik = " << loglikelihood << std::endl;
			if(converged)
				break;
    }
  
		if (iter == maxiter)
			Rcpp::warning("Admixture: did not converge in maximum number of iterations");
  }

  Admixture (const Admixture& rhs, RcppParallel::Split)
    : GL (rhs.GL)
    , K (rhs.K)
    // references that interface with slaves
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
             i2 = i;

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
             i2 = i;

    arma::vec::fixed<3> p = GL.handler(j2,i2);

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

    double msloglik = 0.;
    arma::mat msQ = arma::zeros<arma::mat>(arma::size(Q));
    for (auto ms : msat)
      msloglik += ms.EM(Qmat, msQ);

    RcppParallel::parallelFor(0, GL.sites, *this); 
    //(*this)(0, GL.sites);

    F = arma::sum(Acoef, 1) / (arma::sum(Acoef, 1) + arma::sum(Bcoef, 1));
    Q = arma::sum(Acoef, 2) + arma::sum(Bcoef, 2) + msQ;
    Q.each_col([](arma::vec& x) { x /= arma::accu(x); }); // make sum_k q_{jk} = 1

    return arma::accu(loglik) + msloglik; 
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
