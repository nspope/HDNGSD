
template <>
class sptensor
{
  const arma::uvec dims;
  const arma::umat subs;
  arma::vec vals;

  public:
  sptensor (arma::uvec dims, arma::umat subs, arma::vec vals)
    : dims (dims)
    , subs (subs)
    , vals (subs.n_cols, arma::fill::zeros)
  {
    if (subs.n_rows != dims.n_elem || dims.n_elem != dim)
      Rcpp::stop ("[sptensor] Invalid dimensions");
    for (unsigned i=0; i<dims.n_elem; ++i)
      if (subs.row(i).max() >= dims.at(i))
        Rcpp::stop ("[sptensor] Subscript exceeds size of dimension");
  }

  void operator+= (const sptensor& rhs)
  {
    vals += rhs.vals;
  }

  void zeros (void)
  {
    vals.zeros();
  }

  double& at (const arma::uvec& index)
  {
    return vals.at(arma::find());
  }
}//truly this will be a pain in the ass

namespace Rcpp {
  namespace traits {
    SEXP wrap (const sptensor& obj);


  }
}
