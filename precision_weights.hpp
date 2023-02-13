#ifndef PRECISION_WEIGHTS_H
#define PRECISION_WEIGHTS_H

// [[Rcpp::depends(BH)]]
// #include <boost/multiprecision/cpp_dec_float.hpp>

// namespace mp = boost::multiprecision;

template <class forIt>
double calc_sum_w(forIt first, forIt last, const double& max_sum_w) {
  // mp::cpp_dec_float_100 sum_w(0.0);
  long double sum_w = 0.0;
  for (auto it = first; it != last; ++it) {
    long double exponent = *it - max_sum_w;
    
    //mp::cpp_dec_float_100 exp_prec(exponent);
    long double exp_prec = exponent;
    auto w_prec = std::exp(exp_prec); // other branch with no exponential
  
    sum_w += w_prec;
      *it = w_prec; //.convert_to<double>();
  }
  
    return sum_w; //.convert_to<double>();
}
#endif /* PRECISION_WEIGHTS_H */
