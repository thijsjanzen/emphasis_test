#include <vector>
#include <iostream>
#include <chrono>
#include "emphasis.hpp"
#include "model.hpp"


//' function to perform one step of the E-M algorithm
//' @param brts vector of branching times
//' @param init_pars vector of initial parameter files
//' @param sample_size number of samples
//' @param maxN maximum number of failed trees
//' @param soc number of lineages at the root/crown (1/2)
//' @param max_missing maximum number of species missing
//' @param max_lambda maximum speciation rate
//' @param lower_bound vector of lower bound values for optimization, should
//' be equal in length to the vector of init_pars
//' @param upper_bound vector of upper bound values for optimization, should
//' be equal in length to the vector of init_pars
//' @param xtol_rel relative tolerance for optimization
//' @param num_threads number of threads used.
//' @return a list with the following components: 
//' \itemize{
//'  \item{trees}{list of trees}
//'  \item{rejected}{number of rejected trees}
//'  \item{rejected_overruns}{number of trees rejected due to too large size}
//'  \item{rejected_lambda}{number of trees rejected due to lambda errors}
//'  \item{rejected_zero_weights}{number of trees rejected to zero weight}
//'  \item{time_elapsed}{time used}
//'  \item{weights}{vector of weights}
//'  \item{fhat}{vector of fhat values}
//'  \item{logf}{vector of logf values}
//'  \item{logg}{vector of logg values}
//' }
//' @export
// [[Rcpp::export(name = "e_cpp")]]
void rcpp_mce(const std::vector<double>& brts,
              const std::vector<double>& init_pars,      
              int sample_size,
              int maxN,
              int soc,
              int max_missing,               
              double max_lambda,             
              const std::vector<double>& lower_bound,  
              const std::vector<double>& upper_bound,  
              double xtol_rel,                     
              int num_threads)
{
  auto model = emphasis::Model(lower_bound, upper_bound);

  auto E = emphasis::E_step(sample_size,
                            maxN,
                            init_pars,
                            brts,
                            model,
                            soc,
                            max_missing,
                            max_lambda,
                            num_threads);
  return;
}


std::vector<double> rcpp_mce_nothrow(const std::vector<double>& brts,       
                                     const std::vector<double>& init_pars,      
                                     int sample_size,
                                     int maxN,
                                     int soc,
                                     int max_missing,               
                                     double max_lambda,             
                                     const std::vector<double>& lower_bound,  
                                     const std::vector<double>& upper_bound,  
                                     double xtol_rel,                     
                                     int num_threads)
{
  auto model = emphasis::Model(lower_bound, upper_bound);
  
  try {
  auto E = emphasis::E_step(sample_size,
                            maxN,
                            init_pars,
                            brts,
                            model,
                            soc,
                            max_missing,
                            max_lambda,
                            num_threads);
  std::vector<double> out = {E.fhat, 
                             double(E.rejected_lambda),
                             double(E.rejected_overruns), 
                             double(E.rejected_zero_weights)};
  return out;
  } catch (const emphasis::emphasis_error_E& E) {
    return {42, //E.E_.fhat, 
            double(E.E_.rejected_lambda),
            double(E.E_.rejected_overruns), 
            double(E.E_.rejected_zero_weights)};
  }
}

// [[Rcpp::export]]
std::vector<std::vector<double>> rcpp_mce_grid(const std::vector<std::vector<double>> pars_R,
                                  const std::vector<double>& brts,       
                                  int sample_size,
                                  int maxN,
                                  int soc,
                                  int max_missing,               
                                  double max_lambda,             
                                  const std::vector<double>& lower_bound,  
                                  const std::vector<double>& upper_bound,  
                                  double xtol_rel,                     
                                  int num_threads) {
  
  std::vector< std::vector<double> > results(pars_R.size(), std::vector<double>(4, 0.0));
  
  // copy parameters to C++ object, for multithreading
  std::vector< std::vector<double> > pars(pars_R.size());
  std::vector<double> row_entry(pars_R[0].size());
  for (int i = 0; i < pars_R.size(); ++i) {
    for(int j = 0; j < pars_R[0].size(); ++j) {
      row_entry[j] = pars_R[i][j];
    }
    pars[i] = row_entry;
  }
  
//  const int grainsize = pars.size() / std::max<unsigned>(1, std::min<unsigned>(std::thread::hardware_concurrency(), num_threads));
 // tbb::parallel_for(tbb::blocked_range<unsigned>(0, pars.size(), grainsize), [&](const tbb::blocked_range<unsigned>& r) {
      
    for (unsigned i = 0; i < pars.size(); ++i) {
      ///  std::cout << i << "\n";
      std::vector<double> local_pars = pars[i];
      results[i] = rcpp_mce_nothrow(brts,       
                                     local_pars,      
                                     sample_size,
                                     maxN,
                                     soc,
                                     max_missing,               
                                     max_lambda,             
                                     lower_bound,  
                                     upper_bound,  
                                     xtol_rel,                     
                                     num_threads);
   }
  
  return results;
}


int main() {
    
    std::vector<double> brts = {35.012473, 32.530357, 30.880633, 30.399471, 23.095866, 18.049627, 11.039829, 10.894030, 8.479031,  8.289813,  7.980300, 7.711562,  6.137002,  5.493738,  4.191253,  3.078151,
        3.026430,  2.456891,  1.836070,  1.262732}; // read_brts();
    
    int num_pars = 25;
    
    emphasis::reng_t reng_ = make_random_engine<emphasis::reng_t>();
    std::vector< std::vector< double >> pars;
    std::uniform_real_distribution<double> dist;
    for (int i = 0; i < num_pars; ++i) {
        std::vector<double> row = {dist(reng_),
                                   dist(reng_) * 2,
                                   -0.04 + 0.08 * dist(reng_),
            0.0
        };
        pars.push_back(row);
    }
    
    auto lb = {0.0, 0.0, -0.05, 0.0};
    auto ub = {1.0, 2.0, 0.05, 0.0};
    
    auto clock_start = std::chrono::system_clock::now();
    
    
    rcpp_mce_grid(pars,
                  brts,
                  100, // sample_size
                  100, // maxN
                  2, // soc
                  1000, //max_missing
                  1000, // max_lambda
                  lb, // lower_bound
                  ub, // upper_bound
                  0.1, // xtol_rel
                  1); // num_threads
    
    auto clock_end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = clock_end - clock_start;
    std::cout << "this took: " << elapsed_seconds.count() << "seconds\n";
    
    return 0;
}
