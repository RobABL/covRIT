#include <Rcpp.h>
#include <algorithm>
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

// [[Rcpp::export]]
bool has_inter(List const& interaction, DataFrame const& data, LogicalVector const& isFactor){
  int nrows = as<NumericVector>(data[0]).size();
  int inter_sz = interaction.size();

  for(int i = 0;i<nrows;++i){ // iterate over data rows
    for(int j = 0;j<inter_sz;++j){ // iterate over interaction attributes
      List attr_inter = as<List>(interaction[j]);
      int attr_idx = attr_inter["attr_idx"];

      double instance_val = as<NumericVector>(data[attr_idx])[i];
      std::vector<double> values = as<std::vector<double> >(attr_inter["values"]);

      if(isFactor[attr_idx]){ // Categorical attribute
        if(std::find(values.begin(),values.end(),instance_val) == values.end()){ // doesn't contain the value
          return false;
        }
      }
      else{ // Continuous attribute
        if(instance_val < values[0] || instance_val > values[1]){
          return false;
        }
      }
    }
  }
  
  return true;
}
