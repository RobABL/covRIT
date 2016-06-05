#include <Rcpp.h>
#include "interaction.h"
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <functional>
#include <stack>
#include <unordered_map>
#include <sys/types.h>
#include <unistd.h>

// [[Rcpp::plugins(cpp11)]]

using namespace std;
typedef _Bind<uniform_int_distribution<int>(mersenne_twister_engine<long unsigned int, 32ul, 624ul, 397ul, 31ul, 2567483615ul, 11ul, 4294967295ul, 7ul, 2636928640ul, 15ul, 4022730752ul, 18ul, 1812433253ul>)> Generator;

Generator create_PRNG(Rcpp::DataFrame const& x){
  int nrows = Rcpp::as<vector<double> >(x[0]).size();
  auto seed = chrono::high_resolution_clock::now().time_since_epoch().count()*::getpid();
  Generator int_rand = bind(uniform_int_distribution<int>(0,nrows-1),mt19937(seed));
  return int_rand;
}

void insert_interaction(unordered_map<string,Interaction>& map,Interaction inter,
Rcpp::List const& datas,vector<double> const& theta){
  string repres = inter.as_string();
  
  if(map.count(repres) > 0) // element already in map, do nothing
    return;
    
  if(inter.check_for_map(datas,theta))
    map[repres] = inter;
}

vector<double> random_instance(Rcpp::DataFrame const& x, Generator& prng){
  // Select random row number
  int row = prng();
  
  vector<double> instance(x.size());
  for(int i=0;i<x.size();++i){
    vector<double> col = Rcpp::as<vector<double> >(x[i]);
    instance[i] = col[row];
  }
  
  return instance;
}

// DFS search on tree. Places results in leaves parameter.
void tree(unordered_map<string,Interaction>& leaves,Rcpp::List const& datas, int cls,
vector<double> const& theta, vector<bool> const& factor,double epsilon_cont,double epsilon_cat,
vector<double> const& spans, int depth,int min_inter_sz,int branch,int split_nb,Generator& prng,bool es){
  stack<Interaction, vector<Interaction> > frontier;
  Rcpp::DataFrame class_data = Rcpp::as<Rcpp::DataFrame>(datas[cls]);
  
  // Init root
  int nb_class = datas.size();
  vector<double> root_instance = random_instance(class_data,prng);
  Interaction root(root_instance,factor,spans,epsilon_cont,epsilon_cat,0,nb_class);
  frontier.push(root);

  while(!frontier.empty()){
    Interaction parent = frontier.top();
    frontier.pop();
    int c_depth = parent.get_depth();
    int eff_branch = (c_depth <= 0 || c_depth%split_nb != 0) ? 1 : branch;
    
    int nb_valid = 0; // number of valid children for this parent node
    for(int i=0;i<eff_branch;++i){
      int instance_nb = prng();
      Interaction child(parent);
      child.intersect(instance_nb,class_data);
      child.set_depth(c_depth + 1);
      if(child.size() >= min_inter_sz && child.check_prevs(datas,theta,es)){ // child is valid
        if(child.size() == min_inter_sz || child.get_depth() == depth){ // child is leaf
          insert_interaction(leaves,child,datas,theta);
        }
        else{ // child is valid but not leaf
          frontier.push(child); 
        }
        nb_valid++;
      }
    }
    
    if(nb_valid == 0) // If no child is valid, add parent as a leaf
      insert_interaction(leaves,parent,datas,theta);
  }
}

// [[Rcpp::export]]
Rcpp::List cpp_cov_RIT(Rcpp::List const& datas,Rcpp::NumericVector const& theta,
Rcpp::LogicalVector const& factor,double epsilon_cont,double epsilon_cat, Rcpp::NumericVector const& spans,
int n_trees, int depth,int branch,int split_nb,int min_inter_sz,bool es){
  
  unordered_map<string,Interaction> res;
  
  // Transform types
  vector<bool> cFactor = Rcpp::as<vector<bool> >(factor);
  vector<double> c_theta = Rcpp::as<vector<double> >(theta);
  vector<double> c_spans = Rcpp::as<vector<double> >(spans);
  
  for(int c=0;c<datas.size();++c){// Iterate over classes
    Rcpp::DataFrame class_data = Rcpp::as<Rcpp::DataFrame>(datas[c]);
    Generator gen = create_PRNG(class_data);
    
    for(int t=0;t<n_trees;++t){// Iterate over n_trees
      //Rcpp::Rcout << "Tree " << t << " of class " << c << endl;
      tree(res,datas,c,c_theta,cFactor,epsilon_cont,epsilon_cat,c_spans,depth*split_nb,min_inter_sz,branch,split_nb,gen,es);
      Rcpp::checkUserInterrupt();
    }
  }
  
  // Transform to R type
  Rcpp::List ret(res.size());
  Rcpp::CharacterVector keys(res.size());
  int i = 0;
  for(auto& kv : res){
    keys[i] = kv.first;
    ret[i] = kv.second.as_List(datas);
    i++;
  }
  ret.names() = keys;
  return ret;
}