#include <Rcpp.h>
#include "interaction.h"
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

// [[Rcpp::plugins(cpp11)]]

using namespace std;

vector<double> const& InteractionItem::get_values() const{
  return values;
}

void InteractionItem::intersect(double data){
  if(factor){
    if(find(values.begin(),values.end(),data) == values.end())
      values.push_back(data);
    if(values.size() > max_span)
      remove = true;
  }
  else{
    if(data<values[0])
      values[0] = data;
    
    if(data>values[1])
      values[1] = data;
      
    if(values[1] - values[0] > max_span)
      remove = true;
  }
}

bool InteractionItem::match(double data) const{
  if(factor)
    return find(values.begin(),values.end(),data) != values.end();
  else
    return data >= values[0] && data <= values[1];
}

string InteractionItem::as_string() const{
  stringstream ss;
  ss << attr_idx + 1 << "=";
  if(factor){
    ss << "{";
    for(double e : values)
      ss << e;
    ss << "}";
  }
  else
    ss << "[" << values[0] << "," << values[1] << "]";
    
  return ss.str();
}

int InteractionItem::get_idx() const{
    return attr_idx;
}

bool InteractionItem::to_remove() const{
    return remove;
}

InteractionItem::InteractionItem() : remove(true) {
  
}

InteractionItem::InteractionItem(double p_value,int p_idx,double epsilon,bool p_factor,double p_span) : 
values(1,p_value), attr_idx(p_idx), factor(p_factor), remove(false) {

    if(!factor) // form interval
      values.push_back(p_value);
    
    if(factor)
      max_span = epsilon*(p_span-1) + 1;
    else
      max_span = epsilon*p_span;
}

int Interaction::size() const{
    return items.size();
}

string Interaction::as_string() const{
  stringstream ss;
  
  for(int i=0;i<size();++i){
    ss << items[i].as_string();
    if(i < size()-1) // Not last
      ss << " ";
  }
  
  return ss.str();
}

Rcpp::List Interaction::as_List(Rcpp::List const& datas){
  Rcpp::List inter(size());
  vector<string> info({"attr_idx","values"});
  
  // Interaction
  for(int i=0;i<size();++i){
    Rcpp::List elem(2);
    elem.names() = info;
    elem[0] = items[i].get_idx();
    
    vector<double> const& vec_values = items[i].get_values();
    Rcpp::NumericVector r_values(vec_values.size());
    for(int j = 0;j<vec_values.size();++j)
      r_values[j] = vec_values[j];
      
    elem[1] = r_values;
    
    inter[i] = elem;
  }
  
  // If sims not valid, compute them
  if(!prevs_valid)
    compute_prevs(datas);
    
  // Add prevs
  Rcpp::NumericVector s(prevs.size());
  s.names() = datas.names();
  for(int c=0;c<prevs.size();++c)
    s[c] = prevs[c];
    
  // Result
  vector<string>info2({"interaction","prevalence"});
  Rcpp::List res(2);
  res[0] = inter;
  res[1] = s;
  res.names() = info2;
  
  return res;
}

bool Interaction::check_for_map(Rcpp::List const& datas,vector<double> const& theta){
  bool low = check_prevs(datas,theta,true); // at least one low prevalence
  if(!low) // Optimization: all high
    return false;
    
  bool high = false; // at least one high prevalence
  
  for(int c=0;c<prevs.size();++c){
    if(prevs[c] > theta[c]){
      high = true;
    }
  }
  return low && high;
}

bool Interaction::check_prevs(Rcpp::List const& datas,vector<double> const& theta,bool es){
  if(!es)
    return true;
    
  if(!prevs_valid)
    compute_prevs(datas);
    
  for(int c=0;c<prevs.size();++c){
    if(prevs[c] < theta[c])
      return true;
  }
  return false;
}

void Interaction::compute_prevs(Rcpp::List const& datas){
  for(int c=0;c<prevs.size();++c){ // Iterate over classes
    
    Rcpp::DataFrame class_data = Rcpp::as<Rcpp::DataFrame>(datas[c]);
    int nrows = Rcpp::as<Rcpp::NumericVector>(class_data[0]).size();
    vector<double> diffs(nrows);
    int count = 0;
    
    for(int i=0;i<nrows;++i){ // Interate over instances in class c
      bool matching = true;
      for(int j=0;j<size() && matching;++j){ // Iterate over Interaction items
        int attr_idx = items[j].get_idx();
        double instance_val = Rcpp::as<Rcpp::NumericVector>(class_data[attr_idx])[i];
        
        matching = items[j].match(instance_val);
      }
      
      if(matching)
        ++count;
    }
    
    prevs[c] = double(count)/double(nrows);
  }
  
  prevs_valid = true;
}

void Interaction::intersect(vector<double> const& instance){
    // Update items
    for(int i=0;i<items.size();++i){
        int idx = items[i].get_idx();
        items[i].intersect(instance[idx]);
    }

    // Remove items
    vector<InteractionItem> new_items(items.size());
    auto it = copy_if(items.begin(),items.end(),new_items.begin(),[](InteractionItem i){return !i.to_remove();});
    new_items.resize(distance(new_items.begin(),it));
    
    prevs_valid = false;
      
    items = new_items;
}

void Interaction::intersect(int row,Rcpp::DataFrame const& data){
    // Update items
    for(int i=0;i<items.size();++i){
        int idx = items[i].get_idx();
        double val = Rcpp::as<Rcpp::NumericVector>(data[idx])[row];
        items[i].intersect(val);
    }

    // Remove items
    vector<InteractionItem> new_items(items.size());
    auto it = copy_if(items.begin(),items.end(),new_items.begin(),[](InteractionItem i){return !i.to_remove();});
    new_items.resize(distance(new_items.begin(),it));
    
    prevs_valid = false;
      
    items = new_items;
}

int Interaction::get_depth() const{
  return depth;
}

void Interaction::set_depth(int i){
  depth = i;
}

Interaction::Interaction(){

}

Interaction::Interaction(vector<double> const& instance,vector<bool> const& factor,
vector<double> const& max_spans,double epsilon_cont,double epsilon_cat,int p_depth,int nb_class):
items(instance.size()), depth(p_depth), prevs_valid(false), prevs(nb_class){

    for(int i=0;i<instance.size();++i){
        if(factor[i]){
            InteractionItem item(instance[i],i,epsilon_cat,factor[i],max_spans[i]);
            items[i] = item;
        }
        else{
            InteractionItem item(instance[i],i,epsilon_cont,factor[i],max_spans[i]);
            items[i] = item;
        }
    }
}
