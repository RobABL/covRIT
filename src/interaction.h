#ifndef INTERACT_H
#define INTERACT_H

#include <vector>
#include <string>
#include <Rcpp.h>

using namespace std;

class InteractionItem {
    private:
        vector<double> values;
        int attr_idx;
        double max_span;
        bool factor;
        bool remove;

    public:
        void intersect(double);
        bool match(double) const;
        int get_idx() const;
        vector<double> const& get_values() const;
        bool to_remove() const;
        string as_string() const;
        InteractionItem();
        InteractionItem(double,int,double,bool,double);
};

class Interaction{
    private:
        vector<InteractionItem> items;
        int depth;
        vector<double> prevs;
        bool prevs_valid;
        void compute_prevs(Rcpp::List const&);

    public:
        int size() const;
        void intersect(vector<double> const&);
        void intersect(int,Rcpp::DataFrame const&);
        bool check_for_map(Rcpp::List const&,vector<double> const&);
        bool check_prevs(Rcpp::List const&,vector<double> const&,bool);
        string as_string() const;
        Rcpp::List as_List(Rcpp::List const&);
        int get_depth() const;
        void set_depth(int);
        Interaction();
        Interaction(vector<double> const&,vector<bool> const&,vector<double> const&,double,double,int,int);
};

#endif
