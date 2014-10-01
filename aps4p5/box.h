#ifndef BOX_H
#define BOX_H

#include "containers.h"
#include "goto_tools.h"
#include <time.h>

#define box_exception -1
#define min_pts_per_box 10

class box{

private:
    
    int diagnostic,pts_per_box;
    array_2d<int> tree;
  
    asymm_array_2d<int> box_contents;
    
    array_2d<double> *data;
    
    array_1d<double> tree_values;
    array_2d<double> box_min,box_max;
    array_1d<double> maxs,mins,norm_max,norm_min;
    
    array_1d<int> tree_ct;
    
    void build_tree();
    void organize(array_1d<int>&,int,int,int,
        int,array_1d<double>&,array_1d<double>&);
    int split_box(int,int,int);
    
    double time_add_srch,time_split;
    
    double split_error(array_1d<int>&,int,int*,int*);
    
    double time_search;
    int ct_search;
    
public:
    
    box(array_2d<double>*,int);
    box(array_2d<double>*,int,array_1d<double>&,array_1d<double>&);
    
    void initialize(array_2d<double>*,int,array_1d<double>&,array_1d<double>&);

    ~box();
    
    void set_pts_per(int);
    
    int get_ct_search();
    double get_time_search();
    
    double distance(array_1d<double>&,array_1d<double>&);
    double distance(int,array_1d<double>&);
    double distance(array_1d<double>&,int);
    
    void verify_tree();
    
    int find_box(array_1d<double>&);
    int find_box(array_1d<double>&,int*,int*);
    
    int add_pt();
    int get_nboxes();
    int get_n_small_boxes();
    int get_n_optimal_boxes();
    int get_smallest_box();
    int get_biggest_box();
    double get_mean_box(double*);
    int get_ntree();
    int get_contents(int);
    int get_contents(int,int);
    
    int get_pts();
    int get_dim();
    
    double get_box_max(int,int);
    double get_box_min(int,int);
    
    double get_max(int) const;
    double get_min(int) const;
    
    double get_time_add_srch();
    double get_time_split();
    
    array_1d<double>* get_pt(int);
    double get_pt(int,int) const;
    void get_pt(int,array_1d<double>&);
    
    void nn_srch(array_1d<double>&,array_1d<int>&,array_1d<double>&);
    void nn_srch(int,array_1d<int>&,array_1d<double>&);
    
    void nn_srch(array_1d<double>&,array_1d<int>&,array_1d<double>&,array_1d<int>&);
    void nn_srch(int,array_1d<int>&,array_1d<double>&,array_1d<int>&);
    
    array_1d<int>* get_box(int);
    void refactor();
    
    void add_to_search_time(double);
    
    void get_tree_cts(array_1d<int>&);
    void get_avg_box_bounds(array_1d<double>&,array_1d<double>&,
        array_1d<double>&,array_1d<double>&);
    
    
};

#endif
