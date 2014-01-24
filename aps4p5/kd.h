//this is the 27 February 2012 rewrite

#include "containers.h"

#ifndef KD_H
#define KD_H

class kd_tree{
 
 friend class gp;
 private:
  int diagnostic;
  array_2d<int> tree;
  
  array_2d<double> data;
  array_1d<double> maxs,mins;
  
  int masterparent;

  double tol;

  int nkernel;
  
  void confirm(int,int,int,int);
  void organize(array_1d<int>&,int,int,int,int,int);
  int find_node(array_1d<double>&);
  void neigh_check(array_1d<double>&,int,array_1d<int>&,array_1d<double>&,
  int,int);
  
  void kernel_check(array_1d<double>&,array_1d<double>&,array_1d<int>&,
  int,int);
  
  void radial_check(array_1d<double>&,double,array_1d<int>&,int,int);
  
   void reassign(int);
   void descend(int);
   

  
   
  
 public:
  int xplr,ktests;
  
  
  kd_tree(array_2d<double>&);
  kd_tree(array_2d<double>&,array_1d<double>&,array_1d<double>&);
  ~kd_tree();
 
   void build_tree(array_2d<double>&);
   void build_tree(array_2d<double>&,array_1d<double>&,array_1d<double>&);
 
 void check_tree();
 void check_tree(int);
 double distance(array_1d<double>&,array_1d<double>&);
 double distance(array_1d<double>&,int);
 double distance(int,array_1d<double>&);

 
 void get_pt(int,array_1d<double>&);
 double get_pt(int,int);
 int get_pts();

 void write_tree(char*);
 void add(array_1d<double>&);
 void remove(int);
 void count(int,int*);
 void nn_srch(array_1d<double>&,int,array_1d<int>&,array_1d<double>&);
 void nn_srch(int,int,array_1d<int>&,array_1d<double>&);
 int kernel_srch(array_1d<double>&,array_1d<double>&,array_1d<int>&);
 
 int radial_srch(array_1d<double>&,double,array_1d<int>&);
 int get_dim();
 int get_diagnostic();
 
 double get_max(int);
 double get_min(int);
 
};

#endif
