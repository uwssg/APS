/*
Based on description in

Bentley, J. L.
Communications of the Association for Computint Machinery
volume 18
p 509
(1975)
*/

#include "containers.h"

#ifndef KD_H
#define KD_H

class kd_tree{
    /*
    This class stores points in an n-dimensional parameter space in a
    KD-tree for easy nearest neighbor searching.
    */
 
    /*this friend declaration makes it possible for the Gaussian Process class 
    to access the data stored in the tree*/
    friend class gp;
    
    public:
        int ktests;
        
        /*Build a tree out of the list of points provided (each row is a 
        different point in parameter space; therefore, the number of columns
        is the dimensionality of the parameter space)*/
        kd_tree(array_2d<double>&);
        
        /*Build a tree as above.  The array_1d<double> arguments are minimum
        and maximum values of each parameter in parameter space.  These
        are not bounds.  max-min is used to normalize distances in parameter
        space when searching for nearest neighbors.*/
        kd_tree(array_2d<double>&,array_1d<double>&,array_1d<double>&);
        
        ~kd_tree();
        
        /*
        These routines provide the back end for building the KD tree
        */
        void build_tree(array_2d<double>&);
        void build_tree(array_2d<double>&,array_1d<double>&,array_1d<double>&);
 
        /*
        These routines will set the maximum and minimum values (used for
        normalizing parameter distances; see above)
        */
        void set_max(int,double);
        void set_min(int,double);
 
        void check_tree();
        void check_tree(int);
        double distance(array_1d<double>&,array_1d<double>&);
        double distance(array_1d<double>&,int);
        double distance(int,array_1d<double>&);
        double distance(int,int);

 
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

    private:
        int diagnostic;
        
        /*
        The array_2d<int> tree stores the structure of the KD tree.
        Each row corresponds to a data point stored in the tree.
        There are four columns.
        
        The 0th column stores the dimension that the tree is branching on
        at that point.
        
        The 1st column is the left-hand daughter (points whose value in the
        dimension specified in the 0th column are less than the current
        point).
        
        The 2nd column is the right-hand daughter (points whose value in
        the dimension specified in the 0th column are greater than or
        equal to the current point).
        
        The 3rd column is the parent of the current point.
        
        Columns are set to -1 if they have no meaningful answer (i.e. the parent
        of the original point or the daughter of a terminal node).
        */
        array_2d<int> tree;
        
        /*
        The array_2d<double> data contains the points that are stored in this
        tree.
        */
        array_2d<double> data;
        
        /*
        maxs and mins are maximum and minimum values in each dimension of
        parameter space.  Note: these are not bounds.  The difference max-min is
        used to normalize parameter space distances when trying to find nearest
        neighbors.  Setting all maxs=1 and all mins=0 will result in nearest
        neighbors reckoned by unnormalized parameter space distances.
        */
        array_1d<double> maxs,mins;
  
        int masterparent;

        double tol;

        int nkernel;
  
        void confirm(int,int,int,int);
        void organize(array_1d<int>&,int,int,int,int,int);
        int find_node(array_1d<double>&);
        void neigh_check(array_1d<double>&,
            int,array_1d<int>&,array_1d<double>&,int,int);
  
        void kernel_check(array_1d<double>&,array_1d<double>&,array_1d<int>&,
            int,int);
  
        void radial_check(array_1d<double>&,double,array_1d<int>&,int,int);
  
        void reassign(int);
        void descend(int);

};

#endif
