/*
Based on description in

Bentley, J. L.
Communications of the Association for Computing Machinery
volume 18
p 509
(1975)
*/

#ifndef KD_H
#define KD_H

#include "time.h"
#include "containers.h"

class kd_tree{
    /*
    This class stores points in an n-dimensional parameter space in a
    KD-tree for easy nearest neighbor searching.
    */
 
    /*this friend declaration makes it possible for the Gaussian Process class 
    to access the data stored in the tree*/
    friend class gp;
    friend class box;
    
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
        
        /*check to make sure the tree is properly constructed;
        if it is, set the global variable diagnostic=1;
        if not, set diagnostic=0*/
        void check_tree();
        
        /*check to make sure the part of the tree that descendes from
        the node specified by the int is properly constructed*/
        void check_tree(int);
        
        /*
        The distance routines all return parameter space distances between
        points.  These points do not have to be on the tree, but they can be.
        
        The distances are Euclidean, except that they are normalized by the
        values stored in the private global variables maxs and mins, i.e.
        
        distance = sqrt{\sum_i [(p1_i - p2_i)/(maxs_i - mins_i)]^2}
        */
        
        /*the parameter space distance between arbitrary points*/
        double distance(array_1d<double>&,array_1d<double>&);
 
        /*the parameter space distance between and arbitrary point and a node
        on the tree*/
        double distance(array_1d<double>&,int);
        double distance(int,array_1d<double>&);
        
        /*the parameter space distance between two nodes on the tree*/
        double distance(int,int);

        /*fill the array-1d<double> with the node specified by the int*/
        void get_pt(int,array_1d<double>&);
        
        /*return a node on the tree one dimension at a time; the first
        int specifies the node; the second int specifies the dimension*/
        double get_pt(int,int);
        
        /*return a pointer to a node on the tree*/
        array_1d<double>* get_pt(int);
        
        /*return the number of points stored in the tree*/
        int get_pts();

        void write_tree(char*);
        
        /*add a point to the tree*/
        void add(array_1d<double>&);
        
        /*removes a node from the tree and then rebuilds the tree*/
        void remove(int);
        
        /*counts the number of nodes descended from a given parent*/
        void count(int,int*);
        
        /*
        The nn_srch routins perform the nearest neighbor searches.
        The first argument specifies the point whose nearest neighbors
        you want to find (if it is just an int, you are looking for the
        nearest neighbors of a node onthe tree)
        
        The second argument specifies the number of nearest neighbors to find.
        
        The third argument will store the indices of the nearest neighbors.
        
        The fourth argument will store the (normalized) parameter space
        distances to those nearest neighbors.
        
        Remember: all distances are normalized by maxs-mins (see documentation
        of the distance() routines)
        */
        void nn_srch(array_1d<double>&,int,array_1d<int>&,array_1d<double>&);
        void nn_srch(int,int,array_1d<int>&,array_1d<double>&);
        
        /*
        kernel_srch will do a neighbor search centered on a point (the first
        argument).  It will return all of the points located within a hyperbox 
        whose half-width is specified by the second argument.  The third
        argument stores the indices of the found neighbors.  The routine returns
        the number of found neighbors.
        
        NOTE: THIS IS NOT WELL-TESTED
        */
        int kernel_srch(array_1d<double>&,array_1d<double>&,array_1d<int>&);
        
        /*
        radial_srch performs a neighbor search centered on a point (the first
        argument).  It returns all points within a given (normalized) radius
        (the second argument).  The third argument stores the indices of the
        found neighbors.  The routine returns the number of found neighbors.
        
        NOTE: THIS IS NOT WELL-TESTED
        */
        int radial_srch(array_1d<double>&,double,array_1d<int>&);
        
        /*return the dimensionality of the parameter space*/
        int get_dim();
        
        /*return diagnostic, the global variable that logs whether or not
        the tree was properly constructed*/
        int get_diagnostic();
        
        /*return the maximum and minimum values in a dimension of parameter
        space*/
        double get_max(int);
        double get_min(int);
        
        double get_search_time();
        int get_search_ct();
        
        int get_search_ct_solo();
        double get_search_time_solo();
        
        void set_search_ct(int);
        void set_search_time(double);
        
        void set_search_ct_solo(int);
        void set_search_time_solo(double);
        
    private:
        
        double search_time,search_time_solo;
        int search_ct,search_ct_solo;
        
        /*a global variable logging whether or not the tree was properly
        constructed; diagnostic=1 if the tree is correct; diagnostic=0 if not*/
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
        
        /*masterparent is the point that is the first node of the tree*/
        int masterparent;
        
        /*a global variable to define the tolerance with which dimensions
        are sorted on the tree*/
        double tol;
        
        /*a global variable used by kernel_srch to keep track of how many
        points were found within the kernel
        
        also used by radial_srch
        */
        int nkernel;
        
        /*this provides the backend for check_tree;
        see source code in kd.cpp for more details*/
        void confirm(int,int,int,int);
        
        /*the iterative backend of build_tree*/
        void organize(array_1d<int>&,int,int,int,int,int);
        
        /*find the node where a new point would be added to the tree;
        this is part of the backend for the add() routine*/
        int find_node(array_1d<double>&);
        
        /*neigh_check provides the back end for nn_srch*/
        void neigh_check(array_1d<double>&,
            int,array_1d<int>&,array_1d<double>&,int,int);
        
        /*kernel_check provides the backend for kernel_srch*/
        void kernel_check(array_1d<double>&,array_1d<double>&,array_1d<int>&,
            int,int);
        
        /*radial_check provides the backend for radial_srch*/
        void radial_check(array_1d<double>&,double,array_1d<int>&,int,int);
        
        /*
        reassign() and descend() provide some of the backend for remove()
        */
        void reassign(int);
        void descend(int);

};

#endif
