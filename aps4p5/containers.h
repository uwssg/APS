#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define letters 500

int compare_char(char*,char*);

template <typename T>
class array_1d{

public:

    array_1d();
    ~array_1d();
    
    T* get_ptr();
    
    T get_data(int) const;
    void add(T);
    void set(int,T);
    void add_val(int,T);
    void subtract_val(int,T);
    void divide_val(int,T);
    void multiply_val(int,T);
    
    void zero();
    void remove(int);
    
    void set_dim(int);
    void decrement_dim();
    void increment_dim();
    int get_dim() const;
    void set_name(char*);
    void set_where(char*) const;
    void print_name();
    void die(int) const;
    
    void reset();
    
    void assert_name(char*);
    void assert_where(char*);
    
    void assert_name_null();
    void assert_where_null();
    
    double normalize();
    double get_norm();
    double get_square_norm();
    
    void add_room(int);
    int get_room();

private:
    
    T *data;
    int dim,room,name_set;
    mutable int where_set;
    char *name;
    mutable char *where_am_i;
    
    

};



template <typename T>
class array_2d{

public:
    
    array_2d(int,int);
    array_2d();
    ~array_2d();
    
    void set_dim(int,int);
    void set_cols(int);
    
    T get_data(int,int) const;
    void set_name(char*);
    void set_where(char*) const;
    void print_name();
    void add_row(array_1d<T>&);
    void set_row(int,array_1d<T>&);
    void set(int,int,T);
    
    void zero();
    
    void add_val(int,int,T);
    void subtract_val(int,int,T);
    void multiply_val(int,int,T);
    void divide_val(int,int,T);
    
    void reset();
    void decrement_rows();
    
    int get_rows() const;
    int get_cols() const;
    void die(int,int) const;
    
    void remove_row(int);
    
    array_1d<T>* operator()(int);
    
private:

   int rows,cols,row_room;
   
   array_1d<T> *data;
   
   char *name;
   mutable char *where_am_i;
   
};

template <typename T>
class asymm_array_2d{

public:
    asymm_array_2d();
    ~asymm_array_2d();
    
    void set_name(char*);
    void set_where(char*) const;
    void add_row(const array_1d<T>&);
    void set_row(int, const array_1d<T>&);
    
    void remove_row(int);
    void zero();
    
    void set(int,int,T);
    T get_data(int,int);
    
    void add(int,T);
    
    void add_val(int,int,T);
    void subtract_val(int,int,T);
    void divide_val(int,int,T);
    void multiply_val(int,int,T);
    
    void replace_row(int,array_1d<T>&);
    
    int get_rows();
    int get_cols(int);
    
    void die(int) const;
    
    void reset();
    
    array_1d<T>* operator()(int);
    
private:
    int rows,row_room;
    array_1d<T> *data;
    char *name;
    mutable char *where_am_i;

};


template <typename T>
void merge_sort(const array_1d<T>&,array_1d<int>&,int,int);

template <typename T>
double sort_and_check(const array_1d<T>&, array_1d<T>&, array_1d<int>&);

template <typename T>
int get_dex(const array_1d<T>&, T);

#endif
