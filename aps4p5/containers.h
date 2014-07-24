/*
This file defines several classes that are used for storing vectors and
matrices.  Generally speaking, values in these functions are set using
something like

myArray.set(i,x)

to set the ith element of myArray.  They are returned using

myArray.get_data(i)

to return the ith element of myArray.  If you ask one of these classes
to return a value that doesn't exist, i.e.

myArray.set_dim(6) //setting myArray to a 0-indexed 6-element array
myArray.get_data(6) //try to return the element indexed as 6

the code will throw an exception.

Each array carries with it a member variable 'name' that is an array of char's
and can be set using

myArray.set_name("myArray")

If you choose to set 'name', then, if that array causes an exception, it will
print the name of the array to the screen, letting you know which array
was illegally manipulated.

The arrays also have a variable 'where_am_i' that can be set

myArray.set_where('inMyFunction')

This exists so that you can tell the array where it is in your code at any
given time.  'where_am_i' will also be printed to the screen in the case of an
exception.

myArray.set_where("nowhere")

sets 'where_am_i' to NULL.  Note that once 'where_am_i' is set, it persists until
it is set again, so if the array is passed between nested subroutines, it will
not necessarily know when it leaves the inner subroutine and reverts to the
outer subroutine, unless you explicitly tell it.

*/

#ifndef CONTAINERS_H
#define CONTAINERS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define letters 500

int compare_char(char*,char*);

template <typename T>
class array_1d{

    /*
    This is a class of 1-dimensional arrays with elements of type T.
    
    Presently, the code will compile with T set to double or int
    To change that, change the 'template class' statements in
    containers.cpp
    
    */

public:

    array_1d();
    ~array_1d();
    
    /*return a pointer to the data in this array*/
    T* get_ptr();
    
    /*return the element of the array specified by int*/
    T get_data(int) const;
    
    /*add an element onto the end of the array*/
    void add(T);
    
    /*set the element of the array specified by the int index to the value T.
    If you try to set an index that is beyond the current size of the array,
    zeros will be added to fill in between the current size of the array
    and the element you are setting*/
    void set(int,T);
    
    /*add the value T to the element of the array indexed by int, i.e.
    array[int] = array[int] + T*/
    void add_val(int,T);
    
    /*subtract the value T from the element of the array indexed by int*/
    void subtract_val(int,T);
    
    /*divide the element of the array indexed by int by the value T, i.e.
    array[iint] = array[int]/T */
    void divide_val(int,T);
    
    /*multiply the element of the array indexed by int by the value T*/
    void multiply_val(int,T);
    
    /*set all of the elements of the array to zero*/
    void zero();
    
    /*remove the element indexed by int from the array; shift all of the elements
    with indexes greater than int down to fill in the gap*/
    void remove(int);
    
    /*set the length of the array to int; if the array is already carrying data
    this will not delete that data*/
    void set_dim(int);
    
    /*reduce the length of the array by one without changing the contents of the
    array*/
    void decrement_dim();
    
    /*increase the length of the array by one by adding a zero to the end of the 
    array*/
    void increment_dim();
    
    /*return the length of the array*/
    int get_dim() const;
    
    /*set the name of the array, so that, if it causes an exception to be thrown,
    you will know which array threw the exception*/
    void set_name(char*);
    
    /*set where_am_i so that, if an exception is thrown, you can tell where in the program
    it happened; note that if an array is passed between multiple nested subroutines, 
    this may not be the most useful message, as where_am_i does not automatically revert
    to previous settings upon leaving a subroutine*/
    void set_where(char*) const;
    
    /*print the name of this array to the screen*/
    void print_name();
    
    /*throw an exception; the argument is the element of the array that was asked for when
    the exception is thrown*/
    void die(int) const;
    
    /*clear the contents of the array; name and where_am_i are untouched*/
    void reset();
    
    /*these routines exist so that array_2d can apply its name to the array_1d's that
    comprise its rows*/
    void assert_name(char*);
    void assert_where(char*);
    
    void assert_name_null();
    void assert_where_null();
    
    /*calculate the Euclidean norm of the array and divide all of the elements thereby.
    Return the calculated norm*/
    double normalize();
    
    /*return the Euclidean norm of the array without normalizing the array*/
    double get_norm();
    
    /*return the square of the Euclidean norm of the array without normalizing
    the array*/
    double get_square_norm();
    
    /*add room for int new elements in the array*/
    void add_room(int);
    
    /*return the amount of room allotted for elements in the array.
    This is not the same as get_dim() which returns the number of occupied
    spaces in the array*/
    int get_room();

private:
    
    /*this is where the data is actually stored for the array*/
    T *data;
    
    /*
    dim is the length of the array
    
    room is the length of data, which is >= dim (so that we don't have reallocate dim
    every time we add an element to the array)
    
    name_set is a flag that keeps track of whether name was set for this array as itself,
    or as part of an array_2d
    */
    int dim,room,name_set;
    
    /*where_set is like name_set for where_am_i*/
    mutable int where_set;
    
    /*the name of this array (for diagnostic purposes)*/
    char *name;
    
    /*where this array is in the program (for diagnostic purposes)*/
    mutable char *where_am_i;
    
    

};



template <typename T>
class array_2d{

    /*
    This is a class for 2-dimensional matrices in which each row has
    the same number of columns.
    
    Each row will be an instantiation of array_1d above.
    
    As with array_1d, the code is presently compiled so that one can
    have an array_2d of either ints or doubles
    */
    
public:
    
    /*
    Constructors.  The first constructor takes as arguments the number of
    rows and the number of columns (in that order).  You can always add rows
    to an array_2d.  Trying to add columns without first resetting the matrix
    will throw an exception.
    
    If you use the constructor that takes no arguments, you must specify the
    number of columns, either by using set_cols(int), set_dim(int,int) or by first adding
    a row to your blank array_2d with add_row(array_1d<T>&).
    */
    array_2d(int,int);
    array_2d();
    
    ~array_2d();
    
    /*set the dimensions of the array_2d. rows first, columns second*/
    void set_dim(int,int);
    
    /*set the number of columns of the array_2d. Once this is set, it cannot
    be changed without first calling reset() and deleting the contents of the
    array_2d*/
    void set_cols(int);
    
    /*return the element of the array_2d indexed by the two arguments (rows first, 
    columns second)*/
    T get_data(int,int) const;
    
    /*set the name member variable (for diagnostic purposes)*/
    void set_name(char*);
    
    /*set where_am_i (for diagnostic purposes)*/
    void set_where(char*) const;
    
    /*print name to screen*/
    void print_name();
    
    /*add the array_1d as a row to this array_2d.  If this array_2d is blank,
    then the number of columns in this array_2d will be set to the length
    of the input array_1d*/
    void add_row(array_1d<T>&);
    
    /*set the row indexed by int to the array_1d provided*/
    void set_row(int,array_1d<T>&);
    
    /*set the element indexed by the two ints to the value provided.
    If you try to set a row that is beyond the current size of this
    array_2d, the code will add rows of just zero onto the end of
    the array_2d until there is room.  If you try to set a column
    that is beyond the current size of this array_2d, the code
    will throw an exception.*/
    void set(int,int,T);
    
    /*set all of the elements of this array_2d to zero*/
    void zero();
    
    /*
    add the provided value to the indexed element, i.e.
    array[int1][int2] = array[int1][int2]+T
    */
    void add_val(int,int,T);
    
    /*subtract the provided value from the indexed element*/
    void subtract_val(int,int,T);
    
    /*multiply the indexed element by the provided value*/
    void multiply_val(int,int,T);
    
    /*divide the indexed element by the provided value*/
    void divide_val(int,int,T);
    
    /*reset the contents of the array_2d; name and where_am_i are untouched*/
    void reset();
    
    /*reduce the number of rows by one; the contents of the array_2d
    are untouched*/
    void decrement_rows();
    
    /*return the number of rows*/
    int get_rows() const;
    
    /*return the number of columns*/
    int get_cols() const;
    
    /*throw an exception; the arguments are for indicating which element
    the code tried to access when the exception was thrown*/
    void die(int,int) const;
    
    /*remove the row indexed by the int.  All of the rows with
    indexes greater than the argument are shifted to fill in the gap*/
    void remove_row(int);
    
    /*return a pointer to the indexed array, i.e.
    
    myArray2d(i) is a pointer to the ith row of myArray2d, so that
    
    *myArray2d(i) behaves just like an array_1d
    */
    array_1d<T>* operator()(int);
    
private:
   
   /*
   rows is the number of rows
   
   cols is the number of columns
   
   row_room is the number of array_1d's that have been allotted to store
   rows (so that the code knows how much room it already has when you call
   add_row)
   */
   int rows,cols,row_room;
   
   /*
   This will be allocated an array of array_1d's to store the rows of this
   array_2d
   */
   array_1d<T> *data;
   
   /*the name of this array_2d for purposes of diagnostics*/
   char *name;
   
   /*the current location of this array_2d for purposes of diagnostics*/
   mutable char *where_am_i;

    /*
    name and where_am_i will be passed down to the array_1d's making up the
    rows of this array_2d
    */
   
};

template <typename T>
class asymm_array_2d{

    /*
    This is a class for a 2-dimensional 'matrix' in which each row has a 
    different number of columns.
    
    The rows are stored as array_1d's.
    
    As before, the code will currently compile so that you can have 
    asymm_array_2d's of either ints or doubles.
    */

public:
    asymm_array_2d();
    ~asymm_array_2d();
    
    /*set the name of this asymm_array_2d for diagnostic purposes*/
    void set_name(char*);
    
    /*set the location of this asymm_array_2d for diagnostic purposes*/
    void set_where(char*) const;
    
    /*add a row to the end of this asymm_array_2d*/
    void add_row(const array_1d<T>&);
    
    /*set the indexed row to the provided array_1d;
    If you set a row beyond the current size of this asymm_array_2d,
    then empty rows will be used to fill in the gaps
    */
    void set_row(int, const array_1d<T>&);
    
    /*
    remove the row indexed by int.  Rows with indexes greater than the provided
    index will be shifted down to fill in the gap
    */
    void remove_row(int);
    
    /*set all of the elements in this asymm_array_2d to zero*/
    void zero();
    
    /*
    set the element indexed by the ints to the provided value.
    
    If you specify an element beyond the present size of the asymm_array_2d,
    zeros will be used to fill in the empty space.
    
    Note that there is no longer any restriction on the specified column
    index as each row in asymm_array_2d is allowed to have a different
    number of columns.
    */
    void set(int,int,T);
    
    /*
    Return the indexed element
    */
    T get_data(int,int);
    
    /*add T to the end of the row specified by int*/
    void add(int,T);
    
    /*add the value T onto the indexed element, i.e.
    asymm_array[int1][int2] = asymm_array[int1][int2] + T
    */
    void add_val(int,int,T);
    
    /*subtract the value T from the indexed element*/
    void subtract_val(int,int,T);
    
    /*divide the indexed element by the value T*/
    void divide_val(int,int,T);
    
    /*multiply the indexed elmement by the value T*/
    void multiply_val(int,int,T);
    
    /*replace the indexed row with the provided array_1d*/
    void replace_row(int,array_1d<T>&);
    
    /*return the number of rows*/
    int get_rows();
    
    /*return the number of columns in the row indexed by int*/
    int get_cols(int);
    
    /*throw an exception; the argument indicates the row index being
    called for when the exception was thrown*/
    void die(int) const;
    
    /*reset the contents of this asymm_array_2d*/
    void reset();
    
    /*return a pointer to the row indexed by int, i.e.
    
    myAsymmArray(i) is a pointer to the array_1d in the ith row
    
    so
    
    *myAsymmArray(i) behaves like an array_1d
    
    */
    array_1d<T>* operator()(int);
    
private:
    
    /*
    rows is the number of rows
    
    row_room is the number of array_1d's allocated for rows (so
    the code knows how much room it has when you call add_row)
    */
    int rows,row_room;
    
    /*this will be allocated as an array of array_1d's for storing the
    rows of this asymm_array_2d*/
    array_1d<T> *data;
    
    /*the name of this asymm_array_2d for diagnostic purposes*/
    char *name;
    
    /*the location of this asymm_array_2d in the code, for diagnostic purposes*/
    mutable char *where_am_i;

    /*name and where_am_i will be passed down to the array_1d's making up the rows
    of this assymm_array_2d*/
};


///////////////*BELOW ARE ROUTINES FOR SORTING array_1d's/////////////

/*a driver for a merge sort as described here

http://en.wikipedia.org/wiki/Merge_sort

The user should not call this routine; it is called by sort_and_check, which
the user should call
*/
template <typename T>
void merge_sort(const array_1d<T>&,array_1d<int>&,int,int);



/*
The first array_1d<T> is the input data

The second array_1d<T> will store the sorted output

The array_1d<int> are unique id's associated with each element of the input data.
It will be rearranged to preserve the relationship between these indexes and the
sorted data.

This routine will check to make sure that the output is sorted from lowest
to highest and that the unique id's are preserved relative to the input data.

The routine returns the maximum error fabs(input(unique_id=i)-output(unique_id=i)).

If this maximum error is greater than 10^-12, then the code throws an exception
*/
template <typename T>
double sort_and_check(const array_1d<T>&, array_1d<T>&, array_1d<int>&);



/*return the index of the element of the array_1d that is closest in value to T;

ASSUMES THAT THE array_1d IS SORTED FROM LOWEST TO HIGHEST
*/
template <typename T>
int get_dex(const array_1d<T>&, T);

#endif
