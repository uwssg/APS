#include "containers.h"

int compare_char(char *s1, char *s2){
 //are two character strings the same?
 //if so, return 1
 //if not, return 0

 int i;
 //printf("\ncomparing %s %s\n",s1,s2);
 for(i=0;i<letters && (s1[i]!=0 || s2[i]!=0);i++){
  if(s1[i]!=s2[i])return 0;
 }
 return 1;
 
}


template <typename T>
array_1d<T>::array_1d(){
    room=0;
    dim=0;
    data=NULL;
    name=NULL;
    where_am_i=NULL;
    
    name_set=0;
    where_set=0;
}

template <typename T>
array_1d<T>::~array_1d(){
    if(data!=NULL){
       delete [] data;
    }
    
    if(name!=NULL && name_set==1){
        delete [] name;
    }
    
    if(where_am_i!=NULL && where_set==1){
        delete [] where_am_i;
    }
    
}

template <typename T>
T* array_1d<T>::get_ptr(){
    return data;
}

template <typename T>
void array_1d<T>::die(int ii) const{
    printf("\nWARNING 1d array\n");
    
    if(name!=NULL)printf("in 1_d array %s\n",name);
    if(where_am_i!=NULL)printf("in routine %s\n",where_am_i);
    
    printf("asked for %d but dim %d\n",ii,dim);
    printf("room %d\n",room);
    
    if(data==NULL){
        printf("data is null\n\n");
    }
    
    int ifail=1;
    
    throw ifail;
    
}

template <typename T>
void array_1d<T>::add_val(int dex, T val){
    
    if(dex<0 || dex>=dim){
        printf("dying from add_val\n");
        die(dex);
    }
    
    data[dex]+=val;
    
    
}

template <typename T>
void array_1d<T>::subtract_val(int dex, T val){
 
    
    if(dex<0 || dex>=dim){
        printf("dying from subtract_val");
        die(dex);
    }
    
    data[dex]-=val;
}

template <typename T>
void array_1d<T>::divide_val(int dex, T val){
    
    if(dex<0 || dex>=dim){
        printf("dying from divide_val\n");
        die(dex);
    }
    
    data[dex]=data[dex]/val;
    
    
}

template <typename T>
void array_1d<T>::multiply_val(int dex, T val){

    if(dex<0 || dex>=dim){
        printf("dying from multiply_val\n");
        die(dex);
    }
    
    data[dex]*=val;
    
    
}

template <typename T>
int array_1d<T>::get_room(){
    return room;
}

template <typename T>
void array_1d<T>::add_room(int ii){

    if(data==NULL && dim>0){
       printf("dying from add_room\n");
       die(0);
    }
    
    if(data==NULL && room>0){
        printf("dying from add_room\n");
        die(0);
    }
    
    if(room==0 && data!=NULL){
        printf("dying from add_room\n");
         die(0);
    }
    
    if(dim>room){
        printf("dying from add_room\n");
        die(0);
    }
    
    T *buffer;
    int i,old_room=room;
    
    if(data==NULL){
        room=ii;
        data=new T[room];
        for(i=0;i<room;i++)data[i]=0;
    }
    else{
        buffer=new T[room];
        for(i=0;i<room;i++)buffer[i]=data[i];
        delete [] data;
        room+=ii;
        data=new T[room];
        for(i=0;i<old_room;i++)data[i]=buffer[i];
        delete [] buffer;
    }
    
    
    
}

template <typename T>
void array_1d<T>::add(T in){

    if(data==NULL && dim>0){
       printf("dying from add\n");
       die(0);
    }
    
    if(data==NULL && room>0){
        printf("dying from add\n");
        die(0);
    }
    
    if(room==0 && data!=NULL){
        printf("dying from add\n");
         die(0);
    }


    if(data==NULL){
        room=2;
	data=new T[room];
	dim=0;
    }
    
    T *buffer;
    int i;
    
    if(dim==room){
        buffer=new T[dim];
	for(i=0;i<dim;i++)buffer[i]=data[i];
	delete [] data;
	room+=5;
	data=new T[room];
	for(i=0;i<dim;i++){
	    data[i]=buffer[i];
	}
	delete [] buffer;
    }
    
    data[dim]=in;
    dim++;
}

template <typename T>
void array_1d<T>::set(int dex, T val){
   
    
    int i;
    
    if(dex<0){
        printf("dying from set with negative dex\n");
        die(dex);
    }
    else if(dex>=dim){
        for(i=dim;i<dex+1;i++)add(0);
	set(dex,val);
    }
    else{
        data[dex]=val;
    } 
}

template <typename T>
void array_1d<T>::set_dim(int ii){
    
    
    if(ii<0){
        printf("tried to set dim %d\n",dim);
	die(ii);
    }
    
    if(ii==0){
        reset();
	return;
    }
    
    T *buffer;
    int j,new_room;
    
    if(data!=NULL && room<ii){
        buffer=new T[room];
	for(j=0;j<room;j++)buffer[j]=data[j];
	delete [] data;
	new_room=ii;
	data=new T[new_room];
	for(j=0;j<room;j++)data[j]=buffer[j];
	delete [] buffer;
	dim=ii;
	room=new_room;
	for(;j<room;j++)data[j]=0;
    }
    else if(data!=NULL && room>=ii){
        dim=ii;
    }
    else if(data==NULL){
        data=new T[ii];
	room=ii;
	dim=ii;
	for(j=0;j<room;j++)data[j]=0;
    }
    
    
}

template <typename T>
void array_1d<T>::zero(){
    int i;
    for(i=0;i<room;i++)data[i]=0;
}

template <typename T>
void array_1d<T>::remove(int dex){

    if(dex<0 || dex>=dim){
        return;
    }
    
    int i;
    for(i=dex+1;i<dim;i++){
        data[i-1]=data[i];
    }
    dim--;

}

template <typename T>
void array_1d<T>::decrement_dim(){
    dim--;
    if(dim<0)dim=0;
}

template <typename T>
void array_1d<T>::increment_dim(){
    int i=dim;
    add(0);
    
    if(dim!=i+1){
         printf("WARNING increment_dim did not work %d %d\n",i,dim);
	 die(dim);
    }
}

template <typename T>
void array_1d<T>::reset(){

    if(data==NULL && (room>0 || dim>0)){
        printf("dying from reset\n");
         die(-1);
    }    
    
    if(data!=NULL && room==0){
        printf("dying from reset\n");
        die(-1);
    }
    
    delete [] data;
    data=NULL;
    room=0;
    dim=0;

    
}

template <typename T>
void array_1d<T>::set_where(char *word) const{

    int i,ct=0;

    
    if(where_am_i!=NULL && where_set==1){
        delete [] where_am_i;
    }
    
    if(compare_char(word,"nowhere")==1){
        where_set=0;
	where_am_i=NULL;
	return;
    }
    
    for(i=0;word[i]!=0;i++)ct++;
    ct++;
    
    where_am_i=new char[ct];
    for(i=0;i<ct && word[i]!=0;i++)where_am_i[i]=word[i];
    where_am_i[ct-1]=0;
    
    where_set=1;

}

template <typename T>
void array_1d<T>::set_name(char *word){
    int i,ct=0;

    
    if(name!=NULL && name_set==1)delete [] name;
    
    if(compare_char(word,"nowhere")==1){
        name_set=0;
	name=NULL;
	return;
    }
    
    for(i=0;word[i]!=0;i++)ct++;
    ct++;
    
    name=new char[ct];
    for(i=0;i<ct && word[i]!=0;i++)name[i]=word[i];
    name[ct-1]=0;
    
    name_set=1;
    
}

template <typename T>
void array_1d<T>::assert_name(char *word){
    if(name_set==1){
        printf("cannot assert name; it has been set internally\n");
	die(0);
    }
    
    name=word;
}

template <typename T>
void array_1d<T>::assert_where(char *word){
    if(where_set==1){
        printf("cannot assert where; it has been set internally\n");
	die(0);
    }
    
    where_am_i=word;
}

template <typename T>
void array_1d<T>::assert_name_null(){
    if(name_set==1){
        delete [] name;
	name_set=0;
    }
    
    name=NULL;
}

template <typename T>
void array_1d<T>::assert_where_null(){
    if(where_set==1){
        delete [] where_am_i;
	where_set=0;
    }
    
    where_am_i=NULL;
}


template <typename T>
void array_1d<T>::print_name(){
    if(name!=NULL)printf("%s\n",name);
}

template <typename T>
int array_1d<T>::get_dim() const{
    return dim;
}

template <typename T>
T array_1d<T>::get_data(int dex) const{

    if(data==NULL){
        printf("dying from get_data because data is null\n");
        die(dex);
    }

    if(dex<0 || dex>=dim){
        printf("dying from get_data because request makes no sense\n");
        die(dex);
    }
    else{
        return data[dex];
    }
    
}


template <typename T>
double array_1d<T>::get_square_norm(){
    
    if(dim<0){
        printf("WARNING 1d array has dim %d\n",dim);
	die(-1);
    }
    
    if(dim==0){
        return 0.0;
    }
    
    int i;
    double ans=0.0;
    for(i=0;i<dim;i++){
        ans+=data[i]*data[i];
    }

    return ans;

}


template <typename T>
double array_1d<T>::get_norm(){
    
    if(dim<0){
        printf("WARNING 1d array has dim %d\n",dim);
	die(-1);
    }
    
    if(dim==0){
        return 0.0;
    }
    
    int i;
    double ans=0.0;
    for(i=0;i<dim;i++){
        ans+=data[i]*data[i];
    }
    ans=sqrt(ans);
    return ans;

}

template <typename T>
double array_1d<T>::normalize(){
    
    if(dim<0){
        printf("WARNING 1d array has dim %d\n",dim);
	die(-1);
    }
    
    if(dim==0){
        return 0.0;
    }
     
    double ans;
    int i;
    ans=0.0;
    for(i=0;i<dim;i++){
        ans+=data[i]*data[i];
    }
    
    if(ans<0.0){
        printf("WARNING square of norm %e\n",ans);
	
	die(-1);
    }
    
    if(ans>0.0){
        ans=sqrt(ans);
        for(i=0;i<dim;i++){
            data[i]=data[i]/ans;
        }
    }
    
    return ans;

}

template <typename T>
array_2d<T>::array_2d(){
    rows=0;
    cols=0;
    row_room=0;
    data=NULL;
    name=NULL;
    where_am_i=NULL;
}

template <typename T>
array_2d<T>::array_2d(int r, int c){
    int i,j;
    rows=0;
    cols=c;
    row_room=r;
   
    data=new array_1d<T>[row_room];

    name=NULL;
    where_am_i=NULL;
    
    
}

template <typename T>
array_2d<T>::~array_2d(){
    int i;
    
    for(i=0;i<row_room;i++){
        try{
            data[i].assert_name_null();
	}
	catch(int iex){
	    printf("in 2d destructor\n");
	    die(0,0);
	}
	
	try{
	    data[i].assert_where_null();
	}
	catch(int iex){
	    printf("in 2d destructor\n");
	    die(0,0);
	}
    }
    
    //printf("calling 2d destructor on %s\n",name);
    
    if(data!=NULL){
        delete [] data;
    }
    
    if(name!=NULL)delete [] name;
    
    if(where_am_i!=NULL)delete [] where_am_i;
    
}

template <typename T>
void array_2d<T>::die(int ir, int ic) const{
    printf("\nWARNING 2d array\n");
    
    if(name!=NULL)printf("in 2d_array %s\n",name);
    if(where_am_i!=NULL)printf("in routine %s\n",where_am_i);
    
    printf("asked for %d %d\n",ir,ic);
    printf("but dimensions are %d %d\n",rows,cols);
    printf("row_room %d\n",row_room);
    
    if(data==NULL){
        printf("data is null\n");
    }
    
    int ifail=1;
    
    throw ifail;
}

template <typename T>
void array_2d<T>::add_row(array_1d<T> &in){
    
    
    if(data==NULL && row_room>0){
        printf("dying from add_row\n");
        die(0,0);
    }
    
    if(data==NULL && rows>0){
        printf("dying from add_row\n");
        die(1,0);
    }
    
    if(data==NULL && cols>0){
        printf("dying from add_row\n");
        die(2,0);
    }
    
    if(row_room<rows){
        printf("dying from add_row\n");
        die(3,0);
    }
    
    if((cols<=0 || row_room<=0) && data!=NULL){
        printf("dying from add_row\n");
        die(4,0);
    }
    
    if(cols>0 && cols!=in.get_dim()){
        printf("dying from add_row\n");
        die(-2,0);
    }
    
    if(cols==0){
        cols=in.get_dim();
    }
    
    if(cols==0){
        printf("about add a row but cols zero\n");
	die(-1,-1);
    }
    
    int i,j;
    
    if(data==NULL){
        row_room=2;
        data=new array_1d<T>[row_room];
	for(i=0;i<row_room;i++){
	    try{
	        data[i].assert_name(name);
	    }
	    catch(int iex){
	        printf("in 2d add row\n");
	        die(0,0);
	    }
	    
	    try{
	        data[i].assert_where(where_am_i);
	    }
	    catch(int iex){
	        printf("in 2d add row\n");
	        die(0,0);
	    }
	} 
    }
    
    array_1d<T> *buffer;
    
    if(rows==row_room){
        buffer=new array_1d<T>[rows];
	for(i=0;i<rows;i++){
	    buffer[i].set_dim(data[i].get_dim());
	    for(j=0;j<data[i].get_dim();j++){
	        buffer[i].set(j,data[i].get_data(j));
	    }
	}
	delete [] data;
	
	i=row_room/2;
	if(i<100)i=100;
	
	row_room+=i;
	data=new array_1d<T>[row_room];
	for(i=0;i<rows;i++){
	    data[i].set_dim(buffer[i].get_dim());
	    for(j=0;j<buffer[i].get_dim();j++){
	        data[i].set(j,buffer[i].get_data(j));
	    }
	}
	delete [] buffer;
	
	for(i=0;i<row_room;i++){
	    try{
	        data[i].assert_name(name);
	    }
	    catch(int iex){
	        printf("in 2d add row\n");
	        die(0,0);
	    }
	    
	    try{
	        data[i].assert_where(where_am_i);
	    }
	    catch(int iex){
	        printf("in 2d add row\n");
	        die(0,0);
	    }
	}
	
    }
    
    data[rows].set_dim(cols);
    for(i=0;i<cols;i++){
        try{
	    data[rows].set(i,in.get_data(i));
	}
	catch(int iex){
	    die(rows,i);
	}
    }
    rows++;
}

template <typename T>
void array_2d<T>::zero(){
    int i;
    for(i=0;i<row_room;i++){
        data[i].zero();
    }
}

template <typename T>
void array_2d<T>::remove_row(int dex){

    if(dex<0 || dex>=rows)return;
    
    int i,j;
    for(i=dex+1;i<rows;i++){
        for(j=0;j<data[i].get_dim();j++){
	    data[i-1].set(j,data[i].get_data(j));
	}
    }
    
    data[rows-1].reset();
    rows--;

}

template <typename T>
void array_2d<T>::set_cols(int ii){
    reset();
    
    row_room=2;
    rows=0;
    cols=ii;
    data=new array_1d<T>[row_room];
    int i;
    for(i=0;i<row_room;i++){
        data[i].set_dim(cols);
    }
    
    for(i=0;i<row_room;i++){
        try{
            data[i].assert_name(name);
	}
	catch(int iex){
	    printf("in 2d set dim\n");
	    die(0,0);
	}
	
	try{
	    data[i].assert_where(where_am_i);
	}
	catch(int iex){
	    printf("in 2d set dim\n");
	    die(0,0);
	}
    }
    
}

template <typename T>
void array_2d<T>::set_dim(int ir, int ic){
    
    if(ir<0 || ic<0){
        printf("tried to set dimensions %d %d\n",ir,ic);
	die(ir,ic);
    }
    
    if(data==NULL && (rows>0 || cols>0)){
        printf("WARNING data is null but rows %d cols %d\n",
	rows,cols);
	if(name!=NULL)printf("name %s\n",name);
	if(where_am_i!=NULL)printf("where %s\n",where_am_i);
	exit(1);
    }
    
    if(data!=NULL && cols<=0){
        printf("WARNING data is not null but rows %d cols %d\n",
	rows,cols);
	if(name!=NULL)printf("name %s\n",name);
	if(where_am_i!=NULL)printf("where %s\n",where_am_i);
	exit(1);
    }
    
    if(ir==rows && ic==cols){
        return;
    }
    
    if(ir==0 && ic==0){
        reset();
	return;
    }
    
    if((ir==0 && ic!=0) || (ic==0 && ir!=0)){
        printf("WARNING trying to set dim %d %d\n",ir,ic);
	die(ir,ic);
    }
    
    int i;
    if(data!=NULL){
        delete [] data;
    }
    
    row_room=ir;
    rows=ir;
    cols=ic;
    data=new array_1d<T>[row_room];
    
    int j;
    for(i=0;i<rows;i++){
        data[i].set_dim(cols);
        for(j=0;j<cols;j++){
	    if(i!=j)data[i].set(j,0);
	    else data[i].set(j,1);
	}
    }
    
    for(i=0;i<row_room;i++){
        try{
            data[i].assert_name(name);
	}
	catch(int iex){
	    printf("in 2d set dim\n");
	    die(0,0);
	}
	
	try{
	    data[i].assert_where(where_am_i);
	}
	catch(int iex){
	    printf("in 2d set dim\n");
	    die(0,0);
	}
    }
    
}

template <typename T>
void array_2d<T>::set(int ir, int ic, T val){
    
    if(ir<0){
        printf("tried to set to negative row\n");
	die(ir,ic);
    }
    
    if(ic<0 || ic>=cols){
        printf("dying from set\n");
        die(ir,ic);
    }
    
    if(cols<=0){
        printf("\nYou cannot use set(int,int) on a 2d array if cols are zero\n");
        die(ir,ic);
    }
    
    if(data==NULL){ 
        printf("dying from set\n");
        die(ir,ic);
    }
    
    int i;
    array_1d<T> vector;
    if(ir>=rows){
        for(i=0;i<cols;i++)vector.set(i,0);
	while(rows<=ir)add_row(vector);
	
    }

    data[ir].set(ic,val);
    
    
    
}

template <typename T>
void array_2d<T>::add_val(int ir, int ic, T val){
    
    
    if(ir>=rows || ic>=cols || data==NULL || ir<0 || ic<0){
        printf("dying from add_val\n");
         die(ir,ic);
    }
    
    data[ir].add_val(ic,val);

}

template <typename T>
void array_2d<T>::subtract_val(int ir, int ic, T val){
    
    
    if(ir>=rows || ic>=cols || data==NULL || ir<0 || ic<0){
        printf("dying from subtract_val\n");
         die(ir,ic);
    }
    
    data[ir].subtract_val(ic,val);

}

template <typename T>
void array_2d<T>::multiply_val(int ir, int ic, T val){
    
    
    if(ir>=rows || ic>=cols || data==NULL || ir<0 || ic<0){
        printf("dying from multiply_val\n");
         die(ir,ic);
    }
    
    data[ir].multiply_val(ic,val);

}

template <typename T>
void array_2d<T>::divide_val(int ir, int ic, T val){
    
    
    if(ir>=rows || ic>=cols || data==NULL || ir<0 || ic<0){
        printf("dying from divide_val\n");
         die(ir,ic);
    }
    
    data[ir].divide_val(ic,val);

}

template <typename T>
void array_2d<T>::set_row(int dex, array_1d<T> &in){

    if(dex<0){
        printf("tried to set to negative row\n");
	die(dex,0);
    }
    
    if(data==NULL && row_room>0){
        printf("dying from set_row\n");
        die(0,0);
    }
    
    if(data==NULL && rows>0){
         printf("dying from set_row\n");
        die(1,0);
    }
    
    if(data==NULL && cols>0){
         printf("dying from set_row\n");
        die(2,0);
    }
    
    if(row_room<rows){   
        printf("dying from set_row\n");
        die(3,0);
    }
    
    if((cols<=0 || row_room<=0) && data!=NULL){
         printf("dying from set_row\n");
        die(4,0);
    }
    
    if(cols>0 && cols!=in.get_dim()){
         printf("dying from set_row\n");
        printf("columns do not match\n");
        die(-2,0);
    }
    
    if(dex<0){
         printf("dying from set_row\n");
        die(dex,0);
    }
    
    int i;
    if(dex>=rows){
        for(i=rows;i<dex+1;i++)add_row(in);
    }
    else{
        for(i=0;i<cols;i++){
	    try{
	        data[dex].set(i,in.get_data(i));
	    }
	    catch(int ifail){
	        die(dex,cols);
	    }
	}
    }
    
}

template <typename T>
void array_2d<T>::decrement_rows(){
    
    if(rows==0){
        printf("WARNING trying to decrement rows but rows already zero\n");
	die(0,0);
    }
    
    rows--;
}

template <typename T>
void array_2d<T>::reset(){
    
    //printf("resetting %s\n",name);
    
    //set_name("resetting");
    
    int i;
    
    if(data==NULL && (rows>0 || cols>0 || row_room>0)){
        printf("resetting but data is null and something is wrong\n");
	die(-1,-1);
    }
    
    if(row_room==0 && data!=NULL){
        die(-1,-1);
    }
    
    if(cols==0 && data!=NULL){
        die(-1,-1);
    }
    
    if(row_room<rows){
	die(-2,-2);
    }
    
    if(data!=NULL){
        delete [] data;
	
	data=NULL;
	row_room=0;
	rows=0;
	cols=0;
    }
    

    
}

template <typename T>
void array_2d<T>::print_name(){
    if(name!=NULL)printf("%s\n",name);
}

template <typename T>
void array_2d<T>::set_where(char *word) const {

    int i,ct=0;
    
    for(i=0;word[i]!=0;i++)ct++;
    ct++;
    
    if(where_am_i!=NULL)delete [] where_am_i;
    where_am_i=new char[ct];
    
    for(i=0;i<ct && word[i]!=0;i++)where_am_i[i]=word[i];
    where_am_i[ct-1]=0;
    
    for(i=0;i<row_room;i++){
        try{
            data[i].assert_where(where_am_i);
	}
	catch(int iex){
	    printf("in 2d set where\n");
	    die(0,0);
	}
    }

}

template <typename T>
void array_2d<T>::set_name(char *word){
    int i,ct=0;
    
    for(i=0;word[i]!=0;i++)ct++;
    ct++;
    
    if(name!=NULL){
        delete [] name;
    }
    
    name=new char[ct];
    for(i=0;i<ct && word[i]!=0;i++)name[i]=word[i];
    name[ct-1]=0;
    
    for(i=0;i<row_room;i++){
        try{
            data[i].assert_name(name);
	}
	catch(int iex){
	    printf("in 2d set name\n");
	    die(0,0);
	}
    }
    
}

template <typename T>
int array_2d<T>::get_rows() const{
    return rows;
}

template <typename T>
int array_2d<T>::get_cols() const{
    return cols;
}

template <typename T>
T array_2d<T>::get_data(int ir, int ic) const{
   
    
    if(data==NULL){
        printf("dying from get_data\n");
        die(ir,ic);
    }
    else if(row_room<rows){
       printf("dying from get_data\n");
        die(ir,ic);
    }
    else if(ir>=rows || ir<0 || ic>=cols || ic<0){
       printf("dying from get_data\n");
        die(ir,ic);
    }
    else{
        return data[ir].get_data(ic);
    }
}

template <typename T>
array_1d<T>* array_2d<T>::operator()(int dex){
    
    if(dex<0 || dex>=rows){
        printf("WARNING asked for row %d but only have %d\n",dex,rows);
    }
    
    return &data[dex];

}

template <typename T>
asymm_array_2d<T>::asymm_array_2d(){
    name=NULL;
    where_am_i=NULL;
    data=NULL;
    rows=0;
    row_room=0;
}

template <typename T>
asymm_array_2d<T>::~asymm_array_2d(){
    int i;
    
    for(i=0;i<row_room;i++){
        try{
            data[i].assert_name_null();
	}
	catch(int iex){
	    printf("in asymm 2d destructor\n");
	    die(0);
	}
	
	try{
	    data[i].assert_where_null();
	}
	catch(int iex){
	    printf("in asymm 2d destructor\n");
	    die(0);
	}
    }
    
    //printf("calling 2d destructor on %s\n",name);
    
    if(data!=NULL){
        delete [] data;
    }
    
    if(name!=NULL)delete [] name;
    
    if(where_am_i!=NULL)delete [] where_am_i;
    

}

template <typename T>
void asymm_array_2d<T>::die(int ir) const{
    printf("\nWARNING asymm 2d array\n");
    
    if(name!=NULL)printf("in 2d_array %s\n",name);
    if(where_am_i!=NULL)printf("in routine %s\n",where_am_i);
    
    printf("asked for %d\n",ir);
    printf("but dimensions are %d\n",rows);
    printf("row_room %d\n",row_room);
    
    if(data==NULL){
        printf("data is null\n");
    }
    
    int ifail=1;
    
    throw ifail;
}

template <typename T>
void asymm_array_2d<T>::zero(){
    int i;
    for(i=0;i<row_room;i++)data[i].zero();
}

template <typename T>
void asymm_array_2d<T>::set_where(char *word) const {

    int i,ct=0;
    
    for(i=0;word[i]!=0;i++)ct++;
    ct++;
    
    if(where_am_i!=NULL)delete [] where_am_i;
    where_am_i=new char[ct];
    
    for(i=0;i<ct && word[i]!=0;i++)where_am_i[i]=word[i];
    where_am_i[ct-1]=0;
    
    for(i=0;i<row_room;i++){
        try{
            data[i].assert_where(where_am_i);
	}
	catch(int iex){
	    printf("in asymm 2d set where\n");
	    die(0);
	}
    }

}

template <typename T>
void asymm_array_2d<T>::set_name(char *word){
    int i,ct=0;
    
    for(i=0;word[i]!=0;i++)ct++;
    ct++;
    
    if(name!=NULL){
        delete [] name;
    }
    
    name=new char[ct];
    for(i=0;i<ct && word[i]!=0;i++)name[i]=word[i];
    name[ct-1]=0;
    
    for(i=0;i<row_room;i++){
        try{
            data[i].assert_name(name);
	}
	catch(int iex){
	    printf("in asymm 2d set name\n");
	    die(0);
	}
    }
    
}

template <typename T>
int asymm_array_2d<T>::get_rows(){
    return rows;
}

template <typename T>
int asymm_array_2d<T>::get_cols(int dex){
    
    if(data==NULL){
        return 0;
        
        //printf("WARNING asking for cols in asymm array 2d\n");
	//die(dex);
    }
    
    if(dex<0){
        printf("WARNING asking for cols in asymm array 2d\n");
	die(dex);
    }
    
    if(dex>=rows){
        return 0;
    }
    
    return data[dex].get_dim();
    
}

template <typename T>
void asymm_array_2d<T>::add_row(const array_1d<T> &in){
    
    if(data==NULL){
        row_room=2;
	rows=0;
	data=new array_1d<T>[row_room];
    }
    
    array_1d<T> *buffer;
    int i,j;
    
    if(rows==row_room){
        buffer=new array_1d<T>[rows];
	for(i=0;i<rows;i++){
	    buffer[i].set_dim(data[i].get_dim());
	    for(j=0;j<data[i].get_dim();j++){
	        buffer[i].set(j,data[i].get_data(j));
	    }
	}
	delete [] data;
	
	i=row_room/2;
	if(i<100)i=100;
	
	row_room+=i;
	data=new array_1d<T>[row_room];
	
	for(i=0;i<rows;i++){
	    data[i].set_dim(buffer[i].get_dim());
	    for(j=0;j<buffer[i].get_dim();j++){
	        data[i].set(j,buffer[i].get_data(j));
	    }
	}
	delete [] buffer;
    }
    
    data[rows].set_dim(in.get_dim());
    for(i=0;i<in.get_dim();i++){
        data[rows].set(i,in.get_data(i));
    }
    rows++;
    
    for(i=0;i<row_room;i++){
	    try{
	        data[i].assert_name(name);
	    }
	    catch(int iex){
	        printf("in asymm 2d add row (asserting name)\n");
	        die(0);
	    }
	    
	    try{
	        data[i].assert_where(where_am_i);
	    }
	    catch(int iex){
	        printf("in asymm 2d add row (asserting where)\n");
	        die(0);
	    }
    } 
    
}

template <typename T>
void asymm_array_2d<T>::set(int ir, int ic, T val){
    
    array_1d<T> empty;
    int i;
 
    while(rows<=ir){
        add_row(empty);
    }
    
    
    data[ir].set(ic,val);

}

template <typename T>
void asymm_array_2d<T>::set_row(int ir, const array_1d<T> &vv){
    
    int i;
    array_1d<T> empty;
    
    while(rows<=ir){
        add_row(empty);
    }
    
    data[ir].reset();
    for(i=0;i<vv.get_dim();i++){
        data[ir].set(i,vv.get_data(i));
    }
    

}

template <typename T>
T asymm_array_2d<T>::get_data(int ir, int ic){
    
    if(ir<0 || ir>=rows){
        printf("WARNING asking for asymm 2d data %d %d but rows %d\n",
	ir,ic,rows);
	die(ir);
    }
    
    
    try{
       return data[ir].get_data(ic); 
    }
    catch(int iex){
        printf("tried to get asymm 2d data %d %d\n",ir,ic);
	die(ir);
    }


}

template <typename T>
void asymm_array_2d<T>::remove_row(int dex){
    
    if(dex<0 || dex>=rows){
        printf("WARNING asking to remove %d from asymm\n",dex);
	die(dex);
    }
    
    int i,j;
    for(i=dex;i<rows-1;i++){
        data[i].set_dim(data[i+1].get_dim());
	for(j=0;j<data[i+1].get_dim();j++){
	    data[i].set(j,data[i+1].get_data(j));
	}
    }
    data[rows-1].reset();
    rows--;
    
    
}

template <typename T>
void asymm_array_2d<T>::reset(){
    
    //printf("resetting %s\n",name);
    
    //set_name("resetting");
    
    int i;
    
    if(data==NULL && (rows>0 || row_room>0)){
        printf("resetting but data is null and something is wrong\n");
	die(-1);
    }
    
    if(row_room==0 && data!=NULL){
        die(-1);
    }

    if(row_room<rows){
	die(-2);
    }
    
    if(data!=NULL){
        delete [] data;	
	data=NULL;
	row_room=0;
	rows=0;

    }
    

    
}

template <typename T>
void asymm_array_2d<T>::add_val(int ir, int ic, T val){
    if(ir<0 || ir>=rows){
        printf("in asymm 2d add_val\n");
	die(ir);
    }
    
    data[ir].add_val(ic,val);
}

template <typename T>
void asymm_array_2d<T>::subtract_val(int ir, int ic, T val){
    if(ir<0 || ir>=rows){
        printf("in asymm 2d subtract_val\n");
	die(ir);
    }
    
    data[ir].subtract_val(ic,val);
}

template <typename T>
void asymm_array_2d<T>::multiply_val(int ir, int ic, T val){
    if(ir<0 || ir>=rows){
        printf("in asymm 2d multiply_val\n");
	die(ir);
    }
    
    data[ir].multiply_val(ic,val);
}

template <typename T>
void asymm_array_2d<T>::divide_val(int ir, int ic, T val){
    if(ir<0 || ir>=rows){
        printf("in asymm 2d divide_val\n");
	die(ir);
    }
    
    data[ir].divide_val(ic,val);
}

template <typename T>
void asymm_array_2d<T>::add(int dex, T val){
    if(dex<0 || dex>=rows){
        printf("in asymm 2d add\n");
	die(dex);
    }
    
    data[dex].add(val);
}

template <typename T>
void asymm_array_2d<T>::replace_row(int dex, array_1d<T> &pt){
    if(dex<0 || dex>=rows){
        printf("WARNING trying to replace row %d in asymm, but only have %d\n",
	dex,rows);
	
	die(dex);
    }
    
    data[dex].reset();
    int i;
    for(i=0;i<pt.get_dim();i++){
        data[dex].add(pt.get_data(i));
    }

}

template <typename T>
array_1d<T>* asymm_array_2d<T>::operator()(int dex){
    
    if(dex<0 || dex>=rows){
        printf("WARNING asked for row %d but only have %d\n",dex,rows);
    }
    
    return &data[dex];

}

template class asymm_array_2d<int>;
template class asymm_array_2d<double>;
template class array_2d<int>;
template class array_2d<double>;
template class array_1d<int>;
template class array_1d<double>;


template <typename T>
void merge_sort(array_1d<T> &in, array_1d<int> &dexes,
                int start, int end){
    
    
    if(in.get_dim()==0) return;
    
    in.set_where("merge_sort");
   
    T nn;
    int i1,i2,el;
    
    el=end-start+1;
    
    if(el<2){
        return;
    }
    
    if(el==2){
    
        if(in.get_data(start)>in.get_data(end)){
	    nn=in.get_data(start);
	    in.set(start,in.get_data(end));
	    in.set(end,nn);
	    
	    i1=dexes.get_data(start);
	    dexes.set(start,dexes.get_data(end));
	    dexes.set(end,i1);
	}
	
	return;
    }
    
    int i_mid,*i_use;
    
    i_mid=(start+end)/2;
    
    merge_sort(in,dexes,start,i_mid);
    merge_sort(in,dexes,i_mid+1,end);
    
    array_1d<T> buffer;
    array_1d<int> dex_buffer;
    
    buffer.set_where("merge_sort");
    dex_buffer.set_where("merge_sort");
    
    buffer.set_name("merge_sort_buffer");
    dex_buffer.set_name("merge_sort_dex_buffer");
    
    
    for(i1=start,i2=i_mid+1;i1<=i_mid || i2<=end;){
        
	if(i2>end){
	    i_use=&i1;
	}
	else if(i1>i_mid){
	    i_use=&i2;
	}
	else if(in.get_data(i1)<in.get_data(i2)){
	    i_use=&i1;
        }
	else{
	    i_use=&i2;
	}
	
	//printf("using %d -- %d %d -- %e %e\n",i_use[0],i1,i2,
	//in.get_data(i1),in.get_data(i2));
	
	buffer.add(in.get_data(i_use[0]));
	dex_buffer.add(dexes.get_data(i_use[0]));
	
	i_use[0]++;
    }
    
    //printf("start %d end %d buffer dim %d\n",start,end,buffer.get_dim());
    
    for(i1=0;i1<el;i1++){
        in.set(start+i1,buffer.get_data(i1));
	dexes.set(start+i1,dex_buffer.get_data(i1));
    }
    

}

template <typename T>
double sort_and_check(const array_1d<T> &in, array_1d<T> &sorted, array_1d<int> &dexes){
    

    
    if(in.get_dim()!=dexes.get_dim()){
        printf("WARNING in sort_and_check in.dim %d dexes.dim %d\n",
	in.get_dim(),dexes.get_dim());
	
	exit(1);
    }
    
    if(in.get_dim()==0)return 0.0;
    
    in.set_where("sort_and_check");
    dexes.set_where("sort_and_check");
    in.set_where("sort_and_check");
    
    array_1d<int> dex_buffer;
    int i,j;
    
    dex_buffer.set_where("sort_and_check");
    
    dex_buffer.set_name("sort_and_check_dex_buffer");
    
    sorted.set_dim(in.get_dim());
    dex_buffer.set_dim(in.get_dim());
    
    for(i=0;i<in.get_dim();i++){
        dex_buffer.set(i,dexes.get_data(i));
	sorted.set(i,in.get_data(i));
    }
    
    merge_sort(sorted,dexes,0,in.get_dim()-1);
    
    double err,maxerr,aa,bb;
    
    int ifailure;
    
    for(i=0;i<in.get_dim();i++){
        if(i<in.get_dim()-1){
	    if(sorted.get_data(i+1)<sorted.get_data(i)){
	        printf("WARNING sort failed to get elements in proper order\n");
		
		ifailure=-1;
		
		throw ifailure;
		
	    }
	}
    
    
        for(j=0;j<in.get_dim() && dexes.get_data(i)!=dex_buffer.get_data(j);j++);
	
        if(dexes.get_data(i)!=dex_buffer.get_data(j)){
	    printf("WARNING dexes did not line up %d %d\n",
	    dexes.get_data(i),dex_buffer.get_data(j));
	    
	    ifailure=-1;
	    
	    throw ifailure;
	}
	
	aa=double(sorted.get_data(i));
	bb=double(in.get_data(j));
	
	err=fabs(aa-bb);
	if(fabs(aa)>0.0)err=err/fabs(aa);
	
	if(i==0 || err>maxerr){
	    maxerr=err;
	}
	
    }
   
    if(maxerr>1.0e-12){
        printf("WARNING associative error in merge_sort was %e\n",maxerr);
	
	try{
	  in.die(0);
	}
	catch(int iex){
	   try{
	       sorted.die(0);
	   }
	   catch(int jex){
	       try{
	           dexes.die(0);
	       }
	       catch(int kex){
	       
	       };
	   }
	}
	
	
	ifailure=-1;
	
	throw ifailure;
    }
    
    
    return maxerr;
    
}

template <typename T>
int get_dex(const array_1d<T> &xx, T target){
    int i;
    
    for(i=0;i<xx.get_dim()-1 && xx.get_data(i)<target;i++);
    

    
    if(i>0 && target-xx.get_data(i-1)<xx.get_data(i)-target){
        //printf("decrementing\n");
        i--;
    }
    
    return i;
}

template void merge_sort<double>(array_1d<double>&,array_1d<int>&,int,int);
template void merge_sort<int>(array_1d<int>&,array_1d<int>&,int,int);

template double sort_and_check<double>(const array_1d<double>&,array_1d<double>&,array_1d<int>&);
template double sort_and_check<int>(const array_1d<int>&,array_1d<int>&,array_1d<int>&);

template int get_dex(const array_1d<int>&,int);
template int get_dex(const array_1d<double>&,double);
