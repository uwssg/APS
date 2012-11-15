//this is the 27 February 2012 rewrite

class kd_tree{
 private:
  int **tree;
  int masterparent;
  int room,roomstep;
  double tol;
  int ncube,callcubes;
  
 
  
  void confirm(int,int,int,int);
  void organize(int*,int,int,int,int);
  int find_node(double*);
  void neigh_check(double*,int,int*,double*,int,int);
   void reassign(int);
   void descend(int);
  
 public:
  int dim,pts,diagnostic,xplr;
  
  int **cdex;
  double **cubecenter,*cubevol,**cubemax,**cubemin;
  
  double **data,*maxs,*mins;

  kd_tree(int,int,double**,double*,double*);
  ~kd_tree();
 
 void check_tree(int);
 double distance(double*,double*);
 void find_cubes();

 void write_tree(char*);
 void add(double*);
 void remove(int);
 void count(int,int*);
 void nn_srch(double*,int,int*,double*);

};
