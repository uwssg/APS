#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "chisq.h"

void chisquared::allot_arrays(){
    int ix,iy;
    
    bases.set_dim(dim,dim);
    widths.set_dim(ncenters,dim);
    centers.set_dim(ncenters,dim);     
    nboundary.set_dim(dim*dim);
    boundary_room.set_dim(dim*dim);

    
    for(ix=0;ix<dim*dim;ix++)nboundary.set(ix,0);
    for(ix=0;ix<dim*dim;ix++)boundary_room.set(ix,3);
    
    boundary=new double**[dim*dim];
    for(ix=0;ix<dim*dim;ix++){
        boundary[ix]=new double*[boundary_room.get_data(ix)];
        for(iy=0;iy<boundary_room.get_data(ix);iy++){
	    boundary[ix][iy]=new double[3];
	}
        
    }
    
    mins.set_dim(dim);
    maxs.set_dim(dim);

    for(ix=0;ix<dim;ix++){
        mins.set(ix,2.0*chisq_exception);
	maxs.set(ix,-2.0*chisq_exception);
    } 
    
    bases.set_name("chisq_bases");
    widths.set_name("chisq_widths");
    centers.set_name("chisq_centers");
    nboundary.set_name("chisq_nboundary");
    boundary_room.set_name("chisq_boundary_room");
    mins.set_name("chisq_mins");
    maxs.set_name("chisq_maxs");
    
}

void chisquared::set_max(int dex, double nn){
    maxs.set(dex,nn);
}

void chisquared::set_min(int dex, double nn){
    mins.set(dex,nn);
}

double chisquared::get_min(int dex){
    return mins.get_data(dex);
}

double chisquared::get_max(int dex){
    return maxs.get_data(dex);
}

double chisquared::get_time_spent(){
    return time_spent;
}

void chisquared::reset_boundary(){
    int i;
    for(i=0;i<dim*dim;i++)nboundary.set(i,0);
}

void chisquared::make_bases(int seed){

    int do_random_bases=1;
    
    if(seed==0){
        seed=int(time(NULL));
    }
    else if(seed<0){
        seed*=-1;
        do_random_bases=0;
    }
    
    if(dice==NULL){
        dice=new Ran(seed);
    }
    
    double nn;
    int i,j,ii,jj,goon;
    
    centers.set_where("chisq_make_bases");
    bases.set_where("chisq_make_bases");
    widths.set_where("chisq_make_bases");
    
    if(do_random_bases==1){
        for(ii=0;ii<dim;ii++){
            goon=1;
	    while(goon==1){
                goon=0;
	        for(i=0;i<dim;i++)bases.set(ii,i,dice->doub()-0.5);
	        for(jj=0;jj<ii;jj++){
	            nn=0.0;
		    for(i=0;i<dim;i++)nn+=bases.get_data(ii,i)*bases.get_data(jj,i);
		    for(i=0;i<dim;i++)bases.subtract_val(ii,i,nn*bases.get_data(jj,i));
	        }
	    
	        nn=0.0;
	        for(i=0;i<dim;i++){
		    nn+=power(bases.get_data(ii,i),2);
	        }
	        if(nn<1.0e-20)goon=1;
	        nn=sqrt(nn);
	        for(i=0;i<dim;i++){
	            bases.divide_val(ii,i,nn);
	        }
	    }
        }
    }//if do_random_bases==1
    else{
        for(i=0;i<dim;i++){
            for(jj=0;jj<dim;jj++){
                if(i==jj)bases.set(i,jj,1.0);
                else bases.set(i,jj,0.0);
            }
        }
    }
    
    double normerr,ortherr;
    for(ii=0;ii<dim;ii++){
        nn=0.0;
	for(i=0;i<dim;i++)nn+=power(bases.get_data(ii,i),2);
	nn=fabs(1.0-nn);
	if(ii==0 || nn>normerr)normerr=nn;
	
	for(jj=ii+1;jj<dim;jj++){
	   nn=0.0;
	   for(i=0;i<dim;i++)nn+=bases.get_data(ii,i)*bases.get_data(jj,i);
	   nn=fabs(nn);
	   if((ii==0 && jj==1) || nn>ortherr)ortherr=nn;
	}
    }
    
    printf("normerr %e ortherr %e\n",normerr,ortherr);
    if(normerr>1.0e-3 || ortherr>1.0e-3){
        death_knell("normerr or ortherr too large");
    }
    
    
    /*centers=new double*[ncenters];
    for(i=0;i<ncenters;i++)centers[i]=new double[dim];
    widths=new double*[ncenters];
    for(i=0;i<ncenters;i++)widths[i]=new double[dim];*/
    
    double rr,theta;
    array_1d<double> trial_center,trial_pt;
    int acceptable,iterations=0;;
    
    trial_center.set_dim(dim);
    trial_pt.set_dim(dim);

    goon=0;
    while(goon==0){
        iterations++;
        
        if(iterations>5000){
            printf("WARNING; chisq was unable to construct an acceptable function\n");
            
            exit(1);
        }
        
	for(ii=0;ii<ncenters;ii++){
            for(i=0;i<dim;i++)centers.set(ii,i,1.0e30);
	    for(i=0;i<dim;i++)widths.set(ii,i,1.0e-5);
        }
	
	goon=1;
        for(ii=0;ii<ncenters;ii++){
	    
	    acceptable=1;
	    
	    for(i=0;i<dim;i++)trial_center.set(i,0.0);
	    
	    if(ii>0){
	        rr=normal_deviate(dice,40.0,20.0);
	        theta=dice->doub()*2.0*pi;
	        
		trial_center.set(0,centers.get_data(0,0)+rr*cos(theta));
		trial_center.set(1,centers.get_data(0,1)+rr*sin(theta));
	 
	    } 
	     
	    for(i=0;i<dim;i++){
                trial_center.add_val(i,normal_deviate(dice,30.0,15.0));
	    }
	    
	    for(i=0;i<dim;i++){
		trial_pt.set(i,0.0);
	        for(j=0;j<dim;j++)trial_pt.add_val(i,trial_center.get_data(j)*bases.get_data(j,i));
	    }
	    
	    rr=(*this)(trial_pt);
	    
	    if(rr<100.0)acceptable=0;
	    
	    if(acceptable==1){
	        for(i=0;i<dim;i++){
		    centers.set(ii,i,trial_center.get_data(i));
		    widths.set(ii,i,fabs(normal_deviate(dice,2.0,0.5))+0.1);
	        }
	    }
	    else ii--;
	
        }
	
	
	
	
	
	for(i=0;i<dim && acceptable==1;i++){
	    for(ii=0;ii<ncenters && acceptable==1;ii++){
                if(centers.get_data(ii,i)-4.0*widths.get_data(ii,i)<-100.0)acceptable=0;
                if(centers.get_data(ii,i)+4.0*widths.get_data(ii,i)>100.0)acceptable=0;
                
	        for(jj=ii+1;jj<ncenters && acceptable==1;jj++){
		    nn=fabs(centers.get_data(ii,i)-centers.get_data(jj,i));
		    if(nn<2.0*widths.get_data(ii,i) || nn<2.0*widths.get_data(jj,i)){
                        acceptable=0;
                        printf("centers %d %d -- %d -- %e %e -- %e %e\n",ii,jj,i,centers.get_data(ii,i),widths.get_data(ii,i),
                        centers.get_data(jj,i),widths.get_data(jj,i));
                    }
		}
	    }
	}
        if(acceptable==0)goon=0;
	
    }
    
    centers.set_where("nowhere");
    bases.set_where("nowhere");
    widths.set_where("nowhere");
    
    time_spent=0.0;
    called=0;
    
    printf("set centers and widths %d %d\n",dim,ncenters);
}

void chisquared::add_to_boundary(array_1d<double> &alpha, int ix, int iy,double chitest){
    
    array_2d<double> buffer;
    array_1d<double> pt;
    double rr;
    int ipt,room;
    int i,j;
    
    buffer.set_name("chisq_add_to_boundary_buffer");
    pt.set_name("chisq_add_to_boundary_pt");
    
    //printf("adding %e\n",chitest);
    if(ix>iy){
        i=ix;
	ix=iy;
	iy=i;
    }
    
    ipt=ix*dim+iy;
    room=boundary_room.get_data(ipt);
    
    if(nboundary.get_data(ipt)>=room){
  
	buffer.set_dim(room,3);
	for(i=0;i<room;i++){
	    for(j=0;j<3;j++)buffer.set(i,j,boundary[ipt][i][j]);
	    delete [] boundary[ipt][i];
	}
        delete [] boundary[ipt];
	boundary_room.add_val(ipt,100);
	boundary[ipt]=new double*[boundary_room.get_data(ipt)];
	for(i=0;i<boundary_room.get_data(ipt);i++)boundary[ipt][i]=new double[3];
	
	for(i=0;i<room;i++){
	    for(j=0;j<3;j++)boundary[ipt][i][j]=buffer.get_data(i,j);
	}
	buffer.reset();
    }
    
    boundary[ipt][nboundary.get_data(ipt)][0]=alpha.get_data(ix);
    boundary[ipt][nboundary.get_data(ipt)][1]=alpha.get_data(iy);
    boundary[ipt][nboundary.get_data(ipt)][2]=chitest;
    
    nboundary.add_val(ipt,1);
    
    rr=0.0;
    for(i=0;i<dim;i++)rr+=alpha.get_data(i)*alpha.get_data(i);
    rr=sqrt(rr);
    if(rr>rr_max)rr_max=rr;
    
    pt.set_dim(dim);
    for(i=0;i<dim;i++)pt.set(i,0.0);
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
	    pt.add_val(i,alpha.get_data(j)*bases.get_data(j,i));
	}
    }
    
    double nn=0.0;
    for(i=0;i<dim;i++)nn+=pt.get_data(i)*pt.get_data(i);
    nn=sqrt(nn);
    if(fabs(nn-rr)>1.0e-4){
        printf("WARNING nn %e rr %e\n",nn,rr);
    }
    
    //not sure why this is here
    //I think I was using it for analysis of cartoons
    /*for(i=0;i<dim;i++){
        if(pt.get_data(i)<mins.get_data(i))mins.set(i,pt.get_data(i));
	if(pt.get_data(i)>maxs.get_data(i))maxs.set(i,pt.get_data(i));
    }*/
    
}

void chisquared::print_mins_maxs(){
    int i;
    double nn;
    nn=0.0;
    printf("mins and maxs\n");
    for(i=0;i<dim;i++){
        printf("p%d -- %e %e -- %e\n",i,mins.get_data(i),maxs.get_data(i),maxs.get_data(i)-mins.get_data(i));
	nn+=power(maxs.get_data(i)-mins.get_data(i),2);
    }
    printf("\nfiducial distance %e\n",sqrt(nn));
}

double chisquared::get_rr_max(){
    return rr_max;
}

int chisquared::get_n_boundary(int ix, int iy){
    
    int i;
    if(ix>iy){
        i=ix;
	ix=iy;
	iy=i;
    }
    
    return nboundary.get_data(ix*dim+iy);
}

double chisquared::get_width(int ic, int ix){
  
    return widths.get_data(ic,ix);
   
}

double chisquared::get_center(int ic, int ix){
    return centers.get_data(ic,ix);
}

double chisquared::get_real_center(int ic, int ix){
    if(ic>=ncenters || ix>=dim){
        return chisq_exception;
    }
    
    int i;
    double ans=0.0;
    for(i=0;i<dim;i++){
        ans+=centers.get_data(ic,i)*bases.get_data(i,ix);
    }
    return ans;
    
}

double chisquared::distance_to_center(int ic, array_1d<double> &pt){
    if(ic<0 || ic>=centers.get_rows()){
        printf("WARHNG asked for center %d but max %d\n",
	ic,centers.get_rows());
	
	exit(1);
    }
    
    array_1d<double> projected;
    
    int i;
    
    for(i=0;i<centers.get_cols();i++){
        projected.set(i,project_to_basis(i,pt));
    }
    
    double dd;
    dd=0.0;
    for(i=0;i<centers.get_cols();i++){
        dd+=power(projected.get_data(i)-centers.get_data(ic,i),2);
    }
    dd=sqrt(dd);
    return dd;
}

double chisquared::get_boundary(int ix, int iy, int ipt, int idim){

    int i;
    if(ix>iy){
        i=ix;
	ix=iy;
	iy=i;
    }
    
    if(ix<0 || iy<0 || ix>=dim || iy>=dim){
        printf("WARNING asked for boundary slog %d %d but dim %d\n",ix,iy,dim);
    }
    
    if(idim>=3 || ipt>=nboundary.get_data(ix*dim+iy))return chisq_exception;
    return boundary[ix*dim+iy][ipt][idim];
}

void chisquared::death_knell(char *word)const{
    printf("%s\n",word);
    exit(1);
}

int chisquared::get_dim(){
    return dim;
}

int chisquared::get_ncenters(){
    return ncenters;
}

chisquared::chisquared(){
    death_knell("meaningless constructor");
};

chisquared::chisquared(int id){
    ncenters=1;
    dim=id;
    
    time_spent=0.0;
    
    boundary=NULL;
    dice=NULL;
    
    rr_max=-1.0;
    called=0;
    
    allot_arrays();
};

chisquared::chisquared(int id, int ic){
    dim=id;
    ncenters=ic;
    
    time_spent=0.0;
 
    boundary=NULL;   
    dice=NULL;
     
    rr_max=-1.0;
    called=0;
    
    allot_arrays();
};



chisquared::~chisquared(){
    int i,ix,iy;
    
    if(boundary!=NULL){
        for(ix=0;ix<dim;ix++){
	    for(iy=ix+1;iy<dim;iy++){
	            for(i=0;i<boundary_room.get_data(ix*dim+iy);i++){
		        delete [] boundary[ix*dim+iy][i];
		    }
		    delete [] boundary[ix*dim+iy];
		
	    }
	}
	delete [] boundary;
    }
    
    if(dice!=NULL){
        delete dice;
    }

}

int chisquared::get_called(){
    return called;
}

void chisquared::decrement_called(){
    called--;
}

double chisquared::operator()(array_1d<double> &v)const{
    death_knell("meaningless operator");
}

void chisquared::build_boundary(double rr){
    death_knell("meaningless build_boundary");
}

void chisquared::get_basis(int ix, array_1d<double> &v){
    int i;
    if(ix<dim){
        for(i=0;i<dim;i++)v.set(i,bases.get_data(ix,i));
    }
    else{
        printf("WARNING called get_basis with %d %d\n",ix,dim);
	exit(1);
    }
}

double chisquared::project_to_basis(int ix, array_1d<double> &vv) const{
    int i;
    double nn=1.0e30;
    if(ix<dim){
        nn=0.0;
	for(i=0;i<dim;i++)nn+=vv.get_data(i)*bases.get_data(ix,i);
    }
    else{
        printf("WARNING called project_to_basis with %d %d\n",ix,dim);
    }
    
    return nn;
}

s_curve::~s_curve(){}

s_curve::s_curve() : chisquared(6), trig_factor(10.0){make_bases(22);}

s_curve::s_curve(int id) : chisquared(id), trig_factor(10.0){make_bases(22);}

s_curve::s_curve(int id, int ic) : chisquared(id,ic), trig_factor(10.0){
        make_bases(22);}


ellipses::~ellipses(){}

ellipses::ellipses() : chisquared(22){make_bases(13);}

ellipses::ellipses(int id) : chisquared(id){make_bases(13);}

ellipses::ellipses(int id, int ic) : chisquared(id,ic){
    //printf("constructed dim %d ncenters %d\n",dim,ncenters);
    make_bases(13);
}

ellipses_integrable::~ellipses_integrable(){}

ellipses_integrable::ellipses_integrable() : ellipses(){
    make_bases(-13);
    
    /*int i,j;
    for(i=1;i<ncenters;i++){
        for(j=0;j<dim;j++){
            widths.set(i,j,widths.get_data(0,j));
        }
    }*/

}

ellipses_integrable::ellipses_integrable(int id) : ellipses(id){
    make_bases(-13);

    /*int i,j;
    for(i=1;i<ncenters;i++){
        for(j=0;j<dim;j++){
            widths.set(i,j,widths.get_data(0,j));
        }
    }*/
    
}

ellipses_integrable::ellipses_integrable(int id, int ic) : ellipses(id,ic){
    make_bases(-13);

    /*int i,j;
    for(i=1;i<ncenters;i++){
        for(j=0;j<dim;j++){
            widths.set(i,j,widths.get_data(0,j));
        }
    }*/


}

linear_ellipses::~linear_ellipses(){}

linear_ellipses::linear_ellipses() : chisquared(22){make_bases(17);}

linear_ellipses::linear_ellipses(int id) : chisquared(id){make_bases(17);}

linear_ellipses::linear_ellipses(int id, int ic) : chisquared(id,ic){
    make_bases(17);
}

double s_curve::operator()(array_1d<double> &in_pt) const{
    
    if(dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    
    double before=double(time(NULL));
    
    called++;

    
    array_1d<double> dd,pt;
    dd.set_dim(ncenters);
    pt.set_dim(dim);
    
    pt.set_name("s_curve_operator_pt");
    dd.set_name("s_curve_operator_dd");
    
    centers.set_where("s_curve_operator");
    widths.set_where("s_curve_operator");
    bases.set_where("s_curve_operator");
    
    
    double theta,xth,yth,dth,dthmin;
    int i;
    
    for(i=0;i<dim;i++)pt.set(i,project_to_basis(i,in_pt));
    
    //for(i=0;i<dim;i++)printf("%e %e\n",pt[i],centers[0][i]);
    dd.set(0,0.0);
    for(i=2;i<dim;i++)dd.add_val(0,power((pt.get_data(i)-centers.get_data(0,i))/widths.get_data(0,i),2));
    
    dthmin=-1.0;
    for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
        xth=trig_factor*sin(theta)+centers.get_data(0,0);
	yth=trig_factor*theta*(cos(theta)-1.0)/fabs(theta)+centers.get_data(0,1);
	dth=power((xth-pt.get_data(0))/widths.get_data(0,0),2)+power((yth-pt.get_data(1))/widths.get_data(0,1),2);
	
	if(dthmin<0.0 || dth<dthmin)dthmin=dth;
    }
    dd.add_val(0,dthmin);
    
    int ii;
    for(ii=1;ii<ncenters;ii++){
        dd.set(ii,0.0);
	for(i=0;i<dim;i++){
	    dd.add_val(ii,power((pt.get_data(i)-centers.get_data(ii,i))/widths.get_data(ii,i),2));
	}
    }
    
    double ddmin;
    for(ii=0;ii<ncenters;ii++){
        if(ii==0 || dd.get_data(ii)<ddmin)ddmin=dd.get_data(ii);
    }
    
    centers.set_where("nowhere");
    bases.set_where("nowhere");
    widths.set_where("nowhere");
    
    time_spent+=double(time(NULL))-before;
    
    return ddmin;
}

double s_curve::distance_to_center(int ic, array_1d<double> &in_pt) {
    
    if(dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    if(ic<0 || ic>=centers.get_rows()){
        printf("WARNING asked for center %d but max %d\n",
	ic,centers.get_rows());
	
	exit(1);
    }
    

    
    array_1d<double> pt;
    
    pt.set_dim(dim);
    
    pt.set_name("s_curve_distance_to_center_pt");
   
    centers.set_where("s_curve_operator");
    widths.set_where("s_curve_operator");
    bases.set_where("s_curve_operator");
    
    
    double theta,xth,yth,dth,dthmin,dd;
    int i;
    
    for(i=0;i<dim;i++)pt.set(i,project_to_basis(i,in_pt));
    
    
    //for(i=0;i<dim;i++)printf("%e %e\n",pt[i],centers[0][i]);
    dd=0.0;
    
    if(ic==0){
        for(i=2;i<dim;i++)dd+=power((pt.get_data(i)-centers.get_data(0,i))/widths.get_data(0,i),2);
    
        dthmin=-1.0;
        for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
            xth=trig_factor*sin(theta)+centers.get_data(0,0);
	    yth=trig_factor*theta*(cos(theta)-1.0)/fabs(theta)+centers.get_data(0,1);
	    dth=power((xth-pt.get_data(0))/widths.get_data(0,0),2)+power((yth-pt.get_data(1))/widths.get_data(0,1),2);
	
	    if(dthmin<0.0 || dth<dthmin)dthmin=dth;
        }
        dd+=dth;
    }
    else{

	for(i=0;i<dim;i++){
	    dd+=power((pt.get_data(i)-centers.get_data(ic,i))/widths.get_data(ic,i),2);
	}
    
    }
    
    dd=sqrt(dd);
    
    centers.set_where("nowhere");
    bases.set_where("nowhere");
    widths.set_where("nowhere");
    
    return dd;
}

void s_curve::build_boundary(double br){
    if(dice==NULL){
        death_knell("you called build_boundary before making bases");
    }

    reset_boundary();

    int ix,iy,ic,ir,i,j;
    
    double tol=0.5;
    double theta,dxdth,dydth,x0,y0;
    double grad[2],norm,dfabsdth,ds,chitest;
    
    array_1d<double> pt,alpha;
    
    alpha.set_dim(dim);
    pt.set_dim(dim);
    
    alpha.set_name("s_curve_build_boundary_alpha");
    pt.set_name("s_curve_build_boundary_pt");
    
    centers.set_where("s_curve_build_boundary");
    bases.set_where("s_curve_build_boundary");
    widths.set_where("s_curve_build_boundary");
    
    ///////////////below we will do ix=0, iy=1 for the 0th center (the S curve itself)
    ix=0;
    iy=1;
    
  
    
    if(widths.get_data(0,0)<widths.get_data(0,1))ds=0.1*widths.get_data(0,0);
    else ds=0.1*widths.get_data(0,1);
    
    for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
        if(theta<0.0)dfabsdth=-1.0;
        else dfabsdth=1.0;
    
        x0=centers.get_data(0,0)+trig_factor*sin(theta);
        y0=centers.get_data(0,1)+trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
    
        for(ir=0;ir<2;ir++){

            dxdth=trig_factor*cos(theta);
	    dydth=trig_factor*((cos(theta)-1.0)/fabs(theta)-sin(theta)*theta/fabs(theta)
	               -(cos(theta)-1.0)*dfabsdth/theta);
	
            if(ir==0){
	        grad[0]=dydth;
	        grad[1]=-1.0*dxdth;
	    }
	    else{
	        grad[0]=-1.0*dydth;
	        grad[1]=dxdth;
	    }
            norm=grad[0]*grad[0]+grad[1]*grad[1];
	    norm=sqrt(norm);
	
	    for(i=0;i<dim;i++)alpha.set(i,centers.get_data(0,i));
	

	    alpha.set(0,x0+grad[0]*ds/norm);
	    alpha.set(1,y0+grad[1]*ds/norm);
	    
            alpha.set(0,x0);
	    alpha.set(1,y0);
	    

	    for(i=0;i<dim;i++){
		pt.set(i,0.0);
	        for(j=0;j<dim;j++)pt.add_val(i,alpha.get_data(j)*bases.get_data(j,i));
	    }
	    chitest=(*this)(pt);
	    
	    if(chitest>br)death_knell("started outside the pale\n");
	    
	    while(chitest<br){
	        
	        alpha.add_val(0,ds*grad[0]/norm);
	        alpha.add_val(1,ds*grad[1]/norm);
	    
	        for(i=0;i<dim;i++){

		    pt.set(i,0.0);
		    for(j=0;j<dim;j++)pt.add_val(i,alpha.get_data(j)*bases.get_data(j,i));
	        }
	        chitest=(*this)(pt);
		
		
	    }
	    
	    if(fabs(chitest-br)<tol){
	        add_to_boundary(alpha,ix,iy,chitest);
	    }
	    /*for(i=0;i<2;i++){
	        printf("%e ",alpha[i]);
	    }
	    printf("%e\n",chitest);*/
	    
		
        }
    }
    
    for(i=0;i<dim;i++)alpha.set(i,centers.get_data(0,i));
    
    double th,s2x,s2y,aa,rr,thmin,thmax;
    for(theta=-1.0*pi;theta<=1.5*pi;theta+=2.0*pi){
        for(i=0;i<dim;i++)alpha.set(i,centers.get_data(0,i));
    
        x0=centers.get_data(0,0)+trig_factor*sin(theta);
        y0=centers.get_data(0,1)+trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
    
        s2x=widths.get_data(0,0)*widths.get_data(0,0);
        s2y=widths.get_data(0,1)*widths.get_data(0,1);
        
	if(theta<0.0){
	    thmin=-0.5*pi;
	    thmax=0.5*pi;
	}
	else{
	   thmin=0.5*pi;
	   thmax=1.5*pi;
	}
	
        for(th=thmin;th<thmax;th+=0.01){
        
            aa=cos(th)*cos(th)/s2x+sin(th)*sin(th)/s2y;
	    rr=br/aa;
	    rr=sqrt(rr);
	
	    alpha.set(0,x0+rr*cos(th));
	    alpha.set(1,y0+rr*sin(th));
	
	    for(i=0;i<dim;i++){
		pt.set(i,0.0);
	        for(j=0;j<dim;j++)pt.add_val(i,alpha.get_data(j)*bases.get_data(j,i));
	    }
    
            chitest=(*this)(pt);
	    
	    if(fabs(chitest-br)<5.0){
	        add_to_boundary(alpha,ix,iy,chitest);
	    }
	    /*for(i=0;i<2;i++)printf("%e ",alpha[i]);
	    printf("%e\n",chitest);*/
        }

    }
    
    /////////////////////now do ix=0,1; iy>=2
    double xx,yy,yth,xth;

    for(iy=2;iy<dim;iy++){for(ir=0;ir<2;ir++){
        for(i=0;i<dim;i++)alpha.set(i,centers.get_data(0,i));
        for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
           alpha.set(0,centers.get_data(0,0)+trig_factor*sin(theta));
           alpha.set(1,centers.get_data(0,1)
                   +trig_factor*theta*(cos(theta)-1.0)/fabs(theta));
       
       
           if(ir==0)alpha.set(iy,centers.get_data(0,iy)+sqrt(br)*widths.get_data(0,iy));
           else alpha.set(iy,centers.get_data(0,iy)-sqrt(br)*widths.get_data(0,iy));
       
           for(i=0;i<dim;i++){
	       pt.set(i,0.0);
	       for(j=0;j<dim;j++)pt.add_val(i,alpha.get_data(j)*bases.get_data(j,i));
           } 
           chitest=(*this)(pt);
           
	   if(fabs(chitest-br)<5.0){
	       add_to_boundary(alpha,0,iy,chitest);
	       add_to_boundary(alpha,1,iy,chitest);
	   }
        }
    }}
   

    double thetamin,thetamax,dtheta;

    for(ix=0;ix<2;ix++){
        if(ix==1){
            thetamin=-1.0*pi;
	    thetamax=1.5*pi;
	    dtheta=2.0*pi;
        }
        else{
            thetamin=-0.5*pi;
	    thetamax=0.6*pi;
	    dtheta=1.0*pi;
    
        }


        s2x=widths.get_data(0,ix)*widths.get_data(0,ix);
        for(iy=2;iy<dim;iy++){
            s2y=widths.get_data(0,iy)*widths.get_data(0,iy);
            for(i=0;i<dim;i++)alpha.set(i,centers.get_data(0,i));
            for(theta=thetamin;theta<thetamax;theta+=dtheta){
	     
	         if(ix==0){
	             xth=centers.get_data(0,0)+trig_factor*sin(theta);
	             alpha.set(1,centers.get_data(0,1)
		        +trig_factor*theta*(cos(theta)-1.0)/fabs(theta));
	         }
	         else{
	             xth=centers.get_data(0,1)+trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
	             alpha.set(0,centers.get_data(0,0)+trig_factor*sin(theta));
	         }
	     
	     
	         for(th=0.0;th<2.0*pi;th+=0.01){
	             aa=(cos(th)*cos(th)/s2x+sin(th)*sin(th)/s2y);
		     rr=br/aa;
		     rr=sqrt(rr);
		 
		     alpha.set(iy,centers.get_data(0,iy)+rr*sin(th));
		     alpha.set(ix,xth+rr*cos(th));
	         
		 
		     for(i=0;i<dim;i++){
			 pt.set(i,0.0);
		         for(j=0;j<dim;j++){
		             pt.add_val(i,alpha.get_data(j)*bases.get_data(j,i));
		         }
		     }
		     chitest=(*this)(pt);
		     if(fabs(chitest-br)<5.0){
		         add_to_boundary(alpha,ix,iy,chitest);
		     }
		 
	         }
	     
	    }
        }
    }
    
    for(ic=0;ic<ncenters;ic++){
        for(ix=0;ix<dim;ix++){
	    for(iy=ix+1;iy<dim;iy++){
	        if(ic>0 || ix>1){
		    
		    for(i=0;i<dim;i++)alpha.set(i,centers.get_data(ic,i));
		    s2x=power(widths.get_data(ic,ix),2);
		    s2y=power(widths.get_data(ic,iy),2);
		    for(theta=0.0;theta<=2.0*pi;theta+=0.01){
		        aa=power(cos(theta),2)/s2x+power(sin(theta),2)/s2y;
			rr=br/aa;
			rr=sqrt(rr);
			alpha.set(ix,centers.get_data(ic,ix)+rr*cos(theta));
			alpha.set(iy,centers.get_data(ic,iy)+rr*sin(theta));
		        
			for(i=0;i<dim;i++){
			    pt.set(i,0.0);
			    for(j=0;j<dim;j++)pt.add_val(i,alpha.get_data(j)*bases.get_data(j,i));
			}
			
			chitest=(*this)(pt);
			if(fabs(chitest-br)<5.0){
			    //if(ic==1 && iy==1)printf("adding %d %d %d\n",ix,iy,ic);
			
			    add_to_boundary(alpha,ix,iy,chitest);
			}
		    }
		
		}
	    }
	}
    }
    
    centers.set_where("nowhere");
    widths.set_where("nowhere");
    bases.set_where("nowhere");
  
    
}

double ellipses::operator()(array_1d<double> &in_pt) const{
    
    if(dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    double before=double(time(NULL));
    
    int ii,ix;
    double dd,ddmin,nn;
    
    centers.set_where("ellipse_operator");
    widths.set_where("ellipse_operator");
    bases.set_where("ellipse_operator");
    
    called++;
    
    for(ii=0;ii<ncenters;ii++){
        dd=0.0;
        for(ix=0;ix<dim;ix++){
	    nn=project_to_basis(ix,in_pt);
	    dd+=power((centers.get_data(ii,ix)-nn)/widths.get_data(ii,ix),2);
	    //printf("%e\n",centers[ii][ix]-nn);
	}
	if(ii==0 || dd<ddmin)ddmin=dd;
    }
    
    centers.set_where("nowhere");
    widths.set_where("nowhere");
    bases.set_where("nowhere");
    
    time_spent+=double(time(NULL))-before;
    
    return ddmin;
    
}

void ellipses::build_boundary(double br){

    if(dice==NULL){
         death_knell("you called build_boundary before making bases");
    } 
    
    reset_boundary();
    
    double tol=0.5;
    int ix,iy,ic,i,j;
    double theta,rr,aa,chitest;
    
    /*
    printf("in build_boundary centers are\n");
    for(i=0;i<dim;i++){
        for(j=0;j<ncenters;j++){
            printf("%e ",centers.get_data(j,i));
        }
        printf("\n");
    }
    */
    
    array_1d<double> alpha,pt;
    alpha.set_dim(dim);
    pt.set_dim(dim);
    
    pt.set_name("ellipse_build_boundary_pt");
    alpha.set_name("ellipse_build_boundary_alpha");
    
    centers.set_where("ellipse_build_boundary");
    bases.set_where("ellipse_build_boundary");
    widths.set_where("ellipse_build_boundary");
    
    
    for(ic=0;ic<ncenters;ic++){
        for(ix=0;ix<dim;ix++){
            for(iy=ix+1;iy<dim;iy++){
	        for(i=0;i<dim;i++)alpha.set(i,centers.get_data(ic,i));
	        for(theta=0.0;theta<=2.0*pi;theta+=0.01){
	            aa=power(cos(theta)/widths.get_data(ic,ix),2)+power(sin(theta)/widths.get_data(ic,iy),2);
		    rr=br/aa;
		    rr=sqrt(rr);
		    alpha.set(ix,centers.get_data(ic,ix)+rr*cos(theta));
		    alpha.set(iy,centers.get_data(ic,iy)+rr*sin(theta));
		    
		    
		    for(i=0;i<dim;i++){
		        pt.set(i,0.0);
			for(j=0;j<dim;j++)pt.add_val(i,alpha.get_data(j)*bases.get_data(j,i));
		    }
		    
		    chitest=(*this)(pt);
		    if(fabs(chitest-br)<tol){
		        add_to_boundary(alpha,ix,iy,chitest);
		    }
		    else{
		        printf("failed to add %d %e rr %e\n",ic,chitest,rr);
			exit(1);
		    }
	        }
	    }
        }
    }
    
    centers.set_where("nowhere");
    bases.set_where("nowhere");
    widths.set_where("nowhere");
    
}

double linear_ellipses::operator()(array_1d<double> &in_pt) const{

    if(dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    double before=double(time(NULL));
    
    int ii,ix;
    double dd,ddmin,nn;
    
    centers.set_where("linear_ellipse_operator\n");
    bases.set_where("linear_ellipse_operator\n");
    widths.set_where("linear_ellipse_operator\n");
    
    called++;
    
    for(ii=0;ii<ncenters;ii++){
        dd=0.0;
	for(ix=0;ix<dim;ix++){
	    nn=project_to_basis(ix,in_pt);
	    dd+=power((centers.get_data(ii,ix)-nn)/widths.get_data(ii,ix),2);
	}
	
	if(ii==0 || dd<ddmin)ddmin=dd;
    }
    
    centers.set_where("nowhere");
    bases.set_where("nowhere");
    widths.set_where("nowhere");
    
    time_spent+=double(time(NULL))-before;
    
    return sqrt(ddmin);
}

void linear_ellipses::build_boundary(double br){
    if(dice==NULL){
         death_knell("you called build_boundary before making bases");
    } 
    
    reset_boundary();
    
    double tol=0.5;
    int ic,ix,iy,i,j;
 
    double theta,rr,aa,chitest;
    
    array_1d<double> alpha,pt;
    
    alpha.set_dim(dim);
    pt.set_dim(dim);
    
    alpha.set_name("linear_ellipse_build_boundary_alpha");
    pt.set_name("linear_ellipse_build_boundary_pt");
    
    centers.set_where("linear_ellipse_build_boundary");
    bases.set_where("linear_ellipse_build_boundary");
    widths.set_where("linear_ellipse_build_boundary");
    
 
    
    for(ic=0;ic<ncenters;ic++){
        for(ix=0;ix<dim;ix++){
	    for(iy=ix+1;iy<dim;iy++){
	        for(i=0;i<dim;i++)alpha.set(i,centers.get_data(ic,i));
		for(theta=0.0;theta<=2.0*pi;theta+=0.01){
		    aa=power(cos(theta)/widths.get_data(ic,ix),2)+power(sin(theta)/widths.get_data(ic,iy),2);
		    rr=br*br/aa;
		    rr=sqrt(rr);
		    
		    alpha.set(ix,centers.get_data(ic,ix)+rr*cos(theta));
		    alpha.set(iy,centers.get_data(ic,iy)+rr*sin(theta));
		    for(i=0;i<dim;i++){
		        pt.set(i,0.0);
			for(j=0;j<dim;j++)pt.add_val(i,alpha.get_data(j)*bases.get_data(j,i));
		    }
		    chitest=(*this)(pt);
		    if(fabs(chitest-br)<tol){
		        add_to_boundary(alpha,ix,iy,chitest);
		    }
		    else{
		       printf("failed to add %d %e\n",ic,chitest);
		       exit(1);
		    }
		}
	    }
	}
    }
    
    centers.set_where("nowhere");
    bases.set_where("nowhere");
    widths.set_where("nowhere");
    
}

void ellipses_integrable::integrate_boundary(int ix1, int ix2, double lim, char *filename){
    integrate_boundary(ix1,ix2,lim,filename,5.0);
}

void ellipses_integrable::integrate_boundary(int ix1, int ix2, double lim, char *filename, double spread){
    
    int i,j;
    /*
    printf("in integrate_boundary, centers are\n");
    for(i=0;i<dim;i++){
        for(j=0;j<ncenters;j++){
            printf("%e ",centers.get_data(j,i));
        }
        printf("\n");
    }
    */
    
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            if(i==j){
                if(fabs(bases.get_data(i,j)-1.0)>1.0e-6){
                    printf("WARNING bases %d %d %e\n",i,j,bases.get_data(i,j));
                    exit(1);
                }
            }
            else{
                if(fabs(bases.get_data(i,j))>1.0e-6){
                    printf("WARNING bases %d %d %e\n",i,j,bases.get_data(i,j));
                    exit(1);
                }
            }
        }
    }
    
    /*for(i=0;i<dim;i++){
        for(j=1;j<ncenters;j++){
            if(fabs(widths.get_data(0,i)-widths.get_data(j,i))/widths.get_data(0,i)>1.0e-6){
                printf("WARNING widths do not agree %e %e %d %d\n",
                widths.get_data(0,i),widths.get_data(j,i),j,i);
                exit(1);
            }
        }
    }*/
    
    int icenter,k,imin,ct=0,duplicate;
    array_1d<double> chiarr,xarr,yarr,trial;
    array_1d<int> dexes;
    double xx,yy,xwidth,ywidth,nn,mm,ddmin,ddshld,total=0.0;
    double xbound,ybound,chival,xstep,ystep;
    
    array_1d<double> l_marginalized_volume;
    
    for(icenter=0;icenter<ncenters;icenter++){
        l_marginalized_volume.set(icenter,0.0);
        for(i=0;i<dim;i++){
            if(i!=ix1 && i!=ix2){
                l_marginalized_volume.add_val(icenter,log(widths.get_data(icenter,i)));
            }
        }
    }
    
    /*for(icenter=0;icenter<ncenters;icenter++){
        printf("l_vol %e \n",l_marginalized_volume.get_data(icenter));
    }*/
    
    for(icenter=0;icenter<ncenters;icenter++){
        if(icenter==0 || widths.get_data(icenter,ix1)>xbound)xbound=widths.get_data(icenter,ix1);
        if(icenter==0 || widths.get_data(icenter,ix1)<xwidth)xwidth=widths.get_data(icenter,ix1);
        
        if(icenter==0 || widths.get_data(icenter,ix2)>ybound)ybound=widths.get_data(icenter,ix2);
        if(icenter==0 || widths.get_data(icenter,ix2)<ywidth)ywidth=widths.get_data(icenter,ix2);
    }
    
    dexes.set_dim(100);
    xarr.set_dim(100);
    yarr.set_dim(100);
    chiarr.set_dim(100);
    
    array_1d<int> doubled_up;
    doubled_up.set_name("doubled_up");
    
    xstep=0.1*xwidth;
    ystep=0.1*ywidth;
     
    for(icenter=0;icenter<ncenters;icenter++){
        //xwidth=widths.get_data(icenter,ix1);
        //ywidth=widths.get_data(icenter,ix2);
        
        //xbound=widths.get_data(icenter,ix1);
        //ybound=widths.get_data(icenter,ix2);
        
        for(xx=centers.get_data(icenter,ix1)-spread*xbound;xx<centers.get_data(icenter,ix1)+(spread+0.1)*xbound;xx+=xstep){
            for(yy=centers.get_data(icenter,ix2)-spread*ybound;yy<centers.get_data(icenter,ix2)+(spread+0.1)*ybound;yy+=ystep){
                
                for(k=0;k<ncenters;k++){
                    for(i=0;i<dim;i++)trial.set(i,centers.get_data(k,i));
                    trial.set(ix1,xx);
                    trial.set(ix2,yy);
                    nn=0.0;
                    for(i=0;i<dim;i++){
                        nn+=power((trial.get_data(i)-centers.get_data(k,i))/widths.get_data(k,i),2);
                    }
                    
                    if(k==0 || nn<ddmin){
                        ddmin=nn;
                        imin=k;
                    }
                    if(k==icenter)ddshld=nn;
                }
                
                /*if(imin!=icenter){
                    
                    printf("dim %d %d\n",ix1,ix2);
                    printf("CURIOUS icenter %d imin %d\n",imin,icenter);
                    
                    printf("%e %e -- %e %e %e %e -- %e %e %e %e\n",
                    xx,yy,
                    centers.get_data(icenter,ix1),widths.get_data(icenter,ix1),
                    centers.get_data(icenter,ix2),widths.get_data(icenter,ix2),
                    centers.get_data(imin,ix1),widths.get_data(imin,ix1),
                    centers.get_data(imin,ix2),widths.get_data(imin,ix2));
                    
                    printf("%e %e\n",ddmin,ddshld);
                    
                    
                    exit(1);
                }*/
                
                for(i=0;i<dim;i++)trial.set(i,centers.get_data(icenter,i));
                trial.set(ix1,xx);
                trial.set(ix2,yy);
                
                chival=(*this)(trial)-2.0*l_marginalized_volume.get_data(icenter);
                
                
                if(icenter==imin){
                    duplicate=-1;
                }
                else{
                    duplicate=-1;
                    for(i=0;i<xarr.get_dim() && duplicate<0;i++){
                        if(fabs(xx-xarr.get_data(i))<xstep && fabs(yy-yarr.get_data(i))<ystep){
                            duplicate=i;
                        }
                    }
               
                }
                
                if(duplicate<0){
                    xarr.set(ct,xx);
                    yarr.set(ct,yy);
                    chiarr.set(ct,chival);
                    dexes.set(ct,ct);
                    ct++;
                }
                else{
                    nn=exp(-0.5*chiarr.get_data(duplicate))+exp(-0.5*chival);
                    chiarr.set(duplicate,-2.0*log(nn));
                
                }
                total+=exp(-0.5*chival);
                
                if(dexes.get_room()==ct){
                   //printf("adding room %d %d\n",dexes.get_dim(),ct);
                   dexes.add_room(1000);
                   xarr.add_room(1000);
                   yarr.add_room(1000);
                   chiarr.add_room(1000);
                   //printf("done %d\n",dexes.get_dim());
                   //exit(1);
                }
                
                //if(ct%1000==0)printf("    %d\n",ct);
                
            }//loop over yy
        }//loop overxx
    }
    
    if(dexes.get_dim()!=ct){
        printf("WARNING dexes %d but ct %d\n",dexes.get_dim(),ct);
        exit(1);
    }
    printf("time to sort %d -- %d\n",dexes.get_dim(),ct);
    
    array_1d<double> chisorted;
    sort_and_check(chiarr,chisorted,dexes);
    
    array_2d<double> scatter_data;
    scatter_data.set_cols(2);
    
    double sum=0.0;
    FILE *output;
    
    
    for(i=0;i<dexes.get_dim() && sum<lim*total;i++){
        j=dexes.get_data(i);
        sum+=exp(-0.5*chiarr.get_data(j));
        scatter_data.set(i,0,xarr.get_data(j));
        scatter_data.set(i,1,yarr.get_data(j));
        
    }

    
    //printf("moving on to boundary\n");
    
    array_1d<double> min,max;
    min.set(0,0.0);
    min.set(1,0.0);
    max.set(0,0.1*xwidth);
    max.set(1,0.1*ywidth);
    
    kd_tree scatter_tree(scatter_data,min,max);
    
    array_2d<double> boundary_pts;
    array_1d<int> neigh;
    array_1d<double> ddneigh;
  
    for(i=0;i<scatter_data.get_rows();i++){
        scatter_tree.nn_srch(*scatter_data(i),5,neigh,ddneigh);
        
        if(ddneigh.get_data(4)>1.001){
            boundary_pts.add_row(*scatter_data(i));
        }
    }
    
    kd_tree boundary_tree(boundary_pts,min,max);
    
    array_1d<double> origin,last_pt;
    int last_dex=0;
    
    origin.set(0,boundary_tree.get_pt(0,0));
    origin.set(1,boundary_tree.get_pt(0,1));
    last_pt.set(0,origin.get_data(0));
    last_pt.set(1,origin.get_data(1));
    
    output=fopen(filename,"w");

    while(boundary_tree.get_pts()>1){
        fprintf(output,"%e %e\n",last_pt.get_data(0),last_pt.get_data(1));
        
        boundary_tree.remove(last_dex);
        
        boundary_tree.nn_srch(last_pt,1,neigh,ddneigh);
        last_dex=neigh.get_data(0);
        last_pt.set(0,boundary_tree.get_pt(last_dex,0));
        last_pt.set(1,boundary_tree.get_pt(last_dex,1));
    }
    
    fprintf(output,"%e %e\n",last_pt.get_data(0),last_pt.get_data(1));
    fprintf(output,"%e %e\n",origin.get_data(0),origin.get_data(1));
    
    fclose(output);
    
    called=0;

}
