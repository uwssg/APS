#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "chisq.h"


void chisquared::set_characteristic_error(double nn){
    death_knell("used meaningless set characteristic error");
}

void chisquared::activate_small_m(){
    death_knell("used nonesnese small m");
}

void chisquared::allot_arrays(){
    int ix,iy;
    
    bases=new double*[dim];
    for(ix=0;ix<dim;ix++)bases[ix]=new double[dim];
    widths=new double*[ncenters];
    for(ix=0;ix<ncenters;ix++)widths[ix]=new double[dim];
    centers=new double*[ncenters];
    for(ix=0;ix<ncenters;ix++)centers[ix]=new double[dim];
    nboundary=new int[dim*dim];
    boundary_room=new int[dim*dim];
    
    for(ix=0;ix<dim*dim;ix++)nboundary[ix]=0;
    for(ix=0;ix<dim*dim;ix++)boundary_room[ix]=3;
    
    boundary=new double**[dim*dim];
    for(ix=0;ix<dim*dim;ix++){
        boundary[ix]=new double*[boundary_room[ix]];
        for(iy=0;iy<boundary_room[ix];iy++){
	    boundary[ix][iy]=new double[3];
	}
        
    }
    
    mins=new double[dim];
    maxs=new double[dim];
    for(ix=0;ix<dim;ix++){
        mins[ix]=1.0e30;
	maxs[ix]=-1.0e30;
    } 
    
}

void chisquared::reset_boundary(){
    int i;
    for(i=0;i<dim*dim;i++)nboundary[i]=0;
}

void chisquared::make_bases(int seed){
    if(bases==NULL){
        death_knell("you called make_bases before allotting bases\n");
    }
    
    if(widths==NULL){
        death_knell("you called make_bases before allotting widths\n");
    }
    
    if(centers==NULL){
        death_knell("you called make_bases before allotting centers\n");
    }
    
    if(seed<0)seed=int(time(NULL));
    
    dice=new Ran(seed);
    
    double nn;
    int i,j,ii,jj,goon;
    
    for(ii=0;ii<dim;ii++){
        goon=1;
	while(goon==1){
            goon=0;
	    for(i=0;i<dim;i++)bases[ii][i]=dice->doub()-0.5;
	    for(jj=0;jj<ii;jj++){
	        nn=0.0;
		for(i=0;i<dim;i++)nn+=bases[ii][i]*bases[jj][i];
		for(i=0;i<dim;i++)bases[ii][i]-=nn*bases[jj][i];
	    }
	    
	    nn=0.0;
	    for(i=0;i<dim;i++){
	        nn+=bases[ii][i]*bases[ii][i];
	    }
	    if(nn<1.0e-20)goon=1;
	    nn=sqrt(nn);
	    for(i=0;i<dim;i++){
	        bases[ii][i]=bases[ii][i]/nn;
	    }
	}
    }
    
    double normerr,ortherr;
    for(ii=0;ii<dim;ii++){
        nn=0.0;
	for(i=0;i<dim;i++)nn+=bases[ii][i]*bases[ii][i];
	nn=fabs(1.0-nn);
	if(ii==0 || nn>normerr)normerr=nn;
	
	for(jj=ii+1;jj<dim;jj++){
	   nn=0.0;
	   for(i=0;i<dim;i++)nn+=bases[ii][i]*bases[jj][i];
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
    double *trial_center,*trial_pt;
    int acceptable;
    trial_center=new double[dim];
    trial_pt=new double[dim];
    

    
    goon=1;
    while(goon==1){
        
	for(ii=0;ii<ncenters;ii++){
            for(i=0;i<dim;i++)centers[ii][i]=1.0e30;
	    for(i=0;i<dim;i++)widths[ii][i]=1.0e-5;
        }
	
	goon=0;
        for(ii=0;ii<ncenters;ii++){
	    
	    acceptable=1;
	    
	    for(i=0;i<dim;i++)trial_center[i]=0.0;
	    
	    if(ii>0){
	        rr=normal_deviate(dice,40.0,20.0);
	        theta=dice->doub()*2.0*pi;
	    
	        trial_center[0]=centers[0][0]+rr*cos(theta);
	        trial_center[1]=centers[0][1]+rr*sin(theta);
	    } 
	     
	    for(i=0;i<dim;i++){
	        if(i<2 && ii>0)trial_center[i]+=normal_deviate(dice,3.0,3.0);
		else trial_center[i]+=normal_deviate(dice,30.0,15.0);
	    }
	    
	    for(i=0;i<dim;i++){
	        trial_pt[i]=0.0;
	        for(j=0;j<dim;j++)trial_pt[i]+=trial_center[j]*bases[j][i];
	    }
	    
	    rr=(*this)(trial_pt);
	    
	    if(rr<100.0)acceptable=0;
	    
	    if(acceptable==1){
	        for(i=0;i<dim;i++){
		    centers[ii][i]=trial_center[i];
	            //widths[ii][i]=dice->doub()*1.0+0.1;
		    widths[ii][i]=fabs(normal_deviate(dice,1.0,0.5))+0.1;
	        }
	    }
	    else ii--;
	
        }
	
	
	
	
	
	for(i=0;i<dim && goon==0;i++){
	    for(ii=0;ii<ncenters && goon==0;ii++){
	        for(jj=ii+1;jj<ncenters && goon==0;jj++){
		    nn=fabs(centers[ii][i]-centers[jj][i]);
		    if(nn<2.0*widths[ii][i] || nn<2.0*widths[jj][i])acceptable=0;
		}
	    }
	}
	
    }
    
    delete [] trial_center;
    delete [] trial_pt;
    printf("set centers and widths\n");
}

void chisquared::add_to_boundary(double *alpha,int ix, int iy,double chitest){
    
    double **buffer,rr,*pt;
    int ipt,room;
    int i,j;
    //printf("adding %e\n",chitest);
    if(ix>iy){
        i=ix;
	ix=iy;
	iy=i;
    }
    
    ipt=ix*dim+iy;
    room=boundary_room[ipt];
    
    if(nboundary[ipt]>=room){
        buffer=new double*[room];
	for(i=0;i<room;i++){
	    buffer[i]=new double[3];
	    for(j=0;j<3;j++)buffer[i][j]=boundary[ipt][i][j];
	    delete [] boundary[ipt][i];
	}
        delete [] boundary[ipt];
	boundary_room[ipt]+=100;
	boundary[ipt]=new double*[boundary_room[ipt]];
	for(i=0;i<boundary_room[ipt];i++)boundary[ipt][i]=new double[3];
	
	for(i=0;i<room;i++){
	    for(j=0;j<3;j++)boundary[ipt][i][j]=buffer[i][j];
	    delete [] buffer[i];
	}
	delete [] buffer;
    }
    
    boundary[ipt][nboundary[ipt]][0]=alpha[ix];
    boundary[ipt][nboundary[ipt]][1]=alpha[iy];
    boundary[ipt][nboundary[ipt]][2]=chitest;
    
    nboundary[ipt]++;
    
    rr=0.0;
    for(i=0;i<dim;i++)rr+=alpha[i]*alpha[i];
    rr=sqrt(rr);
    if(rr>rr_max)rr_max=rr;
    
    pt=new double[dim];
    for(i=0;i<dim;i++)pt[i]=0.0;
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
	    pt[i]+=alpha[j]*bases[j][i];
	}
    }
    
    double nn=0.0;
    for(i=0;i<dim;i++)nn+=pt[i]*pt[i];
    nn=sqrt(nn);
    if(fabs(nn-rr)>1.0e-4){
        printf("WARNING nn %e rr %e\n",nn,rr);
    }
    
    for(i=0;i<dim;i++){
        if(pt[i]<mins[i])mins[i]=pt[i];
	if(pt[i]>maxs[i])maxs[i]=pt[i];
    }
    
    delete [] pt;
    
    
}

void chisquared::print_mins_maxs(){
    int i;
    double nn;
    nn=0.0;
    printf("mins and maxs\n");
    for(i=0;i<dim;i++){
        printf("p%d -- %e %e -- %e\n",i,mins[i],maxs[i],maxs[i]-mins[i]);
	nn+=power(maxs[i]-mins[i],2);
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
    
    return nboundary[ix*dim+iy];
}

double chisquared::get_width(int ic, int ix){
    if(widths!=NULL && ic<ncenters && ix<dim){
        return widths[ic][ix];
    }
    else{
        return 1.0e30;
    }
}

double chisquared::get_center(int ic, int ix){
    if(centers!=NULL && ic<ncenters && ix<dim){
        return centers[ic][ix];
    }
    else{
        return 1.0e30;
    }
}

double chisquared::get_real_center(int ic, int ix){
    if(centers==NULL || ic>=ncenters || ix>=dim){
        return 1.e30;
    }
    
    int i;
    double ans=0.0;
    for(i=0;i<dim;i++){
        ans+=centers[ic][i]*bases[i][ix];
    }
    return ans;
    
}

double chisquared::get_boundary(int ix, int iy, int ipt, int idim){

    int i;
    if(ix>iy){
        i=ix;
	ix=iy;
	iy=i;
    }

    if(idim>=3 || ipt>=nboundary[ix*dim+iy])return 1.0e30;
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
    
    bases=NULL;
    widths=NULL;
    centers=NULL;
    boundary=NULL;
    nboundary=NULL;
    boundary_room=NULL;
    maxs=NULL;
    mins=NULL;
    
    rr_max=-1.0;
    called=0;
    time_spent=0.0;
    
    allot_arrays();
};

chisquared::chisquared(int id, int ic){
    dim=id;
    ncenters=ic;
    
    bases=NULL;
    widths=NULL;
    centers=NULL;
    boundary=NULL;
    nboundary=NULL;
    boundary_room=NULL;
    maxs=NULL;
    mins=NULL;
    
    rr_max=-1.0;
    called=0;
    time_spent=0.0;
    
    allot_arrays();
};



chisquared::~chisquared(){
    int i,ix,iy;
    if(bases!=NULL){
        for(i=0;i<dim;i++)delete [] bases[i];
	delete [] bases;
    }
    
    if(widths!=NULL){
        for(i=0;i<ncenters;i++){
	    delete [] widths[i];
	}
	delete [] widths;
    }
    
    if(boundary!=NULL){
        for(ix=0;ix<dim;ix++){
	    for(iy=ix+1;iy<dim;iy++){
	        
	            for(i=0;i<boundary_room[ix*dim+iy];i++){
		        delete [] boundary[ix*dim+iy][i];
		    }
		    delete [] boundary[ix*dim+iy];
		
	    }
	}
	delete [] boundary;
	delete [] nboundary;
    }
    
    if(boundary_room!=NULL){
         delete [] boundary_room;
    }
    
    if(centers!=NULL){
        for(i=0;i<ncenters;i++)delete [] centers[i];
	delete [] centers;
    }
    
    if(dice!=NULL){
        delete dice;
    }
    
    if(maxs!=NULL)delete [] maxs;
    if(mins!=NULL)delete [] mins;
}

void chisquared::set_max_min(
    int dd, double *mn, double *mx){
    
    
    
    int i;
    double nn;
    
    
    if(mins==NULL)mins=new double[dim];
    if(maxs==NULL)maxs=new double[dim];
    
    for(i=0;i<dim;i++){
        mins[i]=mn[i];
	maxs[i]=mx[i];
	if(mins[i]>maxs[i]){
	    nn=mins[i];
	    mins[i]=maxs[i];
	    maxs[i]=nn;
	}
    }    
}

int chisquared::get_called(){
    return called;
}

double chisquared::get_time(){
    return time_spent;
}

double chisquared::operator()(double *v)const{
    death_knell("meaningless operator");
}

void chisquared::build_boundary(double rr){
    death_knell("meaningless build_boundary");
}

void chisquared::get_basis(int ix, double *v){
    int i;
    if(bases!=NULL && ix<dim){
        for(i=0;i<dim;i++)v[i]=bases[ix][i];
    }
    else{
        printf("WARNING called get_basis with %d %d\n",ix,dim);
    }
}

double chisquared::project_to_basis(int ix, double *vv) const{
    int i;
    double nn=1.0e30;
    if(bases!=NULL && ix<dim){
        nn=0.0;
	for(i=0;i<dim;i++)nn+=vv[i]*bases[ix][i];
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

linear_ellipses::~linear_ellipses(){}

linear_ellipses::linear_ellipses() : chisquared(22){make_bases(17);}

linear_ellipses::linear_ellipses(int id) : chisquared(id){make_bases(17);}

linear_ellipses::linear_ellipses(int id, int ic) : chisquared(id,ic){
    make_bases(17);
}

double s_curve::operator()(double *in_pt) const{
    
    if(dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    called++;
    
    double before=double(time(NULL));
    double *dd,*pt;
    dd=new double[ncenters];
    pt=new double[dim];
    
    double theta,xth,yth,dth,dthmin;
    int i;
    
    for(i=0;i<dim;i++)pt[i]=project_to_basis(i,in_pt);
    
    //for(i=0;i<dim;i++)printf("%e %e\n",pt[i],centers[0][i]);
    
    
    dd[0]=0.0;
    for(i=2;i<dim;i++)dd[0]+=power((pt[i]-centers[0][i])/widths[0][i],2);
    
    dthmin=-1.0;
    for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
        xth=trig_factor*sin(theta)+centers[0][0];
	yth=trig_factor*theta*(cos(theta)-1.0)/fabs(theta)+centers[0][1];
	dth=power((xth-pt[0])/widths[0][0],2)+power((yth-pt[1])/widths[0][1],2);
	
	if(dthmin<0.0 || dth<dthmin)dthmin=dth;
    }
    dd[0]+=dthmin;
    
    int ii;
    for(ii=1;ii<ncenters;ii++){
        dd[ii]=0.0;
	for(i=0;i<dim;i++){
	    dd[ii]+=power((pt[i]-centers[ii][i])/widths[ii][i],2);
	}
    }
    
    double ddmin;
    for(ii=0;ii<ncenters;ii++){
        if(ii==0 || dd[ii]<ddmin)ddmin=dd[ii];
    }
    
    delete [] dd;
    delete [] pt;
    
    time_spent+=double(time(NULL))-before;
    
    return ddmin;
}

void s_curve::build_boundary(double br){
    if(dice==NULL){
        death_knell("you called build_boundary before making bases");
    }

    reset_boundary();

    int ix,iy,ic,ir,i,j;
    
    
    double theta,dxdth,dydth,x0,y0;
    double grad[2],norm,*alpha,*pt,dfabsdth,ds,chitest;
    
    alpha=new double[dim];
    pt=new double[dim];
    
    
    ///////////////below we will do ix=0, iy=1 for the 0th center (the S curve itself)
    ix=0;
    iy=1;
    
  
    
    if(widths[0][0]<widths[0][1])ds=0.1*widths[0][0];
    else ds=0.1*widths[0][1];
    
    for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
        if(theta<0.0)dfabsdth=-1.0;
        else dfabsdth=1.0;
    
        x0=centers[0][0]+trig_factor*sin(theta);
        y0=centers[0][1]+trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
    
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
	
	    for(i=0;i<dim;i++)alpha[i]=centers[0][i];
	

	
	    alpha[0]=x0+grad[0]*ds/norm;
	    alpha[1]=y0+grad[1]*ds/norm;
	    
	    alpha[0]=x0;
	    alpha[1]=y0;
	    for(i=0;i<dim;i++){
	        pt[i]=0.0;
	        for(j=0;j<dim;j++)pt[i]+=alpha[j]*bases[j][i];
	    }
	    chitest=(*this)(pt);
	    
	    if(chitest>br)death_knell("started outside the pale\n");
	    
	    while(chitest<br){
	   
	        alpha[0]+=ds*grad[0]/norm;
	        alpha[1]+=ds*grad[1]/norm;
	    
	        for(i=0;i<dim;i++){
	            pt[i]=0.0;
		    for(j=0;j<dim;j++)pt[i]+=alpha[j]*bases[j][i];
	        }
	        chitest=(*this)(pt);
		
		
	    }
	    
	    if(fabs(chitest-br)<5.0){
	        add_to_boundary(alpha,ix,iy,chitest);
	    }
	    /*for(i=0;i<2;i++){
	        printf("%e ",alpha[i]);
	    }
	    printf("%e\n",chitest);*/
	    
		
        }
    }
    
    for(i=0;i<dim;i++)alpha[i]=centers[0][i];
    
    double th,s2x,s2y,aa,rr,thmin,thmax;
    for(theta=-1.0*pi;theta<=1.5*pi;theta+=2.0*pi){
        for(i=0;i<dim;i++)alpha[i]=centers[0][i];
    
        x0=centers[0][0]+trig_factor*sin(theta);
        y0=centers[0][1]+trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
    
        s2x=widths[0][0]*widths[0][0];
        s2y=widths[0][1]*widths[0][1];
        
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
	
	    alpha[0]=x0+rr*cos(th);
	    alpha[1]=y0+rr*sin(th);
	
	    for(i=0;i<dim;i++){
	        pt[i]=0.0;
	        for(j=0;j<dim;j++)pt[i]+=alpha[j]*bases[j][i];
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
        for(i=0;i<dim;i++)alpha[i]=centers[0][i];
        for(theta=-1.0*pi;theta<=1.0*pi;theta+=0.01){
           alpha[0]=centers[0][0]+trig_factor*sin(theta);
           alpha[1]=centers[0][1]
                   +trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
       
       
           if(ir==0)alpha[iy]=centers[0][iy]+sqrt(br)*widths[0][iy];
           else alpha[iy]=centers[0][iy]-sqrt(br)*widths[0][iy];
       
           for(i=0;i<dim;i++){
               pt[i]=0;
	       for(j=0;j<dim;j++)pt[i]+=alpha[j]*bases[j][i];
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


        s2x=widths[0][ix]*widths[0][ix];
        for(iy=2;iy<dim;iy++){
            s2y=widths[0][iy]*widths[0][iy];
            for(i=0;i<dim;i++)alpha[i]=centers[0][i];
            for(theta=thetamin;theta<thetamax;theta+=dtheta){
	     
	         if(ix==0){
	             xth=centers[0][0]+trig_factor*sin(theta);
	             alpha[1]=centers[0][1]
		        +trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
	         }
	         else{
	             xth=centers[0][1]+trig_factor*theta*(cos(theta)-1.0)/fabs(theta);
	             alpha[0]=centers[0][0]+trig_factor*sin(theta);
	         }
	     
	     
	         for(th=0.0;th<2.0*pi;th+=0.01){
	             aa=(cos(th)*cos(th)/s2x+sin(th)*sin(th)/s2y);
		     rr=br/aa;
		     rr=sqrt(rr);
		 
		     alpha[iy]=centers[0][iy]+rr*sin(th);
		     alpha[ix]=xth+rr*cos(th);
	         
		 
		     for(i=0;i<dim;i++){
		         pt[i]=0.0;
		         for(j=0;j<dim;j++){
		             pt[i]+=alpha[j]*bases[j][i];
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
		    
		    for(i=0;i<dim;i++)alpha[i]=centers[ic][i];
		    s2x=power(widths[ic][ix],2);
		    s2y=power(widths[ic][iy],2);
		    for(theta=0.0;theta<=2.0*pi;theta+=0.01){
		        aa=power(cos(theta),2)/s2x+power(sin(theta),2)/s2y;
			rr=br/aa;
			rr=sqrt(rr);
			alpha[ix]=centers[ic][ix]+rr*cos(theta);
			alpha[iy]=centers[ic][iy]+rr*sin(theta);
		        
			for(i=0;i<dim;i++){
			    pt[i]=0.0;
			    for(j=0;j<dim;j++)pt[i]+=alpha[j]*bases[j][i];
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
    
    delete [] pt;
    delete [] alpha;
    
}

double ellipses::operator()(double *in_pt) const{
    
    if(dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    double before=double(time(NULL));
    
    int ii,ix;
    double dd,ddmin,nn;
    
    called++;
    
    for(ii=0;ii<ncenters;ii++){
        dd=0.0;
        for(ix=0;ix<dim;ix++){
	    nn=project_to_basis(ix,in_pt);
	    dd+=power((centers[ii][ix]-nn)/widths[ii][ix],2);
	    //printf("%e\n",centers[ii][ix]-nn);
	}
	if(ii==0 || dd<ddmin)ddmin=dd;
    }
    
    time_spent+=double(time(NULL))-before;
    
    return ddmin;
    
}

void ellipses::build_boundary(double br){

    if(dice==NULL){
         death_knell("you called build_boundary before making bases");
    } 
    
    reset_boundary();
    
    int ix,iy,ic,i,j;
    double theta,rr,aa,*alpha,*pt,chitest;
    
    alpha=new double[dim];
    pt=new double[dim];
    
    for(ic=0;ic<ncenters;ic++){
        for(ix=0;ix<dim;ix++){
            for(iy=ix+1;iy<dim;iy++){
	        for(i=0;i<dim;i++)alpha[i]=centers[ic][i];
	        for(theta=0.0;theta<=2.0*pi;theta+=0.01){
	            aa=power(cos(theta)/widths[ic][ix],2)+power(sin(theta)/widths[ic][iy],2);
		    rr=br/aa;
		    rr=sqrt(rr);
		    alpha[ix]=centers[ic][ix]+rr*cos(theta);
		    alpha[iy]=centers[ic][iy]+rr*sin(theta);
		    
		    
		    for(i=0;i<dim;i++){
		        pt[i]=0.0;
			for(j=0;j<dim;j++)pt[i]+=alpha[j]*bases[j][i];
		    }
		    
		    chitest=(*this)(pt);
		    if(fabs(chitest-br)<5.0){
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
    
    delete [] alpha;
    delete [] pt;
    
}

double linear_ellipses::operator()(double *in_pt) const{

    if(dice==NULL){
         death_knell("you called operator before making bases");
    } 
    
    double before=double(time(NULL));
    
    int ii,ix;
    double dd,ddmin,nn;
    
    called++;
    
    for(ii=0;ii<ncenters;ii++){
        dd=0.0;
	for(ix=0;ix<dim;ix++){
	    nn=project_to_basis(ix,in_pt);
	    dd+=power((centers[ii][ix]-nn)/widths[ii][ix],2);
	}
	
	if(ii==0 || dd<ddmin)ddmin=dd;
    }
    
    time_spent+=double(time(NULL))-before;
    
    return sqrt(ddmin);
}

void linear_ellipses::build_boundary(double br){
    if(dice==NULL){
         death_knell("you called build_boundary before making bases");
    } 
    
    reset_boundary();
    
    int ic,ix,iy,i,j;
    double *pt,*alpha;
    double theta,rr,aa,chitest;
    
    pt=new double[dim];
    alpha=new double[dim];
    
    for(ic=0;ic<ncenters;ic++){
        for(ix=0;ix<dim;ix++){
	    for(iy=ix+1;iy<dim;iy++){
	        for(i=0;i<dim;i++)alpha[i]=centers[ic][i];
		for(theta=0.0;theta<=2.0*pi;theta+=0.01){
		    aa=power(cos(theta)/widths[ic][ix],2)+power(sin(theta)/widths[ic][iy],2);
		    rr=br*br/aa;
		    rr=sqrt(rr);
		    
		    alpha[ix]=centers[ic][ix]+rr*cos(theta);
		    alpha[iy]=centers[ic][iy]+rr*sin(theta);
		    for(i=0;i<dim;i++){
		        pt[i]=0.0;
			for(j=0;j<dim;j++)pt[i]+=alpha[j]*bases[j][i];
		    }
		    chitest=(*this)(pt);
		    if(fabs(chitest-br)<5.0){
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
    
    delete [] pt;
    delete [] alpha;
    
}

udder_likelihood::udder_likelihood() : chisquared(6){
    foundn3=-1;
    foundp3=-1;
}

udder_likelihood::~udder_likelihood(){}

double udder_likelihood::operator()(double *v) const{

     double d1,d2,base1,base2,amp1,amp2,base,chisquared;
   
     int i;
    
    double before=double(time(NULL));
    
    called++;
    d1=power(v[0]-3.0,2);
    d2=power(v[0]+3.0,2);
    
    for(i=1;i<6;i++){
        d1+=power(v[i],2);
	d2+=power(v[i],2);
    }
   
    
     chisquared=1300.0+0.5*d1+0.5*d2-153.0*exp(-2.0*d1)-100.0*exp(-1.0*d2);
    
     time_spent+=double(time(NULL))-before;
    
     if(chisquared<=1280.669){
         if(d1<d2 && foundp3<0)foundp3=called;
	 if(d2<d1 && foundn3<0)foundn3=called;
     }
    
     return chisquared;

}

int udder_likelihood::get_n3(){
    return foundn3;
}

int udder_likelihood::get_p3(){
    return foundp3;
}

int udder_likelihood::get_type(){
    return LK_TYPE_UDDER;
}

#ifdef _WMAP7_

wmap_likelihood::wmap_likelihood() : chisquared(6){}

wmap_likelihood::~wmap_likelihood(){}

int wmap_likelihood::get_type(){
    return LK_TYPE_WMAP;
}

double wmap_likelihood::operator()(double *v) const{
  
  //this is the function that calls the likelihood function
  //to evaluate chi squared
  
  //*v is the point in parameter space you have chosen to evaluate
  
  //as written, it calls the CAMB likelihood function for the CMB anisotropy 
  //spectrum
  
  int i,start,k;
  double params[14],chisquared,omm,base1,base2,amp1,amp2,base;
  double d1,d2,sncc,cc1,cc2,dcc,ccmaxc;
  double cltt[3000],clte[3000],clee[3000],clbb[3000],el[3000];
  
  double *dir;
 
  double before=double(time(NULL));
  called++;
  
  FILE *output;

  for(i=0;i<dim;i++){
      if(v[i]<mins[i] || v[i]>maxs[i]){
          time_spent+=double(time(NULL))-before;
      
          return 2.0*exception;
      }
  }
  
  while((v[0]+v[1])/(v[2]*v[2])>1.0){
    v[0]=0.9*v[0];
    v[1]=0.9*v[1];
    v[2]=1.1*v[2];
    //in the event that total omega_matter>1
    
  }
  
  for(i=0;i<6;i++)params[i]=v[i];
 
  for(i=0;i<3000;i++)el[i]=0;
  params[2]=100.0*params[2];
  params[5]=exp(params[5])*1.0e-10;
  

  camb_wrap_(params,el,cltt,clte,clee,clbb); //a function to call CAMB
  
  for(start=0;el[start]<1.0;start++);

  //printf("cltt start %e %d\n",cltt[start],start);
  if(cltt[start]>=-10.0){
  wmaplikeness_(&cltt[start],&clte[start],&clee[start],\
  &clbb[start],&chisquared); //a function to call the WMAP likelihood code
  }
  else chisquared=2.0*exception;

 //printf("done with likelihood\n");

  if(chisquared<0.01)chisquared=2.0*exception; 
  			//in case the model was so pathological that the
			//likelihood code crashed and returned
			//chisquared=0  (this has been known to happen)
  

 
 
 /////////////////
 
  time_spent+=double(time(NULL))-before;
 
  //printf("got chisquared %e\n",chisquared);
  return chisquared;
  
}

#endif
