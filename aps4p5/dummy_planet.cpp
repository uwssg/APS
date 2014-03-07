#include "dummy_planet.h"

dummy_planet::dummy_planet() : chisquared(15){

    c1.set(0,5191.0); 
    c1.set(1,0.015);
    c1.set(2,-0.4783912);
    
    c1.set(3,259.7); 
    c1.set(4,0.40);
    c1.set(5,-0.1645573);
    
    c1.set(6,44.349); 
    c1.set(7,0.09);
    c1.set(8,0.2048294);
    
    c1.set(9,14.65164); 
    c1.set(10,0.014); 
    c1.set(11,-0.2729630);
    
    c1.set(12,2.81705); 
    c1.set(13,0.066);
    c1.set(14,-0.2742881);
    
     c2.set(0,5205.0); 
     c2.set(1,0.024); 
     c2.set(2,-0.5694518);
     
     c2.set(3,259.8); 
     c2.set(4,0.25);
     c2.set(5,-0.1965394);
     
     c2.set(6,44.342); 
     c2.set(7,0.05);
     c2.set(8,0.326152);
     
     c2.set(9,14.65160); 
     c2.set(10,0.014); 
     c2.set(11,-0.2397775);
     
     c2.set(12,0.736539); 
     c2.set(13,0.17);
     c2.set(14,-0.08269329);
     
     int i;
     for(i=0;i<5;i++){
         ll.set(i*3,1.0);
         ll.set(i*3+1,0.5);
         ll.set(i*3+2,0.5);
     }
     
     ll.set(12,0.001);
     ll.set(0,100.0);
     ll.set(9,0.1);
     

}

dummy_planet::~dummy_planet(){}

double dummy_planet::operator()(array_1d<double> &pp) const{
    double d1=0.0,d2=0.0;
    called++;
    
    int i,j;
    for(i=0;i<15;i++){
        d1+=power((pp.get_data(i)-c1.get_data(i))/ll.get_data(i),2);
        d2+=power((pp.get_data(i)-c2.get_data(i))/ll.get_data(i),2);
    }
    
    double dmin,chimin,nn;
    array_1d<double> dd,center;
    
    if(d1<d2){
        dmin=d1;
        chimin=800.0;
        for(i=0;i<15;i++)center.set(i,c1.get_data(i));
       
    }
    else{
        dmin=d2;
        chimin=520.0;
        for(i=0;i<15;i++)center.set(i,c2.get_data(i));
    }
    
    double ans;
    ans=chimin;
    for(i=0;i<5;i++){
        nn=0.0;
        for(j=0;j<3;j++){
            nn+=power((center.get_data(i*3+j)-pp.get_data(i*3+j))/ll.get_data(i*3+j),2);
        }
        dd.set(i,nn);
        if(nn<0.0){
            printf("WARNING nn %e\n",nn);
        }
    }
    
    nn=100.0*dd.get_data(4);
    if(nn>2000.0)nn=2000.0;
    
    ans+=nn;
    ans+=5.0*sin(2.0*pi*dd.get_data(3)/0.1);
    ans+=10.0*dd.get_data(2);
    ans+=10.0*dd.get_data(1);
    ans+=0.1*dd.get_data(0);
  
    return ans;

}

