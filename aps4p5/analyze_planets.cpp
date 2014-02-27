#include <stdio.h>
#include <stdlib.h>
#include "containers.h"
#include "goto_tools.h"

main(int iargc, char *argv[]){

char inname[500],outname[500];

inname[0]=0;
sprintf(outname,"analysis_output.sav");

int i,j;

for(i=1;i<iargc;i++){
    if(argv[i][0]=='-'){
        switch(argv[i][1]){
            case 'i':
                i++;
                for(j=0;j<499 && argv[i][j]!=0;j++)inname[j]=argv[i][j];
                inname[j]=0;
            break;
            case 'o':
                i++;
                for(j=0;j<499 && argv[i][j]!=0;j++){
                    outname[j]=argv[i][j];
                }
                outname[j]=0;
            break;
        }
    }
}

if(inname[0]==0){
    printf("WARNING: you have to give me an input name\n");
    exit(1);
}

FILE *input;
char word[100];
input=fopen(inname,"r");
word[0]=0;
i=0;
while(compare_char(word,"ling")==0){
    fscanf(input,"%s",word);
    i++;
}
int nparams=i-5;

int npts=0,nplanets=nparams/3;

double nn,chi,chimin=exception;
array_1d<double> min_periods,periods;

FILE *output;
output=fopen(outname,"w");
fprintf(output,"# npts chi periods...\n");
while(fscanf(input,"%le",&nn)>0){
    npts++;
    periods.set(0,nn);
    for(i=1;i<3;i++){
        fscanf(input,"%le",&nn);
    }
    
    for(i=1;i<nplanets;i++){
        fscanf(input,"%le",&nn);
        periods.set(i,nn);
        for(j=1;j<3;j++)fscanf(input,"%le",&nn);
    }
    
    fscanf(input,"%le",&chi);
    for(i=0;i<3;i++)fscanf(input,"%le",&nn);
    
    if(chi<chimin){
        chimin=chi;
        
        fprintf(output,"%d %le ",npts,chimin);
        
        for(i=0;i<nplanets;i++){
            min_periods.set(i,periods.get_data(i));
            fprintf(output,"%le ",periods.get_data(i));
        }
        fprintf(output,"\n");
    }
}
fclose(output);
fclose(input);
printf("nplanets %d npts %d\n",nplanets,npts);


}
