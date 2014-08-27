/*
This code gives an example of how to use the aps_extractor class defined
in aps_extractor.h and aps_extractor.cpp to take the output from APS and
return both Frequentist confidence limits and Bayesian credible limits.

This file can be compiled using

make aps_extract

the executable aps_extract will then run this file
*/

#include "aps_extractor.h"

main(){

char sizeflag[100];
char outname[letters];

aps_extractor apsExtractorObj;

/*set the name of the file where the APS data is stored

note that the data file referred to here is not a part of the repository
*/
apsExtractorObj.set_filename("aps_output/WMAP/apsWMAP_unitSphere_GaussLearn_output.sav");


/*
set the delta chisquared so that 

chisquared_lim = chisquared_min + delta chisquared

12.6 is the delta chisquared corresponding to the 95% Frequentist 
confidence limit for a chisquared distribution with 6 degrees of freedom
*/
apsExtractorObj.set_delta_chi(12.6);

/*one could alternative set chisquared_lim by hand using*/
//apsExtractorObj.set_target(1283.3)

/*optionally set the maximum number of points from the data file to use*/
//apsExtractorObj.set_cutoff(47000);

/*select out only those points with chisquared < chisquared_lim and print them to the output file*/
apsExtractorObj.write_good_points("aps_processed/WMAP/frequentist/gaussLearn_good_points.sav");

/*plot chisquared_min discovered as a function of the number of samples*/
apsExtractorObj.plot_chimin("aps_chi_min_test.sav");

/*
The code below draws Frequentist and Bayesian limits in 2-dimensional sub-spaces of the
parameter space
*/
int i,j;
for(i=0;i<6;i++){
    for(j=i+1;j<6;j++){
        if(i!=3 && j!=3){
            sprintf(outname,"aps_processed/WMAP/bayesian/gaussLearn_%d_%d_bayes.sav",i,j);
            /*
            Below:
            the first argument is the name of the output file.
            the second and third arguments denote which dimensions to plot.
            the fourth argument is the confidence limit to plot
            */
            apsExtractorObj.draw_bayesian_bounds(outname,i,j,0.95);
            
            sprintf(outname,"aps_processed/WMAP/frequentist/gaussLearn_%d_%d_frequentist.sav",i,j);
            apsExtractorObj.write_good_points(outname,i,j);
        }
    }
}


}
