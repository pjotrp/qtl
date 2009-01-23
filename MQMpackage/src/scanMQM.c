/**********************************************************************
 * 
 * scanMQM.c
 *
 * copyright (c) 2009
 *
 * last modified Feb 2009
 * first written Feb 2009
 *
 * C functions for the R/qtl package
 * Contains: R_scanMQM, scanMQM
 *
 **********************************************************************/
using namespace std;

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "MQMdata.h"
#include "MQMsupport.h"

double neglect=1000; // eliminate unlikely genotype configurations
int maxNaug=10000; // maximum size of augmented dataset
int imaxNaug=1000; // maximum size of augmented data for individual i
int maxdimX=50; // maximum size of design matrix in regression
int em=1000; // maximum number of em iterations
double alfa=0.02; // alfa used in selection procedure
double windowsize=25.0; // used in mapQTL procedure
double stepsize=5; // size of steps when moving QTL along chromosomes (for output)
double stepmin=-20; // start moving QTL at position stepmin cM (for output)
double stepmax=220; // move QTL up to stepmax (for output)
long *idum; // for monte carlo simulation or permutation

vector mapdistance, r;
matrix XtWX;
cmatrix Xt;
vector XtWY;
ivector chr;
matrix Frun;
vector informationcontent;
int Nfam, Nrun=0;
int run=-1;
char REMLorML;
char fitQTL='n';
char dominance='n';
char perm_simu;

char ok, defset;

void OK()
{    ok='0';
     if (defset=='n') {cout << "OK (y/n)?"; cin >> ok; if (ok=='n') exit(1); }
}

double Lnormal(double residual, double variance)
{      double Likelihood;
       Likelihood=exp(-pow(residual/sqrt(variance),2.0)/2.0 - log(sqrt(2.0*acos(-1.0)*variance)));
       return Likelihood;
}

double absdouble(double x)
{      double z; z= (x<0 ? -x : x); return z;}

int mod(int a, int b)
{      int c;
       c= a/b;
       return a-b*c;
}

/**********************************************************************
 * 
 * R_scanMQM
 * 
 **********************************************************************/

void R_scanMQM(){

} /* end of function R_scanMQM */


/**********************************************************************
 * 
 * scanMQM
 *
 * 
 **********************************************************************/

void scanMQM(){

} /* end of function scanMQM */

int main(){
 return (1);   
}

/* end of scanMQM.c */
