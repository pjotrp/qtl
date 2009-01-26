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


#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Utils.h>
#include <R_ext/Lapack.h>
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

void reorg_geno(int n_ind, int n_pos, int *geno, int ***Geno)
{
  int i;

  *Geno = (int **)R_alloc(n_ind, sizeof(int *));

  (*Geno)[0] = geno;
  for(i=1; i< n_ind; i++) 
    (*Geno)[i] = (*Geno)[i-1] + n_pos;

}

void reorg_pheno(int n_ind, int n_mar, double *pheno, double ***Pheno)
{
  int i;

  *Pheno = (double **)R_alloc(n_mar, sizeof(double *));

  (*Pheno)[0] = pheno;
  for(i=1; i< n_mar; i++) 
    (*Pheno)[i] = (*Pheno)[i-1] + n_ind;
}

/**********************************************************************
 * 
 * scanMQM
 *
 * 
 **********************************************************************/

void scanMQM(int Nind, int Nmark,int Npheno, int Nfam,int **Geno, double **Pheno){
   ivector f1genotype;
   f1genotype = newivector(Nmark);

   //Rprintf("Printing Genotype matrix\n");
   for(int i=0; i< Nmark; i++){
     for(int j=0; j< Nind; j++){ 
       Rprintf("%d ",Geno[i][j]);
       f1genotype[j] = 12;         
     }
     Rprintf("\n");
   }

   //Rprintf("Printing Phenotype matrix\n");
   for(int i=0; i< Npheno; i++){
     for(int j=0; j< Nind; j++){
  //     Rprintf("%f ",Pheno[i][j]);         
     }
    // Rprintf("\n");
   }
   Rprintf("We got all the needed information, so lets start with the MQM\n");
   idum = (long *)R_alloc(1, sizeof(long *));
   idum[0]=-1;
   cvector cofactor;
   char real_simu;
   real_simu = '0';
   cmatrix markername;
   
   markername= newcmatrix(Nmark,20);
   cofactor= newcvector(Nmark);   
   mapdistance= newvector(Nmark);

   chr= newivector(Nmark);        

   if (real_simu=='1'){
     simuF2(Nind, Nmark, cofactor, Geno, Pheno[0]);
   }
   analyseF2(Nind, Nmark, cofactor, Geno, Pheno[0], f1genotype);
   return;
}  /* end of function scanMQM */



/**********************************************************************
 * 
 * R_scanMQM
 * 
 **********************************************************************/

void R_scanMQM(int *Nind,int *Nmark,int *Npheno, int *Nfam,int *geno, double *pheno){
   int **Geno;
   double **Pheno;   
   reorg_geno(*Nind,*Nmark,geno,&Geno);
   reorg_pheno(*Nind,*Npheno,pheno,&Pheno);
   //Done with reorganising lets start executing the main loop
   
   scanMQM(*Nind,*Nmark,*Npheno, *Nfam,Geno,Pheno);
} /* end of function R_scanMQM */


