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
#include "scanMQM.h"
#include "MQMdata.h"

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

  *Geno = (int **)R_alloc(n_pos, sizeof(int *));

  (*Geno)[0] = geno;
  for(i=1; i< n_pos; i++) 
    (*Geno)[i] = (*Geno)[i-1] + n_ind;

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


/**********************************************************************
 * 
 * scanMQM
 *
 * 
 **********************************************************************/

void scanMQM(int Nind, int Nmark,int Npheno, int Nfam,int **Geno, double **Pheno){
   //ivector f1genotype;
  // f1genotype = newivector(Nmark);

   Rprintf("Printing Genotype matrix\n");
   for(int i=0; i< Nind; i++){
     for(int j=0; j< Nmark; j++){ 
       Rprintf("%d ",Geno[i][j]);
     //  f1genotype[j] = 12;         
     }
     Rprintf("\n");
   }

   Rprintf("Printing Phenotype matrix\n");
   for(int i=0; i< Npheno; i++){
     for(int j=0; j< Nind; j++){
       Rprintf("%f ",Pheno[i][j]);         
     }
     Rprintf("\n");
   }
     Rprintf("We got all the needed information, so lets start with the MQM\n"); 
}  /* end of function scanMQM */

int main(){
    int Nmark, saveNmark, Nind, Nsteps, saveNind, fam=0, f2genotype, i, ii, j;
  //  idum= new long[1];
  //  idum[0]=-1;

    vector y; //cvariance
    cvector cross, datafilef1marker, datafilef2marker, datafilephenotype, skipstring, cofactor;
    ivector f1genotype;
    cmatrix markername, marker;
    char ch, real_simu;
    double readf;

    run= -1;
    real_simu = '0';
    Nind=100;
    Nmark=100;

   // markername= newcmatrix(Nmark,20);
  //  cofactor= newcvector(Nmark);       // iscofactor?
  //  mapdistance= newvector(Nmark);
  //  chr= newivector(Nmark);            // chromosome number
  //  f1genotype= newivector(Nmark);     // parent genotype
  
  //  if (real_simu=='1') simuF2(Nind, Nmark, cofactor, marker, y);  // simulate cvariance
  //  analyseF2(Nind, Nmark, cofactor, marker, y, f1genotype);
  //  cout << "Analysis of data of family " << fam << " finished" << endl;
//
  // ---- Write output file
 // ofstream fff("mqm_out.txt", ios::out | ios::app);
 // fff << endl;
 // double moveQTL= stepmin;
//  int chrnumber=1;
 // cout << "-1- " << Nsteps << " " << Nrun << " " << Nfam << endl;
 
  // chr pos Frun    information 
  // 1  -20  97.4561 0.677204
  // 1  -15 103.29   0.723067
  // 1  -10 108.759  0.777696
  // 1   -5 113.737  0.842778
  // 1    0 118.112  0.920356
  // 1    5 120.051  0.928594
  // 1   10 114.469  0.959548

//  for (ii=0; ii<Nsteps; ii++)
//  {   fff << chrnumber << " " << moveQTL << " " << Frun[ii][Nrun]
//          << " " << ((informationcontent[ii]/Nfam)/(Nrun+1)) << endl;
//      if (moveQTL+stepsize<=stepmax) moveQTL+= stepsize;
//      else { moveQTL= stepmin; chrnumber++; }
//  }
//  fff << ":" << endl;
//  cout << "-2- " << endl;
//  if (Nrun>0) // ((Nfam>1)&&(Nrun>0))  // multiple permutations or simulations
//  {  for (ii=0; ii<Nsteps; ii++)
//       for (i=0; i<Nrun; i++)
 //        Frun[0][i]= (Frun[0][i]<Frun[ii][i] ? Frun[ii][i] : Frun[0][i]);
 //    cout << endl;
//     cout << "Cumulative distribution of maximum test statistic value in "
//          << Nrun << " permutations or simulations" << endl;
 //    fff  << "Cumulative distribution of maximum test statistic value in "
//          << Nrun << " permutations or simulations" << endl;
//     sort1(Nrun,Frun[0]);
//     if (Nrun>1)
 //      for (i=1; /* (((double)run/( (double)Nrun+1.0))<0.1)*/ i<Nrun+1; i++)
 //      {   cout << setprecision(8) << ( (double)i/( (double)Nrun+1.0) )
 //               << " " << setprecision(8) << Frun[0][i-1] << endl;
 //          fff << setprecision(8) << ( (double)i/( (double)Nrun+1.0) )
 //              << " " << setprecision(8) << Frun[0][i-1] << endl;
 //      }
 //    fff.close();
//  }
 // delmatrix(Frun,Nsteps);
  return 0;   
}

/* end of scanMQM.c */
