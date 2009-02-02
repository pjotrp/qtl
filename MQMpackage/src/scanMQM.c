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
//#include <alloc.h> // for alloc,free & coreleft()
#include "MQMdata.h"
#include "MQMsupport.h"

double neglect=100; // eliminate unlikely genotype configurations
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
int Nfam=1, Nrun=0;
int run=-1;
char REMLorML='0';
char fitQTL='n';
char dominance='n';
char perm_simu='1';

char ok='0', defset='1';

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

void reorg_int(int n_ind, int n_mar, int *pheno, int ***Pheno)
{
  int i;

  *Pheno = (int **)R_alloc(n_mar, sizeof(int *));

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

void scanMQM(int Nind, int Nmark,int Npheno, int Nfam,int **Geno,int **Chromo, 
			 double **Dist, double **Pheno, int **Cofactors, int Backwards, int RMLorML,double Alfa,int Emiter,
			 double Windowsize,double Steps,
			 double Stepmi,double Stepma, double **QTL){
   ivector f1genotype;
   f1genotype = newivector(Nmark);
   cmatrix markers;
   markers= newcmatrix(Nmark,Nind);
   //Rprintf("Printing Genotype matrix\n");
   for(int i=0; i< Nmark; i++){
     f1genotype[i] = 12;
     for(int j=0; j< Nind; j++){ 
	    markers[i][j] = '9';
	    if(Geno[i][j] == 1){
		 markers[i][j] = '0';
		}
	    if(Geno[i][j] == 2){
		 markers[i][j] = '2';
		}
	    if(Geno[i][j] == 3){
		 markers[i][j] = '1';
		}
		
	//	Rprintf("%d ",markers[i][j]);
     }
   //  Rprintf("\n");
   }

   //Rprintf("Printing Chromosome matrix\n");
   for(int i=0; i< Nmark; i++){
     // Rprintf("%d\n",Chromo[0][i]);         
   }

   //Rprintf("Printing Phenotype matrix\n");
   for(int i=0; i< Npheno; i++){
     for(int j=0; j< Nind; j++){
      // Rprintf("%f\n",Pheno[i][j]);         
     }
    // Rprintf("\n");
   }   

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
   int cnt=0;   
   for(int i=0; i< Nmark; i++){
      //Filling chromosome information in the MQM style
      chr[i] = Chromo[0][i];
	  cofactor[i] = '0';
	  if(Cofactors[0][i] == 1){
        cofactor[i] = '1';
		cnt++;
	  }
	  if(Cofactors[0][i] == 2){
	    cofactor[i] = '2';
		cnt++;
	  }
	  if(cnt > (Nmark/2)){
	    Rprintf("ERROR: More than half of the markers are to be cofactors, this is not allowed\n");
		return;
	  }
	//  Rprintf("Cofactor %d : %c\n",i,cofactor[i]);
      mapdistance[i]=999.0;
	  mapdistance[i]=Dist[0][i];
	//  Rprintf("Mapdist %d: %f <-> %f\n",i,mapdistance[i],Dist[0][i]);
   }
   if (real_simu=='1'){
     simuF2(Nind, Nmark, cofactor, markers, Pheno[0]);
   }
   //SETTING UP GLOBAL VARIABLES
   windowsize=Windowsize; // used in mapQTL procedure
   stepsize=Steps; // size of steps when moving QTL along chromosomes (for output)
   stepmin=Stepmi; // start moving QTL at position stepmin cM (for output)
   stepmax=Stepma; // move QTL up to stepmax (for output)
   REMLorML='0';
   if(RMLorML == 1){
		REMLorML='1';
   }
   em=Emiter; // maximum number of em iterations
   alfa=Alfa; // alfa used in selection procedure
   Rprintf("We got all the needed information, so lets start with the MQM\n");   
   analyseF2(Nind, Nmark, cofactor, markers, Pheno[0], f1genotype, Backwards,QTL);
   return;
}  /* end of function scanMQM */



/**********************************************************************
 * 
 * R_scanMQM
 * 
 **********************************************************************/

void R_scanMQM(int *Nind,int *Nmark,int *Npheno, int *Nfam,
			   int *geno,int *chromo, double *dist, double *pheno, 
			   int *cofactors, int *backwards, int *RMLorML,double *alfa,int *emiter,
			   double *windowsize,double *steps,
			   double *stepmi,double *stepma, double *qtl){
   int **Geno;
   int **Chromo;
   double **Dist;  
   double **Pheno;   
   double **QTL;   
   int **Cofactors;
   
   //Reorganise the pointers into arrays, ginletons are just cast into the function
   reorg_geno(*Nind,*Nmark,geno,&Geno);
   reorg_int(*Nmark,1,chromo,&Chromo);   
   reorg_pheno(*Nmark,1,dist,&Dist);
   reorg_pheno((*chromo) * (( (*stepmi)+(*stepma))/ (*steps)),1,qtl,&QTL);
   reorg_pheno(*Nind,*Npheno,pheno,&Pheno);
   reorg_int(*Nmark,1,cofactors,&Cofactors);  
   //Done with reorganising lets start executing the main loop
   
   scanMQM(*Nind,*Nmark,*Npheno, *Nfam,Geno,Chromo,Dist,Pheno,Cofactors,*backwards,*RMLorML,*alfa,*emiter,*windowsize,*steps,*stepmi,*stepma,QTL);
} /* end of function R_scanMQM */


