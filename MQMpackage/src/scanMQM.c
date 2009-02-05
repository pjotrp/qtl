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
#include <R_ext/PrtUtil.h>
#include <R_ext/RS.h> /* for Calloc, Realloc */
#include <R_ext/Utils.h>
#include "MQMdata.h"
#include "MQMsupport.h"

//long *idum; // for monte carlo simulation or permutation

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

void scanMQM(int Nind, int Nmark,int Npheno,int **Geno,int **Chromo, 
			 double **Dist, double **Pheno, int **Cofactors, int Backwards, int RMLorML,double Alfa,int Emiter,
			 double Windowsize,double Steps,
			 double Stepmi,double Stepma,int NRUN, double **QTL){
	
	ivector f1genotype;
	cmatrix markers;
	cvector cofactor;
	vector mapdistance;
	
	markers= newcmatrix(Nmark,Nind);
	f1genotype = newivector(Nmark);
	cofactor= newcvector(Nmark);  
	mapdistance= newvector(Nmark);
	
	int cnt=0;  	
   
	//Rprintf("Converting Genotype matrix\n");
	for(int i=0; i< Nmark; i++){
		f1genotype[i] = 12;
		//receiving mapdistances
		mapdistance[i]=999.0;
		mapdistance[i]=Dist[0][i];
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
		for(int j=0; j< Nind; j++){ 
			markers[i][j] = '9';
			if(Geno[i][j] == 1){				//AA
				markers[i][j] = '0';
			}
			if(Geno[i][j] == 2){				//AB
				markers[i][j] = '1';
			}
			if(Geno[i][j] == 3){				//BB
				markers[i][j] = '2';
			}
			if(Geno[i][j] == 4){				//AA of AB
				markers[i][j] = '4';
			}
			if(Geno[i][j] == 5){				//BB of AB
				markers[i][j] = '3';
			}
		}
	}

	Rprintf("We got all the needed information, so lets start with the MQM\n");   
	analyseF2(Nind, Nmark, cofactor, markers, Pheno[0], f1genotype, Backwards,QTL,&mapdistance,Chromo,NRUN,RMLorML,Windowsize,Steps,Stepmi,Stepma,Alfa,Emiter);
	//Rprintf("Starting Cleanup\n");
	delcmatrix(markers,Nmark);
	Free(f1genotype);
	Free(cofactor);
	Free(mapdistance);
	return;
}  /* end of function scanMQM */



/**********************************************************************
 * 
 * R_scanMQM
 * 
 **********************************************************************/

void R_scanMQM(int *Nind,int *Nmark,int *Npheno,
			   int *geno,int *chromo, double *dist, double *pheno, 
			   int *cofactors, int *backwards, int *RMLorML,double *alfa,int *emiter,
			   double *windowsize,double *steps,
			   double *stepmi,double *stepma, int *nRun, double *qtl){
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
   
   scanMQM(*Nind,*Nmark,*Npheno,Geno,Chromo,Dist,Pheno,Cofactors,*backwards,*RMLorML,*alfa,*emiter,*windowsize,*steps,*stepmi,*stepma,*nRun,QTL);
} /* end of function R_scanMQM */


