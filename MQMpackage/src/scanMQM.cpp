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

extern "C"
{

#include <R.h>
#include <Rdefines.h>
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
			 double Stepmi,double Stepma,int NRUN,int out_Naug,int **INDlist, double **QTL, int re_estimate,int crosstype,int domi){
	
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
	for(int i=0; i< Nmark; i++){
		for(int j=0; j< Nind; j++){
			//Some lame ass checks to see if the cross really is the cross we got (So BC can't contain 3's (BB) and RILS can't contain 2's (AB)
			if(Geno[i][j] > 3 && crosstype != 1){
				Rprintf("INFO: Stange genotype pattern, switching to F2\n");
				crosstype = 1;
				break;
			}
			if(Geno[i][j] == 3 && crosstype == 2){
				Rprintf("INFO: Stange genotype pattern, switching from BC to F2\n");
				crosstype = 1;
				break;
			}
			//IF we have a RIL and find AB then Rqtl messed up, so we have a BC genotype
			if(Geno[i][j] == 2 && crosstype == 3){
				Rprintf("INFO: Stange genotype pattern, switching from RISELF to BC\n");
				crosstype = 2;
				break;
			}
			
		}
		//Rprintf("\n");
	}
	
	char reestimate = 'y';
	if(re_estimate == 0){
		reestimate = 'n';
	}
	char cross = 'F';
	if(crosstype == 1){
		cross = 'F';	
	}
	if(crosstype == 2){
		cross = 'B';	
	}
	if(crosstype == 3){
		cross = 'R';	
	}	
	char dominance='n';
	if(domi != 0){
		dominance='y';
	}	
	Rprintf("We got all the needed information, so lets start with the MQM\n");   
	analyseF2(Nind, Nmark, cofactor, markers, Pheno[(Npheno-1)], f1genotype, Backwards,QTL,&mapdistance,Chromo,NRUN,RMLorML,Windowsize,Steps,Stepmi,Stepma,Alfa,Emiter,out_Naug,INDlist,reestimate,cross,dominance);
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
			   double *stepmi,double *stepma, int *nRun,int *out_Naug,int *indlist,  double *qtl,int *reestimate,int *crosstype,int *domi){
   int **Geno;
   int **Chromo;
   double **Dist;  
   double **Pheno;   
   double **QTL;   
   int **Cofactors;
   int **INDlist;
   
   //Reorganise the pointers into arrays, ginletons are just cast into the function
   reorg_geno(*Nind,*Nmark,geno,&Geno);
   reorg_int(*Nmark,1,chromo,&Chromo);   
   reorg_pheno(*Nmark,1,dist,&Dist);
   //Here we have  the assumption that step.min is negative this needs to be split in 2
   reorg_pheno((*chromo) * (((*stepma)-(*stepmi))/ (*steps)),1,qtl,&QTL);
   reorg_pheno(*Nind,*Npheno,pheno,&Pheno);
   reorg_int(*Nmark,1,cofactors,&Cofactors);  
   reorg_int(*out_Naug,1,indlist,&INDlist);  
   //Done with reorganising lets start executing
   
   scanMQM(*Nind,*Nmark,*Npheno,Geno,Chromo,Dist,Pheno,Cofactors,*backwards,*RMLorML,*alfa,*emiter,*windowsize,*steps,*stepmi,*stepma,*nRun,*out_Naug,INDlist,QTL, *reestimate,*crosstype,*domi);
} /* end of function R_scanMQM */

}
