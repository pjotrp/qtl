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
	
	int cof_cnt=0;
 	
   	//Change all the markers from Karl format to our own
	change_coding(&Nmark,&Nind,Geno,markers);

	for(int i=0; i< Nmark; i++){
		f1genotype[i] = 12;
		//receiving mapdistances
		mapdistance[i]=999.0;
		mapdistance[i]=Dist[0][i];
	 	cofactor[i] = '0';
		if(Cofactors[0][i] == 1){
			cofactor[i] = '1';
			cof_cnt++;
		}
		if(Cofactors[0][i] == 2){
			cofactor[i] = '2';
			cof_cnt++;
		}
		if(cof_cnt+10 > Nind){
			Rprintf("ERROR: Setting this many cofactors would leave less than 10 degrees of freedom.\n");
			return;
		}
	}

	char reestimate = 'y';
	if(re_estimate == 0){
		reestimate = 'n';
	}
	//determine what kind of cross we have
	char cross = determin_cross(&Nmark,&Nind,Geno,&crosstype);
	//set dominance accordingly
	if(cross != 'F'){
		Rprintf("INFO: Dominance setting ignored (dominance=0)\n");   
		domi = 0;
	}else{
		domi= domi;
	}

	char dominance='n';
	if(domi != 0){
		dominance='y';
	}	
	cout << "INFO: All the needed information, so lets start with the MQM\n";   
	analyseF2(Nind, Nmark, &cofactor, markers, Pheno[(Npheno-1)], f1genotype, Backwards,QTL,&mapdistance,Chromo,NRUN,RMLorML,Windowsize,Steps,Stepmi,Stepma,Alfa,Emiter,out_Naug,INDlist,reestimate,cross,dominance);
	if(re_estimate){
		Rprintf("INFO: Sending back the reestimated map used during analysis\n");
		for(int i=0; i< Nmark; i++){
			Dist[0][i]=mapdistance[i];
		}
	}
	if(Backwards){
		Rprintf("INFO: Sending back the model\n");
		for(int i=0; i< Nmark; i++){
			Cofactors[0][i]=cofactor[i];
		}
	}	
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
   reorg_pheno(2*(*chromo) * (((*stepma)-(*stepmi))/ (*steps)),1,qtl,&QTL);
   reorg_pheno(*Nind,*Npheno,pheno,&Pheno);
   reorg_int(*Nmark,1,cofactors,&Cofactors);  
   reorg_int(*out_Naug,1,indlist,&INDlist);  
   //Done with reorganising lets start executing
   
   scanMQM(*Nind,*Nmark,*Npheno,Geno,Chromo,Dist,Pheno,Cofactors,*backwards,*RMLorML,*alfa,*emiter,*windowsize,*steps,*stepmi,*stepma,*nRun,*out_Naug,INDlist,QTL, *reestimate,*crosstype,*domi);
} /* end of function R_scanMQM */

int count_lines(char *file){
	//NUM: number of elements on 1 line
	int cnt=0;
	char line[100];
	ifstream file_stream(file, ios::in);
	while (!file_stream.eof()){
        file_stream >> line;
		cnt++;
	}
	file_stream.close();
	return cnt;
}


int main(){

	char *genofile = "geno.dat";
	char *phenofile = "pheno.dat";
	char *mposfile = "markerpos.txt";
	char *chrfile = "chrid.dat";
    double **QTL;  
	ivector f1genotype;
	ivector chr;
	cvector cofactor;
	vector mapdistance;
	cmatrix markers;
	ivector INDlist;
       
	int cnt=0;
	int cInd=0; //Curretn individual
	int nInd;
	int nMark;
	
	nInd=count_lines(phenofile);
	printf("# of individuals: %d\n",nInd);
	nMark=count_lines(chrfile);
	printf("# of markers: %d\n",nMark);	
    f1genotype = newivector(nMark);	
	cofactor= newcvector(nMark);  
	mapdistance= newvector(nMark);
	markers= newcmatrix(nMark,nInd);
	double pheno_value[nInd];
	chr = newivector(nMark);
	INDlist= newivector(nInd);
	double pos[nMark];

	char peek_c;

	ifstream geno(genofile, ios::in);
	while (!geno.eof()){
        if(cnt < nMark){
          	geno >> markers[cnt][cInd];
		//	printf(" %c",markers[cnt][cInd]);
			cnt++;
        }else{
		//	printf("\n");
			cnt = 0;
			cInd++;
		}	
	}
	geno.close();
	printf("Genotypes done %d %d\n",cInd,cnt);
	cnt = 0;
	ifstream pheno(phenofile, ios::in);
	while (!pheno.eof()){
		pheno >> pheno_value[cnt];
	//	printf("%f\n",pheno_value[cnt]);
		cnt++;
	}
	pheno.close();
	printf("Phenotype done %d\n",cnt);
	cnt = 0;
	ifstream mpos(mposfile, ios::in);
	while (!mpos.eof()){
		peek_c=mpos.peek();
    	if(peek_c=='\t' or peek_c == ' '){
           	mpos >> pos[cnt];
   //         printf("%f\n",pos[cnt]);
            cnt++;
		}else{
            mpos >> peek_c;
        }
	}	
	mpos.close();

    printf("Positions done %d\n",cnt);	
	cnt = 0;	
	ifstream chrstr(chrfile, ios::in);
	while (!chrstr.eof()){
		chrstr >> chr[cnt];
//		printf("%d\n",chr[cnt]);
		cnt++;
	}
	chrstr.close();
	printf("Chromosomes done %d\n",cnt);
	for(int i=0; i< nMark; i++){
    	cofactor[i] = '0';
    	f1genotype[nMark] = 12;
    	mapdistance[i]=999.0;
		mapdistance[i]=pos[i];
    }
	for(int i=0; i< nInd; i++){
    	INDlist[i] = i;
    }

 	Rprintf("Cofactor done, starting analyseF2\n",cnt);
	//ALL information is read in or calculated, so we gonna start MQM, however Rprintf crashes MQM
   	analyseF2(nInd, nMark, &cofactor, markers, pheno_value, f1genotype, 0,QTL, &mapdistance,&chr,0,0,5,5,0,220,0.05,1000,0,&INDlist,0,1,0);
	return 1;
}


}
