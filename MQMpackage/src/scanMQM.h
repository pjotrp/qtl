/**********************************************************************
 * 
 * scanMQM.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
 * first written Feb, 2009
 *
 * C functions for the R/qtl package
 * Contains: R_scanMQM, scanMQM
 *
 **********************************************************************/

/**********************************************************************
 * 
 * R_scanMQM
 * Wrapper for call from R;
 * 
 **********************************************************************/

void R_scanMQM(int *Nind,int *Nmark,int *Npheno, int *Nfam,
			   int *geno,int *chromo, double *dist, double *pheno, 
			   int *cofactors, int *backwards, int *RMLorML,double *alfa,int *emiter,
			   double *windowsize,double *steps,
			   double *stepmi,double *stepma, double *qtl);


/**********************************************************************
 * 
 * scanMQM
 *
 * 
 **********************************************************************/

void scanMQM(int Nind, int Nmark,int Npheno, int Nfam,int **Geno,int **Chromo, 
			 double **Dist, double **Pheno, int **Cofactors, int Backwards, int RMLorML,double Alfa,int Emiter,
			 double Windowsize,double Steps,
			 double Stepmi,double Stepma, double **QTL);

/**********************************************************************
 * 
 * Helper functions
 *
 *
 **********************************************************************/


double Lnormal(double residual, double variance);
double absdouble(double x);
int mod(int a, int b);

/* end of scanMQM.h */
