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

void R_scanMQM(int *Nind,int *Nmark,int *Npheno, int *Nfam,int *geno, double *pheno);


/**********************************************************************
 * 
 * scanMQM
 *
 * 
 **********************************************************************/

void scanMQM(int Nind, int Nmark,int Npheno, int Nfam,int **Geno, double **Pheno);

/**********************************************************************
 * 
 * Helper functions
 *
 *
 **********************************************************************/


void OK();
double Lnormal(double residual, double variance);
double absdouble(double x);
int mod(int a, int b);

/* end of scanMQM.h */
