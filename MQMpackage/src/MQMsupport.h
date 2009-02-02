/**********************************************************************
 * 
 * MQMsupport.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
 * first written Feb, 2009
 *
 * C external functions used by the MQM algorithm
 * Contains: 
 *
 **********************************************************************/

 /*
 * simuF2 for every individual calculate a random cvariance (y). Next the * markers are walked and depending on type the cvariance is adjusted by +/- 1
 */
void simuF2(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y);

/* analyseF2 - analyse one F2 family */

void analyseF2(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, ivector f1genotype, int Backwards, double **QTL);
double probleft(char c, int jloc, cvector imarker, vector r, cvector position);
double probright(char c, int jloc, cvector imarker, vector r, cvector position);

/* 
 * augmentdata inserts missing data - does not use the phenotypes - augments the
 * marker data. Phenotype checking is done in the EM step.
 */

void augmentdata(cmatrix marker, vector y, cmatrix *augmarker, vector *augy, ivector *augind, int *Nind, int *Naug, int Nmark, cvector position, vector r);

/* ML estimation of recombination frequencies via EM;  calculation of multilocus genotype probabilities;  ignorance of unlikely genotypes*/
void rmixture(cmatrix marker, vector weight, vector r, cvector position, ivector ind, int Nind, int Naug, int Nmark);

/* ML estimation of parameters in mixture model via EM;
*/
double QTLmixture(cmatrix loci, cvector cofactor, vector r, cvector position,
              vector y, ivector ind, int Nind, int Naug, int Nloci, double *variance, int em, vector *weight);

/* regression of trait on multiple cofactors
   y=xb+e with weight w   (xtwx)b=(xtw)y   b=inv(xtwx)(xtw)y */

double regression(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, vector* weight, ivector ind, int Naug, double *variance, vector Fy, char biasadj);


/* backward elimination in regression of trait on multiple cofactors
   routine subX haalt uit matrices voor volledige model de submatrices voor submodellen;
   matrices XtWX en Xt van volledig model worden genoemd fullxtwx en fullxt;
   analoog vector XtWY wordt full xtwy genoemd;
*/
double backward(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, vector weight, int* ind, int Naug, double logLfull, double variance, double F1, double F2, cvector* newcofactor, vector r, cvector position);

/* mapQTL */
double mapQTL(int Nind, int Nmark, cvector cofactor, cvector selcofactor, cmatrix marker, cvector position, vector mapdistance, vector y, vector r, ivector ind, int Naug, double variance, char printoutput);


/*
-----------------------------------------------------------------------
subroutines from book 'Numerical Recipees in C' for calculating F-probabilities and 
for generating randomly permuted trait data for other tasks
-----------------------------------------------------------------------*/

void ludcmp(matrix m, int dim, ivector ndx, int *d);
void lusolve(matrix lu, int dim, ivector ndx, vector b);
double gammln(double xx);
double betai(double a, double b, double x);
double betacf(double a, double b, double x);
double inverseF(int df1, int df2, double alfa);
double ran2(long *idum);
double randomnormal(long *idum);
void sort1(int n, vector ra);
void sort2(int n, double *ra, ivector rb);

/* end of MQMsupport.h */
