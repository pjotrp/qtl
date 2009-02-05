/**********************************************************************
 * 
 * MQMdata.h
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
 * first written Feb, 2009
 *
 * Several datastructures needed by the MQM algorithm are defined here
 * Contains: 
 *
 **********************************************************************/


/*------------------------------------------------------------------------
Datastructures for matrix and vector calculus
------------------------------------------------------------------------ */
typedef double** matrix;
typedef double*  vector;
typedef char**   cmatrix;
typedef char*    cvector;
typedef int*  ivector;

/*
Data augmentation routing
*/

void R_augdata(int *geno,double *dist,double *pheno,int *auggeno,double *augPheno,int *Nind,int *Naug,int *Nmark, int *Npheno, int *maxaug, int *maxiaug,double *neglect);

int augdata(cmatrix marker, vector y, cmatrix* augmarker, vector *augy, ivector* augind, int *Nind, int* Naug, int Nmark, cvector position, vector r,int maxNaug,int imaxNaug,double neglect);

/*------------------------------------------------------------------------
Basic routines for matrix and vector calculus
------------------------------------------------------------------------ */
vector newvector(int dim);
ivector newivector(int dim);
cvector newcvector(int dim);
matrix newmatrix(int rows, int cols);
void   printmatrix(matrix m, int rows, int cols);
void   printcmatrix(cmatrix m, int rows, int cols);
cmatrix newcmatrix(int rows, int cols);
void delmatrix(matrix m, int rows);
void delcmatrix(cmatrix m, int rows);
void copyvector(vector vsource, vector vdestination, int dim);
/* end of MQMdata.h */
