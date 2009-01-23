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

void R_scanMQM();

/**********************************************************************
 * 
 * scanMQM
 * the workhorse function
 *
 **********************************************************************/

void scanMQM();

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
