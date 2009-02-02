/**********************************************************************
 * 
 * MQMdata.c
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
 * first written Feb, 2009
 *
 * Several basic routines needed by the MQM algorithm are defined here
 * Contains: 
 *
 **********************************************************************/
#include <R.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/RS.h> /* for Calloc, Realloc */
#include <R_ext/Utils.h>
#include "MQMdata.h"

vector newvector(int dim)
{      vector v;
       v = (double *)Calloc(dim, double);
       if (v==NULL) { 
         warning("Not enough memory for new vector of dimension %d",(dim+1)); 
       }
       return v;
}

ivector newivector(int dim)
{      ivector v;
       v = (int *)Calloc(dim, int);
       if (v==NULL) { 
         warning("Not enough memory for new vector of dimension %d",(dim+1));
       }
       return v;
}

cvector newcvector(int dim)
{      cvector v;
       v = (char *)Calloc(dim, char);
       if (v==NULL) { 
         warning("Not enough memory for new vector of dimension %d",(dim+1));
       }
       return v;
}

matrix newmatrix(int rows, int cols){
    matrix m;
    m = (double **)Calloc(rows, double*);
    if (m==NULL) { 
        warning("Not enough memory for new double matrix");
    }
    for (int i=0; i<rows; i++){
		m[i]= newvector(cols);
	}
    return m;
}

void printmatrix(matrix m, int rows, int cols){
      
	for (int r=0; r<rows; r++){
		for (int c=0; c<cols; c++){
			Rprintf("%d",m[r][c]);
		}
        Rprintf("\n");
	}
}

void printcmatrix(cmatrix m, int rows, int cols){
      
	for (int r=0; r<rows; r++){   
		for (int c=0; c<cols; c++){
			Rprintf("%c",m[r][c]);
		}
        Rprintf("\n");
	}
}

cmatrix newcmatrix(int rows, int cols){
    cmatrix m;
    m = (char **)Calloc(rows, char*);
    if (m==NULL) {
		warning("Not enough memory for new char matrix");
	}
    for (int i=0; i<rows; i++){
		m[i]= newcvector(cols);
	}
	return m;
}

void delmatrix(matrix m, int rows){
      
	Free(m);
}

void delcmatrix(cmatrix m, int rows){
     
     Free(m);
}

void copyvector(vector vsource, vector vdestination, int dim){
    
	for (int i=0; i<dim; i++){
		vdestination[i]= vsource[i];
	}
}



/* end of MQMdata.h */
