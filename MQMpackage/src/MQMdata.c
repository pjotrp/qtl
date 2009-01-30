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
#include "MQMdata.h"

//#include <alloc.h> // for coreleft()


vector newvector(int dim)
{      vector v;
       v = (double *)R_alloc(dim, sizeof(double));
       if (v==NULL) { 
         Rprintf("Not enough memory for new vector of dimension %d",(dim+1));
         exit(1); 
       }
       return v;
}

ivector newivector(int dim)
{      ivector v;
       v = (int *)R_alloc(dim, sizeof(int));
       if (v==NULL) { 
         Rprintf("Not enough memory for new vector of dimension %d",(dim+1));
         exit(1); 
       }
       return v;
}

cvector newcvector(int dim)
{      cvector v;
       v = (char *)R_alloc(dim, sizeof(char));
       if (v==NULL) { 
         Rprintf("Not enough memory for new vector of dimension %d",(dim+1));
         exit(1); 
       }
       return v;
}

matrix newmatrix(int rows, int cols)
{      matrix m;
       m = (double **)R_alloc(rows, sizeof(double *));
       if (m==NULL) { 
         Rprintf("Not enough memory for new double matrix");
         exit(1); 
       }
       for (int i=0; i<rows; i++) m[i]= newvector(cols);
       return m;
}

void   printmatrix(matrix m, int rows, int cols)
{      for (int r=0; r<rows; r++)
       {   for (int c=0; c<cols; c++) Rprintf("%d",m[r][c]);
           Rprintf("\n");
       }
}

void   printcmatrix(cmatrix m, int rows, int cols)
{      for (int r=0; r<rows; r++)
       {   for (int c=0; c<cols; c++) Rprintf("%c",m[r][c]);
           Rprintf("\n");
       }
}

cmatrix newcmatrix(int rows, int cols)
{      cmatrix m;
       m = (char **)R_alloc(rows, sizeof(char *));
       if (m==NULL) { 
         Rprintf("Not enough memory for new char matrix");
         exit(1); 
       }
       for (int i=0; i<rows; i++) m[i]= newcvector(cols);
       return m;
}

void delmatrix(matrix m, int rows)
{      
	  for (int i=0; i<rows; i++) free((void*)m[i]);
      free((void*)m);
}

void delcmatrix(cmatrix m, int rows)
{     
	  for (int i=0; i<rows; i++) free((void*)m[i]);
      free((void*)m);
}

void copyvector(vector vsource, vector vdestination, int dim)
{    for (int i=0; i<dim; i++) vdestination[i]= vsource[i];
}



/* end of MQMdata.h */
