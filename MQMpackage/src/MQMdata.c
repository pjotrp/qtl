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
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include "MQMdata.h"

vector newvector(int dim)
{      vector v;
       v= new double[dim];
       if (v==NULL) { cout << " not enough memory for new vector of dimension " << (dim+1); exit(1); }
       return v;
}

ivector newivector(int dim)
{      ivector v;
       v= new int[dim];
       if (v==NULL) { cout << " not enough memory for new vector of dimension " << (dim+1); exit(1); }
       return v;
}

cvector newcvector(int dim)
{      cvector v;
       v= new char[dim];
       if (v==NULL) { cout << " not enough memory for new char vector of dimension " << (dim+1); exit(1); }
       return v;
}

matrix newmatrix(int rows, int cols)
{      matrix m;
       m=new double*[rows];
       if (m==NULL) { cout << " not enough memory for new double matrix"; exit(1); }
       for (int i=0; i<rows; i++) m[i]= newvector(cols);
       return m;
}

void   printmatrix(matrix m, int rows, int cols)
{      for (int r=0; r<rows; r++)
       {   for (int c=0; c<cols; c++) cout << setw(5) << m[r][c] << " ";
           cout << endl;
       }
}

void   printcmatrix(cmatrix m, int rows, int cols)
{      for (int r=0; r<rows; r++)
       {   for (int c=0; c<cols; c++) cout << setw(2) << m[r][c];
           cout << endl;
       }
}

cmatrix newcmatrix(int rows, int cols)
{      cmatrix m;
       m=new char*[rows];
       if (m==NULL) { cout << " not enough memory for new char matrix"; exit(1); }
       for (int i=0; i<rows; i++) m[i]= newcvector(cols);
       return m;
}

void delmatrix(matrix m, int rows)
{      for (int i=0; i<rows; i++) delete[] m[i];
       delete[] m;
}

void delcmatrix(cmatrix m, int rows)
{      for (int i=0; i<rows; i++) delete[] m[i];
       delete[] m;
}

void copyvector(vector vsource, vector vdestination, int dim)
{    for (int i=0; i<dim; i++) vdestination[i]= vsource[i];
}



/* end of MQMdata.h */
