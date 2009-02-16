/**********************************************************************
 * 
 * MQMprob.c
 *
 * copyright (c) 2009 Danny Arends
 * last modified Feb, 2009
 * first written Feb, 2009
 *
 * C external functions used by the MQM algorithm
 * Contains: 
 *
 **********************************************************************/
#include <R.h>
#include <math.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/RS.h> /* for Calloc, Realloc */
#include <R_ext/Utils.h>
#include "scanMQM.h"
#include "MQMdata.h"
#include "MQMsupport.h"
#include "Regression.h"
#include "MQMmapQTL.h"
#include "MQMmixture.h"

//extern long *idum; // for monte carlo simulation or permutation


double probleft(char c, int jloc, cvector imarker, vector r, cvector position)
{    double nrecom;
     if ((position[jloc]=='L')||(position[jloc]=='U')) return (c=='1' ? 0.50 : 0.25);
     else if ((c=='1')&&(imarker[jloc-1]=='1')) return r[jloc-1]*r[jloc-1]+(1.0-r[jloc-1])*(1.0-r[jloc-1]);
     else
     {  nrecom= absdouble(c-imarker[jloc-1]);
        if      (nrecom==0) return (1.0-r[jloc-1])*(1.0-r[jloc-1]);
        else if (nrecom==1) return (c=='1' ? 2.0*r[jloc-1]*(1.0-r[jloc-1]) : r[jloc-1]*(1.0-r[jloc-1]));
        else return r[jloc-1]*r[jloc-1];
     }
}

double probright(char c, int jloc, cvector imarker, vector r, cvector position)
{    double nrecom, prob0, prob1, prob2;
     if ((position[jloc]=='R')||(position[jloc]=='U')) return 1.0;
     else if ((imarker[jloc+1]=='0')||(imarker[jloc+1]=='1')||(imarker[jloc+1]=='2'))
     {   if ((c=='1')&&(imarker[jloc+1]=='1'))
            return r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]);
         else
         {  nrecom= absdouble(c-imarker[jloc+1]);
            if      (nrecom==0) return (1.0-r[jloc])*(1.0-r[jloc]);
            else if (nrecom==1) return (imarker[jloc+1]=='1' ? 2.0*r[jloc]*(1.0-r[jloc]) : r[jloc]*(1.0-r[jloc]));
            else return r[jloc]*r[jloc];
         }
     }
     else if (imarker[jloc+1]=='3')
     {  if      (c=='0') { prob1= 2.0*r[jloc]*(1.0-r[jloc]);
                           prob2= r[jloc]*r[jloc]; }
        else if (c=='1') { prob1= r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]);
                           prob2= r[jloc]*(1.0-r[jloc]); }
        else             { prob1= 2.0*r[jloc]*(1.0-r[jloc]);
                           prob2= (1.0-r[jloc])*(1-r[jloc]); }
        return prob1*probright('1',jloc+1,imarker,r,position) +
               prob2*probright('2',jloc+1,imarker,r,position);
     }
     else if (imarker[jloc+1]=='4')
     {  if      (c=='0') { prob0= (1.0-r[jloc])*(1.0-r[jloc]);
                           prob1= 2.0*r[jloc]*(1.0-r[jloc]); }
        else if (c=='1') { prob0= r[jloc]*(1.0-r[jloc]);
                           prob1= r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]); }
        else             { prob0= r[jloc]*r[jloc];
                           prob1= 2.0*r[jloc]*(1.0-r[jloc]); }
        return prob0*probright('0',jloc+1,imarker,r,position) +
               prob1*probright('1',jloc+1,imarker,r,position);
     }
     else // (imarker[j+1]=='9')
     {  if      (c=='0') { prob0= (1.0-r[jloc])*(1.0-r[jloc]);
                           prob1= 2.0*r[jloc]*(1.0-r[jloc]);
                           prob2= r[jloc]*r[jloc]; }
        else if (c=='1') { prob0= r[jloc]*(1.0-r[jloc]);
                           prob1= r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]);
                           prob2= r[jloc]*(1.0-r[jloc]); }
        else             { prob0= r[jloc]*r[jloc];
                           prob1= 2.0*r[jloc]*(1.0-r[jloc]);
                           prob2= (1.0-r[jloc])*(1.0-r[jloc]); }
        return prob0*probright('0',jloc+1,imarker,r,position) +
               prob1*probright('1',jloc+1,imarker,r,position) +
               prob2*probright('2',jloc+1,imarker,r,position);
     }
}

/* end of MQMprob.c */
