/**********************************************************************
 * 
 * MQMsupport.c
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
#include "MQMprob.h"
#include "Regression.h"

/* regression of trait on multiple cofactors  y=xb+e with weight w
*							(xtwx)b=(xtw)y  
*							b=inv(xtwx)(xtw)y 
*/

double regression(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y,
                vector *weight, ivector ind, int Naug,
                double *variance, vector Fy, char biasadj,char fitQTL,char dominance)
{    // cout << "regression IN" << endl;
     /*
     cofactor[j] at locus j:
     '0': no cofactor at locus j
     '1': cofactor at locus j
     '2': QTL at locus j, but QTL effect is not included in the model
     '3': QTL at locu j and QTL effect is included in the model
     */
	//for (int j=0; j<Naug; j++){
	//   Rprintf("J:%d,COF:%d,VAR:%f,WEIGHT:%f,Trait:%f,IND[j]:%d\n",j,cofactor[j],*variance,(*weight)[j],y[j],ind[j]);
    //}

	matrix XtWX;
	cmatrix Xt;
	vector XtWY;
	int dimx=1, j, jj;
    for (int j=0; j<Nmark; j++)
    if (cofactor[j]=='1') dimx+= (dominance=='n' ? 1 : 2);  // per QTL only additivity !!
    else if (cofactor[j]=='2') { dimx+=1;}
	
	XtWX= newmatrix(dimx+2,dimx+2);
    Xt= newcmatrix(dimx+2,Naug);
	XtWY= newvector(dimx+2);
    
	dimx=1;	
	for (j=0; j<Nmark; j++)
     if ((cofactor[j]=='1')||(cofactor[j]=='3')) dimx+= (dominance=='y' ? 2 : 1);
     cvector xtQTL; // '0'=mu; '1'=cofactor; '2'=QTL (additive); '3'= QTL (dominance);
     xtQTL= newcvector(dimx);
     int jx=0;
     for (int i=0; i<Naug; i++) Xt[jx][i]= '1';
     xtQTL[jx]= '0';

     for (j=0; j<Nmark; j++)
     if (cofactor[j]=='1') // cofactor (not a QTL moving along the chromosome)
     {  jx++;
        xtQTL[jx]= '1';
        if (dominance=='y')
        {  for (int i=0; i<Naug; i++)
           if      (marker[j][i]=='1') { Xt[jx][i]=48; Xt[jx+1][i]=49; }  //ASCII code 47, 48 en 49 voor -1,0,1;
           else if (marker[j][i]=='0') { Xt[jx][i]=47; Xt[jx+1][i]=48; } // '/' stands for -1
           else                        { Xt[jx][i]=49; Xt[jx+1][i]=48; }
           jx++;
           xtQTL[jx]= '1';
        }
        else
        {  for (int i=0; i<Naug; i++)
           if      (marker[j][i]=='1') { Xt[jx][i]=48; }  //ASCII code 47, 48 en 49 voor -1,0,1;
           else if (marker[j][i]=='0') { Xt[jx][i]=47; } // '/' stands for -1
           else                        { Xt[jx][i]=49; }
        }
     }
     else if (cofactor[j]=='3') // QTL
     {  jx++;
        xtQTL[jx]= '2';
        if (dominance=='y')
        {  jx++;
           xtQTL[jx]= '3';
        }
     }

    //Rprintf("calculate xtwx and xtwy\n");
     /* calculate xtwx and xtwy */
     double xtwj, yi, wi, calc_i;
     for (j=0; j<dimx; j++)
     {   XtWY[j]= 0.0;
         for (jj=0; jj<dimx; jj++) XtWX[j][jj]= 0.0;
     }


     if (fitQTL=='n')
     for (int i=0; i<Naug; i++)
     {  yi= y[i];
        wi= (*weight)[i];
        for (j=0; j<dimx; j++)
        {   xtwj= ((double)Xt[j][i]-48.0)*wi;
            XtWY[j]+= xtwj*yi;
            for (jj=0; jj<=j; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
        }
     }
     else // QTL is moving along the chromosomes
     for (int i=0; i<Naug; i++)
     {   wi= (*weight)[i]+ (*weight)[i+Naug]+ (*weight)[i+2*Naug];
         yi= y[i];
         for (j=0; j<=dimx; j++)
         if (xtQTL[j]<='1')
         {  xtwj= ((double)Xt[j][i]-48.0)*wi;
            XtWY[j]+= xtwj*yi;
            for (jj=0; jj<=j; jj++)
            if (xtQTL[jj]<='1') XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
            else if (xtQTL[jj]=='2') // QTL: additive effect if QTL='0' or '2'
            {  // QTL=='0'
               XtWX[j][jj]+= ((double)(Xt[j][i]-48.0))*(*weight)[i]*(47.0-48.0);
               // QTL=='2'
               XtWX[j][jj]+= ((double)(Xt[j][i]-48.0))*(*weight)[i+2*Naug]*(49.0-48.0);
            }
            else // (xtQTL[jj]=='3')  QTL: dominance effect only if QTL='1'
            {  // QTL=='1'
               XtWX[j][jj]+= ((double)(Xt[j][i]-48.0))*(*weight)[i+Naug]*(49.0-48.0);
            }
         }
         else if (xtQTL[j]=='2') // QTL: additive effect if QTL='0' or '2'
         {  xtwj= -1.0*(*weight)[i]; // QTL=='0'
            XtWY[j]+= xtwj*yi;
            for (jj=0; jj<j; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
            XtWX[j][j]+= xtwj*-1.0;
            xtwj= 1.0*(*weight)[i+2*Naug]; // QTL=='2'
            XtWY[j]+= xtwj*yi;
            for (jj=0; jj<j; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
            XtWX[j][j]+= xtwj*1.0;
         }
         else // (xtQTL[j]=='3') QTL: dominance effect only if QTL='1'
         {  xtwj= 1.0*(*weight)[i+Naug]; // QTL=='1'
            XtWY[j]+= xtwj*yi;
            // j-1 is for additive effect, which is orthogonal to dominance effect
            for (jj=0; jj<j-1; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
            XtWX[j][j]+= xtwj*1.0;
         }
     }
     for (j=0; j<dimx; j++)
     for (jj=j+1; jj<dimx; jj++) XtWX[j][jj]= XtWX[jj][j];

     /* solve equations */
    // Rprintf("solve equations\n");
     // printmatrix(XtWX,dimx,dimx);
     int d;
     ivector indx;
     indx= newivector(dimx);
     ludcmp(XtWX,dimx,indx,&d);
     lusolve(XtWX,dimx,indx,XtWY);
     // luinvert(xtwx, inv, dimx, indx);
     // cout << "Parameter estimates" << endl;
     // for (jj=0; jj<dimx; jj++) cout << jj << " " << XtWY[jj] << endl;

     long double *indL;
    // long double indL[Nind];
	 indL = (long double *)Calloc(Nind, long double);
     int newNaug;
     newNaug= (fitQTL=='n' ? Naug : 3*Naug);
     vector fit, resi;
     fit= newvector(newNaug);
     resi= newvector(newNaug);
     // cout << "Calculate residuals" << endl;
     if (*variance<0)
     {    *variance= 0.0;
          if (fitQTL=='n')
          for (int i=0; i<Naug; i++)
          {   fit[i]= 0.0;
              for (j=0; j<dimx; j++)
              fit[i]+=((double)Xt[j][i]-48.0)*XtWY[j];
              resi[i]= y[i]-fit[i];
              *variance += (*weight)[i]*pow(resi[i],2.0);
          }
          else
          for (int i=0; i<Naug; i++)
          {   fit[i]= 0.0; fit[i+Naug]= 0.0; fit[i+2*Naug]= 0.0;
              for (j=0; j<dimx; j++)
              if (xtQTL[j]<='1')
              {   calc_i =((double)Xt[j][i]-48.0)*XtWY[j];
                  fit[i]+= calc_i; fit[i+Naug]+= calc_i; fit[i+2*Naug]+= calc_i;
              }
              else if (xtQTL[j]=='2')
              {   fit[i]+=-1.0*XtWY[j];
                  fit[i+2*Naug]+=1.0*XtWY[j];
              }
              else
              fit[i+Naug]+=1.0*XtWY[j];
              resi[i]= y[i]-fit[i];
              resi[i+Naug]= y[i]-fit[i+Naug];
              resi[i+2*Naug]= y[i]-fit[i+2*Naug];
              *variance +=(*weight)[i]*pow(resi[i],2.0);
              *variance +=(*weight)[i+Naug]*pow(resi[i+Naug],2.0);
              *variance +=(*weight)[i+2*Naug]*pow(resi[i+2*Naug],2.0);
          }
          *variance/= (biasadj=='n' ? Nind : Nind-dimx); // to compare results with Johan; variance/=Nind;
          if (fitQTL=='n')
          for (int i=0; i<Naug; i++) Fy[i]= Lnormal(resi[i],*variance);
          else
          for (int i=0; i<Naug; i++)
          {   Fy[i]       = Lnormal(resi[i],*variance);
              Fy[i+Naug]  = Lnormal(resi[i+Naug],*variance);
              Fy[i+2*Naug]= Lnormal(resi[i+2*Naug],*variance);
          }
     }
     else
     {    if (fitQTL=='n')
          for (int i=0; i<Naug; i++)
          {   fit[i]= 0.0;
              for (j=0; j<dimx; j++)
              fit[i]+=((double)Xt[j][i]-48.0)*XtWY[j];
              resi[i]= y[i]-fit[i];
              Fy[i]  = Lnormal(resi[i],*variance); // ????
          }
          else
          for (int i=0; i<Naug; i++)
          {   fit[i]= 0.0; fit[i+Naug]= 0.0; fit[i+2*Naug]= 0.0;
              for (j=0; j<dimx; j++)
              if (xtQTL[j]<='1')
              {   calc_i =((double)Xt[j][i]-48.0)*XtWY[j];
                  fit[i]+= calc_i; fit[i+Naug]+= calc_i; fit[i+2*Naug]+= calc_i;
              }
              else if (xtQTL[j]=='2')
              {   fit[i]+=-1.0*XtWY[j];
                  fit[i+2*Naug]+=1.0*XtWY[j];
              }
              else
                  fit[i+Naug]+=1.0*XtWY[j];
              resi[i]= y[i]-fit[i];
              resi[i+Naug]= y[i]-fit[i+Naug];
              resi[i+2*Naug]= y[i]-fit[i+2*Naug];
              Fy[i]       = Lnormal(resi[i],*variance);
              Fy[i+Naug]  = Lnormal(resi[i+Naug],*variance);
              Fy[i+2*Naug]= Lnormal(resi[i+2*Naug],*variance);
          }
     }
	Free(fit);
    Free(resi);

    /* calculation of logL */
    // cout << "calculate logL" << endl;
    long double logL=0.0;
    for (int i=0; i<Nind; i++) indL[i]= 0.0;
    if (fitQTL=='n')
    for (int i=0; i<Naug; i++) indL[ind[i]]+=(*weight)[i]*Fy[i];
    else
    for (int i=0; i<Naug; i++)
    {   indL[ind[i]]+=(*weight)[i]*       Fy[i];
        indL[ind[i]]+=(*weight)[i+Naug]*  Fy[i+Naug];
        indL[ind[i]]+=(*weight)[i+2*Naug]*Fy[i+2*Naug];
    }
    for (int i=0; i<Nind; i++){
		//Sum up log likelyhoods for each individual
		logL+= log(indL[i]);
	}

	Free(indL);
    Free(indx);
    Free(xtQTL);
	delmatrix(XtWX,dimx+2);
	delcmatrix(Xt,dimx+2);
	Free(XtWY);    
	return (double)logL;
}

 /* LU decomposition (from Johan via Numerical Recipes in C) Given an n x n matrix a[1..n][1..n], this routine replaces it by the LU
   decomposition of a rowwise permutation of itself.   A and n are input.  a is output. indx[1..n] is an output vector which records the row
   permutation effected by the partial pivoting; d is output as +-1 depending on whether the number of row interchanges was even or odd,
   respectively. This routine is used in combination with lusolve to solve
   linear equations or to invert a matrix.
*/
void ludcmp(matrix m, int dim, ivector ndx, int *d)
{   int r,c,rowmax,i;
    double max,temp,sum;
    vector scale, swap;
    scale= newvector(dim);
    *d=1;
  //  Rprintf("dim: %d, d: %d\n",dim,*d);
    for (r=0; r<dim; r++)
    {   for (max=0.0, c=0; c<dim; c++) if ((temp=fabs(m[r][c])) > max) max=temp;
        if (max==0.0) {warning("Singular matrix.");}
        scale[r]=1.0/max;
    }
    for (c=0; c<dim; c++)
    {   for (r=0; r<c; r++)
        {   for (sum=m[r][c], i=0; i<r; i++) sum-= m[r][i]*m[i][c];
            m[r][c]=sum;
        }
        for (max=0.0, rowmax=c, r=c; r<dim; r++)
        {   for (sum=m[r][c], i=0; i<c; i++) sum-= m[r][i]*m[i][c];
            m[r][c]=sum;
            if ((temp=scale[r]*fabs(sum)) > max) { max=temp; rowmax=r; }
        }
        if (max==0.0) {warning("singular matrix"); }
        if (rowmax!=c)
        {  swap=m[rowmax]; m[rowmax]=m[c]; m[c]=swap;
           scale[rowmax]=scale[c]; (*d)= -(*d);
        }
        ndx[c]=rowmax;
        temp=1.0/m[c][c];
        for (r=c+1; r<dim; r++) m[r][c]*=temp;
    }
    Free(scale);
}

/* Solve the set of n linear equations AX=B.
Here a[1..n][1..n] is input as the LU decomposition of A.
b[1..n] is input as the right hand side vector B, and returns
with the solution vector X.
a, n and indx are not modified by this routine and can be left
for successive calls with different right-hand sides b.
*/
void lusolve(matrix lu, int dim, ivector ndx, vector b)
{    int r,c;
     double sum;
     for (r=0; r<dim; r++)
     {   sum=b[ndx[r]];
         b[ndx[r]]=b[r];
         for (c=0; c<r; c++) sum-= lu[r][c]*b[c];
         b[r]=sum;
     }
     for (r=dim-1; r>-1; r--)
     {   sum=b[r];
         for (c=r+1; c<dim; c++) sum-= lu[r][c]*b[c];
         b[r]=sum/lu[r][r];
     }
}
 
 
 
/* functions gammln, betacf, betai necessary to calculate F(P,df1,df2) */
double gammln(double xx)
{     double x,tmp,ser;
      static double cof[6]={76.18009173, -86.50532033, 24.01409822,
                           -1.231739516, 0.120858003e-2, -0.536382e-5};
      // if (xx<1) cout << "warning: full accuracy only for xx>1; xx= " << xx << endl;
      x=xx-1.0;
      tmp=x+5.5;
      tmp-= (x+0.5)*log(tmp);
      ser=1.0;
      for (int j=0; j<=5; j++) { x+=1.0; ser += cof[j]/x; }
      // delete[] cof;
      return -tmp+log(2.50662827465*ser);
}

double betacf(double a, double b, double x)
{     double qap,qam,qab,em,tem,d,bz,bm=1.0,bp,bpp,az=1.0,am=1.0,ap,app,aold;
      int m;
      qab=a+b;
      qap=a+1.0;
      qam=a-1.0;
      bz=1.0-qab*x/qap;
      for (m=1; m<=100; m++)
      {   em=(double)m;
          tem=em+em;
          d=em*(b-em)*x/((qam+tem)*(a+tem));
          ap=az+d*am;
          bp=bz+d*bm;
          d= -(a+em)*(qab+em)*x/((qap+tem)*(a+tem));
          app=ap+d*az;
          bpp=bp+d*bz;
          aold=az;
          am=ap/bpp;
          bm=bp/bpp;
          az=app/bpp;
          bz=1.0;
          if ( absdouble((az-aold)/az)  < 3.0e-7) return az;
      }
      warning("a or b too big or max number of iterations too small\n");
      return 0.0;
}

double betai(double a, double b, double x)
{     double bt;
      if (x<0.0 || x>1.0) { warning("x not between 0 and 1\n");}
      if (x==0.0 || x==1.0) bt=0.0;
      else bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
      if (x<(a+1.0)/(a+b+2.0)) return bt*betacf(a,b,x)/a;
      else return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double inverseF(int df1, int df2, double alfa)
{      double prob=0.0, minF=0.0, maxF=100.0, halfway=50.0, absdiff=1.0;
       int count=0;
       while ((absdiff>0.001)&&(count<100))
       {     count++;
             halfway= (maxF+minF)/2.0;
             prob= betai(df2/2.0,df1/2.0,df2/(df2+df1*halfway));
             if (prob<alfa) maxF= halfway;
             else minF= halfway;
             absdiff= fabs(prob-alfa);
       }
       Rprintf("prob=%f alfa=%f\n",prob,alfa);
       return halfway;
}

double randomnormal(long *idum)
{     static int iset=0;
      static double gset;
      double fac,rsq,v1,v2;
      if (iset==0)
      {  do
         { v1= 2.0*ran2(idum)-1.0;
           v2= 2.0*ran2(idum)-1.0;
           rsq= v1*v1 + v2*v2;
         } while (rsq>=1.0 || rsq==0.0);
         fac= sqrt(-2.0*log(rsq)/rsq);
         gset= v1*fac;
         iset=1;
         return v2*fac;
      }
      else
      {  iset=0;
         return gset;
      }
}

/* routine generates random numbers (from Numerical Recipes in C) */
double ran2(long *idum)
{     static long idum2=123456789;
      static long iy=0;
      static long iv[32]; //iv[NTAB]
      long IM1=2147483563, IM2=2147483399, IMM1, IA1=40014, IA2=40692;
      long IQ1=53668, IQ2=52774, IR1=12211, IR2=3791, NTAB=32;
      double AM, NDIV, EPS=1.2E-7, RNMX;
      AM= (1.0/IM1);
      IMM1= IM1-1;
      NDIV= (1.0+IMM1)/NTAB;
      RNMX= 1.0-EPS;
      int j;
      long k;
      double temp;
      if (*idum<=0)
      {  if (-(*idum)<1) *idum=1;
         else *idum= -(*idum);
         idum2= (*idum);
         for (j=NTAB+7;j>=0;j--)
         {   k= (*idum)/IQ1;
            *idum= IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum<0) *idum += IM1;
            if (j<NTAB) iv[j] = *idum;
         }
         iy= iv[0];
      }
      k= (*idum)/IQ1;
      *idum= IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum<0) *idum += IM1;
      k= idum2/IQ2;
      idum2= IA2*(idum2-k*IQ2)-k*IR2;
      if (idum2<0) idum2 += IM2;
      j=iy/NDIV;
      iy= iv[j]-idum2;
      iv[j]= *idum;
      if (iy<1) iy+= IMM1;
      // delete[] iv; iv is on the stack
      if ((temp=AM*iy)>RNMX) return RNMX;
      else return temp;
}

void sort2(int n, double *ra, ivector rb)
{    int l,ir,i,j;
     double rra;
     int rrb;
     l=(n>>1)+1;
     ir=n;
     for (;;)
     {   if (l>1) { rra=ra[--l-1]; rrb=rb[l-1]; }
         else
         { rra= ra[ir-1];
           rrb= rb[ir-1];
           ra[ir-1]=ra[0];
           rb[ir-1]=rb[0];
           if (--ir==1)
           { ra[0]=rra; rb[0]=rrb; return; }
         }
         i=l;
         j=l << 1;
         while (j<=ir)
         {   if (j<ir && ra[j-1]<ra[j+1-1]) ++j;
             if (rra < ra[j-1])
             {  ra[i-1]= ra[j-1];
                rb[i-1]= rb[j-1];
                j+= (i=j);
             }
             else j= ir+1;
         }
         ra[i-1] =rra;
         rb[i-1]= rrb;
     }
}

void sort1(int n, double *ra)
{    int l,ir,i,j;
     double rra;
     l=(n>>1)+1;
     ir=n;
     for (;;)
     {   if (l>1) rra=ra[--l-1];
         else
         { rra= ra[ir-1];
           ra[ir-1]=ra[0];
           if (--ir==1) { ra[0]=rra; return; }
         }
         i=l;
         j=l << 1;
         while (j<=ir)
         {     if (j<ir && ra[j-1]<ra[j+1-1]) ++j;
               if (rra < ra[j-1])
               {  ra[i-1]= ra[j-1];
                  j += (i=j);
               }
               else j= ir+1;
         }
         ra[i-1] =rra;
     }
}
 
 
/* end of MQMsupport.c */
