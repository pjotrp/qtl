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
#include "MQMsupport.h"
#include "Regression.h"

//extern long *idum; // for monte carlo simulation or permutation



/*
 * analyseF2 - analyse one F2 family
 *
 */
void analyseF2(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, ivector f1genotype, int Backwards, 
			   double **QTL,vector *mapdistance,int **Chromo,int Nrun,int RMLorML, double windowsize,double stepsize,
			   double stepmin,double stepmax,double alfa,int em,int out_Naug,int **INDlist)
{    
    int Naug;
	int run=0;
    cvector position;
	vector informationcontent;
	char dominance='n';
	char perm_simu='1';
	ivector chr;
	matrix Frun;

	vector r;
	r= newvector(Nmark);
    position= newcvector(Nmark);
	char REMLorML='0';
    char fitQTL='n';
	
	chr= newivector(Nmark);
	Rprintf("Starting MQM analysis\n\n");
    Rprintf("Filling the chromosome matrix\n");
	for(int i=0; i< Nmark; i++){
		chr[i] = Chromo[0][i];
	}

	if(RMLorML == 1){
		REMLorML='1';
	}

	Rprintf("Calculating relative genomepositions of the markers\n");
	for (int j=0; j<Nmark; j++){
        if (j==0)
        { if (chr[j]==chr[j+1]) position[j]='L'; else position[j]='U'; }
        else if (j==(Nmark-1))
        { if (chr[j]==chr[j-1]) position[j]='R'; else position[j]='U'; }
        else if (chr[j]==chr[j-1])
        { if (chr[j]==chr[j+1]) position[j]='M'; else position[j]='R'; }
        else
        { if (chr[j]==chr[j+1]) position[j]='L'; else position[j]='U'; }
    }
    
	Rprintf("Estimating recombinant frequencies\n");	
    for (int j=0; j<Nmark; j++){   
		r[j]= 999.0;
		if ((position[j]=='L')||(position[j]=='M')){
			r[j]= 0.5*(1.0-exp(-0.02*((*mapdistance)[j+1]-(*mapdistance)[j])));
		}
    }
	
	Rprintf("Initialize Frun and informationcontent to 0.0\n");	// ---- Initialize Frun and informationcontent to 0.0
	int Nsteps;
	Nsteps= chr[Nmark-1]*((stepmax-stepmin)/stepsize+1);	
    Frun= newmatrix(Nsteps,Nrun+1);
    informationcontent= newvector(Nsteps);
    for (int i=0; i<Nrun+1; i++){
		for (int ii=0; ii<Nsteps; ii++){
			Frun[ii][i]= 0.0;
		}
	}
    for (int ii=0; ii<Nsteps; ii++){
		informationcontent[ii]= 0.0;
	}

    char dropj='y';
    int jj=0;
   //  Rprintf("any triple of non-segregating markers is considered to be the result of:\n");
   //  Rprintf("identity-by-descent (IBD) instead of identity-by-state (IBS)\n");
  //   Rprintf("no (segregating!) cofactors are fitted in such non-segregating IBD regions\n");
    for (int j=0; j<Nmark; j++){
		if (mod(f1genotype[j],11)!=0){
			dropj='n';
		}else if (cofactor[j]=='0'){
			dropj='y';
		}else if (position[j]=='L'){
			// (cofactor[j]!='0') cofactor at non-segregating marker
			// test whether next segregating marker is nearby (<20cM)
			dropj='y';
            if ((((*mapdistance)[j+1]-(*mapdistance)[j])<20)&&(mod(f1genotype[j+1],11)!=0)) dropj='n';
            else if (position[j+1]!='R')
            if ((((*mapdistance)[j+2]-(*mapdistance)[j])<20)&&(mod(f1genotype[j+2],11)!=0)) dropj='n';
        }else if (position[j]=='M'){
			dropj='y';
            if ((((*mapdistance)[j]-(*mapdistance)[j-1])<20)&&(mod(f1genotype[j-1],11)!=0)) dropj='n';
            else if ((((*mapdistance)[j+1]-(*mapdistance)[j])<20)&&(mod(f1genotype[j+1],11)!=0)) dropj='n';
        }else if (position[j]=='R'){
			dropj='y';
            if ((((*mapdistance)[j]-(*mapdistance)[j-1])<20)&&(mod(f1genotype[j-1],11)!=0)) dropj='n';
            else if (position[j-1]!='L')
            if ((((*mapdistance)[j]-(*mapdistance)[j-2])<20)&&(mod(f1genotype[j-2],11)!=0)) dropj='n';
        }
		if (dropj=='n'){  
            marker[jj]= marker[j];
            cofactor[jj]= cofactor[j];
            (*mapdistance)[jj]= (*mapdistance)[j];
            chr[jj]= chr[j];
            r[jj]= r[j];
            position[jj]= position[j];
            jj++;
        }else if (cofactor[j]=='1'){  
            Rprintf("cofactor at chr %d is dropped\n",chr[j]);
        }
    }
    Nmark= jj;
    for (int j=0; j<Nmark; j++){
		r[j]= 999.0;
        if (j==0)
        { if (chr[j]==chr[j+1]) position[j]='L'; else position[j]='U'; }
        else if (j==(Nmark-1))
        { if (chr[j]==chr[j-1]) position[j]='R'; else position[j]='U'; }
        else if (chr[j]==chr[j-1])
        { if (chr[j]==chr[j+1]) position[j]='M'; else position[j]='R'; }
        else
        { if (chr[j]==chr[j+1]) position[j]='L'; else position[j]='U'; }
    }
    for (int j=0; j<Nmark; j++){
		if ((position[j]=='L')||(position[j]=='M')){
			r[j]= 0.5*(1.0-exp(-0.02*((*mapdistance)[j+1]-(*mapdistance)[j])));
			if (r[j]<0){
				Rprintf("error: recombination frequency is negative\n");
				Rprintf("chr=%d mapdistance=%d\n",chr[j],(*mapdistance)[j]); 
				Rprintf("position=%d r[j]=%d\n",position[j], r[j]);
				return;
			}
		}
    }

    ivector newind;
    vector newy;
    cmatrix newmarker;
    double ymean=0.0, yvari=0.0;
    for (int i=0; i<Nind; i++) ymean += y[i];
    ymean/= Nind;
    for (int i=0; i<Nind; i++) yvari += pow(y[i]-ymean,2);
    yvari/= (Nind-1);
    Rprintf("ymean=%f yvari=%f\n",ymean,yvari);
	
	//Fix for not doing dataaugmentation, we just copy the current as the augmented and set Naug to Nind
	Naug=Nind;
	Nind=out_Naug;
	newind= newivector(Naug);
	newy= newvector(Naug);
	newmarker= newcmatrix(Nmark,Naug);
    for (int i=0; i<Naug; i++){
		newy[i]= y[i];
        newind[i]= INDlist[0][i];
        for (int j=0; j<Nmark; j++){
			newmarker[j][i]= marker[j][i];
		}
    }
//	 if(augdata(marker,y,&newmarker,&newy,&newind,&Nind,&Naug,Nmark,position,r,maxNaug,imaxNaug,neglect)==1){
//		Rprintf("Data augmentation finished succesfull\n");
//		Rprintf("# Unique individuals before augmentation:%d\n",prior);
//		Rprintf("# Unique selected individuals:%d\n",Nind);
//		Rprintf("# Marker p individual:%d\n",Nmark);
//		Rprintf("# Individuals after augmentation:%d\n",Naug);
//	}else{
//		Rprintf("Data augmentation failed\n");
//		return;
//	}
//    Rprintf("Naug:%d\n",Naug);

    vector newweight;
    newweight= newvector(Naug);
    rmixture(newmarker, newweight, r, position, newind,Nind, Naug, Nmark);
	//for (int j=0; j<Nmark; j++){
	//	Rprintf("r(%d)= %f\n",j,r[j]);
	//}
    /* eliminate individuals with missing trait values */
    //We can skip this part iirc because R throws out missing phenotypes beforehand
	int oldNind=Nind;
    for (int i=0; i<oldNind; i++){
		Nind-= ((y[i]==999.0) ? 1 : 0);
	}
   
    int oldNaug=Naug;
    for (int i=0; i<oldNaug; i++){
		Naug-= ((newy[i]==999.0) ? 1 : 0);
	}
    
	vector weight;
    ivector ind;
    marker= newcmatrix(Nmark,Naug);
    y= newvector(Naug);
    ind= newivector(Naug);
    weight= newvector(Naug);
    int newi=0;
    for (int i=0; i<oldNaug; i++)
    if (newy[i]!=999.0){
		y[newi]= newy[i];
        ind[newi]= newind[i];
        weight[newi]= newweight[i];
        for (int j=0; j<Nmark; j++) marker[j][newi]= newmarker[j][i];
        newi++;
    }
    int diff;
    for (int i=0; i<Naug-1; i++){
		diff=ind[i+1]-ind[i];
        if  (diff>1)
        for (int ii=i+1; ii<Naug; ii++) ind[ii]=ind[ii]-diff+1;
    }
    delcmatrix(newmarker,Nmark);
    Free(newy);
    Free(newind);
    Free(newweight);

 //    vector Fy;
 //    Fy= newvector(Naug);
    double variance=-1.0;
	double logLfull;
    cvector selcofactor;
    selcofactor= newcvector(Nmark); /* selected cofactors */

    int dimx=1;
    for (int j=0; j<Nmark; j++){
		if (cofactor[j]=='1'){
			dimx+= (dominance=='n' ? 1 : 2);  // per QTL only additivity !!
		}else if (cofactor[j]=='2'){
			dimx+=1;  /* sex of the mouse */
		}
	}
     double F1, F2;
     F1= inverseF(1,Nind-dimx,alfa);
     F2= inverseF(2,Nind-dimx,alfa);
     Rprintf("F(%f,1,%d)=%f\n",F1,(Nind-dimx),alfa);
     Rprintf("F(%f,2,%d)=%f\n",F2,(Nind-dimx),alfa);
     F2= 2.0* F2; // 9-6-1998 using threshold x*F(x,df,alfa)

     weight[0]= -1.0;
     logLfull= QTLmixture(marker,cofactor,r,position,y,ind,Nind,Naug,Nmark,&variance,em,&weight,REMLorML,fitQTL,dominance);
     Rprintf("log-likelihood of full model= %f\n",logLfull);
     Rprintf("residual variance= %f\n",variance);
     Rprintf("Trait mean= %f \nTrait variation= %f\n",ymean,yvari);

     if (Backwards==1)    // use only selected cofactors
         logLfull= backward(Nind, Nmark, cofactor, marker, y, weight, ind, Naug, logLfull,
                    variance, F1, F2, &selcofactor, r, position, &informationcontent, mapdistance,&Frun,run,REMLorML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax);
     if (Backwards==0) // use all cofactors
         logLfull= mapQTL(Nind, Nmark, cofactor, cofactor, marker, position,
                  (*mapdistance), y, r, ind, Naug, variance, 'n', &informationcontent,&Frun,run,REMLorML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax); // printout=='n'

	long *idum;
    idum = (long *)Calloc(1, long*);
    idum[0]=-1;

	 double savevariance= variance;
     double *urand;
     vector maxF;
	 maxF= newvector(Nrun);
     ivector indorder; // individu 0...Nind-1; order will be permuted
     vector yoriginal;
     indorder= newivector(Nind);
     yoriginal= newvector(Nind);
	 urand= newvector(Naug);
	// printf("Gonna start bootstrapping??? Nrun:%d\n",Nrun);
	//This is gonna be removed and put into R
     if (Nrun>0){  
        urand[0]= ran2(idum);
        for (int i=0; i<Naug; i++) yoriginal[ind[i]]= y[i];

        for (run=1; run<Nrun; run++)
        {   
		    R_CheckUserInterrupt(); /* check for ^C */
			Rprintf("Run = %d\n",run);
            if (perm_simu=='0')
            {  for (int i=0; i<Nind; i++) indorder[i]= i;
               for (int i=0; i<Nind; i++) urand[i]= ran2(idum);
               sort2(Nind,urand,indorder);  // y[individu 0...Nind-1] are permuted
               for (int i=0; i<Naug; i++) // indorder[ind[i]]== new individu number
                   y[i]= yoriginal[indorder[ind[i]]];
            }
            else
            {  // cout << "Parametric bootstrap" << endl;
               for (int i=0; i<Nind; i++) yoriginal[i]= pow(savevariance,0.5)*randomnormal(idum);
               for (int i=0; i<Naug; i++) y[i]= yoriginal[ind[i]];
            }
            variance= -1.0;
            // logLfull= QTLmixture(marker,cofactor,r,position,y,ind,Nind,Naug,Nmark,variance,em,weight);
			if (Backwards==1){ 
                  maxF[run]= backward(Nind, Nmark, cofactor, marker, y, weight, ind, Naug, logLfull,
                    variance, F1, F2, &selcofactor, r, position,&informationcontent, mapdistance,&Frun,run,REMLorML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax);
			}
			if (Backwards==0){
                  maxF[run]= mapQTL(Nind, Nmark, cofactor, cofactor, marker, position,
                  (*mapdistance), y, r, ind, Naug, variance, 'n',&informationcontent,&Frun,run,REMLorML, fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax);
			}
            // cout << "run " << run <<" ready; maxF= " << maxF[run] << endl;
        }
	}
	//administration and such (should be also)... however we do need the Cumulative distribution of maximum test statistic value
    if (Nrun > 0){
		sort1(Nrun,maxF);
		Rprintf("Cumulative distribution of maximum test statistic value in %d permutations or simulations\n",Nrun);
		for (int i=1; i<Nrun+1; i++){   
			Rprintf(" %f %f\n",( (double)i/( (double)Nrun+1.0) ),maxF[i-1]);
			R_ProcessEvents();
			R_FlushConsole();
        }
    }	
	Rprintf("Analysis of data finished\n");
  // ---- Write output / send it back to R
	double moveQTL= stepmin;
	int chrnumber=1;
  
	Rprintf("-1- %d %d\n",Nsteps,Nrun);
	// chr pos Frun    information 
	// 1  -20  97.4561 0.677204
	// 1  -15 103.29   0.723067
	// 1  -10 108.759  0.777696
	// 1   -5 113.737  0.842778
	// 1    0 118.112  0.920356
	// 1    5 120.051  0.928594
	// 1   10 114.469  0.959548
	
	//Printout output to QTL for usage in R
	//we want the first run we did
	for (int ii=0; ii<Nsteps; ii++){   
		QTL[0][ii] = Frun[ii][0];
		if (moveQTL+stepsize<=stepmax){
			moveQTL+= stepsize;
		} else { 
			moveQTL= stepmin; 
			chrnumber++; 
		}
    }
	Free(urand);
	Free(indorder);
	Free(yoriginal);
	Free(maxF);
	Free(position);
	Free(weight);
	Free(ind);
	delcmatrix(marker,Nmark);
	Free(y);
	Free(selcofactor);
	Free(idum);
	return;
}

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

/* ML estimation of recombination frequencies via EM;
   calculation of multilocus genotype probabilities;
   ignorance of unlikely genotypes*/
void rmixture(cmatrix marker, vector weight, vector r,
              cvector position, ivector ind,
              int Nind, int Naug, int Nmark){   
	int i,j;
    int iem= 0;
    double Nrecom, oldr=0.0, newr, rdelta=1.0;
    vector indweight;
    indweight = newvector(Nind);
    char rknown='n';
	//this should change because we should re-estimate the Recombination frequency for all r[j] != 999.0
	//this then effects to mapQTL and ends because of r[j] errors... perhaps r[j] should come from R/QTL so we don't have to worry about those calculations
    for (j=0; j<Nmark; j++){
		//Rprintf("Recombination frequency: %f at marker %d\n",r[j],j);
		if (r[j]!=999.0){
			rknown='y';
		}
	}
    if (rknown=='y'){
		Rprintf("recombination parameters are not re-estimated\n");
		//rknown='n'; //HAX to make it update recombination frequenties
    }
	//We should use other values then set ones.... em.iter here ?
     while ((iem<100)&&(rdelta>0.001))
     {     R_CheckUserInterrupt(); /* check for ^C */
		   R_ProcessEvents(); /* do some windows/C stuff so R doesn't look so unresponsive */
		  // R_FlushConsole();
		   iem+=1;
           rdelta= 0.0;
           /* calculate weights = conditional genotype probabilities */
           for (i=0; i<Naug; i++) weight[i]=1.0;
           for (j=0; j<Nmark; j++)
           {   if ((position[j]=='L')||(position[j]=='U'))
               for (i=0; i<Naug; i++)
               if (marker[j][i]=='1') weight[i]*= 0.5;
               else weight[i]*= 0.25;
               if ((position[j]=='L')||(position[j]=='M'))
               for (i=0; i<Naug; i++)
               {   Nrecom= absdouble((double)marker[j][i]-marker[j+1][i]);
                   if ((marker[j][i]=='1')&&(marker[j+1][i]=='1'))
                      weight[i]*= (r[j]*r[j]+(1.0-r[j])*(1.0-r[j])); // /2.0;
                   else if (Nrecom==0) weight[i]*= (1.0-r[j])*(1.0-r[j]);
                   else if (Nrecom==1) weight[i]*= ((marker[j+1][i]=='1') ? 2.0*r[j]*(1.0-r[j]) :
                                                                                r[j]*(1.0-r[j]));
                   else                weight[i]*=      r[j] *     r[j];
               }
           }
           for (i=0; i<Nind; i++){ 
               indweight[i]= 0.0;
           }   
           for (i=0; i<Naug; i++){
               indweight[ind[i]]+=weight[i];
           }       
           for (i=0; i<Naug; i++){ 
               weight[i]/=indweight[ind[i]];
           }
           for (j=0; j<Nmark; j++)
           {   if ((position[j]=='L')||(position[j]=='M'))
               {  newr= 0.0;
                  for (i=0; i<Naug; i++)
                  {   Nrecom= absdouble((double)marker[j][i]-marker[j+1][i]);
                      if ((marker[j][i]=='1')&&(marker[j+1][i]=='1'))
                         Nrecom= 2.0*r[j]*r[j]/(r[j]*r[j]+(1-r[j])*(1-r[j]));
                      newr+= Nrecom*weight[i];
                  }
                  if (rknown=='n' && position[j]!='R') //only update if it isn't the last marker of a chromosome ;)
                  {  oldr=r[j];
                     r[j]= newr/(2.0*Nind);
                     rdelta+=pow(r[j]-oldr,2.0);
                  }
                  else rdelta+=0.0;
               }
            }
     }
     
/*   print new estimates of recombination frequencies */
    Rprintf("iem= %d rdelta= %f\n",iem,rdelta);
   // if (rknown=='n'){  
   //     for (j=0; j<Nmark; j++){
	//		Rprintf("r(%d)= %f\n",j,r[j]);
	//	}
   // }
	Free(indweight);
}


/* ML estimation of parameters in mixture model via EM;
*/
double QTLmixture(cmatrix loci, cvector cofactor, vector r, cvector position,
              vector y, ivector ind, int Nind, int Naug,
              int Nloci,
              double *variance, int em, vector *weight,char REMLorML,char fitQTL,char dominance)
{  //  Rprintf("QTLmixture called\n");
    int iem= 0, newNaug, i, j;
    char varknown, biasadj='n';
	double Nrecom, oldlogL=-10000, delta=1.0, calc_i, logP=0.0, Pscale=1.75;
    
	vector indweight, Ploci, Fy;
    
	indweight= newvector(Nind);
    newNaug= (fitQTL=='n' ? Naug : 3*Naug);
    Ploci= newvector(newNaug);
    Fy= newvector(newNaug);
    logP= Nloci*log(Pscale); // only for computational accuracy
	varknown= (((*variance)==-1.0) ? 'n' : 'y' );
	
    if ((REMLorML=='0')&&(varknown=='n')){ 
		Rprintf("variance is being estimated and bias adjusted\n");
	}
    if (REMLorML=='1') { 
		varknown='n'; biasadj='n'; 
	}
    for (i=0; i<newNaug; i++){ 
		Ploci[i]= 1.0;
	}
    if (fitQTL=='n'){
	//Rprintf("FitQTL=N\n");	
		for (j=0; j<Nloci; j++){
		    for (i=0; i<Naug; i++) 
			Ploci[i]*= Pscale;
		    if ((position[j]=='L')||(position[j]=='U')){
				for (i=0; i<Naug; i++) Ploci[i]*= (loci[j][i]=='1' ? 0.5 : 0.25);
			}
		    if ((position[j]=='L')||(position[j]=='M')){
				for (i=0; i<Naug; i++){
					Nrecom= absdouble((double)loci[j][i]-(double)loci[j+1][i]);
					if ((loci[j][i]=='1')&&(loci[j+1][i]=='1')){
						calc_i= (r[j]*r[j]+(1.0-r[j])*(1.0-r[j]));}
					else if (Nrecom==0) {
						calc_i= (1.0-r[j])*(1.0-r[j]);
					}else if (Nrecom==1) {
						calc_i= ((loci[j+1][i]=='1') ? 2.0*r[j]*(1.0-r[j]) : r[j]*(1.0-r[j]));
					}else {
						calc_i= r[j]*r[j];
					}
					Ploci[i]*= calc_i;
				}
			}
		}
	//Rprintf("FitQTL=N Done\n");
	//for (j=0; j<Nloci; j++){
	//    for (i=0; i<Naug; i++){
	//		Rprintf("%c ",loci[j][i]);
	//	}
	//	Rprintf("\n");
	//}
	}else{
	//Rprintf("FitQTL=Y\n");	
     for (j=0; j<Nloci; j++)
     {    for (i=0; i<Naug; i++)
          {   Ploci[i]*= Pscale; Ploci[i+Naug]*= Pscale; Ploci[i+2*Naug]*= Pscale;
              // only for computational accuracy; see use of logP
          }
          if ((position[j]=='L')||(position[j]=='U'))
          {  if (cofactor[j]<='1')
             for (i=0; i<Naug; i++)
             {   calc_i= (loci[j][i]=='1' ? 0.5 : 0.25);
                 Ploci[i]*= calc_i; Ploci[i+Naug]*= calc_i; Ploci[i+2*Naug]*= calc_i;
             }
             else
             for (i=0; i<Naug; i++)
             {   Ploci[i]*= 0.25; Ploci[i+Naug]*= 0.5; Ploci[i+2*Naug] *= 0.25; }
                 // QTL='0', '1' or'2'
          }
          if ((position[j]=='L')||(position[j]=='M'))
          {  if ((cofactor[j]<='1')&&(cofactor[j+1]<='1'))
             for (i=0; i<Naug; i++)
             {  Nrecom= absdouble((double)loci[j][i]-(double)loci[j+1][i]);
                if ((loci[j][i]=='1')&&(loci[j+1][i]=='1'))
                   calc_i= (r[j]*r[j]+(1.0-r[j])*(1.0-r[j]));
                else if (Nrecom==0) calc_i= (1.0-r[j])*(1.0-r[j]);
                else if (Nrecom==1) calc_i= ((loci[j+1][i]=='1')
                        ? 2.0*r[j]*(1.0-r[j]) : r[j]*(1.0-r[j]));
                else calc_i= r[j]*r[j];
                Ploci[i]*= calc_i; Ploci[i+Naug]*= calc_i; Ploci[i+2*Naug]*= calc_i;
             }
             else if (cofactor[j]<='1') // locus j+1 == QTL
             for (i=0; i<Naug; i++)
             {  // QTL=='0'
                Nrecom= absdouble((double)loci[j][i]-(double)'0');
                if (Nrecom==0) calc_i= (1.0-r[j])*(1.0-r[j]);
                else if (Nrecom==1) calc_i= r[j]*(1.0-r[j]);
                else calc_i= r[j]*r[j];
                Ploci[i]*= calc_i;
                // QTL=='1'
                Nrecom= absdouble((double)loci[j][i]-(double)'1');
                if (loci[j][i]=='1')
                   calc_i= (r[j]*r[j]+(1.0-r[j])*(1.0-r[j]));
                else if (Nrecom==0) calc_i= (1.0-r[j])*(1.0-r[j]);
                else if (Nrecom==1) calc_i= 2.0*r[j]*(1.0-r[j]);
                else calc_i= r[j]*r[j];
                Ploci[i+Naug]*= calc_i;
                // QTL=='2'
                Nrecom= absdouble((double)loci[j][i]-(double)'2');
                if (Nrecom==0) calc_i= (1.0-r[j])*(1.0-r[j]);
                else if (Nrecom==1) calc_i= r[j]*(1.0-r[j]);
                else calc_i= r[j]*r[j];
                Ploci[i+2*Naug]*= calc_i;
             }
             else // locus j == QTL
             for (i=0; i<Naug; i++)
             {  // QTL=='0'
                Nrecom= absdouble((double)loci[j+1][i]-(double)'0');
                if (Nrecom==0) calc_i= (1.0-r[j])*(1.0-r[j]);
                else if (Nrecom==1) calc_i= ((loci[j+1][i]=='1')
                        ? 2.0*r[j]*(1.0-r[j]) : r[j]*(1.0-r[j]));
                else calc_i= r[j]*r[j];
                Ploci[i]*= calc_i;
                //test= 0;
                //test+= (calc_i==0 ? 1 : 0);
                // QTL=='1'
                Nrecom= absdouble((double)loci[j+1][i]-(double)'1');
                if (loci[j+1][i]=='1')
                   calc_i= (r[j]*r[j]+(1.0-r[j])*(1.0-r[j]));
                else if (Nrecom==0) calc_i= (1.0-r[j])*(1.0-r[j]);
                else if (Nrecom==1) calc_i= ((loci[j+1][i]=='1')
                        ? 2.0*r[j]*(1.0-r[j]) : r[j]*(1.0-r[j]));
                else calc_i= r[j]*r[j];
                Ploci[i+Naug]*= calc_i;
                //test+= (calc_i==0 ? 1 : 0);
                // QTL=='2'
                Nrecom= absdouble((double)loci[j+1][i]-(double)'2');
                if (Nrecom==0) calc_i= (1.0-r[j])*(1.0-r[j]);
                else if (Nrecom==1) calc_i= ((loci[j+1][i]=='1')
                        ? 2.0*r[j]*(1.0-r[j]) : r[j]*(1.0-r[j]));
                else calc_i= r[j]*r[j];
                Ploci[i+2*Naug]*= calc_i;
             }
          }
	 }
	 }
	// Rprintf("fitQTL's done\n");
     if ((*weight)[0]== -1.0)
     {  for (i=0; i<Nind; i++) indweight[i]= 0.0;
		if (fitQTL=='n')
        {  for (i=0; i<Naug; i++) indweight[ind[i]]+=Ploci[i];
           for (i=0; i<Naug; i++) (*weight)[i]= Ploci[i]/indweight[ind[i]];
        }
        else
        {  for (i=0; i<Naug; i++) indweight[ind[i]]+=Ploci[i]+Ploci[i+Naug]+Ploci[i+2*Naug];
           for (i=0; i<Naug; i++)
           {   (*weight)[i]       = Ploci[i]/indweight[ind[i]];
               (*weight)[i+Naug]  = Ploci[i+Naug]/indweight[ind[i]];
               (*weight)[i+2*Naug]= Ploci[i+2*Naug]/indweight[ind[i]];
           }
        }
     }
	// Rprintf("Weights done\n");
     //Rprintf("Individual->trait->cofactor->weight\n");
    // for (int j=0; j<Nind; j++){
	//    Rprintf("%d->%f,%d,%f\n",j,y[j],cofactor[j],(*weight)[j]);
	 //}	
     double logL=0;
     vector indL;
     indL= newvector(Nind);
     while ((iem<em)&&(delta>1.0e-5))
     {  
		R_CheckUserInterrupt(); /* check for ^C */
		R_ProcessEvents();
		//R_FlushConsole();
           iem+=1;
           if (varknown=='n') *variance=-1.0;
        //   Rprintf("Checkpoint_b\n");           
           logL= regression(Nind, Nloci, cofactor, loci, y,
                 weight, ind, Naug, variance, Fy, biasadj,fitQTL,dominance);
           logL=0.0;
        //   Rprintf("regression ready\n");
           for (i=0; i<Nind; i++) indL[i]= 0.0;
           if (fitQTL=='n') // no QTL fitted
           for (i=0; i<Naug; i++)
           {   (*weight)[i]= Ploci[i]*Fy[i];
               indL[ind[i]]= indL[ind[i]] + (*weight)[i];
           }
           else // QTL moved along the chromosomes
           for (i=0; i<Naug; i++)
           {  (*weight)[i]= Ploci[i]*Fy[i];
              (*weight)[i+Naug]  = Ploci[i+Naug]*  Fy[i+Naug];
              (*weight)[i+2*Naug]= Ploci[i+2*Naug]*Fy[i+2*Naug];
              indL[ind[i]]+=(*weight)[i]+(*weight)[i+Naug]+(*weight)[i+2*Naug];
           }
           for (i=0; i<Nind; i++) logL+=log(indL[i])-logP;
           for (i=0; i<Nind; i++) indweight[i]= 0.0;
           if (fitQTL=='n')
           {  for (i=0; i<Naug; i++) indweight[ind[i]]+=(*weight)[i];
              for (i=0; i<Naug; i++) (*weight)[i]/=indweight[ind[i]];
           }
           else
           {  for (i=0; i<Naug; i++)
                  indweight[ind[i]]+=(*weight)[i]+(*weight)[i+Naug]+(*weight)[i+2*Naug];
              for (i=0; i<Naug; i++)
              {   (*weight)[i]       /=indweight[ind[i]];
                  (*weight)[i+Naug]  /=indweight[ind[i]];
                  (*weight)[i+2*Naug]/=indweight[ind[i]];
              }
           }
           delta= absdouble(logL-oldlogL);
           oldlogL= logL;
     }
  //   Rprintf("EM Finished\n");
     // bias adjustment after finished ML estimation via EM
     if ((REMLorML=='0')&&(varknown=='n'))
     {  
       // Rprintf("Checkpoint_c\n");
        *variance=-1.0;
        biasadj='y';
        logL= regression(Nind, Nloci, cofactor, loci, y,
              weight, ind, Naug, variance, Fy, biasadj,fitQTL,dominance);
        logL=0.0;
        for (int _i=0; _i<Nind; _i++) indL[_i]= 0.0;
        if (fitQTL=='n')
        for (i=0; i<Naug; i++)
        {   (*weight)[i]= Ploci[i]*Fy[i];
            indL[ind[i]]+=(*weight)[i];
        }
        else
        for (i=0; i<Naug; i++)
        {   (*weight)[i]= Ploci[i]*Fy[i];
            (*weight)[i+Naug]= Ploci[i+Naug]*Fy[i+Naug];
            (*weight)[i+2*Naug]= Ploci[i+2*Naug]*Fy[i+2*Naug];
            indL[ind[i]]+=(*weight)[i];
            indL[ind[i]]+=(*weight)[i+Naug];
            indL[ind[i]]+=(*weight)[i+2*Naug];
        }
        for (i=0; i<Nind; i++) logL+=log(indL[i])-logP;
        for (i=0; i<Nind; i++) indweight[i]= 0.0;
        if (fitQTL=='n')
        {  for (i=0; i<Naug; i++) indweight[ind[i]]+=(*weight)[i];
           for (i=0; i<Naug; i++) (*weight)[i]/=indweight[ind[i]];
        }
        else
        {  for (i=0; i<Naug; i++)
           {   indweight[ind[i]]+=(*weight)[i];
               indweight[ind[i]]+=(*weight)[i+Naug];
               indweight[ind[i]]+=(*weight)[i+2*Naug];
           }
           for (i=0; i<Naug; i++)
           {   (*weight)[i]       /=indweight[ind[i]];
               (*weight)[i+Naug]  /=indweight[ind[i]];
               (*weight)[i+2*Naug]/=indweight[ind[i]];
           }
        }
     }
	// for (i=0; i<Nind; i++) Rprintf("IND %d Ploci: %f Fy: %f UNLOG:%f LogL:%f LogL-LogP: %f\n",i,Ploci[i],Fy[i],indL[i],log(indL[i]),log(indL[i])-logP);
	Free(Fy);
	Free(Ploci);
	Free(indweight);
	Free(indL);
    return logL;
}

/* backward elimination in regression of trait on multiple cofactors
   routine subX haalt uit matrices voor volledige model de submatrices voor submodellen;
   matrices XtWX en Xt van volledig model worden genoemd fullxtwx en fullxt;
   analoog vector XtWY wordt full xtwy genoemd;
*/
double backward(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, vector weight, int* ind, int Naug, double logLfull, double variance, double F1, double F2, cvector* newcofactor, vector r, cvector position,vector *informationcontent,vector *mapdistance,matrix *Frun,int run,char REMLorML,char fitQTL,char dominance,int em, double windowsize,double stepsize,
			  double stepmin,double stepmax)
{    int dropj=0, Ncof=0;
     double maxlogL, savelogL, maxF=0.0; //, minlogL=logLfull, maxFtest=0.0;
     char finished='n'; //, biasadj='n';
     vector logL;
     logL = newvector(Nmark);
     savelogL= logLfull;
     maxlogL= logLfull-10000;
	 //Rprintf("Backward started\n");
     for (int j=0; j<Nmark; j++)
     {   (*newcofactor)[j]= cofactor[j];
         Ncof+=(cofactor[j]!='0');
     }
     while ((Ncof>0)&&(finished=='n'))
     {     for (int j=0; j<Nmark; j++)
           {   if ((*newcofactor)[j]=='1')
               {  //Rprintf("Drop marker %d\n",j);
                  (*newcofactor)[j]='0';
                  if (REMLorML=='1') variance= -1.0;
                  logL[j]= QTLmixture(marker,(*newcofactor),r,position,y,ind,Nind,Naug,Nmark,&variance,em,&weight,REMLorML,fitQTL,dominance);
                  (*newcofactor)[j]='1';
               }
               else if ((*newcofactor)[j]=='2')
               {  //Rprintf("Drop marker %d\n",j);
                  (*newcofactor)[j]='0';
                  if (REMLorML=='1') variance= -1.0;
                  logL[j]=  QTLmixture(marker,(*newcofactor),r,position,y,ind,Nind,Naug,Nmark,&variance,em,&weight,REMLorML,fitQTL,dominance);
                  (*newcofactor)[j]='2';
               }
               else if ((*newcofactor)[j]!='0') Rprintf("Something is wrong");
           }
           /* nu bepalen welke cofactor 0 kan worden (=verwijderd) */
           maxlogL= logLfull-10000.0;
           for (int j=0; j<Nmark; j++)
           {   if ((*newcofactor)[j]!='0')
               if (logL[j]>maxlogL) { maxlogL= logL[j]; dropj = j; }
           }
           if  ( ((*newcofactor)[dropj]=='1') && ( F2> 2.0*(savelogL-maxlogL)) )
           {   savelogL= maxlogL;
               (*newcofactor)[dropj]= '0'; Ncof-=1;
               Rprintf("Marker %d is dropped, resulting in logL of reduced model = %f\n",(dropj+1),savelogL);
           }
           else if  ( ((*newcofactor)[dropj]=='2') && (F1> 2.0*(savelogL-maxlogL)) )
           {   savelogL= maxlogL;
               (*newcofactor)[dropj]= '0'; Ncof-=1;
               Rprintf("marker %d is dropped, resulting in logL of reduced model = %f\n",(dropj+1),savelogL);
           }
           else /* ready */
           {   finished='y';
               for (int j=0; j<Nmark; j++)
               if ((*newcofactor)[j]=='1') Rprintf("Marker %d is selected\n",(j+1));
           }
     }
     for (int j=0; j<Nmark; j++)
     if ((*newcofactor)[j]!='0') Rprintf("Marker %d is in final model\n",(j+1));

     maxF= mapQTL(Nind, Nmark, cofactor, (*newcofactor), marker, position,
           (*mapdistance), y, r, ind, Naug, variance, 'n', informationcontent,Frun,run,REMLorML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax); // printoutput='n'
     //Rprintf("Backward selection finished\n");
     Free(logL);
     return maxF;
}

/* mapQTL */
double mapQTL(int Nind, int Nmark, cvector cofactor, cvector selcofactor, cmatrix marker, cvector position, vector mapdistance, vector y, 
			  vector r, ivector ind, int Naug, double variance, char printoutput,vector *informationcontent,matrix *Frun,int run,char REMLorML,char fitQTL,char dominance,int em, double windowsize,double stepsize,
			  double stepmin,double stepmax)
{      
       Rprintf("mapQTL function called\n");
       int Nloci, j, jj, jjj=0;
       vector Fy;
       Fy= newvector(Naug);
       cvector QTLcofactor, saveQTLcofactor;
       QTLcofactor= newcvector(Nmark+1);
       saveQTLcofactor= newcvector(Nmark+1);
       double infocontent;
       vector info0, info1, info2, weight;
       info0= newvector(Nind);
       info1= newvector(Nind);
       info2= newvector(Nind);
       weight= newvector(Naug);
       weight[0]= -1.0;

       /* fit QTL on top of markers (but: should also be done with routine QTLmixture()
       for exact ML) */

       cvector newcofactor;
       newcofactor= newcvector(Nmark);
       vector cumdistance;
       double QTLlikelihood=0.0;
       cumdistance= newvector(Nmark+1);
       for (j=0; j<Nmark; j++)
       {   if (position[j]=='L')
              cumdistance[j]= -50*log(1-2.0*r[j]);
           else if (position[j]=='M')
              cumdistance[j]= cumdistance[j-1]-50*log(1-2.0*r[j]);
       }
       double savelogL=999.0; // log-likelihood of model with all selected cofactors


      // ofstream fff("mqm_out.txt", ios::out | ios::app);
      // cout << endl << endl;
	//   Rprintf("DEBUG testing_1");
       /* fit QTL on top of markers (full ML)   fit QTL between markers (full ML) */
       // cout << "please wait (mixture calculus may take quite a lot of time)" << endl;
       /* estimate variance in mixture model with all marker cofactors */
       // cout << "estimate variance in mixture model with all cofactors" << endl;
       variance= -1.0;
       savelogL= 2.0*QTLmixture(marker,cofactor,r,position, y,ind,Nind,Naug,Nmark,&variance,em,&weight,REMLorML,fitQTL,dominance);
	   Rprintf("log-likelihood of full model= %f\n",savelogL/2);
       Nloci= Nmark+1;
       // augment data for missing QTL observations (x 3)
       fitQTL='y';
       int newNaug;
       newNaug= 3*Naug;
	   Free(weight);
       weight= newvector(newNaug);
       weight[0]= 1.0;
       vector weight0;
       weight0= newvector(newNaug);
       weight0[0]= -1.0;
    //   Rprintf("DEBUG testing_2");
//       augmentdataforQTL(marker);
       vector QTLr, QTLmapdistance;
       QTLr= newvector(Nloci);
       QTLmapdistance= newvector(Nloci);
       cvector QTLposition;
       QTLposition= newcvector(Nloci);
       cmatrix QTLloci;
	   QTLloci = (char **)R_alloc(Nloci, sizeof(char *));
       double moveQTL= stepmin;
       char nextinterval= 'n', firsttime='y';
       double maxF=0.0, savebaseNoQTLModel=0.0;
       int baseNoQTLModel=0, step=0;
	 //  Rprintf("DEBUG testing_3");
       for (j=0; j<Nmark; j++){   
	    /* 	fit a QTL in two steps:
			1. move QTL along marker interval j -> j+1 with steps of stepsize=20 cM, starting from -20 cM up to 220 cM
			2. all marker-cofactors in the neighborhood of the QTL are dropped by using cM='windows' as criterium
		*/
         nextinterval= 'n';
         while (nextinterval=='n')
         { // step 1:
		//   Rprintf("DEBUG testing STEP 1");
           if (position[j]=='L')
           {  if (moveQTL<=mapdistance[j])
              {  QTLposition[j]= position[j];
                 QTLposition[j+1]= 'M';
                 QTLr[j]= 0.5*(1.0-exp(-0.02*(mapdistance[j]-moveQTL)));
                 QTLr[j+1]= r[j];
                 QTLloci[j+1]= marker[j];
                 QTLloci[j]= marker[Nloci-1];
                 QTLmapdistance[j]= moveQTL;
                 QTLmapdistance[j+1]= mapdistance[j];
                 if (firsttime=='y') weight[0]= -1.0;
                 moveQTL+= stepsize;
              }
              else if (moveQTL<=mapdistance[j+1])
              {  QTLposition[j]= position[j];
                 QTLposition[j+1]= 'M';
				 //HERE there is an error, a slicht buggyness which should be solved the we can re-estimate rec frec. by EM.... QTLr should only depend on r[j]
                 QTLr[j]= 0.5*(1.0-exp(-0.02*(moveQTL-mapdistance[j])));
                 QTLr[j+1]= 0.5*(1.0-exp(-0.02*(mapdistance[j+1]-moveQTL))); //r[j];
                 QTLloci[j]= marker[j];
                 QTLloci[j+1]= marker[Nloci-1];
                 QTLmapdistance[j]= mapdistance[j];
                 QTLmapdistance[j+1]= moveQTL;
                 moveQTL+= stepsize;
              }
              else nextinterval= 'y';
           }
           else if (position[j]=='M')
           {  if (moveQTL<=mapdistance[j+1])
              {  QTLposition[j]= position[j];
                 QTLposition[j+1]= 'M';
				 //HERE there is an error, a slicht buggyness which should be solved the we can re-estimate rec frec. by EM.... QTLr should only depend on r[j]
                 QTLr[j]= 0.5*(1.0-exp(-0.02*(moveQTL-mapdistance[j]))); //0.0;
                 QTLr[j+1]= 0.5*(1.0-exp(-0.02*(mapdistance[j+1]-moveQTL))); //r[j];
                 QTLloci[j]= marker[j];
                 QTLloci[j+1]= marker[Nloci-1];
                 QTLmapdistance[j]= mapdistance[j];
                 QTLmapdistance[j+1]= moveQTL;
                 moveQTL+= stepsize;
              }
              else nextinterval= 'y';
           }
           else if (position[j]=='R')
           {  if (moveQTL<=stepmax)
              {  QTLposition[j]= 'M';
                 QTLposition[j+1]= 'R';
				 //HERE there is an error, a slicht buggyness which should be soled the we can re-estimate rec frec. by EM.... QTLr should only depend on r[j]
                 QTLr[j]= 0.5*(1.0-exp(-0.02*(moveQTL-mapdistance[j]))); //0.0;
                 QTLr[j+1]= r[j]; // note r[j]=999.0
                 QTLloci[j]= marker[j];
                 QTLloci[j+1]= marker[Nloci-1];
                 QTLmapdistance[j]= mapdistance[j];
                 QTLmapdistance[j+1]= moveQTL;
                 moveQTL+= stepsize;
              }
              else
              { nextinterval= 'y';
                moveQTL= stepmin;
              }
           }
           else if (position[j]=='U')
           {  QTLposition[j]= 'L';
              QTLposition[j+1]= 'R'; //position[j] ?? 'R' ?
              QTLr[j]= 0.0;
              QTLr[j+1]= r[j];
              QTLloci[j+1]= marker[j];
              QTLloci[j]= marker[Nloci-1];
              QTLmapdistance[j]= mapdistance[j];
              QTLmapdistance[j+1]= mapdistance[j];
              if (firsttime=='y') weight[0]= -1.0;
              nextinterval= 'y';
              moveQTL= stepmin;
           }
           if (nextinterval=='n')
           {  // QTLcofactor[j]= '0';
              // QTLcofactor[j+1]= '0';
              for (jj=0; jj<j; jj++)
              {   QTLposition[jj]= position[jj];
                  QTLr[jj]= r[jj];
                  QTLloci[jj]= marker[jj];
                  QTLmapdistance[jj]= mapdistance[jj];
                  QTLcofactor[jj]= selcofactor[jj];
              }
              for (jj=j+1; jj<Nmark; jj++)
              {   QTLposition[jj+1]= position[jj];
                  QTLr[jj+1]= r[jj];
                  QTLloci[jj+1]= marker[jj];
                  QTLcofactor[jj+1]= selcofactor[jj];
                  QTLmapdistance[jj+1]= mapdistance[jj];
                  QTLcofactor[jj+1]= selcofactor[jj];
              }
              // step 2:
			//  Rprintf("DEBUG testing STEP 2");
              if ((position[j]=='L')&&((moveQTL-stepsize)<=mapdistance[j]))
              {  QTLcofactor[j]= '0';
                 QTLcofactor[j+1]=
                       (((QTLmapdistance[j+1]-QTLmapdistance[j])<windowsize) ? '0' : selcofactor[j]);
              }
              else
              {  QTLcofactor[j+1]= '0';
                 QTLcofactor[j]=
                       (((QTLmapdistance[j+1]-QTLmapdistance[j])<windowsize) ? '0' : selcofactor[j]);
              }
              if ((position[j]=='L')||(position[j]=='M'))
              {   jjj=j+2;
                  while (QTLposition[jjj]=='M')
                  { if ((position[j]=='L')&&((moveQTL-stepsize)<=mapdistance[j]))
                       QTLcofactor[jjj]=
                       (((QTLmapdistance[jjj]-QTLmapdistance[j])<windowsize) ? '0' : QTLcofactor[jjj]);
                    else
                       QTLcofactor[jjj]=
                       (((QTLmapdistance[jjj]-QTLmapdistance[j+1])<windowsize) ? '0' : QTLcofactor[jjj]);
                    jjj++;
                  }
                  QTLcofactor[jjj]=
                  (((QTLmapdistance[jjj]-QTLmapdistance[j+1])<windowsize) ? '0' : QTLcofactor[jjj]);
              }
              if ((position[j]=='M')||(position[j]=='R'))
              {   jjj=j-1;
                  while (QTLposition[jjj]=='M')
                  { QTLcofactor[jjj]= (((QTLmapdistance[j+1]-QTLmapdistance[jjj])<windowsize) ? '0' : QTLcofactor[jjj]);
                    jjj--;
                  }
                  QTLcofactor[jjj]= (((QTLmapdistance[j+1]-QTLmapdistance[jjj])<windowsize) ? '0' : QTLcofactor[jjj]);
              }

              // fit no-QTL model at current map position (cofactors only)
              if (firsttime=='y')
              {  for (jj=0; jj<Nloci; jj++) saveQTLcofactor[jj]= QTLcofactor[jj];
                 baseNoQTLModel=1;
                 firsttime='n';
              }
              else
              {  baseNoQTLModel=0;
                 for (jj=0; jj<Nloci; jj++) baseNoQTLModel+= (saveQTLcofactor[jj]==QTLcofactor[jj] ? 0 : 1);
              }
//              cout << "fit NoQTL model(1=y; 0=n)= " << baseNoQTLModel << endl;
              if (baseNoQTLModel!=0) // new base no-QTL model
              {  if ((position[j]=='L')&&((moveQTL-stepsize)<=mapdistance[j])) QTLcofactor[j]= '2';
                 else QTLcofactor[j+1]= '2';
                 QTLlikelihood= -2.0*QTLmixture(QTLloci,QTLcofactor,QTLr,QTLposition,y,ind,Nind,Naug,Nloci,&variance,em,&weight0,REMLorML,fitQTL,dominance);
				 Rprintf("log-likelihood of NO QTL model= %f\n",QTLlikelihood/-2);
				 weight0[0]= -1.0;
                 savebaseNoQTLModel= QTLlikelihood;
                 if ((position[j]=='L')&&((moveQTL-stepsize)<=mapdistance[j])) QTLcofactor[j]= '0';
                 else QTLcofactor[j+1]= '0';
                 for (jj=0; jj<Nloci; jj++) saveQTLcofactor[jj]= QTLcofactor[jj];
              }
              else
                 QTLlikelihood= savebaseNoQTLModel;

              // fit QTL-model (plus cofactors) at current map position
              // '3'= QTL
              if ((position[j]=='L')&&((moveQTL-stepsize)<=mapdistance[j])) QTLcofactor[j]= '3';
              else QTLcofactor[j+1]= '3';
              if (REMLorML=='1') weight[0]= -1.0;
              QTLlikelihood+=2.0*QTLmixture(QTLloci,QTLcofactor,QTLr,QTLposition,y,ind,Nind,Naug,Nloci,&variance,em,&weight,REMLorML,fitQTL,dominance);
			  //this is the place we error at, because the likelyhood is not correct.
			  if (QTLlikelihood<-0.05) { Rprintf("Error Negative QTLlikelihood=%f  versus BASE MODEL:%f QTL at %d\n",QTLlikelihood,savebaseNoQTLModel,j);} //return 0;}
              maxF= (maxF<QTLlikelihood ? QTLlikelihood : maxF);
              if (run>0) (*Frun)[step][run]+= QTLlikelihood;
              else (*Frun)[step][0]+= QTLlikelihood;

            /* 	Each individual has condition multilocus probabilities for being 0, 1 or 2 at the QTL.
				Calculate the maximum per individu. Calculate the mean of this maximum, averaging over all individuals
				This is the information content plotted.
			*/
              infocontent= 0.0;
              for (int i=0; i<Nind; i++)
              {   info0[i]= 0.0; // qq
                  info1[i]= 0.0; // Qq
                  info2[i]= 0.0; // QQ
              }
              for (int i=0; i<Naug; i++)
              {   info0[ind[i]]+= weight[i];
                  info1[ind[i]]+= weight[i+Naug];
                  info2[ind[i]]+= weight[i+2*Naug];
              }
              for (int i=0; i<Nind; i++)
              if (info0[i]<info1[i]) infocontent+= (info1[i]<info2[i] ? info2[i] : info1[i]);
              else infocontent+= (info0[i]<info2[i] ? info2[i] : info0[i]);
              (*informationcontent)[step]+=infocontent/Nind;
              step++;
           }
         }
       }

    fitQTL='n';
	
	Free(info0);
    Free(info1);
    Free(info2);
    Free(weight);
    Free(weight0);
	Free(QTLr);
    Free(QTLposition);
	Free(Fy);
	Free(newcofactor);
    Free(QTLcofactor);
	Free(cumdistance);
	Free(QTLmapdistance);
    //Rprintf("MapQTL finished\n");
    return maxF; //QTLlikelihood;
}
 
/* end of MQMsupport.c */
