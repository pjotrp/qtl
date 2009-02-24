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
 
 extern "C"
{

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



/*
 * analyseF2 - analyse one F2 family
 *
 */
void analyseF2(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, ivector f1genotype, int Backwards, 
			   double **QTL,vector *mapdistance,int **Chromo,int Nrun,int RMLorML, double windowsize,double stepsize,
			   double stepmin,double stepmax,double alfa,int em,int out_Naug,int **INDlist,char reestimate,char crosstype,char dominance)
{    
    int Naug;
	int run=0;
    cvector position;
	vector informationcontent;
	//char dominance='n';
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
    
	//Re-estimation of recombinant frequencies
	rmixture(newmarker, newweight, r, position, newind,Nind, Naug, Nmark, mapdistance,reestimate,crosstype);
	
	//Check if everything still is correct
	for (int j=0; j<Nmark; j++){
		if ((position[j]=='L')||(position[j]=='M')){
			r[j]= 0.5*(1.0-exp(-0.02*((*mapdistance)[j+1]-(*mapdistance)[j])));
			if (r[j]<0){
				Rprintf("error: recombination frequency is negative\n");
				Rprintf("j=%d chr=%d mapdistance=%f mapdistance=%f\n",j,chr[j],(*mapdistance)[j+1],(*mapdistance)[j]); 
				Rprintf("position=%d r[j]=%f\n",position[j], r[j]);
				return;
			}
		}
    }
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
     logLfull= QTLmixture(marker,cofactor,r,position,y,ind,Nind,Naug,Nmark,&variance,em,&weight,REMLorML,fitQTL,dominance,crosstype);
     Rprintf("log-likelihood of full model= %f\n",logLfull);
     Rprintf("residual variance= %f\n",variance);
     Rprintf("Trait mean= %f \nTrait variation= %f\n",ymean,yvari);

     if (Backwards==1)    // use only selected cofactors
         logLfull= backward(Nind, Nmark, cofactor, marker, y, weight, ind, Naug, logLfull,
                    variance, F1, F2, &selcofactor, r, position, &informationcontent, mapdistance,&Frun,run,REMLorML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype);
     if (Backwards==0) // use all cofactors
         logLfull= mapQTL(Nind, Nmark, cofactor, cofactor, marker, position,
                  (*mapdistance), y, r, ind, Naug, variance, 'n', &informationcontent,&Frun,run,REMLorML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype); // printout=='n'

	long *idum;
    idum = (long *)Calloc(1, long*);
    idum[0]=-1;

	 double savevariance= variance;
	 Rprintf("Saved Variance:%f\n",savevariance);
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
                    variance, F1, F2, &selcofactor, r, position,&informationcontent, mapdistance,&Frun,run,REMLorML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype);
			}
			if (Backwards==0){
                  maxF[run]= mapQTL(Nind, Nmark, cofactor, cofactor, marker, position,
                  (*mapdistance), y, r, ind, Naug, variance, 'n',&informationcontent,&Frun,run,REMLorML, fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype);
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

/* backward elimination in regression of trait on multiple cofactors
   routine subX haalt uit matrices voor volledige model de submatrices voor submodellen;
   matrices XtWX en Xt van volledig model worden genoemd fullxtwx en fullxt;
   analoog vector XtWY wordt full xtwy genoemd;
*/
double backward(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, vector weight, int* ind, int Naug, double logLfull, double variance, double F1, double F2, cvector* newcofactor, vector r, cvector position,vector *informationcontent,vector *mapdistance,matrix *Frun,int run,char REMLorML,char fitQTL,char dominance,int em, double windowsize,double stepsize,
			  double stepmin,double stepmax,char crosstype)
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
                  logL[j]= QTLmixture(marker,(*newcofactor),r,position,y,ind,Nind,Naug,Nmark,&variance,em,&weight,REMLorML,fitQTL,dominance,crosstype);
                  (*newcofactor)[j]='1';
               }
               else if ((*newcofactor)[j]=='2')
               {  //Rprintf("Drop marker %d\n",j);
                  (*newcofactor)[j]='0';
                  if (REMLorML=='1') variance= -1.0;
                  logL[j]=  QTLmixture(marker,(*newcofactor),r,position,y,ind,Nind,Naug,Nmark,&variance,em,&weight,REMLorML,fitQTL,dominance,crosstype);
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
           (*mapdistance), y, r, ind, Naug, variance, 'n', informationcontent,Frun,run,REMLorML,fitQTL,dominance, em, windowsize, stepsize, stepmin, stepmax,crosstype); // printoutput='n'
     //Rprintf("Backward selection finished\n");
     Free(logL);
     return maxF;
}

 }
 
/* end of MQMsupport.c */