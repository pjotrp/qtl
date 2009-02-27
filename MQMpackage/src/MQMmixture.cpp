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
#include "MQMprob.h"
#include "Regression.h"

/* ML estimation of recombination frequencies via EM;
    calculation of multilocus genotype probabilities;
    ignorance of unlikely genotypes*/
void rmixture(cmatrix marker, vector weight, vector r,
              cvector position, ivector ind,
              int Nind, int Naug, int Nmark,vector *mapdistance, char reestimate,char crosstype){   
	int i,j;
    int iem= 0;
    double Nrecom, oldr=0.0, newr, rdelta=1.0;
    vector indweight;
    indweight = newvector(Nind);
	vector distance;
    distance= newvector(Nmark+1);

    if (reestimate=='n'){
		Rprintf("INFO: recombination parameters are not re-estimated\n");
    }else{
		Rprintf("INFO: recombination parameters are re-estimated\n");
	}
	//Reestimation of map now works
     while ((iem<1000)&&(rdelta>0.0001))
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
               {   
				   double calc_i = prob(marker,r,i,j,marker[j+1][i],crosstype,0,0,0);
				   weight[i]*=calc_i;
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
                  if (reestimate=='y' && position[j]!='R') //only update if it isn't the last marker of a chromosome ;)
                  {  oldr=r[j];
                     r[j]= newr/(2.0*Nind);
                     rdelta+=pow(r[j]-oldr,2.0);
                  }
                  else rdelta+=0.0;
               }
            }
     }
     
/*   print new estimates of recombination frequencies */

	float last_step = 0.0;
    if (reestimate=='y'){  
        for (j=0; j<Nmark; j++){
			if(position[j+1]=='R'){
				last_step = (*mapdistance)[j+1]-(*mapdistance)[j];
			}
			if(position[j]!='L'){
				if(position[j]!='R'){
					(*mapdistance)[j]= -50*log(1-2.0*r[j])+(*mapdistance)[j-1];
				}else{
					(*mapdistance)[j]= (*mapdistance)[j-1]+last_step;
				}
			}else{
				(*mapdistance)[j]= -50*log(1-2.0*r[j]);
			}
			//Rprintf("r(%d)= %f -> %f\n",j,r[j],(*mapdistance)[j]);
		}
	}
	Rprintf("INFO: Re-estimation of the genetic map took %d iterations, to reach a rdelta of %f\n",iem,rdelta);
	Free(indweight);
}


/* ML estimation of parameters in mixture model via EM;
*/
double QTLmixture(cmatrix loci, cvector cofactor, vector r, cvector position,
              vector y, ivector ind, int Nind, int Naug,
              int Nloci,
              double *variance, int em, vector *weight,char REMLorML,char fitQTL,char dominance,char crosstype){
	//Rprintf("QTLmixture called\n");
    int iem= 0, newNaug, i, j;
    char varknown, biasadj='n';
	double oldlogL=-10000, delta=1.0, calc_i, logP=0.0, Pscale=1.75;
    
	vector indweight, Ploci, Fy;
    
	indweight= newvector(Nind);
    newNaug= (fitQTL=='n' ? Naug : 3*Naug);
    Ploci= newvector(newNaug);
    Fy= newvector(newNaug);
    logP= Nloci*log(Pscale); // only for computational accuracy
	varknown= (((*variance)==-1.0) ? 'n' : 'y' );
    Ploci= newvector(newNaug);	
    if ((REMLorML=='0')&&(varknown=='n')){ 
		Rprintf("INFO: Variance is being estimated and bias adjusted\n");
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
			//Here we have ProbLeft
		    if ((position[j]=='L')||(position[j]=='U')){
				for (i=0; i<Naug; i++){
					calc_i= prob(loci,r,i,j,'1',crosstype,1,0,1);
					Ploci[i]*= calc_i;
				}
			}
		    if ((position[j]=='L')||(position[j]=='M')){
				for (i=0; i<Naug; i++){
					calc_i = prob(loci,r,i,j,loci[j+1][i],'F',0,0,0);
					Ploci[i]*= calc_i;
				}
			}
		}
	}else{
	//Rprintf("FitQTL=Y\n");	
     for (j=0; j<Nloci; j++)
     {    for (i=0; i<Naug; i++)
          {   Ploci[i]*= Pscale; Ploci[i+Naug]*= Pscale; Ploci[i+2*Naug]*= Pscale;
              // only for computational accuracy; see use of logP
          }
          if ((position[j]=='L')||(position[j]=='U'))
          {  
			//Here we don't have any f2 dependancies anymore by using the prob function
			if (cofactor[j]<='1')
             for (i=0; i<Naug; i++)
             {   calc_i= prob(loci,r,i,j,'1',crosstype,1,0,1);
                 Ploci[i]*= calc_i; Ploci[i+Naug]*= calc_i; Ploci[i+2*Naug]*= calc_i;
             }
             else
             for (i=0; i<Naug; i++)
             {
				//startvalues for each new chromosome
				Ploci[i]*= start_prob(crosstype,'0');
				Ploci[i+Naug]*= start_prob(crosstype,'1'); 
				Ploci[i+2*Naug] *= start_prob(crosstype,'2');
			 }
                 // QTL='0', '1' or'2'
          }
          if ((position[j]=='L')||(position[j]=='M'))
          {  if ((cofactor[j]<='1')&&(cofactor[j+1]<='1'))
             for (i=0; i<Naug; i++){  
				calc_i = prob(loci,r,i,j,loci[j+1][i],crosstype,0,0,0);
                Ploci[i]*= calc_i; Ploci[i+Naug]*= calc_i; Ploci[i+2*Naug]*= calc_i;
             }
             else if (cofactor[j]<='1') // locus j+1 == QTL
             for (i=0; i<Naug; i++)
             {  // QTL=='0' What is the prob of finding an '0' at J=1
				calc_i = prob(loci,r,i,j,'0',crosstype,1,0,0);
                Ploci[i]*= calc_i;
                // QTL=='1'
                calc_i = prob(loci,r,i,j,'1',crosstype,1,0,0);
                Ploci[i+Naug]*= calc_i;
                // QTL=='2'
                calc_i = prob(loci,r,i,j,'2',crosstype,1,0,0);
                Ploci[i+2*Naug]*= calc_i;
             }
             else // locus j == QTL
             for (i=0; i<Naug; i++)
             {  // QTL=='0'
                calc_i = prob(loci,r,i,j+1,'0',crosstype,1,-1,0);
                Ploci[i]*= calc_i;
                // QTL=='1'
				calc_i = prob(loci,r,i,j+1,'1',crosstype,1,-1,0);
                Ploci[i+Naug]*= calc_i;
                // QTL=='2'
                calc_i = prob(loci,r,i,j+1,'2',crosstype,1,-1,0);
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
 
 }
/* end of MQMmixture.c */
