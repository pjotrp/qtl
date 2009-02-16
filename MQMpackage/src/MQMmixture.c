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

/* ML estimation of recombination frequencies via EM;
    calculation of multilocus genotype probabilities;
    ignorance of unlikely genotypes*/
void rmixture(cmatrix marker, vector weight, vector r,
              cvector position, ivector ind,
              int Nind, int Naug, int Nmark,vector *mapdistance){   
	int i,j;
    int iem= 0;
    double Nrecom, oldr=0.0, newr, rdelta=1.0;
    vector indweight;
    indweight = newvector(Nind);
	vector distance;
    distance= newvector(Nmark+1);
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
    //    for (j=0; j<Nmark; j++){
	//		if ((position[j]=='L')||(position[j]=='U')){
	//			(*mapdistance)[j]= -50*log(1-2.0*r[j]);
	//			Rprintf("r(%d)= %f -> %f\n",j,r[j],(*mapdistance)[j]);
	//		}
	//	}
	Free(indweight);
}


/* ML estimation of parameters in mixture model via EM;
*/
double QTLmixture(cmatrix loci, cvector cofactor, vector r, cvector position,
              vector y, ivector ind, int Nind, int Naug,
              int Nloci,
              double *variance, int em, vector *weight,char REMLorML,char fitQTL,char dominance){
	//Rprintf("QTLmixture called\n");
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
 
/* end of MQMmixture.c */
