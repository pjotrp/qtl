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
#include "scanMQM.h"
#include "MQMdata.h"
#include "MQMprob.h"


void R_augdata(int *geno,double *dist,double *pheno,int *auggeno,double *augPheno,int *augIND,int *Nind,int *Naug,int *Nmark, int *Npheno, int *maxaug, int *maxiaug,double *neglect,int *chromo){
	int **Geno;
	double **Pheno;
	double **Dist;
	int **NEW;
	int **Chromo;
	double **NEWPheno;
	int **NEWIND;
	int prior = *Nind;
	Rprintf("Starting C-part of the dataaugmentation routine\n");
	ivector new_ind;
    vector new_y,r,mapdistance;
	cvector position;
    cmatrix markers,new_markers;
	ivector chr;
	
	markers= newcmatrix(*Nmark,*Nind);
	new_markers= newcmatrix(*Nmark,*maxaug);
	r = newvector(*Nmark);
	mapdistance = newvector(*Nmark);
	position= newcvector(*Nmark);
	chr= newivector(*Nmark);
	
	//Reorganise the pointers into arrays, Singletons are just cast into the function
	reorg_geno(*Nind,*Nmark,geno,&Geno);
	reorg_int(*Nmark,1,chromo,&Chromo); 
	reorg_pheno(*Nind,*Npheno,pheno,&Pheno);
	reorg_pheno(*Nmark,1,dist,&Dist);
   
    reorg_int(*maxaug,*Nmark,auggeno,&NEW);	   
	reorg_int((*maxiaug)*(*Nind),1,augIND,&NEWIND);	 
	reorg_pheno(*maxaug,1,augPheno,&NEWPheno);	 
	
	//Change all the markers from Karl format to our own
	for(int i=0; i< *Nmark; i++){
		for(int j=0; j< *Nind; j++){ 
			markers[i][j] = '9';
			if(Geno[i][j] == 1){				//AA
				markers[i][j] = '0';
			}
			if(Geno[i][j] == 2){				//AB
				markers[i][j] = '1';
			}
			if(Geno[i][j] == 3){				//BB
				markers[i][j] = '2';
			}
			if(Geno[i][j] == 4){				//AA of AB
				markers[i][j] = '4';
			}
			if(Geno[i][j] == 5){				//BB of AB
				markers[i][j] = '3';
			}
		}
		mapdistance[i]=999.0;
	    mapdistance[i]=Dist[0][i];
	}
	Rprintf("Filling the chromosome matrix\n");

	for(int i=0; i<(*Nmark); i++){
		chr[i] = Chromo[0][i];
	}
	Rprintf("Calculating relative genomepositions of the markers\n");

    for (int j=0; j<(*Nmark); j++)
    {   
        if (j==0)
        { if (chr[j]==chr[j+1]) position[j]='L'; else position[j]='U'; }
        else if (j==(*Nmark-1))
        { if (chr[j]==chr[j-1]) position[j]='R'; else position[j]='U'; }
        else if (chr[j]==chr[j-1])
        { if (chr[j]==chr[j+1]) position[j]='M'; else position[j]='R'; }
        else
        { if (chr[j]==chr[j+1]) position[j]='L'; else position[j]='U'; }
    }
	Rprintf("Estimating recombinant frequencies\n");

	for (int j=0; j<(*Nmark); j++){   
		r[j]= 999.0;
		if ((position[j]=='L')||(position[j]=='M')){
			r[j]= 0.5*(1.0-exp(-0.02*(mapdistance[j+1]-mapdistance[j])));
			if (r[j]<0){
				Rprintf("error: recombination frequency is negative\n");
				Rprintf("position=%d r[j]=%d\n",position[j], r[j]);
				return;
			}
		}
		//Rprintf("recomfreq:%d,%f\n",j,r[j]);
    }

	if(augdata(markers, Pheno[(*Npheno-1)], &new_markers, &new_y, &new_ind, Nind, Naug, *Nmark, position, r,*maxaug,*maxiaug,*neglect)==1){
		//Data augmentation finished succesfully
		//Push it back into RQTL format
		for (int i=0; i<(*Nmark); i++){   
			for (int j=0; j<(*Naug); j++){
				NEWPheno[0][j] = new_y[j];
				NEWIND[0][j] = new_ind[j];
				NEW[i][j] = 9;
				if(new_markers[i][j] == '0'){
					NEW[i][j] = 1;
				}
				if(new_markers[i][j] == '1'){
					NEW[i][j] = 2;
				}
				if(new_markers[i][j] == '2'){
					NEW[i][j] = 3;
				}
				if(new_markers[i][j] == '3'){
					NEW[i][j] = 5;
				}
				if(new_markers[i][j] == '4'){
					NEW[i][j] = 4;
				}				
			}
		}
		delcmatrix(new_markers,(*Nmark));
		delcmatrix(markers,(*Nmark));
		Free(mapdistance);
		Free(position);
		Free(r);
		Free(chr);
		Rprintf("Data augmentation finished succesfull\n");
		Rprintf("# Unique individuals before augmentation:%d\n",prior);
		Rprintf("# Unique selected individuals:%d\n",*Nind);
		Rprintf("# Marker p individual:%d\n",*Nmark);
		Rprintf("# Individuals after augmentation:%d\n",*Naug);
	}else{
		//Unsuccessfull data augmentation exit
		*Naug = *Nind;
		for (int i=0; i<(*Nmark); i++){   
			for (int j=0; j<(*Naug); j++){
				NEWPheno[0][j] = Pheno[0][j];
				NEW[i][j] = 9;
				if(markers[i][j] == '0'){
					NEW[i][j] = 1;
				}
				if(markers[i][j] == '1'){
					NEW[i][j] = 2;
				}
				if(markers[i][j] == '2'){
					NEW[i][j] = 3;
				}
				if(markers[i][j] == '3'){
					NEW[i][j] = 5;
				}
				if(markers[i][j] == '4'){
					NEW[i][j] = 4;
				}					
			}
		}
		delcmatrix(new_markers,(*Nmark));
		delcmatrix(markers,(*Nmark));
		Free(mapdistance);
		Free(position);
		Free(r);
		Free(chr);		
		Rprintf("Data augmentation failed\n");
	}
    return;
}

int augdata(cmatrix marker, vector y, cmatrix* augmarker, vector *augy, ivector* augind, int *Nind, int *Naug, int Nmark, cvector position, vector r,int maxNaug,int imaxNaug,double neglect){

	int jj;
    int newNind=(*Nind);
    (*Naug)= maxNaug; /* maximum size of augmented dataset */
	cmatrix newmarker;
    vector newy;
    cvector imarker;
    ivector newind;

    newmarker= newcmatrix(Nmark+1,*Naug);
    newy= newvector(*Naug);
    newind= newivector(*Naug);
    imarker= newcvector(Nmark);

    int iaug=0;      // iaug keeps track of current augmented individual
    int maxiaug=0;   // highest reached(?)
    int saveiaug=0;  // previous iaug
    double prob0, prob1, prob2, sumprob,
           prob0left, prob1left, prob2left,
           prob0right, prob1right, prob2right;
    double probmax;
    vector newprob, newprobmax;
    newprob= newvector(*Naug);
    newprobmax= newvector(*Naug);
	Rprintf("Parameters: MAXaug=%d,MAXindaug=%d,Neglect=%f\n",maxNaug, imaxNaug, neglect);
    // ---- foreach individual create one in the newmarker matrix
    for (int i=0; i<(*Nind); i++)
    {   newind[iaug]=i-((*Nind)-newNind);  // index of individuals
         newy[iaug]= y[i];               // cvariance
         newprob[iaug]= 1.0;
         probmax= 1.0;
         for (int j=0; j<Nmark; j++) newmarker[j][iaug]=marker[j][i];
         for (int j=0; j<Nmark; j++)
         {   maxiaug=iaug;
             if ((maxiaug-saveiaug)<=imaxNaug)  // within bounds for individual?
               for (int ii=saveiaug; ii<=maxiaug; ii++)
               {   if (newmarker[j][ii]=='3')
                   {  for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];
                      prob1left= probleft('1',j,imarker,r,position);
                      prob2left= probleft('2',j,imarker,r,position);
                      prob1right= probright('1',j,imarker,r,position);
                      prob2right= probright('2',j,imarker,r,position);
                      prob1= prob1left*prob1right;
                      prob2= prob2left*prob2right;
                      if (ii==saveiaug) probmax= (prob2>prob1 ? newprob[ii]*prob2 : newprob[ii]*prob1);
                      if (prob1>prob2)
                      {  if (probmax/(newprob[ii]*prob2)<neglect)
                         {  iaug++;
                            newmarker[j][iaug]= '2';
                            newprob[iaug]= newprob[ii]*prob2left;
                            newprobmax[iaug]= newprob[iaug]*prob2right;
                            for (jj=0; jj<Nmark; jj++)
                            { if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii]; }
                            newind[iaug]=i-((*Nind)-newNind);
                            newy[iaug]=y[i];
                         }
                         newmarker[j][ii]= '1';
                         newprobmax[ii]= newprob[ii]*prob1;
                         newprob[ii]= newprob[ii]*prob1left;
                      }
                      else
                      {  if (probmax/(newprob[ii]*prob1)<neglect)
                         {  iaug++;
                            newmarker[j][iaug]= '1';
                            newprob[iaug]= newprob[ii]*prob1left;
                            newprobmax[iaug]= newprob[iaug]*prob1right;
                            for (jj=0; jj<Nmark; jj++)
                            {   if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii]; }
                            newind[iaug]=i-((*Nind)-newNind);
                            newy[iaug]=y[i];
                         }
                         newmarker[j][ii]= '2';
                         newprobmax[ii]= newprob[ii]*prob2;
                         newprob[ii]*= prob2left;
                      }
                      probmax= (probmax>newprobmax[ii] ? probmax : newprobmax[ii]);
                   }
                   else if (newmarker[j][ii]=='4')
                     {  for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];
                        prob0left= probleft('0',j,imarker,r,position);
                        prob1left= probleft('1',j,imarker,r,position);
                        prob0right= probright('0',j,imarker,r,position);
                        prob1right= probright('1',j,imarker,r,position);
                        prob0= prob0left*prob0right;
                        prob1= prob1left*prob1right;
                        if (ii==saveiaug) probmax= (prob0>prob1 ? newprob[ii]*prob0 : newprob[ii]*prob1);
                        if (prob1>prob0)
                        {  if (probmax/(newprob[ii]*prob0)<neglect)
                           {  iaug++;
                              newmarker[j][iaug]= '0';
                              newprob[iaug]= newprob[ii]*prob0left;
                              newprobmax[iaug]= newprob[iaug]*prob0right;
                              for (jj=0; jj<Nmark; jj++)
                              {   if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii]; }
                              newind[iaug]=i-((*Nind)-newNind);
                              newy[iaug]=y[i];
                           }
                           newmarker[j][ii]= '1';
                           newprobmax[ii]= newprob[ii]*prob1;
                           newprob[ii]*= prob1left;
                        }
                        else
                        {  if (probmax/(newprob[ii]*prob1)<neglect)
                           {  iaug++;
                              newmarker[j][iaug]= '1';
                              newprob[iaug]= newprob[ii]*prob1left;
                              newprobmax[iaug]= newprob[iaug]*prob1right;
                              for (jj=0; jj<Nmark; jj++)
                              {   if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                              }
                              newind[iaug]=i-((*Nind)-newNind);
                              newy[iaug]=y[i];
                           }
                           newmarker[j][ii]= '0';
                           newprobmax[ii]= newprob[ii]*prob0;
                           newprob[ii]*= prob0left;
                        }
                        probmax= (probmax>newprobmax[ii] ? probmax : newprobmax[ii]);
                     }
                   else if (newmarker[j][ii]=='9')
                     {  for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];
                        prob0left= probleft('0',j,imarker,r,position);
                        prob1left= probleft('1',j,imarker,r,position);
                        prob2left= probleft('2',j,imarker,r,position);
                        prob0right= probright('0',j,imarker,r,position);
                        prob1right= probright('1',j,imarker,r,position);
                        prob2right= probright('2',j,imarker,r,position);
                        prob0= prob0left*prob0right;
                        prob1= prob1left*prob1right;
                        prob2= prob2left*prob2right;
                        if (ii==saveiaug)
                        {  if ((prob2>prob1)&&(prob2>prob0)) probmax= newprob[ii]*prob2;
                           else if ((prob1>prob0)&&(prob1>prob2)) probmax= newprob[ii]*prob1;
                           else probmax= newprob[ii]*prob0;
                        }
                        if ((prob2>prob1)&&(prob2>prob0))
                        {  if (probmax/(newprob[ii]*prob1)<neglect)
                           {  iaug++;
                              newmarker[j][iaug]= '1';
                              newprob[iaug]= newprob[ii]*prob1left;
                              newprobmax[iaug]= newprob[iaug]*prob1right;
                              for (jj=0; jj<Nmark; jj++)
                              {   if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                              }
                              newind[iaug]=i-((*Nind)-newNind);
                              newy[iaug]=y[i];
                           }
                           if (probmax/(newprob[ii]*prob0)<neglect)
                           {  iaug++;
                              newmarker[j][iaug]= '0';
                              newprob[iaug]= newprob[ii]*prob0left;
                              newprobmax[iaug]= newprob[iaug]*prob0right;
                              for (jj=0; jj<Nmark; jj++)
                              {   if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                              }
                              newind[iaug]=i-((*Nind)-newNind);
                              newy[iaug]=y[i];
                           }
                           newmarker[j][ii]= '2';
                           newprobmax[ii]= newprob[ii]*prob2;
                           newprob[ii]*= prob2left;

                        }
                        else if ((prob1>prob2)&&(prob1>prob0))
                        {  if (probmax/(newprob[ii]*prob2)<neglect)
                           {  iaug++;
                              newmarker[j][iaug]= '2';
                              newprob[iaug]= newprob[ii]*prob2left;
                              newprobmax[iaug]= newprob[iaug]*prob2right;
                              for (jj=0; jj<Nmark; jj++)
                              {   if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                              }
                              newind[iaug]=i-((*Nind)-newNind);
                              newy[iaug]=y[i];
                           }
                           if (probmax/(newprob[ii]*prob0)<neglect)
                           {  iaug++;
                              newmarker[j][iaug]= '0';
                              newprob[iaug]= newprob[ii]*prob0left;
                              newprobmax[iaug]= newprob[iaug]*prob0right;
                              for (jj=0; jj<Nmark; jj++)
                              {   if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                              }
                              newind[iaug]=i-((*Nind)-newNind);
                              newy[iaug]=y[i];
                           }
                           newmarker[j][ii]= '1';
                           newprobmax[ii]= newprob[ii]*prob1;
                           newprob[ii]*= prob1left;
                        }
                        else
                        {  if (probmax/(newprob[ii]*prob1)<neglect)
                           {  iaug++;
                              newmarker[j][iaug]= '1';
                              newprob[iaug]= newprob[ii]*prob1left;
                              newprobmax[iaug]= newprob[iaug]*prob1right;
                              for (jj=0; jj<Nmark; jj++)
                              {   if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                              }
                              newind[iaug]=i-((*Nind)-newNind);
                              newy[iaug]=y[i];
                           }
                           if (probmax/(newprob[ii]*prob2)<neglect)
                           {  iaug++;
                              newmarker[j][iaug]= '2';
                              newprob[iaug]= newprob[ii]*prob2left;
                              newprobmax[iaug]= newprob[iaug]*prob2right;
                              for (jj=0; jj<Nmark; jj++)
                              {   if (jj!=j) newmarker[jj][iaug]=newmarker[jj][ii];
                              }
                              newind[iaug]=i-((*Nind)-newNind);
                              newy[iaug]=y[i];
                           }
                           newmarker[j][ii]= '0';
                           newprobmax[ii]= newprob[ii]*prob0;
                           newprob[ii]*= prob0left;
                        }
                        probmax= (probmax>newprobmax[ii] ? probmax : newprobmax[ii]);
                     }
                   else // newmarker[j][ii] is observed
                   {  for (jj=0; jj<Nmark; jj++) imarker[jj]= newmarker[jj][ii];
                      newprob[ii]*= probleft(newmarker[j][ii],j,imarker,r,position);
                   }

                   if (iaug+3>maxNaug)
                   {       
                    Rprintf("ERROR in augmentation routine: dataset too large after augmentation\n");
                    Rprintf("Recall procedure with larger value for parameter maxaug or lower for the parameter neglect\n");
					// Better not free them, we don't know if the arrays already contain something, perhaps not... then we would segfault in R
					//Free(newy);
					//Free(newmarker);
					//Free(newind);
					//Free(newprob);
					//Free(newprobmax);
					//Free(imarker);
                    return 0;
                   }
               }
             if ((iaug-saveiaug+1)>imaxNaug)
             {  newNind-= 1;
                iaug= saveiaug-1;
				Rprintf("Individual %d is eliminated\n",i);
             }
             sumprob= 0.0;
             for (int ii=saveiaug; ii<=iaug; ii++) sumprob+= newprob[ii];
             for (int ii=saveiaug; ii<=iaug; ii++) newprob[ii]/= sumprob;
         }
         iaug++;
         saveiaug=iaug;
     }
     *Naug= iaug;
     *Nind= newNind;
     *augmarker= newcmatrix(Nmark,*Naug);
     *augy= newvector(*Naug);
     *augind = newivector(*Naug);
     for (int i=0; i<(*Naug); i++)
     {   (*augy)[i]= newy[i];
         (*augind)[i]= newind[i];
         for (int j=0; j<Nmark; j++) (*augmarker)[j][i]= newmarker[j][i];
     }
	Free(newy);
	Free(newmarker);
	Free(newind);
	Free(newprob);
	Free(newprobmax);
	Free(imarker);
	return 1;
}


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



/* end of MQMdata.c */
