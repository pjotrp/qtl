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

double start_prob(char crosstype,char c){
	switch(crosstype){
		case 'F':
			return (c=='1' ? 0.5 : 0.25);
		break;
		case 'R':
			return (c=='1' ? 0.0 : 0.5);
		break;
		case 'B':
			return (c=='2' ? 0.0 : 0.5);
		break;
	}
	return 0.0;
}

double prob(cmatrix loci, vector r, int i, int j,char c,char crosstype,int JorC,int ADJ,int start){
	//Compares loci[j][i] versus loci[j+1][i]
	//OR if JorC is set to 1 loci[j][i] versus compareto
	//Specify an ADJ to adjust loci[j][i] to a specific location in the r[j+ADJ]
	
	double calc_i=0.0;
	double Nrecom;
	char compareto;
	
	if(JorC==1){
		//printf("C %d %d\n",i,j);
		compareto = c;
	}else{
		//printf("loci[j+1][i] %d\n",j);
		compareto = loci[j+1][i];
	}
	switch(crosstype){
		case 'F':
				if(start){
					return (loci[j][i]=='1' ? 0.5 : 0.25);
				}
				//printf("before Nrecom\n",j);				
				Nrecom= absdouble((double)loci[j][i]-(double)compareto);
				if ((loci[j][i]=='1')&&(compareto=='1')){
					//printf("SCase %c <-> %c:\n",compareto,loci[j][i]);
					calc_i= (r[j+ADJ]*r[j+ADJ]+(1.0-r[j+ADJ])*(1.0-r[j+ADJ]));}
				else if (Nrecom==0) {
					//printf("Nrecom=0 %c <-> %c:\n",compareto,loci[j][i]);
					calc_i= (1.0-r[j+ADJ])*(1.0-r[j+ADJ]);
				}else if (Nrecom==1) {
					//printf("Nrecom=1 %c <-> %c:\n",compareto,loci[j][i]);
					if(ADJ!=0){
						calc_i= ((loci[j][i]=='1') ? 2.0*r[j+ADJ]*(1.0-r[j+ADJ]) : r[j+ADJ]*(1.0-r[j+ADJ]));
					}else{
						calc_i= ((compareto=='1') ? 2.0*r[j+ADJ]*(1.0-r[j+ADJ]) : r[j+ADJ]*(1.0-r[j+ADJ]));
					}
				}else {
					//printf("Nrecom=2 %c <-> %c:\n",compareto,loci[j][i]);
					calc_i= r[j+ADJ]*r[j+ADJ];
				}
				//printf("after IF\n",j);
			break;
		case 'R':
				if(start){
					return 0.5;
				}
				if(compareto=='1' && JorC){
					return 0.0; // No chance in hell finding a 1 in an RIL
				}
				Nrecom= absdouble((double)loci[j][i]-(double)compareto);
				if(Nrecom==0){
					//No recombination has a chance of r[j]
					calc_i = 1.0-r[j+ADJ];
				}else{
					// Recombination between markers has a chance of r[j-1]
					calc_i = r[j+ADJ];
				}
			break;
		case 'B':
				if(start){
					return 0.5;
				}
				if(compareto=='2' && JorC){
					return 0.0; // No chance in hell finding a 2 in a BC
				}
				Nrecom= absdouble((double)loci[j][i]-(double)compareto);
				if(Nrecom==0){
					//No recombination has a chance of r[j]
					calc_i =  (1.0-r[j+ADJ]);
				}else{
					// Recombination between markers has a chance of r[j-1]
					calc_i = r[j+ADJ];
				}
			break;			
	}
	return calc_i;
}

double probright(char c, int jloc, cvector imarker, vector r, cvector position,char crosstype){
	//This is for an F2 population, where 'c'==1 stands for H (so it has two times higher chance than A or B
	double nrecom, prob0, prob1, prob2;
    if ((position[jloc]=='R')||(position[jloc]=='U')){
		//We're at the end of a chromosome or an unknown marker
		return 1.0;
	}
	switch(crosstype){
		case 'F':	
		if ((imarker[jloc+1]=='0')||(imarker[jloc+1]=='1')||(imarker[jloc+1]=='2')){
			//NEXT marker is known 
			if ((c=='1')&&(imarker[jloc+1]=='1')){
				//special case in which we observe a H after an H then we can't know if we recombinated or not
				return r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]);
			}else{
				//The number of recombinations between observed marker and the next marker
				nrecom = absdouble(c-imarker[jloc+1]);
				if(nrecom==0){
					//No recombination			
					return (1.0-r[jloc])*(1.0-r[jloc]);
				}else if (nrecom==1){
					if(imarker[jloc+1]=='1'){
						//the chances of having a H after 1 recombination are 2 times the chance of being either A or B
						return 2.0*r[jloc]*(1.0-r[jloc]);
					}else{ 
						//Chance of 1 recombination
						return r[jloc]*(1.0-r[jloc]);
					}
				}else{
					//Both markers could have recombinated which has a very low chance
					return r[jloc]*r[jloc];
				}
			}
		}else if (imarker[jloc+1]=='3'){
			//SEMI unknown next marker known is it is not an A
			if(c=='0'){
				//Observed marker is an A
				prob1= 2.0*r[jloc]*(1.0-r[jloc]);
				prob2= r[jloc]*r[jloc]; 
			}else if (c=='1') { 
				//Observed marker is an H
				prob1= r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]);
				prob2= r[jloc]*(1.0-r[jloc]); 
			}else{
				//Observed marker is an B
				prob1= 2.0*r[jloc]*(1.0-r[jloc]);
				prob2= (1.0-r[jloc])*(1-r[jloc]); 
			}
			return prob1*probright('1',jloc+1,imarker,r,position,crosstype) + prob2*probright('2',jloc+1,imarker,r,position,crosstype);
		}else if (imarker[jloc+1]=='4'){
			//SEMI unknown next marker known is it is not a B
			if(c=='0'){
				//Observed marker is an A
				prob0= (1.0-r[jloc])*(1.0-r[jloc]);
				prob1= 2.0*r[jloc]*(1.0-r[jloc]); 
			}else if (c=='1') { 
				//Observed marker is an H
				prob0= r[jloc]*(1.0-r[jloc]);
				prob1= r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]); 
			}else{
				//Observed marker is an B
				prob0= r[jloc]*r[jloc];
				prob1= 2.0*r[jloc]*(1.0-r[jloc]); 
			}
			return prob0*probright('0',jloc+1,imarker,r,position,crosstype) + prob1*probright('1',jloc+1,imarker,r,position,crosstype);
		}else{
		// Unknown next marker so estimate all posibilities (imarker[j+1]=='9')
			if(c=='0'){
				//Observed marker is an A
				prob0= (1.0-r[jloc])*(1.0-r[jloc]);
				prob1= 2.0*r[jloc]*(1.0-r[jloc]);
				prob2= r[jloc]*r[jloc]; 
			}else if (c=='1') { 
				//Observed marker is an H
				prob0= r[jloc]*(1.0-r[jloc]);
				prob1= r[jloc]*r[jloc]+(1.0-r[jloc])*(1.0-r[jloc]);
				prob2= r[jloc]*(1.0-r[jloc]); 
			}else{
				//Observed marker is an B
				prob0= r[jloc]*r[jloc];
				prob1= 2.0*r[jloc]*(1.0-r[jloc]);
				prob2= (1.0-r[jloc])*(1.0-r[jloc]); 
			}
			return prob0*probright('0',jloc+1,imarker,r,position,crosstype) + prob1*probright('1',jloc+1,imarker,r,position,crosstype) + prob2*probright('2',jloc+1,imarker,r,position,crosstype);
		}
		break;
		case 'R':
				if(c=='1'){
					return 0.0;
				}
				if ((imarker[jloc+1]=='0')||(imarker[jloc+1]=='2')){
					nrecom = absdouble(c-imarker[jloc+1]);
					if(nrecom==0){
						return (1.0-r[jloc]);
					}else{
						return r[jloc];
					}
				}else{
					if(c=='0'){//Both markers could have recombinated which has a very low chance
						prob0= (1.0-r[jloc]);
						prob2= r[jloc]; 
					}else{
						prob0= r[jloc];
						prob2= (1.0-r[jloc]);
					}
					return prob0*probright('0',jloc+1,imarker,r,position,crosstype) + prob2*probright('2',jloc+1,imarker,r,position,crosstype);
				}
			break;
		case 'B':
				if(c=='2'){
					return 0.0;
				}
				if ((imarker[jloc+1]=='0')||(imarker[jloc+1]=='1')){
					nrecom = absdouble(c-imarker[jloc+1]);
					if(nrecom==0){
						return 1.0-r[jloc];
					}else{
						return r[jloc];
					}
				}else{
					if(c=='0'){//Both markers could have recombinated which has a very low chance
						prob0= 1.0-r[jloc];
						prob2= r[jloc]; 
					}else{
						prob0= r[jloc];
						prob2= 1.0-r[jloc];
					}
					return prob0*probright('0',jloc+1,imarker,r,position,crosstype) + prob2*probright('1',jloc+1,imarker,r,position,crosstype);
				}
			break;
	}
	return 0.0;
}

}

/* end of MQMprob.c */
