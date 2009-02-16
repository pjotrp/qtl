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


double probleft(char c, int jloc, cvector imarker, vector r, cvector position){
	//This is for an F2 population, where 'c'==1 stands for H (so it has two times higher chance than A or B
	double nrecom;
	//IF position = L (start of chromosome) position = U (unknown position)
    if ((position[jloc]=='L')||(position[jloc]=='U')){
		if(c=='1'){
			return 0.50;
		}else{
			return 0.25;
		}
	}else if ((c=='1')&&(imarker[jloc-1]=='1')){
		//special case in which we observe a H after an H then we can't know if we recombinated or not
		return r[jloc-1]*r[jloc-1]+(1.0-r[jloc-1])*(1.0-r[jloc-1]);
	}else{
		//The number of recombinations between observed marker and the previous marker
		nrecom= absdouble(c-imarker[jloc-1]);
        if(nrecom==0){
			//No recombination
			return (1.0-r[jloc-1])*(1.0-r[jloc-1]);
		}else if (nrecom==1){
			//1 recombination
			if(c=='1'){
				//the chances of having a H after 1 recombination are 2 times the chance of being either A or B
				return 2.0*r[jloc-1]*(1.0-r[jloc-1]);
			}else{
				//Chance of 1 recombination
				return r[jloc-1]*(1.0-r[jloc-1]);
			}
		}else{
			//Both markers could have recombinated which has a very low chance
			return r[jloc-1]*r[jloc-1];
		}
    }
}

double probleft_RIL(char c, int jloc, cvector imarker, vector r, cvector position){
	//This is for an RIL population, 0=A , 2=B... there are no H types XD
	double nrecom;
	//IF position = L (start of chromosome) position = U (unknown position) we have 50 % chance being either
    if ((position[jloc]=='L')||(position[jloc]=='U')){
		return 0.50;
	}else{
		//The number of recombinations between observed marker and the previous marker
		nrecom= absdouble(c-imarker[jloc-1]);
        if(nrecom==0){
			//No recombination has a chance of r[j]
			return (1.0-r[jloc-1]);
		}else{
			// Recombination between markers has a chance of r[j-1]
			return r[jloc-1];
		}
    }
}

double probright(char c, int jloc, cvector imarker, vector r, cvector position){
	//This is for an F2 population, where 'c'==1 stands for H (so it has two times higher chance than A or B
	double nrecom, prob0, prob1, prob2;
    if ((position[jloc]=='R')||(position[jloc]=='U')){
		//We're at the end of a chromosome or an unknown marker
		return 1.0;
	}else if ((imarker[jloc+1]=='0')||(imarker[jloc+1]=='1')||(imarker[jloc+1]=='2')){
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
		return prob1*probright('1',jloc+1,imarker,r,position) + prob2*probright('2',jloc+1,imarker,r,position);
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
		return prob0*probright('0',jloc+1,imarker,r,position) + prob1*probright('1',jloc+1,imarker,r,position);
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
        return prob0*probright('0',jloc+1,imarker,r,position) + prob1*probright('1',jloc+1,imarker,r,position) + prob2*probright('2',jloc+1,imarker,r,position);
    }
}

double probright_RIL(char c, int jloc, cvector imarker, vector r, cvector position){
	//This is for an RIL  population
	double nrecom, prob0, prob1;
    if ((position[jloc]=='R')||(position[jloc]=='U')){
		//We're at the end of a chromosome or an unknown marker
		return 1.0;
	}else if ((imarker[jloc+1]=='0')||(imarker[jloc+1]=='2')){
		//NEXT marker is known
		//The number of recombinations between observed marker and the next marker
		nrecom = absdouble(c-imarker[jloc+1]);
        if(nrecom==0){
		//No recombination			
			return (1.0-r[jloc]);
		}else{
			return r[jloc];
			
		}
	}else if (imarker[jloc+1]=='1'||imarker[jloc+1]=='3'||imarker[jloc+1]=='4'){
		Rprintf("Unwanted situation: Marker at %d = %d. Illegal RIL marker",jloc+1,imarker[jloc+1]);
		return 0.0;
	}else{
	// Unknown next marker so estimate all posibilities (imarker[j+1]=='9')
		if(c=='0'){
			//Observed marker is an A
			prob0= (1.0-r[jloc]);
			prob1= r[jloc]; 
		}else if (c=='2') {
			//Observed marker is an B
			prob0= r[jloc];
            prob1= (1.0-r[jloc]);
		}else{
			//Observed marker is an H (or a C or a D)
			Rprintf("Unwanted situation: Marker at %d = %d. Illegal RIL marker",jloc,imarker[jloc]);
			return 0.0;
		}
        return prob0*probright_RIL('0',jloc+1,imarker,r,position) + prob1*probright_RIL('1',jloc+1,imarker,r,position);
    }
}

/* end of MQMprob.c */
