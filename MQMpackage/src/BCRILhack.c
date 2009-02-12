#include <R.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/RS.h> /* for Calloc, Realloc */
#include <R_ext/Utils.h>
#include "scanMQM.h"
#include "MQMdata.h"
#include "MQMsupport.h"

vector MARKmixture(char cross, cmatrix loci, cvector cofactor, vector r, char fitQTL, int NewNaug, cvector position){
	int Nloci,newNaug,Naug,j,i;
	char varknown, biasadj='n';
    double Nrecom, oldlogL=-10000, delta=1.0, calc_i, logP=0.0, Pscale=1.75;
    vector indweight, Ploci, Fy;	
	
	if(cross =='F'){
		newNaug= (fitQTL=='n' ? Naug : 3*Naug);
	}else{
		newNaug= (fitQTL=='n' ? Naug : 2*Naug);
	}
	for (j=0; j<Nloci; j++){
		for (i=0; i<Naug; i++)
            Ploci[i]*= Pscale; // only for computational accuracy; see use of logP
        if ((position[j]=='L')||(position[j]=='U'))
        for (i=0; i<Naug; i++) Ploci[i]*= (loci[j][i]=='1' ? 0.5 : 0.25);  //HERE we should edit for BC/RISELF
        if ((position[j]=='L')||(position[j]=='M'))
        {   for (i=0; i<Naug; i++)
            {   Nrecom= absdouble((double)loci[j][i]-(double)loci[j+1][i]);
                if ((loci[j][i]=='1')&&(loci[j+1][i]=='1'))
					calc_i= (r[j]*r[j]+(1.0-r[j])*(1.0-r[j]));
                else if (Nrecom==0) calc_i= (1.0-r[j])*(1.0-r[j]);
                else if (Nrecom==1) calc_i= ((loci[j+1][i]=='1')
                        ? 2.0*r[j]*(1.0-r[j]) : r[j]*(1.0-r[j]));
                else calc_i= r[j]*r[j];
                Ploci[i]*= calc_i;
             }
        }
    }
	return Ploci;
}