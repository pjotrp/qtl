/*
Routine MQM:
Regression on multiple cofactors
Selection of important cofactors by backward elimination
*/

using namespace std;

#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
//#include <alloc.h> // for coreleft()

/*
Routine MQM for analysing multiple families - should also work with one family
R.C. Jansen
------------------------------------------------------------------------
Symbol definitions for instructions-file mqm.in
------------------------------------------------------------------------
Nind          number of individuals per family
Nfam          number of families
Nmark         number of markers
cross         type of cross (F2 only in current version)
filem         file with molecular marker data
filey         file with quantitative trait data
dominance     y/n (assumed, yes or no)
Plan for Ritsert Jansen's MQM Code
==================================

- Standardized ML/EM and linear regression solutions
- Simplify algorithms and reduce use of 'special' values and globals
- Introduce unit tests
- Embedding in R
- Increased performance and multi-threading
- Flexible input model
- Flexibility for number of tested QTLs at the same time
- Variable number of absorbing surrounding markers for LR
- Flexibility for tweaking parameters and models (MC EM)
REMLorML      0/1 (REML or ML analysis)
permutations  number of permutations (->significance level)
------------------------------------------------------------------------
Format for specifying map in mqm.in (marker-by-marker; each marker at a new row)
------------------------------------------------------------------------
chromosome    number of chromosome
markername    any name
mapposition   in cM, measured from beginning of chromosome
cofactor      indicator (blank, or * for cofactor) - depends on size of
              population
------------------------------------------------------------------------
symbol definitions of some default parameters in routine
------------------------------------------------------------------------
em=100        maximum number of em iterations
alfa=0.02     alfa used in backward selection procedure (threshold for dropping)
maxNaug=10000 maximum size of augmented dataset
imaxNaug=1000 maximum size of augmented data for individual i (1000 genotypes)
neglect=100   eliminate genotypes 100 times less likely
              than the most likely configuration
maxdimX=50    maximum size of design matrix X in regression (max number of cofactors?)
windowsize=20 drop cofactor if less than windowsize away from QTL (20 cM)
------------------------------------------------------------------------
symbol definitions of some data structures in routine
------------------------------------------------------------------------
marker        markervalues for Nmark markers, recoded to
              0 homozyote parent A
              1 heterozygote
              2 homozygote parent B
              3 not parent A
              4 not parent B
              7 homozygote parent A, not segregating within the family - not used
              8 homozygote parent B, not segregating within the family - not used
              9 missing
y             trait values
chr           chromosome-numbers for Nmark markers
Ncof          number of markers used as cofactors
cofactor      markers used as cofactors (0 or 1 per marker)
mapdistance   measured in cM from beginning of chromosome (dimension Nmark)
position      position of marker on it's chromosome
              L left
              M middle
              R right
              U unlinked
r             recombination frequencies between flanking markers
ind           individual's number (numbered 0 up to Nind-1)
Naug          size of augmented dataset
Nrun          permutations (or simulations if set)
------------------------------------------------------------------------
short description of subroutines in MQM
------------------------------------------------------------------------
analyseF2     you call this routine after having read all the marker and
              trait data. It will do a single-family F2 analysis with
              cofactors as specified. All other routines will be called
              within routine analyseF2
augmentdata   genotypic marker information is incomplete when
              - marker scores are missing or dominant, or
              - markers are non-segregating but you want to assume a
              segregating QTL on top of this marker (~very closely linked),
              and all likely possible configurations are generated
augmentdata   if you study a certain map position for the presence of a QTL,
forQTL        then QTL genotypic scores are missing. All likely
              configurations are generated
rmixture      maximum-likelihood estimation of recombination frequencies
              via the EM algorithm, using multilocus information (default:
              the recombination frequencies are not estimated but taken from
              mqm.in)
QTLmixture    maximum-likelihood estimation of parameters in the mixture model
              via the EM algorithm, using multilocus information, but assuming
              known recombination frequencies
regression    performs weighted regression of trait on genotype (QTL and
              cofactors) for augmented data
backward      selects "important" cofactors in weighted regression of trait
              on genotype (cofactors) using the augmented data
mapQTL        moves a QTL along the chromosome and calculated at each map
              position the QTL likelihood. Uses either all cofactors, or
              selected cofactors only

Interactions are not (yet) included, but could be added to the model. There
is the question of multiple testing and overfitting, though. Maybe if we do
one at a time and drop the uninteresting ones. Or add cofactors which don't
interact.
------------------------------------------------------------------------
*/
double neglect=1000; // eliminate unlikely genotype configurations
int maxNaug=10000; // maximum size of augmented dataset
int imaxNaug=1000; // maximum size of augmented data for individual i
int maxdimX=50; // maximum size of design matrix in regression
int em=1000; // maximum number of em iterations
double alfa=0.02; // alfa used in selection procedure
double windowsize=25.0; // used in mapQTL procedure
double stepsize=5; // size of steps when moving QTL along chromosomes (for output)
double stepmin=-20; // start moving QTL at position stepmin cM (for output)
double stepmax=220; // move QTL up to stepmax (for output)
long *idum; // for monte carlo simulation or permutation
/*
------------------------------------------------------------------------
datastructures and routines for matrix and vector calculus
------------------------------------------------------------------------ */
typedef double** matrix;
typedef double*  vector;
typedef char**   cmatrix;
typedef char*    cvector;
typedef int*  ivector;

vector newvector(int dim)
{      vector v;
       v= new double[dim];
       if (v==NULL) { cout << " not enough memory for new vector of dimension " << (dim+1); exit(1); }
       return v;
}

ivector newivector(int dim)
{      ivector v;
       v= new int[dim];
       if (v==NULL) { cout << " not enough memory for new vector of dimension " << (dim+1); exit(1); }
       return v;
}

cvector newcvector(int dim)
{      cvector v;
       v= new char[dim];
       if (v==NULL) { cout << " not enough memory for new char vector of dimension " << (dim+1); exit(1); }
       return v;
}

matrix newmatrix(int rows, int cols)
{      matrix m;
       m=new double*[rows];
       if (m==NULL) { cout << " not enough memory for new double matrix"; exit(1); }
       for (int i=0; i<rows; i++) m[i]= newvector(cols);
       return m;
}

void   printmatrix(matrix m, int rows, int cols)
{      for (int r=0; r<rows; r++)
       {   for (int c=0; c<cols; c++) cout << setw(5) << m[r][c] << " ";
           cout << endl;
       }
}

void   printcmatrix(cmatrix m, int rows, int cols)
{      for (int r=0; r<rows; r++)
       {   for (int c=0; c<cols; c++) cout << setw(2) << m[r][c];
           cout << endl;
       }
}

cmatrix newcmatrix(int rows, int cols)
{      cmatrix m;
       m=new char*[rows];
       if (m==NULL) { cout << " not enough memory for new char matrix"; exit(1); }
       for (int i=0; i<rows; i++) m[i]= newcvector(cols);
       return m;
}

void delmatrix(matrix m, int rows)
{      for (int i=0; i<rows; i++) delete[] m[i];
       delete[] m;
}

void delcmatrix(cmatrix m, int rows)
{      for (int i=0; i<rows; i++) delete[] m[i];
       delete[] m;
}

void copyvector(vector vsource, vector vdestination, int dim)
{    for (int i=0; i<dim; i++) vdestination[i]= vsource[i];
}

// void luinvert(matrix a, matrix inv, int n, int *indx);
void ludcmp(matrix m, int dim, ivector ndx, int &d);
void lusolve(matrix lu, int dim, ivector ndx, vector b);
/*------------------------------------------------------------------------
declaration of subroutines
--------------------------------------------------------------------------*/
char dominance='n';
void analyseF2(int Nind, int Nmark, cvector cofactor,
       cmatrix marker, vector y, ivector f1genotype);
void simuF2(int Nind, int Nmark, cvector cofactor,
       cmatrix marker, vector y);
double probleft(char c, int j, cvector imarker, vector r, cvector position);
double probright(char c, int j, vector imarker, vector r, cvector position);
void augmentdata(cmatrix marker, vector y, cmatrix &augmarker, vector &augy,
     ivector &augind, int &Nind, int &Naug, int Nmark,
     cvector position, vector r);
void rmixture(cmatrix marker, vector weight, vector r,
     cvector position, ivector ind,
     int Nind, int Naug, int Nmark);
double regression(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y,
       vector weight, ivector ind, int Naug,
       double &variance, vector Fy, char biasadj);
double backward(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, vector weight,
       ivector ind, int Naug, double logLfull, double variance,
       double F1, double F2, cvector &newcofactor, vector r, cvector position);
double mapQTL(int Nind, int Nmark, cvector cofactor,
     cvector selcofactor, cmatrix marker, cvector position, vector mapdistance,
     vector y, vector r, ivector ind, int Naug,
     double variance, char printoutput);
double QTLmixture(cmatrix loci, cvector cofactor, vector r, cvector position,
     vector y, ivector ind, int Nind, int Naug, int Nloci,
     double &variance, int em, vector &weight);
/*
-----------------------------------------------------------------------
subroutines from book 'Numerical Recipees in C'
for calculating F-probabilities and
for generating randomly permuted trait data
for other tasks
-----------------------------------------------------------------------*/
double gammln(double xx);
double betai(double a, double b, double x);
double betacf(double a, double b, double x);
double inverseF(int df1, int df2, double alfa);
double ran2(long *idum);
double randomnormal(long *idum);
void sort1(int n, vector ra);
void sort2(int n, double *ra, ivector rb);

/* void select()
{ if (strcmp(cross,"F2")==0)  printf(cross);
  else printf("BC");
} */
char ok, defset;
void OK()
{    ok='0';
     if (defset=='n') {cout << "OK (y/n)?"; cin >> ok; if (ok=='n') exit(1); }
}

double Lnormal(double residual, double variance)
{      double Likelihood;
       Likelihood=exp(-pow(residual/sqrt(variance),2.0)/2.0 - log(sqrt(2.0*acos(-1.0)*variance)));
       return Likelihood;
}

double absdouble(double x)
{      double z; z= (x<0 ? -x : x); return z;}

int mod(int a, int b)
{      int c;
       c= a/b;
       return a-b*c;
}
/*
-----------------------------------------------------------------------
beginning of the MAIN part of the MQM routine
(reads data and calls routine analyseF2)
-----------------------------------------------------------------------*/
vector mapdistance, r;
matrix XtWX;
cmatrix Xt;
vector XtWY;
ivector chr;
matrix Frun;
vector informationcontent;
int Nfam, Nrun=0;
int run=-1;
char REMLorML;
char fitQTL='n';
char perm_simu;

int main()
{ cout << "Welcome to routine MQM" << endl;
  //------------------ data structures -----------------------------------
  int Nmark, saveNmark, Nind, saveNind, fam=0, f2genotype, i, ii, j;
  idum= new long[1];
  idum[0]=-1;


  vector y; //cvariance
  cvector cross, datafilef1marker, datafilef2marker, datafilephenotype,
          skipstring, cofactor;
  ivector f1genotype;
  cmatrix markername, marker;
  char ch, real_simu;
  double readf;

  //------------------ allocate memory for vectors and matrices ----------
  cross= newcvector(4);  // char[4]
  datafilef1marker= newcvector(50);
  datafilef2marker= newcvector(50);
  datafilephenotype= newcvector(50);
  skipstring= newcvector(50);
  ifstream ff("mqm_in.txt", ios::in);

// Nind= 10
// Nfam= 60
// Nmark= 200

  ff >> skipstring >> Nind >> skipstring >> Nfam >> skipstring >> Nmark;
  saveNmark= Nmark;
  saveNind= Nind;
//  ff.close();

  markername= newcmatrix(Nmark,20);
  cofactor= newcvector(Nmark);       // iscofactor?
  mapdistance= newvector(Nmark);
  chr= newivector(Nmark);            // chromosome number
  f1genotype= newivector(Nmark);     // parent genotype

// chr name distance cofactor
// 1 L116 32
// 1 L129 58
// 1 L152 104
// 1 L154 108 *
// 1 L155 110
// 1 L156 112

  for (j=0; j<Nmark; j++)
  {   mapdistance[j]=999.0;
      ff >> chr[j];
      ff >> markername[j];
      do ff.get(ch); while (ch==' ');
      if (ch=='*') cofactor[j]='1';
      else if (ch=='s') cofactor[j]='2'; /* sexe of the mouse */
      else if (ch=='\n' || ch=='\r') cofactor[j]='0';
      else
      {   ff.putback(ch);
          ff >> mapdistance[j];
          do ff.get(ch); while (ch==' ');
          if (ch=='*') cofactor[j]='1';
          else if (ch=='s') cofactor[j]='2'; /* sexe of the mouse */
          else if (ch=='\n' || ch=='\r') cofactor[j]='0';
      }
  }


// cross= F2
// filef1m= f1c.mar
// filef2m= f2c.mar
// filey= f2c.qua
// dominance= n
// REMLorML= 0
// defset= y
// real_simu= 0
// permutation= 0

  ff >> skipstring >> cross
     >> skipstring >> datafilef1marker  >> skipstring >> datafilef2marker
     >> skipstring >> datafilephenotype
     >> skipstring >> dominance >> skipstring >> REMLorML
     >> skipstring >> defset    >> skipstring >> real_simu;
  if (ff.get(ch))
  {  ff.putback(ch);
     ff >> skipstring >> perm_simu >> Nrun;
  }
  cout << "Nind= " << Nind << endl;
  cout << "Nmark= " << Nmark << endl;
  for (j=0; j<Nmark; j++)
  cout << setw(3) << j << setw(3) << chr[j] << "  "
       << markername[j] << " " << setw(2) << cofactor[j] << " " << mapdistance[j] << endl;
  cout << "cross= " << cross << endl;
  cout << "datafile for markers (F1)= " << datafilef1marker << endl;
  cout << "datafile for markers (F2)= " << datafilef2marker << endl;
  cout << "datafile for phenotypes= " << datafilephenotype << endl;
  cout << "dominance= " << dominance << endl;
  cout << "REMLorML= " << REMLorML << endl;
  if (Nrun>0) (perm_simu=='0' ?
               cout << "number of permutations=" << Nrun << endl :
               cout << "number of simulations=" << Nrun << endl);
  int Nsteps;
  Nsteps= chr[Nmark-1]*((stepmax-stepmin)/stepsize+1);
  // out: 10 -20 220 5 490
  cout << "chr " << chr[Nmark-1] << ", stepmin " << stepmin << ", stepmax "
       << stepmax << ", stepsize " << stepsize << ", Nsteps " << Nsteps << endl;
  // exit(1);

	// ---- Initialize Frun and informationcontent to 0.0
  Frun= newmatrix(Nsteps,Nrun+1);
  informationcontent= newvector(Nsteps);
  for (i=0; i<Nrun+1; i++)
    for (ii=0; ii<Nsteps; ii++) Frun[ii][i]= 0.0;
  for (ii=0; ii<Nsteps; ii++) informationcontent[ii]= 0.0;
  OK();
  ff.close();

  ifstream fpheno(datafilephenotype, ios::in);
  if (!fpheno.is_open()) {
    printf("Failed to open %s",datafilephenotype);
    exit(2);
  }
  fpheno.close();

//  if (datafilef1marker[0]!='-')
//  {  ifstream f1(datafilef1marker, ios::in);
//     f1.close();
//  }

  ifstream f2(datafilef2marker, ios::in);
   if (!f2.is_open()) {
       printf("Failed to open %s",datafilef2marker);
       exit(2);
   }
  f2.close();
  ofstream fff("mqm_out.txt", ios::out | ios::app);
  fff << endl << endl;
  fff.close();

  //------------------ analysis of data of families ---------------------
  // repeat for number of families (e.g. 60) families...
  for (fam=0; fam<Nfam; fam++)
  {  cout << "Analysis of data of family " << fam << endl;
     Nmark= saveNmark;
     Nind= saveNind;
     y = newvector(Nind);  // cvariance
     marker = newcmatrix(Nmark,Nind);

  //------------------ read instructions-file mqm.in ---------------------
  //  rewind(ff);
  // FIXME: reading file again for every family...
  ff.open("mqm_in.txt", ios::in);

  // Nind= 10
  // Nfam= 60
  // Nmark= 200
  // 1 L116 32   chr name distance cofactor
  // 1 L129 58 
  // 1 L152 104
  // 1 L154 108 *
  // 1 L155 110
  // 1 L156 112

  ff >> skipstring >> Nind >> skipstring >> Nfam >> skipstring >> Nmark;
  for (j=0; j<Nmark; j++)
  {   mapdistance[j]=999.0;
      ff >> chr[j];
      ff >> markername[j];
      do ff.get(ch); while (ch==' ');
      if (ch=='*') cofactor[j]='1';
      else if (ch=='s') cofactor[j]='2'; /* sexe of the mouse */
      else if (ch=='\n' || ch=='\r') cofactor[j]='0';
      else
      {   ff.putback(ch);
          ff >> mapdistance[j];
          do ff.get(ch); while (ch==' ');
          if (ch=='*') cofactor[j]='1';
          else if (ch=='s') cofactor[j]='2'; /* sexe of the mouse */
          else if (ch=='\n' || ch=='\r') cofactor[j]='0';
      }
  }
  ff.close();

   //------------------ analyse ONE family --------------------------------
   //fff.open("mqm.out", ios::out | ios::app);
   //fff << endl << endl << fam << endl;
   //fff.close();

   //------------------ read quantitative trait data for family under study
   //------------------ open data files for reading -----------------------
 	 //  f2c.qua 
   // 118.3  125.8  130.8  142.6  111.8  122.1  124.9  117.2  144.5  120.2
   // 124.7  100.3  125.6  106.8  125.1  117.8  109.4  117.4  116.8  102.7  124.5
   // 124.2  129.6  142.1  134.3  116.1  145.0  137.1  139.6  128.0  131.7  113.2

   fpheno.open(datafilephenotype, ios::in);
   if (!fpheno.is_open()) {
       printf("Failed to open %s",datafilephenotype);
       exit(2);
   }
   for (ii=0; ii<fam*Nind; ii++) fpheno >> readf;
   for (i=0; i<Nind; i++)
   {   fpheno >> y[i];
       // y[i]= log(y[i]+1.0);
   }
   cout << "IMPORTANT: trait data are log-transformed during analysis !!!" << endl;
   fpheno.close();
   //------------------ read f1 and f2 marker data ------------------------
//   f1.open("f1.mar", ios::in);
   if (datafilef1marker[0]!='-')
   {  ifstream f1(datafilef1marker, ios::in);
      // f1c.mar (60 F1's or families)
      // *L116        11 11 11 11 21 22 21 12 11 11 11 11 11 21 12 21 11 21 12 12 12 22 21 12 12 11 11 21 22 21 11 12 11 21 22 12 12 12 21 21 22 21 11 12 11 12 22 11 12 11 11 21 21 21 21 22 22 21 11 12
      // *L129        21 22 22 11 12 21 21 12 21 22 12 21 22 12 11 12 12 21 21 12 12 22 11 21 12 11 12 11 11 12 11 12 21 12 11 22 21 21 21 21 21 12 12 22 12 22 21 21 12 21 11 11 11 21 22 22 21 12 21 21


       if (!f1.is_open()) {
         printf("Failed to open %s",datafilef1marker);
         exit(2);
       }
       // f2c.mar (10 individuals per family - 600 F2's)
       // *L116        11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 (etc.)
      f2.open(datafilef2marker, ios::in);
       if (!f2.is_open()) {
         printf("Failed to open %s",datafilef2marker);
         exit(2);
       }
      // ---- For each marker j read (rewound) F1 and F2 input files
      for (j=0; j<Nmark; j++)
      {   f1 >> skipstring;
          // ---- skip fam F1 genotypes and read marker F1 genotype
          for (i=0; i<(fam+1); i++) f1 >> f1genotype[j];
          f2 >> skipstring;
          // if ((f1genotype[j]==11)||(f1genotype[j]==22))    // geen uitsplitsing in F2 !!
          if (mod(f1genotype[j],11)==0)    // F1 is 11 or 22: geen uitsplitsing in F2 !!
          {  for (i=0; i<(Nind*Nfam); i++) f2 >> f2genotype;  // skip all F2
             for (i=0; i<Nind; i++) marker[j][i]= '9'; // (f1genotype[j]==11 ? '9' : '9');
          }
          else // F1 parent is heterozygous - set marker to 0,1,2
          {  for (i=0; i<(fam*Nind); i++) f2 >> f2genotype;   // skip to file pos
             for (i=0; i<Nind; i++) // read F2 genotype
             {   f2 >> f2genotype;
                 // F2 individual is heterozygous?
                 if (mod(f2genotype,11)!=0) marker[j][i]= '1';
                 // F2 individual is homozygous for 2nd allele of F1 parent?
                 else if (mod(f2genotype,10)==mod(f1genotype[j],10)) marker[j][i]= '2';
                 // F2 individual is homozygous for 1st allele of F1 parent?
                 else marker[j][i]= '0';
                 /*if (f2genotype==11) marker[j][i]= (f1genotype[j]==12 ? '0' : '2');
                 else if ((f2genotype==12)||(f2genotype==21)) marker[j][i]= '1';
                 else if (f2genotype==22) marker[j][i]= (f1genotype[j]==12 ? '2' : '0');
                 else { cout << "error reading marker data"; exit(1); }*/
             }
             for (i=(fam+1)*Nind; i<(Nind*Nfam); i++) f2 >> f2genotype; // skip all F2
          }
          for (i=fam+1; i<Nfam; i++) f1 >> f2genotype; // skip all F1
      }
      f1.close();
      f2.close();
   }
   else // Different routing for datafilef1marker[0]=='-'
   {  f2.open(datafilef2marker, ios::in);
       if (!f2.is_open()) {
         printf("Failed to open %s",datafilef2marker);
         exit(2);
       }

      for (j=0; j<Nmark; j++)
      {   f1genotype[j]= 12;
          f2 >> skipstring;
          for (i=0; i<Nind; i++)
          {   do f2.get(ch); while ((ch==' ')||(ch=='\n')||ch=='\r');
              if (ch=='A') marker[j][i]='0';
              else if (ch=='H') marker[j][i]= '1';
              else if (ch=='B') marker[j][i]= '2';
              else if (ch=='C') marker[j][i]= '3';
              else if (ch=='D') marker[j][i]= '4';
              else if (ch=='-') marker[j][i]= '9';
              else if (ch=='U') marker[j][i]= '9';
              else { cout << "error reading marker data; ch=" << ch; exit(1); }
          }
      }
      f2.close();
   }

    // ---- At this point all input files have been parsed for one family
 
    //     delete[] datafilemarker, datafilephenotype, skipstring;
    //     delete[] mapdistance;
    run= -1;
    if (real_simu=='1') simuF2(Nind, Nmark, cofactor, marker, y);  // simulate cvariance
    analyseF2(Nind, Nmark, cofactor, marker, y, f1genotype);
    cout << "Analysis of data of family " << fam << " finished" << endl;
  } // foreach family

  // ---- Write output file
  fff.open("mqm_out.txt", ios::out | ios::app);
  fff << endl;
  double moveQTL= stepmin;
  int chrnumber=1;
  cout << "-1- " << Nsteps << " " << Nrun << " " << Nfam << endl;
 
  // chr pos Frun    information 
  // 1  -20  97.4561 0.677204
  // 1  -15 103.29   0.723067
  // 1  -10 108.759  0.777696
  // 1   -5 113.737  0.842778
  // 1    0 118.112  0.920356
  // 1    5 120.051  0.928594
  // 1   10 114.469  0.959548

  for (ii=0; ii<Nsteps; ii++)
  {   fff << chrnumber << " " << moveQTL << " " << Frun[ii][Nrun]
          << " " << ((informationcontent[ii]/Nfam)/(Nrun+1)) << endl;
      if (moveQTL+stepsize<=stepmax) moveQTL+= stepsize;
      else { moveQTL= stepmin; chrnumber++; }
  }
  fff << ":" << endl;
  cout << "-2- " << endl;
  if (Nrun>0) // ((Nfam>1)&&(Nrun>0))  // multiple permutations or simulations
  {  for (ii=0; ii<Nsteps; ii++)
       for (i=0; i<Nrun; i++)
         Frun[0][i]= (Frun[0][i]<Frun[ii][i] ? Frun[ii][i] : Frun[0][i]);
     cout << endl;
     cout << "Cumulative distribution of maximum test statistic value in "
          << Nrun << " permutations or simulations" << endl;
     fff  << "Cumulative distribution of maximum test statistic value in "
          << Nrun << " permutations or simulations" << endl;
     sort1(Nrun,Frun[0]);
     if (Nrun>1)
       for (i=1; /* (((double)run/( (double)Nrun+1.0))<0.1)*/ i<Nrun+1; i++)
       {   cout << setprecision(8) << ( (double)i/( (double)Nrun+1.0) )
                << " " << setprecision(8) << Frun[0][i-1] << endl;
           fff << setprecision(8) << ( (double)i/( (double)Nrun+1.0) )
               << " " << setprecision(8) << Frun[0][i-1] << endl;
       }
     fff.close();
  }
  delmatrix(Frun,Nsteps);
  return 0;
}
/*
-----------------------------------------------------------------------
end of the MAIN part of the MQM routine
below follows code for all (sub)routines
-----------------------------------------------------------------------*/

/*
 * simuF2 for every individual calculate a random cvariance (y). Next the
 * markers are walked and depending on type the cvariance is adjusted by +/- 1
 */
void simuF2(int Nind, int Nmark, cvector cofactor,
       cmatrix marker, vector y)
{    //long *idum;
     //idum= new long[1];
     //idum[0]=-1;
     for (int i=0; i<Nind; i++)
     {   y[i]= pow(5,0.5)*randomnormal(idum);
         for (int j=0; j<Nmark; j++)
           if (cofactor[j]=='1')
             if (marker[j][i]=='0') y[i] -= 1.0;
              else if (marker[j][i]=='2') y[i] += 1.0;
         /*else if (marker[j][i]!='1')
         {  cout << "Marker score is missing..." << marker[j][i]; exit(1); }
          family is not segregating for marker[j], i.e. QTL has only a
          family effect, which can be ignored in analysis per family */
     }
     //delete[] idum;
}

/*
 * analyseF2 - analyse one F2 family
 *
 */

void analyseF2(int Nind, int Nmark, cvector cofactor,
       cmatrix marker, vector y, ivector f1genotype)
{    int Naug;
     cvector position;
     r= newvector(Nmark);
     position= newcvector(Nmark);
     for (int j=0; j<Nmark; j++)
     {   r[j]= 999.0;
         if (j==0)
            { if (chr[j]==chr[j+1]) position[j]='L'; else position[j]='U'; }
         else if (j==(Nmark-1))
            { if (chr[j]==chr[j-1]) position[j]='R'; else position[j]='U'; }
         else if (chr[j]==chr[j-1])
            { if (chr[j]==chr[j+1]) position[j]='M'; else position[j]='R'; }
         else
            { if (chr[j]==chr[j+1]) position[j]='L'; else position[j]='U'; }
     }
     for (int j=0; j<Nmark; j++)
     {   if ((position[j]=='L')||(position[j]=='M'))
         r[j]= 0.5*(1.0-exp(-0.02*(mapdistance[j+1]-mapdistance[j])));
         // cout << "r(" << setw(2) << j << ")=" << r[j] << endl;
     }

     char dropj='y';
     int jj=0;
     cout << "any triple of non-segregating markers is considered to be the result of: " << endl;
     cout << "identity-by-descent (IBD) instead of identity-by-state (IBS)" << endl;
     cout << "no (segregating!) cofactors are fitted in such non-segregating IBD regions" << endl;
     OK();
     for (int j=0; j<Nmark; j++)
     {   // if ((f1genotype[j]==12)||(f1genotype[j]==21)) dropj='n';
         if (mod(f1genotype[j],11)!=0) dropj='n';
         else if (cofactor[j]=='0') dropj='y';
         else if (position[j]=='L') // (cofactor[j]!='0') cofactor at non-segregating marker
         // test whether next segregating marker is nearby (<20cM)
         {  dropj='y';
            if (((mapdistance[j+1]-mapdistance[j])<20)&&(mod(f1genotype[j+1],11)!=0)) dropj='n';
            else if (position[j+1]!='R')
            if (((mapdistance[j+2]-mapdistance[j])<20)&&(mod(f1genotype[j+2],11)!=0)) dropj='n';
         }
         else if (position[j]=='M')
         {  dropj='y';
            if (((mapdistance[j]-mapdistance[j-1])<20)&&(mod(f1genotype[j-1],11)!=0)) dropj='n';
            else if (((mapdistance[j+1]-mapdistance[j])<20)&&(mod(f1genotype[j+1],11)!=0)) dropj='n';
         }
         else if (position[j]=='R')
         {  dropj='y';
            if (((mapdistance[j]-mapdistance[j-1])<20)&&(mod(f1genotype[j-1],11)!=0)) dropj='n';
            else if (position[j-1]!='L')
            if (((mapdistance[j]-mapdistance[j-2])<20)&&(mod(f1genotype[j-2],11)!=0)) dropj='n';
         }
         if (dropj=='n')
         {  marker[jj]= marker[j];
            cofactor[jj]= cofactor[j];
            mapdistance[jj]= mapdistance[j];
            chr[jj]= chr[j];
            r[jj]= r[j];
            position[jj]= position[j];
            jj++;
         }
         else if (cofactor[j]=='1')
         {  cout << "cofactor at chr " << chr[j] << " is dropped" << endl;
            OK();
         }
     }

     Nmark= jj;

     for (int j=0; j<Nmark; j++)
     {   r[j]= 999.0;
         if (j==0)
            { if (chr[j]==chr[j+1]) position[j]='L'; else position[j]='U'; }
         else if (j==(Nmark-1))
            { if (chr[j]==chr[j-1]) position[j]='R'; else position[j]='U'; }
         else if (chr[j]==chr[j-1])
            { if (chr[j]==chr[j+1]) position[j]='M'; else position[j]='R'; }
         else
            { if (chr[j]==chr[j+1]) position[j]='L'; else position[j]='U'; }
     }
     for (int j=0; j<Nmark; j++)
     {   if ((position[j]=='L')||(position[j]=='M'))
         r[j]= 0.5*(1.0-exp(-0.02*(mapdistance[j+1]-mapdistance[j])));
         // cout << "r(" << setw(2) << j << ")=" << r[j] << endl;
         if (r[j]<0)
         {  cout << "error: recombination frequency is negative" << endl;
            cout << "chr=" << chr[j] << " mapdistance=" << mapdistance[j] 
                 << " position=" << position[j] << " r[j]=" << r[j] << endl;
            exit(1);
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
     cout << "ymean=" << ymean << " yvari=" << yvari << endl;

     augmentdata(marker,y,newmarker,newy,newind,Nind,Naug,Nmark,position,r);
     cout << "Nind= " << Nind << "; Naug= " << Naug << endl;
     OK();

     delcmatrix(marker,Nmark);
     vector newweight;
     newweight= newvector(Naug);
     rmixture(newmarker, newweight, r, position, newind,Nind, Naug, Nmark);

     /* eliminate individuals with missing trait values */
     int oldNind=Nind;
     for (int i=0; i<oldNind; i++) Nind-= ((y[i]==999.0) ? 1 : 0);
     delete[] y;
     int oldNaug=Naug;
     for (int i=0; i<oldNaug; i++) Naug-= ((newy[i]==999.0) ? 1 : 0);
     vector weight;
     ivector ind;
     marker= newcmatrix(Nmark,Naug);
     y= newvector(Naug);
     ind= newivector(Naug);
     weight= newvector(Naug);
     int newi=0;
     for (int i=0; i<oldNaug; i++)
     if (newy[i]!=999.0)
     {  y[newi]= newy[i];
        ind[newi]= newind[i];
        weight[newi]= newweight[i];
        for (int j=0; j<Nmark; j++) marker[j][newi]= newmarker[j][i];
        newi++;
     }
     int diff;
     for (int i=0; i<Naug-1; i++)
     {   diff=ind[i+1]-ind[i];
         if  (diff>1)
         for (int ii=i+1; ii<Naug; ii++) ind[ii]=ind[ii]-diff+1;
     }
     delcmatrix(newmarker,Nmark);
     delete[] newy;
     delete[] newind;
     delete[] newweight;

 //    vector Fy;
 //    Fy= newvector(Naug);
     double variance=-1.0, logLfull;
     cvector selcofactor;
     selcofactor= newcvector(Nmark); /* selected cofactors */

     int dimx=1;
     for (int j=0; j<Nmark; j++)
     if (cofactor[j]=='1') dimx+= (dominance=='n' ? 1 : 2);  // per QTL only additivity !!
     else if (cofactor[j]=='2') { dimx+=1; cout << "sex of mouse" << endl; } /* sex of the mouse */
     double F1, F2;
     F1= inverseF(1,Nind-dimx,alfa);
     F2= inverseF(2,Nind-dimx,alfa);
     cout << "F(" << F1 << ",1," << (Nind-dimx) << ")=0.02" << endl;
     cout << "F(" << F2 << ",2," << (Nind-dimx) << ")=0.02" << endl;
     F2= 2.0* F2; // 9-6-1998 using threshold x*F(x,df,alfa)
     OK();

     XtWX= newmatrix(dimx+2,dimx+2);
     Xt= newcmatrix(dimx+2,Naug); // 3*Naug);  --- there is note it should be 3x
     XtWY= newvector(dimx+2);
     weight[0]= -1.0;
     logLfull= QTLmixture(marker,cofactor,r,position,y,ind,Nind,Naug,Nmark,variance,em,weight);
     cout << "log-likelihood of full model= " << logLfull << endl;
     cout << "residual variance= " << variance << "OK ?" << endl;
     cout << "ymean= " << ymean << " yvari= " << yvari << endl;
     OK();
     cout << endl;
     cout << endl << "Use all cofactors or selected cofactors in mapQTL ?" << endl;
     cout << "Answer= (1=selected; 0=all)" << endl;
     OK();
     cout << endl;

     if (ok=='1')    // use only selected cofactors
         logLfull= backward(Nind, Nmark, cofactor, marker, y, weight, ind, Naug, logLfull,
                    variance, F1, F2, selcofactor, r, position);
     else if (ok=='0') // use all cofactors
         logLfull= mapQTL(Nind, Nmark, cofactor, cofactor, marker, position,
                  mapdistance, y, r, ind, Naug, variance, 'n'); // printout=='n'
     OK();

     //long *idum;
     //idum= new long[1];
     //idum[0]=-1;
     double savevariance= variance;
     double *urand;
     vector maxF;
     ivector indorder; // individu 0...Nind-1; order will be permuted
     vector yoriginal;
     indorder= newivector(Nind);
     yoriginal= newvector(Nind);
     if (Nrun>0)
     {  urand= new double[Naug];
        urand[0]= ran2(idum);
        maxF= newvector(Nrun);
        for (int i=0; i<Naug; i++) yoriginal[ind[i]]= y[i];

        for (run=0; run<Nrun; run++)
        {   cout << "run= " << run << endl;
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
            if (ok=='1')    // use only selected cofactors
                  maxF[run]= backward(Nind, Nmark, cofactor, marker, y, weight, ind, Naug, logLfull,
                    variance, F1, F2, selcofactor, r, position);
            else // if (ok=='0') // use all cofactors
                  maxF[run]= mapQTL(Nind, Nmark, cofactor, cofactor, marker, position,
                  mapdistance, y, r, ind, Naug, variance, 'n');
            // cout << "run " << run <<" ready; maxF= " << maxF[run] << endl;
        }
        /* if (Nfam==1)
        {  sort1(Nrun,maxF);
           cout << endl;
           cout << "Cumulative distribution of maximum test statistic value in "
                << Nrun << " permutations or simulations" << endl;
           ofstream fff("mqm.out", ios::out | ios::app);
           fff << endl;
           fff << "Cumulative distribution of maximum test statistic value in "
               << Nrun << " permutations or simulations" << endl;
           for (i=1; i<Nrun+1; i++)
           {   cout << setprecision(8) << ( (double)i/( (double)Nrun+1.0) )
                    << " " << setprecision(8) << maxF[i-1] << endl;
               fff << setprecision(8) << ( (double)i/( (double)Nrun+1.0) )
                   << " " << setprecision(8) << maxF[i-1] << endl;
           }
           fff.close();
        } */
        delete[] urand;
        delete[] indorder;
        delete[] yoriginal;
        delete[] maxF;
        run= -1;
     }
     delete[] position;
     delete[] weight;
     delete[] ind;
     delcmatrix(marker,Nmark);
     delete[] y;
     delete[] selcofactor;
     delmatrix(XtWX,dimx+2);
     delcmatrix(Xt,dimx+2);
     delete[] XtWY;
     //delete[] idum;
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

/* 
 * augmentdata inserts missing data - does not use the phenotypes - augments the
 * marker data. Phenotype checking is done in the EM step.
 */

void augmentdata(cmatrix marker, vector y, cmatrix &augmarker, vector &augy,
     ivector &augind, int &Nind, int &Naug, int Nmark,
     cvector position, vector r)
{    int jj;
     int newNind=Nind;
     Naug= maxNaug; /* maximum size of augmented dataset */
     cmatrix newmarker;
     vector newy;
     cvector imarker;
     ivector newind;
     newmarker= newcmatrix(Nmark+1,Naug);
     newy= newvector(Naug);
     newind= newivector(Naug);
     imarker= newcvector(Nmark);
     int iaug=0;      // iaug keeps track of current augmented individual
     int maxiaug=0;   // highest reached(?)
     int saveiaug=0;  // previous iaug
     double prob0, prob1, prob2, sumprob,
            prob0left, prob1left, prob2left,
            prob0right, prob1right, prob2right;
     double probmax;
     vector newprob, newprobmax;
     newprob= newvector(Naug);
     newprobmax= newvector(Naug);
     cout << "maximum Naug=" << Naug << endl;
     OK();
     // ---- foreach individual create one in the newmarker matrix
     for (int i=0; i<Nind; i++)
     {   newind[iaug]=i-(Nind-newNind);  // index of individuals
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
                            newind[iaug]=i-(Nind-newNind);
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
                            newind[iaug]=i-(Nind-newNind);
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
                              newind[iaug]=i-(Nind-newNind);
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
                              newind[iaug]=i-(Nind-newNind);
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
                              newind[iaug]=i-(Nind-newNind);
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
                              newind[iaug]=i-(Nind-newNind);
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
                              newind[iaug]=i-(Nind-newNind);
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
                              newind[iaug]=i-(Nind-newNind);
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
                              newind[iaug]=i-(Nind-newNind);
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
                              newind[iaug]=i-(Nind-newNind);
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
                   {  cout << "warning in augmentation routine: dataset too large after augmentation" << endl;
                      cout << "recall procedure with larger value for parameter neglect or maxNaug";
                      exit(1);
                   }
               }
             if ((iaug-saveiaug+1)>imaxNaug)
             {  newNind-= 1;
                iaug= saveiaug-1;
                cout << "individual " << i << " is eliminated, because it is not informative enough" << endl;
                ofstream fff("mqm_out.txt", ios::out | ios::app);
                fff << "individual " << i << " is eliminated, because it is not informative enough" << endl;
                fff.close();
             }

             sumprob= 0.0;
             for (int ii=saveiaug; ii<=iaug; ii++) sumprob+= newprob[ii];
             for (int ii=saveiaug; ii<=iaug; ii++) newprob[ii]/= sumprob;
         }
         iaug++;
         saveiaug=iaug;
     }
     cout << "data augmentation ready" << endl;
     Naug= iaug;
     Nind= newNind;
     augmarker= newcmatrix(Nmark,Naug);
     augy= newvector(Naug);
     augind= newivector(Naug);
     for (int i=0; i<Naug; i++)
     {   augy[i]= newy[i];
         augind[i]= newind[i];
         for (int j=0; j<Nmark; j++) augmarker[j][i]= newmarker[j][i];
     }
     delcmatrix(newmarker,Nmark);
     delete[] newy;
     delete[] newind;
     delete[] newprob;
     delete[] newprobmax;
     delete[] imarker;
}


/* ML estimation of recombination frequencies via EM;
   calculation of multilocus genotype probabilities;
   ignorance of unlikely genotypes*/
void rmixture(cmatrix marker, vector weight, vector r,
              cvector position, ivector ind,
              int Nind, int Naug, int Nmark)
{    int i,j;
     int iem= 0;
     double Nrecom, oldr=0.0, newr, rdelta=1.0;
     vector indweight;
     indweight = newvector(Nind);
     char rknown='n';
     for (j=0; j<Nmark; j++)
     if (r[j]!=999.0) rknown='y';
     if (rknown=='y')
     {  cout << "recombination parameters are not re-estimated" << endl;
        OK();
        if (ok=='1') rknown='n';
     }
     while ((iem<100)&&(rdelta>0.001))
     {     iem+=1;
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
           for (i=0; i<Nind; i++) indweight[i]= 0.0;
           for (i=0; i<Naug; i++) indweight[ind[i]]+=weight[i];
           for (i=0; i<Naug; i++) weight[i]/=indweight[ind[i]];

           for (j=0; j<Nmark; j++)
           {   if ((position[j]=='L')||(position[j]=='M'))
               {  newr= 0.0;
                  for (i=0; i<Naug; i++)
                  {   Nrecom= absdouble((double)marker[j][i]-marker[j+1][i]);
                      if ((marker[j][i]=='1')&&(marker[j+1][i]=='1'))
                         Nrecom= 2.0*r[j]*r[j]/(r[j]*r[j]+(1-r[j])*(1-r[j]));
                      newr+= Nrecom*weight[i];
                  }
                  if (rknown=='n')
                  {  oldr=r[j];
                     r[j]= newr/(2.0*Nind);
                     rdelta+=pow(r[j]-oldr,2.0);
                  }
                  else rdelta=0.0;
               }
            }
     }
/*   print new estimates of recombination frequencies */
     cout << "iem= " << iem << " rdelta= " << rdelta << endl;
     if (rknown=='n')
     {  for (j=0; j<Nmark; j++)
        if ((position[j]=='L')||(position[j]=='M'))
        cout << "r(" << setw(2) << j << ")= " << setw(15) << r[j] << endl;
        OK();
     }
     delete[] indweight;
}


/* ML estimation of parameters in mixture model via EM;
*/
double QTLmixture(cmatrix loci, cvector cofactor, vector r, cvector position,
              vector y, ivector ind, int Nind, int Naug,
              int Nloci,
              double &variance, int em, vector &weight)
{    // cout << "QTLmixture IN" << endl;

     int iem= 0, newNaug, i, j;
     char varknown, biasadj='n';
     double Nrecom, oldlogL=-10000, delta=1.0, calc_i, logP=0.0, Pscale=1.75;
     vector indweight, Ploci, Fy;
     indweight= newvector(Nind);
     newNaug= (fitQTL=='n' ? Naug : 3*Naug);
     Ploci= newvector(newNaug);
     Fy= newvector(newNaug);
     logP= Nloci*log(Pscale); // only for computational accuracy
     varknown= ((variance==-1.0) ? 'n' : 'y' );
//     if ((REMLorML=='0')&&(varknown=='n')) cout << "variance is being estimated and bias adjusted" << endl;
     if (REMLorML=='1') { varknown='n'; biasadj='n'; }

     /* calculate weights (= conditional genotype probabilities) given
     marker data only */
     for (i=0; i<newNaug; i++) Ploci[i]= 1.0;
     if (fitQTL=='n')
     for (j=0; j<Nloci; j++)
     {    for (i=0; i<Naug; i++)
              Ploci[i]*= Pscale; // only for computational accuracy; see use of logP
          if ((position[j]=='L')||(position[j]=='U'))
          for (i=0; i<Naug; i++) Ploci[i]*= (loci[j][i]=='1' ? 0.5 : 0.25);
          if ((position[j]=='L')||(position[j]=='M'))
          {  for (i=0; i<Naug; i++)
             {  Nrecom= absdouble((double)loci[j][i]-(double)loci[j+1][i]);
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
     else
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
     if (weight[0]== -1.0)
     {  for (i=0; i<Nind; i++) indweight[i]= 0.0;
        if (fitQTL=='n')
        {  for (i=0; i<Naug; i++) indweight[ind[i]]+=Ploci[i];
           for (i=0; i<Naug; i++) weight[i]= Ploci[i]/indweight[ind[i]];
        }
        else
        {  for (i=0; i<Naug; i++) indweight[ind[i]]+=Ploci[i]+Ploci[i+Naug]+Ploci[i+2*Naug];
           for (i=0; i<Naug; i++)
           {   weight[i]       = Ploci[i]/indweight[ind[i]];
               weight[i+Naug]  = Ploci[i+Naug]/indweight[ind[i]];
               weight[i+2*Naug]= Ploci[i+2*Naug]/indweight[ind[i]];
           }
        }
     }
     double logL=0;
     vector indL;

     indL= newvector(Nind);
     while ((iem<em)&&(delta>1.0e-5))
     {     // cout << "EM algorithm, cycle " << iem << endl;
           iem+=1;
           if (varknown=='n') variance=-1.0;
           logL= regression(Nind, Nloci, cofactor, loci, y,
                 weight, ind, Naug, variance, Fy, biasadj);
           logL=0.0;
           // cout << "regression ready" << endl;
           for (i=0; i<Nind; i++) indL[i]= 0.0;
           if (fitQTL=='n') // no QTL fitted
           for (i=0; i<Naug; i++)
           {   weight[i]= Ploci[i]*Fy[i];
               indL[ind[i]]+=weight[i];
           }
           else // QTL moved along the chromosomes
           for (i=0; i<Naug; i++)
           {  weight[i]= Ploci[i]*Fy[i];
              weight[i+Naug]  = Ploci[i+Naug]*  Fy[i+Naug];
              weight[i+2*Naug]= Ploci[i+2*Naug]*Fy[i+2*Naug];
              indL[ind[i]]+=weight[i]+weight[i+Naug]+weight[i+2*Naug];
           }
           for (i=0; i<Nind; i++) logL+=log(indL[i])-logP;
           for (i=0; i<Nind; i++) indweight[i]= 0.0;
           if (fitQTL=='n')
           {  for (i=0; i<Naug; i++) indweight[ind[i]]+=weight[i];
              for (i=0; i<Naug; i++) weight[i]/=indweight[ind[i]];
           }
           else
           {  for (i=0; i<Naug; i++)
                  indweight[ind[i]]+=weight[i]+weight[i+Naug]+weight[i+2*Naug];
              for (i=0; i<Naug; i++)
              {   weight[i]       /=indweight[ind[i]];
                  weight[i+Naug]  /=indweight[ind[i]];
                  weight[i+2*Naug]/=indweight[ind[i]];
              }
           }
           delta= absdouble(logL-oldlogL);
           oldlogL= logL;
     }
     // cout << "EM finished" << endl;
     // bias adjustment after finished ML estimation via EM
     if ((REMLorML=='0')&&(varknown=='n'))
     {  variance=-1.0;
        biasadj='y';
        logL= regression(Nind, Nloci, cofactor, loci, y,
              weight, ind, Naug, variance, Fy, biasadj);
        logL=0.0;
        for (int _i=0; _i<Nind; _i++) indL[_i]= 0.0;
        if (fitQTL=='n')
        for (i=0; i<Naug; i++)
        {   weight[i]= Ploci[i]*Fy[i];
            indL[ind[i]]+=weight[i];
        }
        else
        for (i=0; i<Naug; i++)
        {   weight[i]= Ploci[i]*Fy[i];
            weight[i+Naug]= Ploci[i+Naug]*Fy[i+Naug];
            weight[i+2*Naug]= Ploci[i+2*Naug]*Fy[i+2*Naug];
            indL[ind[i]]+=weight[i];
            indL[ind[i]]+=weight[i+Naug];
            indL[ind[i]]+=weight[i+2*Naug];
        }
        for (i=0; i<Nind; i++) logL+=log(indL[i])-logP;
        for (i=0; i<Nind; i++) indweight[i]= 0.0;
        if (fitQTL=='n')
        {  for (i=0; i<Naug; i++) indweight[ind[i]]+=weight[i];
           for (i=0; i<Naug; i++) weight[i]/=indweight[ind[i]];
        }
        else
        {  for (i=0; i<Naug; i++)
           {   indweight[ind[i]]+=weight[i];
               indweight[ind[i]]+=weight[i+Naug];
               indweight[ind[i]]+=weight[i+2*Naug];
           }
           for (i=0; i<Naug; i++)
           {   weight[i]       /=indweight[ind[i]];
               weight[i+Naug]  /=indweight[ind[i]];
               weight[i+2*Naug]/=indweight[ind[i]];
           }
        }
     }
     // cout << "; iem=" << iem << "; delta=" << delta << "; variance=" << variance;
     // cout << "; logL=" << setprecision(8) << logL << endl;
     // cout << "QTLmixture OUT" << endl;
     delete[] Fy;
     delete[] Ploci;
     delete[] indweight;
     delete[] indL;
     return logL;
}

/* regression of trait on multiple cofactors
   y=xb+e with weight w
   (xtwx)b=(xtw)y
   b=inv(xtwx)(xtw)y */

double regression(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y,
                vector weight, ivector ind, int Naug,
                double &variance, vector Fy, char biasadj)
{    // cout << "regression IN" << endl;
     /*
     cofactor[j] at locus j:
     '0': no cofactor at locus j
     '1': cofactor at locus j
     '2': QTL at locus j, but QTL effect is not included in the model
     '3': QTL at locu j and QTL effect is included in the model
     */
     int dimx=1, j, jj;
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

     // cout << "calculate xtwx and xtwy" << endl;
     /* calculate xtwx and xtwy */
     double xtwj, yi, wi, calc_i;
     for (j=0; j<dimx; j++)
     {   XtWY[j]= 0.0;
         for (jj=0; jj<dimx; jj++) XtWX[j][jj]= 0.0;
     }


     if (fitQTL=='n')
     for (int i=0; i<Naug; i++)
     {  yi= y[i];
        wi= weight[i];
        for (j=0; j<dimx; j++)
        {   xtwj= ((double)Xt[j][i]-48.0)*wi;
            XtWY[j]+= xtwj*yi;
            for (jj=0; jj<=j; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
        }
     }
     else // QTL is moving along the chromosomes
     for (int i=0; i<Naug; i++)
     {   wi= weight[i]+ weight[i+Naug]+ weight[i+2*Naug];
         yi= y[i];
         for (j=0; j<dimx; j++)
         if (xtQTL[j]<='1')
         {  xtwj= ((double)Xt[j][i]-48.0)*wi;
            XtWY[j]+= xtwj*yi;
            for (jj=0; jj<=j; jj++)
            if (xtQTL[jj]<='1') XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
            else if (xtQTL[jj]=='2') // QTL: additive effect if QTL='0' or '2'
            {  // QTL=='0'
               XtWX[j][jj]+= ((double)(Xt[j][i]-48.0))*weight[i]*(47.0-48.0);
               // QTL=='2'
               XtWX[j][jj]+= ((double)(Xt[j][i]-48.0))*weight[i+2*Naug]*(49.0-48.0);
            }
            else // (xtQTL[jj]=='3')  QTL: dominance effect only if QTL='1'
            {  // QTL=='1'
               XtWX[j][jj]+= ((double)(Xt[j][i]-48.0))*weight[i+Naug]*(49.0-48.0);
            }
         }
         else if (xtQTL[j]=='2') // QTL: additive effect if QTL='0' or '2'
         {  xtwj= -1.0*weight[i]; // QTL=='0'
            XtWY[j]+= xtwj*yi;
            for (jj=0; jj<j; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
            XtWX[j][j]+= xtwj*-1.0;
            xtwj= 1.0*weight[i+2*Naug]; // QTL=='2'
            XtWY[j]+= xtwj*yi;
            for (jj=0; jj<j; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
            XtWX[j][j]+= xtwj*1.0;
         }
         else // (xtQTL[j]=='3') QTL: dominance effect only if QTL='1'
         {  xtwj= 1.0*weight[i+Naug]; // QTL=='1'
            XtWY[j]+= xtwj*yi;
            // j-1 is for additive effect, which is orthogonal to dominance effect
            for (jj=0; jj<j-1; jj++) XtWX[j][jj]+= xtwj*((double)Xt[jj][i]-48.0);
            XtWX[j][j]+= xtwj*1.0;
         }
     }
     for (j=0; j<dimx; j++)
     for (jj=j+1; jj<dimx; jj++) XtWX[j][jj]= XtWX[jj][j];

     /* solve equations */
     // cout << "solve equations" << endl;
     // printmatrix(XtWX,dimx,dimx);
     int d;
     ivector indx;
     indx= newivector(dimx);
     ludcmp(XtWX,dimx,indx,d);
     lusolve(XtWX,dimx,indx,XtWY);
     // luinvert(xtwx, inv, dimx, indx);
     // cout << "Parameter estimates" << endl;
     // for (jj=0; jj<dimx; jj++) cout << jj << " " << XtWY[jj] << endl;

     long double *indL;
     indL= new long double[Nind];
     int newNaug;
     newNaug= (fitQTL=='n' ? Naug : 3*Naug);
     vector fit, resi;
     fit= newvector(newNaug);
     resi= newvector(newNaug);
     // cout << "Calculate residuals" << endl;
     if (variance<0)
     {    variance= 0.0;
          if (fitQTL=='n')
          for (int i=0; i<Naug; i++)
          {   fit[i]= 0.0;
              for (j=0; j<dimx; j++)
              fit[i]+=((double)Xt[j][i]-48.0)*XtWY[j];
              resi[i]= y[i]-fit[i];
              variance+=weight[i]*pow(resi[i],2.0);
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
              variance+=weight[i]*pow(resi[i],2.0);
              variance+=weight[i+Naug]*pow(resi[i+Naug],2.0);
              variance+=weight[i+2*Naug]*pow(resi[i+2*Naug],2.0);
          }
          variance/= (biasadj=='n' ? Nind : Nind-dimx); // to compare results with Johan; variance/=Nind;
          if (fitQTL=='n')
          for (int i=0; i<Naug; i++) Fy[i]= Lnormal(resi[i],variance);
          else
          for (int i=0; i<Naug; i++)
          {   Fy[i]       = Lnormal(resi[i],variance);
              Fy[i+Naug]  = Lnormal(resi[i+Naug],variance);
              Fy[i+2*Naug]= Lnormal(resi[i+2*Naug],variance);
          }
     }
     else
     {    if (fitQTL=='n')
          for (int i=0; i<Naug; i++)
          {   fit[i]= 0.0;
              for (j=0; j<dimx; j++)
              fit[i]+=((double)Xt[j][i]-48.0)*XtWY[j];
              resi[i]= y[i]-fit[i];
              Fy[i]  = Lnormal(resi[i],variance); // ????
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
              Fy[i]       = Lnormal(resi[i],variance);
              Fy[i+Naug]  = Lnormal(resi[i+Naug],variance);
              Fy[i+2*Naug]= Lnormal(resi[i+2*Naug],variance);
          }
     }
     delete[] fit;
     delete[] resi;

     /* calculation of logL */
     // cout << "calculate logL" << endl;
     long double logL=0.0;
     for (int i=0; i<Nind; i++) indL[i]= 0.0;
     if (fitQTL=='n')
     for (int i=0; i<Naug; i++) indL[ind[i]]+=weight[i]*Fy[i];
     else
     for (int i=0; i<Naug; i++)
     {   indL[ind[i]]+=weight[i]*       Fy[i];
         indL[ind[i]]+=weight[i+Naug]*  Fy[i+Naug];
         indL[ind[i]]+=weight[i+2*Naug]*Fy[i+2*Naug];
     }
     for (int i=0; i<Nind; i++) logL+= log(indL[i]);

     // if (biasadj=='y') cout << "Nind= " << Nind << " Degrees of Freedom (df)= " << (Nind-dimx) << endl;
     // cout << "regression OUT" << endl;
     delete[] indL;
     delete[] indx;
     delete[] xtQTL;
     return (double)logL;
}


/* backward elimination in regression of trait on multiple cofactors
   routine subX haalt uit matrices voor volledige model de submatrices voor submodellen;
   matrices XtWX en Xt van volledig model worden genoemd fullxtwx en fullxt;
   analoog vector XtWY wordt full xtwy genoemd;
*/
double backward(int Nind, int Nmark, cvector cofactor, cmatrix marker, vector y, vector weight,
              int* ind, int Naug, double logLfull, double variance,
              double F1, double F2, cvector &newcofactor, vector r, cvector position)
{    int dropj=0, Ncof=0;
     double maxlogL, savelogL, maxF=0.0; //, minlogL=logLfull, maxFtest=0.0;
     char finished='n'; //, biasadj='n';
     vector logL;
     logL = newvector(Nmark);
     savelogL= logLfull;
     maxlogL= logLfull-10000;
     for (int j=0; j<Nmark; j++)
     {   newcofactor[j]= cofactor[j];
         Ncof+=(cofactor[j]!='0');
     }
     while ((Ncof>0)&&(finished=='n'))
     {     for (int j=0; j<Nmark; j++)
           {   if (newcofactor[j]=='1')
               {  cout << "drop marker " << j;
                  newcofactor[j]='0';
                  if (REMLorML=='1') variance= -1.0;
                  logL[j]= QTLmixture(marker,newcofactor,r,position,y,ind,Nind,Naug,Nmark,variance,em,weight);
                  newcofactor[j]='1';
               }
               else if (newcofactor[j]=='2')
               {  cout << "drop marker " << j;
                  newcofactor[j]='0';
                  if (REMLorML=='1') variance= -1.0;
                  logL[j]=  QTLmixture(marker,newcofactor,r,position,y,ind,Nind,Naug,Nmark,variance,em,weight);
                  newcofactor[j]='2';
               }
               else if (newcofactor[j]!='0') cout << " something is wrong ";
           }
           /* nu bepalen welke cofactor 0 kan worden (=verwijderd) */
           maxlogL= logLfull-10000.0;
           for (int j=0; j<Nmark; j++)
           {   if (newcofactor[j]!='0')
               if (logL[j]>maxlogL) { maxlogL= logL[j]; dropj = j; }
           }
           if  ( (newcofactor[dropj]=='1') && ( F2> 2.0*(savelogL-maxlogL)) )
           {   savelogL= maxlogL;
               newcofactor[dropj]= '0'; Ncof-=1;
               cout << "marker "<< dropj << " is dropped; logL of reduced model = " << savelogL << endl;
           }
           else if  ( (newcofactor[dropj]=='2') && (F1> 2.0*(savelogL-maxlogL)) )
           {   savelogL= maxlogL;
               newcofactor[dropj]= '0'; Ncof-=1;
               cout << "marker "<< dropj << " is dropped; logL of reduced model = " << savelogL << endl;
           }
           else /* ready */
           {   finished='y';
               for (int j=0; j<Nmark; j++)
               if (newcofactor[j]=='1') cout << "marker " << j << " is selected" << endl;
               OK();
           }
     }
     for (int j=0; j<Nmark; j++)
     if (newcofactor[j]!='0') cout << "marker " << j << " is in final model" << endl;

     maxF= mapQTL(Nind, Nmark, cofactor, newcofactor, marker, position,
           mapdistance, y, r, ind, Naug, variance, 'n'); // printoutput='n'
     cout << "backward selection finished" << endl;
     delete[] logL;
     return maxF;
}

/* mapQTL */
double mapQTL(int Nind, int Nmark, cvector cofactor,
       cvector selcofactor, cmatrix marker, cvector position, vector mapdistance,
       vector y, vector r, ivector ind, int Naug, double variance,
       char printoutput)
{      int Nloci, j, jj, jjj=0;
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
       newcofactor= new char[Nmark];
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


       ofstream fff("mqm_out.txt", ios::out | ios::app);
       cout << endl << endl;

       /* fit QTL on top of markers (full ML)
          fit QTL between markers (full ML) */
       // cout << "please wait (mixture calculus may take quite a lot of time)" << endl;
       /* estimate variance in mixture model with all marker cofactors */
       // cout << "estimate variance in mixture model with all cofactors" << endl;
       variance= -1.0;
       savelogL= 2.0*QTLmixture(marker,cofactor,r,position,
                     y,ind,Nind,Naug,Nmark,variance,em,weight);
       Nloci= Nmark+1;
       // augment data for missing QTL observations (x 3)
       fitQTL='y';
       int newNaug;
       newNaug= 3*Naug;
       delete[] weight;
       weight= newvector(newNaug);
       weight[0]= 1.0;
       vector weight0;
       weight0= newvector(newNaug);
       weight0[0]= -1.0;

//       augmentdataforQTL(marker);
       vector QTLr, QTLmapdistance;
       QTLr= newvector(Nloci);
       QTLmapdistance= newvector(Nloci);
       cvector QTLposition;
       QTLposition= newcvector(Nloci);
       cmatrix QTLloci;
       QTLloci= new char*[Nloci];
       double moveQTL= stepmin;
       char nextinterval= 'n', firsttime='y';
       double maxF=0.0, savebaseNoQTLModel=0.0;
       int baseNoQTLModel=0, step=0;

       for (j=0; j<Nmark; j++)
       {   /* fit a QTL in two steps:
              1. move QTL along marker interval j -> j+1 with steps
                 of stepsize=20 cM, starting from -20 cM up to 220 cM
              2. all marker-cofactors in the neighborhood of the QTL
                 are dropped by using cM='windows' as criterium

           */
         nextinterval= 'n';
         while (nextinterval=='n')
         { // step 1:
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
                 QTLlikelihood= -2.0*
                 QTLmixture(QTLloci,QTLcofactor,QTLr,QTLposition,y,ind,Nind,Naug,Nloci,variance,em,weight0);
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
              QTLlikelihood+=2.0*
              QTLmixture(QTLloci,QTLcofactor,QTLr,QTLposition,y,ind,Nind,Naug,Nloci,variance,em,weight);
              if (QTLlikelihood<-0.05) { cout << "error QTLlikelihood=" << QTLlikelihood <<endl; exit(1); }
              maxF= (maxF<QTLlikelihood ? QTLlikelihood : maxF);
              if (run>=0) Frun[step][run]+= QTLlikelihood;
              else Frun[step][Nrun]+= QTLlikelihood;

              /* Each individual has condition multilocus probabilities
                 for being 0, 1 or 2 at the QTL.
                 Calculate the maximum per individu.
                 Calculate the mean of this maximum, averaging over all individuals
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
              informationcontent[step]+=infocontent/Nind;
              step++;
              if (printoutput=='y')
              {  cout << j << " " << chr[j] << " " << setprecision(0) << (moveQTL-stepsize) << " "
                   << setprecision(5) << QTLlikelihood << " " << (infocontent/Nind) << endl;
                 fff << chr[j] << " " << setprecision(0) << (moveQTL-stepsize) << " " <<  setprecision(5)
                   << QTLlikelihood << " " << setprecision(2) << (infocontent/Nind) << endl;
              }
           }
         }
       }
       if (printoutput=='y')
       fff << ":" << endl; // genstat code for end of data
       fff.close();
       OK();
       fitQTL='n';
       delete[] info0;
       delete[] info1;
       delete[] info2;
       delete[] weight;
       delete[] weight0;
       delete[] QTLr;
       delete[] QTLposition;
       delete[] Fy;
       delete[] newcofactor;
       delete[] QTLcofactor;
       delete[] cumdistance;
       delete[] QTLmapdistance;
       return maxF; //QTLlikelihood;
}




/* LU decomposition (from Johan via Numerical Recipes in C)
   Given an n x n matrix a[1..n][1..n], this routine replaces it by the LU
   decomposition of a rowwise permutation of itself.
   A and n are input.
   a is output. indx[1..n] is an output vector which records the row
   permutation effected by the partial pivoting; d is output as +-1
   depending on whether the number of row interchanges was even or odd,
   respectively. This routine is used in combination with lusolve to solve
   linear equations or to invert a matrix.
*/
void ludcmp(matrix m, int dim, ivector ndx, int &d)
{   int r,c,rowmax,i;
    double max,temp,sum;
    vector scale, swap;
    scale= newvector(dim);

    d=1;
    for (r=0; r<dim; r++)
    {   for (max=0.0, c=0; c<dim; c++) if ((temp=fabs(m[r][c])) > max) max=temp;
        if (max==0.0) {cout << "singular matrix"; exit(1);}
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
        if (max==0.0) {cout << "singular matrix"; exit(1);}
        if (rowmax!=c)
        {  swap=m[rowmax]; m[rowmax]=m[c]; m[c]=swap;
           scale[rowmax]=scale[c]; d= -d;
        }
        ndx[c]=rowmax;
        temp=1.0/m[c][c];
        for (r=c+1; r<dim; r++) m[r][c]*=temp;
    }
    delete[] scale;
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
      cout << "a or b too big or max number of iterations too small";
      exit(1); return 0.0;
}

double betai(double a, double b, double x)
{     double bt;
      if (x<0.0 || x>1.0) { cout << "x not between 0 and 1"; exit(1); }
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
       cout << "prob=" << prob << "; alfa=" << alfa << endl;
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
/*
void nrerror(char *error_text)
{    cout << stderr << "Numerical Recipes run-time error" << endl;
     cout << stderr << error_text;
     cout << stderr << "now exiting to system" << endl;
     exit(1);
}

void luinvert(double **lu, double **inv, int dim, int *ndx)
{    int r,c;
     double *b;
     b= new double[dim];
     for (c=0; c<dim; c++)
     {   for (r=0; r<dim; r++) b[r]=0.0;
         b[c]=1.0;
         lusolve(lu,dim,ndx,b);
         for (r=0; r<dim; r++) inv[r][c]= b[r];
     }
     delete[] b;
} */
