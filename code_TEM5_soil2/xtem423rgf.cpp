
/************************************************************************************************************
************************************************************************************************************
XTEM423EDF.CPP - Terrestrial Ecosystem Model Version 4.2
		MPI version for GloBays
************************************************************************************************************

Modifications:

19991028 - DWK added bug fixes found and implemented by Gregg
           Christopher and Jim Long
20000616 - DWK eliminates the global variables tempred[][][] and atmspred[][]
20000616 - DWK eliminates the task of writing TEM output files from
           extrapolate()
20000616 - DWK adds the global constant double ZERO
20010314 - JSC adds transient solar radiation as input possibility
20010418- Q.Z. addition of soil thermal model
20020202 - DWK changed include from xtem423e.cpp to xtem423e1.cpp
200612   - J.Tang removed deprecated included files
20080305 - J.Tang changed it into mpi version
20080324 - J.Tang changed the code to read name of parafile and log file as argument
************************************************************************************************************
********************************************************************************************************* */

#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <fcntl.h>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cctype>
#include <cstring>
#include <ctime>
//#include <conio.h>
#include <values.h>


using namespace std;


const int CYCLE = 12;
const int MAXRTIME = 600;

const int MAXPRED = 92+4+10; // changed from 85 to 92 by DWK on 20000202,
							// QZ, plus 4+10
const int MAXNPAR = 80;	//Maximum number of parameters to be tested
const int MAXGRID = 1;
// const int MAXNUMROW = 360; // commented out by DWK on 20000202

const double MISSING = -999999.9;
const double ZERO = 0.000000;  // added by DWK on 20000615

long int kswitch; //added for soil thermal model
int stmflg;   // added for soil thermal model

ofstream flog1;
ifstream fpara;

#if !defined(ELMNT423_C)
 #include "elmnt423.cpp"    //Elmnt Class
#endif

#if !defined(TELM423E1_C) 
 #include "telm423e1.cpp"     //TEMelmnt Class
#endif 

struct Temflags{//claimed for MPI version
	int stmflg;   // added for soil thermal model
	int equil;
	int RTIME;

	int spinflag;
	int numspin;
	int spintime;
	int ispinout;  //added by J. Tang to control output from spin up period
	int totsptime;
	int transtime;

//	int atmstotyr;//not used;  commented out by J.Tang Oct.2006
	int atmsflag;
	int atmsoutfg;
	int temflag;
	int stateflag;

	int cldflag;
	int natmspred;
	int ntempred;
	int totpred;
	char predmap[MAXPRED][9];
	char mdir[60];
	int lhflag[MAXNPAR];
	size_t whichcmnt;//added by J. Tang to do latin hypercube analysis
	int kdinflg;
	int kdoutflg;
	int f1stnd;
//	int glob_count;  commented out by J. Tang Oct. 2006

}tmflgs;


struct Temstac//added for mpi-tem
{
	int avlnflag;
	int nfeed;
	int initbase;
	int baseline;
	int moistlim;
	int strteq;
	int maxyears;
	int runsize;
	int maxnrun;
	int rheqflag;
	double wtol;
	double ctol;
	double ntol;
	int startyr;
	int endyr;
	int diffyr;
	
	double inittol;
	int maxit;
	long maxitmon;
	
}temstac;

struct Fnminout{
   char fnmlonlat[60];
   char fnmstxt[60];
   char fnmelev[60];
   char fnmtveg[60];
   char fnmclds[60];
   char fnmtair[60];
   char fnmprec[60];
   char fnmnirr[60];
   char fnmpar[60];
   char fnmlulc[60];
   char fnmnpp[60];
   char fnmstate[MAXSTATE][60];
   char fnmclmpred[NUMATMS][60];
   char fnmtempred[MAXPRED][60];
   char fnmkdin[60];
   char fnmkdout[60];
   
}fnminout;




Elmnt elmnt;
TEMelmnt telmnt[MAXGRID];


/* *************************************************************
		 Function Declarations
************************************************************* */

int askpred(char pvarname[MAXPRED][9], int spred);
void read_dr(void);	//read driving data
void open_dr(int rank);
void cpflags(const int opt);
void extrapolate(void);//do simulation
void startclm(void);
void starttem(void);

//int preparavar(ofstream& outfl, const vector<int>& flgs);//added by J. Tang 20061208
// *************************************************************

int equil;
int RTIME;

int spinflag;
int numspin;
int spintime;
int totsptime;
int transtime;
int ispinout;

int atmstotyr;
int atmsflag;
int atmsoutfg;
int temflag;
int stateflag;

int cldflag;
int natmspred;
int ntempred;
int totpred;
char predmap[MAXPRED][9];
char mdir[60];	//added for mpi-lhtem
int kdinflg;
int kdoutflg;
int lhflag[MAXNPAR] = {0};
unsigned int whichcmnt;//added by J. Tang to do latin hypercube analysis

int fatalerr;
int end1;
int icount;
int glob_count;
int f1stnd;


char fnmclmpred[NUMATMS][60];
char fnmtempred[MAXPRED][60];
char fnmkdout[60];

FILE* flonlat;
FILE* fstxt;
FILE* felev;
FILE* ftveg;
FILE* fclds;
FILE* ftair;
FILE* fprec;
FILE* fnirr;
FILE* fpar;
FILE* flulc;
FILE* fnpp;
FILE* fstate[MAXSTATE];
ofstream fclmpred[NUMATMS];
ofstream ftempred[MAXPRED];

FILE* fkdin;
ofstream fkdout;
KDdata kdparam;
ofstream otemls;
/* *************************************************************
**********************START MAIN PROGRAM************************
************************************************************* */

int main(int argc,char* argv[])
{
	  
    int i;

    int root = 0;
    char strtmp[160];
    MPI::Init(argc,argv);
  	
    int tasks = MPI::COMM_WORLD.Get_size();
    int myrank0 = MPI::COMM_WORLD.Get_rank();
    int myrank;
    if(argc != 3)
    {
         cout <<"Please run tem as ./xtem tem4.para tem4.log"<<endl;

         MPI::Finalize();
         exit(-1);
    }
    if(myrank0 == root)
    {
         //proc 0, read necessary specifications for the model run
		
         flog1.open(argv[2]);
         fpara.open(argv[1]);

         if(! fpara.is_open())
         {
              cout <<"Cannot open file "<<argv[1]<<" to specify data"<<endl;		
              exit(-1);
         }
         if (! flog1.is_open())
         {
              cout << endl << "Cannot open log file "<<argv[2]<<endl;		
              exit(-1);
         }
	
/* *************************************************************
  Run equilibrium simulation or transient simulation ?
		************************************************************* */
		//
         cout <<endl<<" Enter the generic name of directory where the simulation output would be located. ";
         flog1 <<endl<<" Enter the generic name of directory where the simulation output would be located. ";
         fpara >>mdir;
         cout <<mdir<<endl;
         flog1 <<mdir<<endl;

         cout <<endl << "Start from which subfolder u want to run tem: ";
         flog1 <<endl << "Start from which subfolder u want to run tem: ";
         fpara >>f1stnd;
         cout <<f1stnd<<endl;
         flog1 <<f1stnd<<endl;

         equil = 1;
         cout << endl << "Do you want to run the model only for steady"
               << " state conditions ? " << endl;
         cout << " Enter 0 for transient simulation" << endl;
         cout << " Enter 1 for steady state simulation" << endl;

         fpara >> equil;

         flog1 << endl << "Do you want to run the model only for"
               << " steady state conditions ? " << endl;
         flog1 << " Enter 0 for transient simulation" << endl;
         flog1 << " Enter 1 for steady state simulation" << endl;
         flog1 << "equil = " << equil << endl << endl;

         RTIME = 1;
         ispinout = 0;
         if (equil == 0)
         {
    // Start transient TEM run from equilibrium conditions (i.e. spinflag == 0)
	// or with a "spin up" period to remove potential artifacts associated with
	// changing from equilibrium conditions to transient conditions from model
    // results

              cout << "Do you want to start the transient with a spin"
                   << " up period? " << endl;
              cout << "Enter 0 for no:" << endl;
              cout << "Enter 1 for yes: ";

              fpara >> spinflag;

              flog1 << "Do you want to start the transient with"
                   << " a spin up period? " << endl;
              flog1 << "Enter 0 for no:" << endl;
              flog1 << "Enter 1 for yes: " << endl;
              flog1 << "spinflag = " << spinflag << endl << endl;

              totsptime = 0;
              if (spinflag == 1)
              {
	  // Specify conditions for initializing TEM with
	  // a transient "spin up" period
      //   for a grid cell

                   cout << "How many spins do you want in the spin up period? ";

                   fpara >> numspin;
                   flog1 << "How many spins do you want in the spin up period? " << endl;
                   flog1 << "numspin = " << numspin << endl << endl;

                   cout << "How many years per spin? ";

                   fpara >> spintime;
                   flog1 << "How many years per spin? " << endl;
                   flog1 << "spintime = " << spintime << endl << endl;

                   cout << "Do you want to output the spin up tem data? 1: yes/0: no"<<endl;
                   fpara >> ispinout;
                   flog1 <<"Do you want to output the spin up tem data? 1: yes/0: no  ";
                   flog1 <<ispinout<<endl;
      
                   totsptime = spintime * numspin;
                   flog1 << "totsptime = " << totsptime << endl << endl;
                   RTIME += totsptime;
              }

    // Specify conditions for the "non-spin up" part of the transient TEM run

              cout << endl << "How many years do you run for"
                   << " transient simulations ? " << endl;

              fpara >> transtime;

              flog1 << endl << "How many years do you run for"
	                << " transient simulations ? " << endl;
              flog1 << "transtime = " << transtime << endl << endl;

              RTIME += transtime;
              flog1 << "RTIME = " << RTIME << endl << endl;
         }

/* *************************************************************
What models do you want to run?  - Open associated files
************************************************************* */
 //  STM (soil thermal model)?

         cout << endl << "Do you want to run the SOIL THERMAL"
              << " MODEL (STM) model for soil temperatures?" << endl;
         cout << "  Enter 0 for No" << endl;
         cout << "  Enter 1 for Yes" << endl;

         fpara >> stmflg;

         flog1 << endl << "Do you want to run the SOIL THERMAL"
               << " MODEL (STM) model for soil temperatures?"<< endl;
         flog1 << "  Enter 0 for No" << endl;
         flog1 << "  Enter 1 for Yes" << endl;
         flog1 << " stmflg = " << stmflg << endl << endl;
         if (stmflg == 1) 
         {
// Get snow and soil layer dependent parameters,
// added for soil thermal model integration
              telmnt[0].tem.sthermal.getsnowecd(flog1);   //qtsp44a2.ecd
              telmnt[0].tem.sthermal.getsoillecd(flog1);  // qtla44a2.ecd
              telmnt[0].tem.sthermal.getsoiltecd(flog1);  // qtst44a2.ecd

         }



// POTSCLM model ?

         cout << endl << "Do you want to run the POTSCLM model"
              << " for solar radiation variables?" << endl;
         cout << "  Enter 0 for No" << endl;
         cout << "  Enter 1 for Yes" << endl;

         fpara >> atmsflag;

         flog1 << endl << "Do you want to run the POTSCLM model"
               << " for solar radiation variables?" << endl;
         flog1 << "  Enter 0 for No" << endl;
         flog1 << "  Enter 1 for Yes" << endl;
         flog1 << " atmsflag = " << atmsflag << endl << endl;

         totpred = natmspred = 0;
         if (atmsflag == 1) { startclm(); }

// TEM ?

         cout << endl << "Do you want to run the terrestrial"
              << " ecosystem model (TEM)?" << endl;
         cout << "  Enter 0 for No" << endl;
         cout << "  Enter 1 for Yes" << endl;

         fpara >> temflag;

         flog1 << endl << "Do you want to run the terrestrial"
                << " ecosystem model (TEM)?" << endl;
         flog1 << "  Enter 0 for No" << endl;
         flog1 << "  Enter 1 for Yes" << endl;
         flog1 << "temflag = " << temflag << endl << endl;

         ntempred = 0;
         if (temflag == 1) { starttem(); }

  // Use all records in GIS data sets (i.e. desired coverage)?

         elmnt.ask(flog1);
         cout <<"Do latin hypercube sensitivity analysis for which community? Community";

         fpara >>whichcmnt;
         flog1 <<"Do latin hypercube sensitivity analysis for which community?  Community ";
         flog1 <<whichcmnt<<endl;
         cout <<whichcmnt<<endl; 
  // Run Potsclm and TTEM modules over desired region (modules could potentially
  // run over several grid cells simultaneously with the use of extrapolate()


//-----------------------------------------------
//
         size_t nvar=0;  //number of variables
         size_t ii;
         size_t jj;
         char  lhflname[40];
         ifstream lhtest; 
         cout <<"Enter name of the file where the parameters of interest are documented? "<<endl;
         flog1 <<"Enter name of the file where the parameters of interest are documented? "<<endl;
         fpara >>lhflname;
  
         cout <<lhflname<<endl;
         flog1 <<lhflname<<endl;
         sprintf(strtmp,"%s/%s",mdir, lhflname);
         lhtest.open(strtmp, ifstream::in);
         if(!lhtest.is_open())
         {
              cout <<"Failed to open file "<<strtmp<<endl;
              flog1 <<"Failed to open file "<<strtmp<<endl;
              exit(-1);	
         }
         char line[100];
		//read file head
         for(ii = 0; ii < 7; ii++)
              lhtest.getline(line,100);

		//read section head
         lhtest.getline(line,100);
         int lhlabel;
         char lhvarname[15];
         int lflg;
         double lhval1,lhval2;
		//initialize parameter counter
         jj = 0;
          //soil texture head
         for(ii = 0; ii < 6; ii ++)
         { 
              lhtest.getline(line,100);
              sscanf(line,"%d%s%d",&lhlabel,lhvarname,&lflg);
              lhflag[jj] = lflg;
              jj ++;
         }  
         //read section tail
         lhtest.getline(line,100);
         //read section head
         lhtest.getline(line,100);
         //root parameters
         for (ii = 0; ii < 4; ii ++)
         {
              lhtest.getline(line,100);
              sscanf(line,"%d%s%d",&lhlabel,lhvarname,&lflg);
              lhflag[jj] = lflg;

              jj ++;
         }

         //read section tail
         lhtest.getline(line,100);  
         //read section head
         lhtest.getline(line,100);
         //vegetation parameters
         for(ii = 0; ii < 18; ii ++)
         {
              lhtest.getline(line,100);
              sscanf(line,"%d%s%d",&lhlabel,lhvarname,&lflg);
              lhflag[jj] = lflg;

              jj ++;
         }

         //read section tail
         lhtest.getline(line,100);
         //read section head
         lhtest.getline(line,100);
         //leaf parameters
         for(ii = 0; ii < 4; ii ++)
         {
              lhtest.getline(line,100);
              sscanf(line,"%d%s%d",&lhlabel,lhvarname,&lflg);
              lhflag[jj] = lflg;
              jj ++;
         }

         //read section tail
         lhtest.getline(line,100);
         //read section head
         lhtest.getline(line,100);
         //microbial parameters
         for(ii = 0; ii < 5; ii ++)
         {
              lhtest.getline(line,100);
              sscanf(line,"%d%s%d",&lhlabel,lhvarname,&lflg);

              lhflag[jj] = lflg;
              jj ++;

         } 

         //read section tail
         lhtest.getline(line,100);
         //read section head
         lhtest.getline(line,100);
         //community parameters
         for(ii = 0; ii < 42; ii ++)
         {
              lhtest.getline(line,100);
              sscanf(line,"%d%s%d",&lhlabel,lhvarname,&lflg);
              lhflag[jj] = lflg;	
              jj ++;
			
         }
  
         lhtest.close();  
         for(ii = 0; ii < MAXNPAR; ii ++)
         {
              nvar += lhflag[ii];
         }
         flog1 <<endl<<nvar<<" variables are involved in the latin hypercube sensitivity test"<<endl;



//open and read sample data
	
         FILE*  fiprsmp;
         fiprsmp = fopen("para.smp", "r");
         if(!fiprsmp)
         {
              cout <<"Failed to open file "<<strtmp <<" for parameter input"<<endl;
              flog1 <<"Failed to open file "<<strtmp <<" for parameter input"<<endl;
              MPI::Finalize();
              exit(-1);
         }
  
         char line1[600];

         fgets(line1,600,fiprsmp);//read file head
         string instr;
         istringstream ins;//Declare an input string stream
	 //do regional extrapolation

	
         //read parameters
         strcpy(line1, "");
         if(fgets(line1,600,fiprsmp)==NULL)break;
         line1[strlen(line1)]='\0';
         instr = line1;
         ins.str(instr);
         if(lhflag[0] == 1)  ins >>telmnt[0].tem.soil.pctpora; 
         if(lhflag[1] == 1)  ins >>telmnt[0].tem.soil.pctporb;
         if(lhflag[2] == 1)  ins >>telmnt[0].tem.soil.fldcapa;
         if(lhflag[3] == 1)  ins >>telmnt[0].tem.soil.fldcapb;
         if(lhflag[4] == 1)  ins >>telmnt[0].tem.soil.wiltpta;
         if(lhflag[5] == 1)  ins >>telmnt[0].tem.soil.wiltptb;
		
         if(lhflag[6] == 1)  ins >>telmnt[0].tem.soil.rootza[whichcmnt]; 
         if(lhflag[7] == 1)  ins >>telmnt[0].tem.soil.rootzb[whichcmnt];
         if(lhflag[8] == 1)  ins >>telmnt[0].tem.soil.rootzc[whichcmnt]; 
         if(lhflag[9] == 1)  ins >>telmnt[0].tem.soil.minrootz[whichcmnt]; 
		
         if(lhflag[10] == 1) ins >>telmnt[0].tem.veg.kc[whichcmnt];
         if(lhflag[11] == 1) ins >>telmnt[0].tem.veg.ki[whichcmnt];
         if(lhflag[12] == 1) ins >>telmnt[0].tem.veg.gva[whichcmnt];
         if(lhflag[13] == 1) ins >>telmnt[0].tem.veg.tmin[whichcmnt];
         if(lhflag[14] == 1) ins >>telmnt[0].tem.veg.toptmin[whichcmnt];
         if(lhflag[15] == 1) ins >>telmnt[0].tem.veg.toptmax[whichcmnt];	
         if(lhflag[16] == 1) ins >>telmnt[0].tem.veg.tmax[whichcmnt];
         if(lhflag[17] == 1) ins >>telmnt[0].tem.veg.raq10a0[whichcmnt];
         if(lhflag[18] == 1) ins >>telmnt[0].tem.veg.raq10a1[whichcmnt];
         if(lhflag[19] == 1) ins >>telmnt[0].tem.veg.raq10a2[whichcmnt];
         if(lhflag[20] == 1) ins >>telmnt[0].tem.veg.raq10a3[whichcmnt];
         if(lhflag[21] == 1) ins >>telmnt[0].tem.veg.kn1[whichcmnt];
         if(lhflag[22] == 1) ins >>telmnt[0].tem.veg.labncon[whichcmnt];	
         if(lhflag[23] == 1) ins >>telmnt[0].tem.veg.leafmxc[whichcmnt];
         if(lhflag[24] == 1) ins >>telmnt[0].tem.veg.kleafc[whichcmnt];
         if(lhflag[25] == 1) ins >>telmnt[0].tem.veg.sla[whichcmnt];
         if(lhflag[26] == 1) ins >>telmnt[0].tem.veg.cov[whichcmnt];
         if(lhflag[27] == 1) ins >>telmnt[0].tem.veg.fpcmax[whichcmnt];
		
         if(lhflag[28] == 1) ins >>telmnt[0].tem.veg.minleaf[whichcmnt];
         if(lhflag[29] == 1) ins >>telmnt[0].tem.veg.aleaf[whichcmnt];
         if(lhflag[30] == 1) ins >>telmnt[0].tem.veg.bleaf[whichcmnt];
         if(lhflag[31] == 1) ins >>telmnt[0].tem.veg.cleaf[whichcmnt];
		
         if(lhflag[32] == 1) ins >>telmnt[0].tem.microbe.rhq10[whichcmnt];
         if(lhflag[33] == 1) ins >>telmnt[0].tem.microbe.kn2[whichcmnt];
         if(lhflag[34] == 1) ins >>telmnt[0].tem.microbe.moistmin[whichcmnt];
         if(lhflag[35] == 1) ins >>telmnt[0].tem.microbe.moistopt[whichcmnt];
         if(lhflag[36] == 1) ins >>telmnt[0].tem.microbe.moistmax[whichcmnt];
	 
         if(lhflag[37] == 1) ins >>telmnt[0].tem.vegca[whichcmnt];
         if(lhflag[38] == 1) ins >>telmnt[0].tem.vegcb[whichcmnt];
         if(lhflag[39] == 1) ins >>telmnt[0].tem.strna[whichcmnt];
         if(lhflag[40] == 1) ins >>telmnt[0].tem.strnb[whichcmnt];
         if(lhflag[41] == 1) ins >>telmnt[0].tem.solca[whichcmnt];
         if(lhflag[42] == 1) ins >>telmnt[0].tem.solcb[whichcmnt];
         if(lhflag[43] == 1) ins >>telmnt[0].tem.solna[whichcmnt];
         if(lhflag[44] == 1) ins >>telmnt[0].tem.solnb[whichcmnt];
         if(lhflag[45] == 1) ins >>telmnt[0].tem.avlna[whichcmnt];
         if(lhflag[46] == 1) ins >>telmnt[0].tem.avlnb[whichcmnt];
         if(lhflag[47] == 1) ins >>telmnt[0].tem.stona[whichcmnt];
         if(lhflag[48] == 1) 
         {
              ins >>telmnt[0].tem.stonb[whichcmnt];

              telmnt[0].tem.strnb[whichcmnt] = telmnt[0].tem.stonb[whichcmnt] * 0.9716;
              telmnt[0].tem.stonb[whichcmnt] = telmnt[0].tem.stonb[whichcmnt] * 0.02839;
         }
         if(lhflag[49] == 1) ins >>telmnt[0].tem.veg.unleaf12[whichcmnt];
         if(lhflag[50] == 1) ins >>telmnt[0].tem.veg.prvleafmx[whichcmnt];
         if(lhflag[51] == 1) ins >>telmnt[0].tem.veg.cmaxcut[whichcmnt];
         if(lhflag[52] == 1) ins >>telmnt[0].tem.veg.cmax1a[whichcmnt];
         if(lhflag[53] == 1) ins >>telmnt[0].tem.veg.cmax1b[whichcmnt];
         if(lhflag[54] == 1) ins >>telmnt[0].tem.veg.cmax2a[whichcmnt];
         if(lhflag[55] == 1) ins >>telmnt[0].tem.veg.cmax2b[whichcmnt];
         if(lhflag[56] == 1) ins >>telmnt[0].tem.veg.cfall[whichcmnt];
         if(lhflag[57] == 1) ins >>telmnt[0].tem.veg.kra[whichcmnt];
         if(lhflag[58] == 1) ins >>telmnt[0].tem.veg.krb[whichcmnt];
         if(lhflag[59] == 1) ins >>telmnt[0].tem.microbe.kda[whichcmnt];
         if(lhflag[60] == 1) ins >>telmnt[0].tem.microbe.kdb[whichcmnt];
         if(lhflag[61] == 1) ins >>telmnt[0].tem.microbe.lcclnc[whichcmnt];
         if(lhflag[62] == 1) ins >>telmnt[0].tem.microbe.propftos[whichcmnt];
         if(lhflag[63] == 1) ins >>telmnt[0].tem.veg.nmaxcut[whichcmnt];
         if(lhflag[64] == 1) ins >>telmnt[0].tem.veg.nmax1a[whichcmnt];
         if(lhflag[65] == 1) ins >>telmnt[0].tem.veg.nmax1b[whichcmnt];
         if(lhflag[66] == 1) ins >>telmnt[0].tem.veg.nmax2a[whichcmnt];
         if(lhflag[67] == 1) ins >>telmnt[0].tem.veg.nmax2b[whichcmnt];
         if(lhflag[68] == 1) ins >>telmnt[0].tem.veg.nfall[whichcmnt];
         if(lhflag[69] == 1) ins >>telmnt[0].tem.microbe.nupa[whichcmnt];
         if(lhflag[70] == 1) ins >>telmnt[0].tem.microbe.nupb[whichcmnt];
         if(lhflag[71] == 1) ins >>telmnt[0].tem.soil.nloss[whichcmnt];
         if(lhflag[72] == 1) ins >>telmnt[0].tem.microbe.nfixpar[whichcmnt];
         if(lhflag[73] == 1) ins >>telmnt[0].tem.veg.cneven[whichcmnt];
         if(lhflag[74] == 1) ins >>telmnt[0].tem.veg.cnmin[whichcmnt];
         if(lhflag[75] == 1) ins >>telmnt[0].tem.veg.c2na[whichcmnt];
         if(lhflag[76] == 1) ins >>telmnt[0].tem.veg.c2nb[whichcmnt];
         if(lhflag[77] == 1) ins >>telmnt[0].tem.veg.c2nmin[whichcmnt];
         if(lhflag[78] == 1) ins >>telmnt[0].tem.microbe.cnsoil[whichcmnt];	
         fclose(fiprsmp);

		//read parameter specification finished
//----------------------------------------------------------------------------------

//       cout <<"write flags"<<endl;
         cpflags(1);
         flog1.close();
    }
//-------------------------------------------------------------
//
	//communicating with other procs
    MPI::COMM_WORLD.Bcast(&tmflgs,sizeof(Temflags),MPI::BYTE,root);
    MPI::COMM_WORLD.Bcast(&(telmnt[0]),sizeof(TEMelmnt),MPI::BYTE,root);
    MPI::COMM_WORLD.Bcast(&elmnt,sizeof(Elmnt),MPI::BYTE,root);
    MPI::COMM_WORLD.Bcast(&fnminout,sizeof(Fnminout),MPI::BYTE,root);
    MPI::COMM_WORLD.Bcast(&temstac,sizeof(Temstac),MPI::BYTE,root);
    MPI::COMM_WORLD.Bcast(&basicf,sizeof(Soilthermal),MPI::BYTE,root);  
//------------------------------------------------  
//Initialize non-root nodes
    if(myrank0 != root)
    {		
         cpflags(0);
    }

//------------------------------------------------
//Open files for output

    myrank = myrank0 + f1stnd;
	
//Open log file
    sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, argv[2], myrank);
//	cout <<strtmp<<endl;
    flog1.open(strtmp);
    if (!flog1.is_open())
    {
         cout << endl << "Cannot open log file "<<strtmp<<endl;
         MPI::Finalize();
         exit(-1);
    }




//open input file
    sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmclds, myrank);
    fclds=fopen(strtmp,"r");
    if(fclds==NULL)
    {
         cout <<"Can not open file "<<strtmp<<endl;
         MPI::Finalize();
         exit(-1);
    }	
  
    if(telmnt[0].lonlatflag==0)
    {
         sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmlonlat, myrank);
         flonlat=fopen(strtmp,"r");
         if(flonlat==NULL)
         {
              cout <<"Can not open file "<<strtmp<<endl;
              MPI::Finalize();
              exit(-1);
         }		
    }
    sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmstxt, myrank);
    fstxt=fopen(strtmp,"r");
    if(fstxt==NULL)
    {
         cout <<"Can not open file "<<strtmp<<endl;
         MPI::Finalize();
         exit(-1);
    }	
    sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmtveg, myrank);
    ftveg=fopen(strtmp,"r");
    if(ftveg==NULL)
    {
         cout <<"Can not open file "<<strtmp<<endl;
         MPI::Finalize();
         exit(-1);
    }	
    sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmelev, myrank);

    felev=fopen(strtmp,"r");
    if(felev==NULL)
    {
         cout <<"Can not open file "<<strtmp<<endl;
         MPI::Finalize();
         exit(-1);
    }
    sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmtair, myrank);
    ftair=fopen(strtmp,"r");
    if(ftair==NULL)
    {
         cout <<"Can not open file "<<strtmp<<endl;
         MPI::Finalize();
         exit(-1);
    }	
    sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmprec, myrank);
	
    fprec=fopen(strtmp,"r");
    if(fprec==NULL)
    {
         cout <<"Can not open file "<<strtmp<<endl;
         MPI::Finalize();
         exit(-1);
    }	
    if(tmflgs.atmsflag==0)
    {
         sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmnirr, myrank);

         fnirr=fopen(strtmp,"r");
         if(fnirr==NULL)
	 {
              cout <<"Can not open file "<<strtmp<<endl;
              MPI::Finalize();
              exit(-1);

         }	
         sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmpar, myrank);

         fpar=fopen(strtmp,"r");
         if(fpar==NULL)
         {
              cout <<"Can not open file "<<strtmp<<endl;
              MPI::Finalize();
              exit(-1);

         }	
    }	
    if(tmflgs.equil==0)
    {
         if(telmnt[0].tem.ag.tlulcflag==1)
         {
              sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmlulc, myrank);
              flulc=fopen(strtmp,"r");
              if(flulc==NULL)
              {
                   cout <<"Can not open file "<<strtmp<<endl;
                   MPI::Finalize();
                   exit(-1);

              }
              if(telmnt[0].tem.ag.RAP0flag==0)
              {
                   sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmnpp, myrank);
                   fnpp=fopen(strtmp,"r");
                   if(fnpp==NULL)
                   {
                        cout <<"Can not open file "<<strtmp<<endl;
                        MPI::Finalize();
                        exit(-1);

                   }	
              }		
         }	
    }	
    if(tmflgs.kdinflg==1)
    {
         sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmkdin, myrank);              
         fkdin=fopen(strtmp,"r");
         if(fkdin==NULL)
         {
              cout <<"Can not open file "<<strtmp<<endl;
              MPI::Finalize();
              exit(-1);

         }	
    }
    if(tmflgs.kdoutflg==1)
    {
         sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmkdout, myrank);              
         fkdout.open(strtmp,ofstream::out);
         if(!fkdout.is_open())
         {
               cout <<"Can not open file "<<strtmp<<endl;
               MPI::Finalize();
               exit(-1);

         }	
    }



//open output file
    if(atmsoutfg == 1)
    {
         for (i = 0; i < natmspred; i++) 
         { 
              sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmclmpred[i], myrank);
              fclmpred[i].open(strtmp, ios_base::out); 
         }
    }
//	cout <<"rank "<<myrank<<" :"<<ntempred<<" "<<tmflgs.natmspred<<endl;
    for(i = 0; i < ntempred; i ++)
    {

         sprintf(strtmp, "%s/part-%d/%s-%d", mdir, myrank, fnmtempred[i], myrank);
//		flog1 <<strtmp<<endl;
         ftempred[i].open(strtmp, ios_base::out);
    }
//------------------------------------------------
// Do simulations 

    extrapolate();

	

//end

// Finished processing all elements - close open files

    if (fatalerr != 0)
    {
         if (elmnt.grdcnt != -99 && elmnt.count <= elmnt.grdcnt)
         {
              cout << "FATAL ERROR! Program Terminated" << endl;
         }
         flog1 << "FATAL ERROR! Program Terminated" << endl;
    }
    else
    {
         cout << "Extrapolation successfully completed - Congratulations!" << endl;
         flog1 << "Extrapolation successfully completed - Congratulations!" << endl;
    }

    if (atmsflag == 1)
    {
         fclose(fclds);
         if (telmnt[0].lonlatflag == 0) { fclose(flonlat); }
         if (atmsoutfg == 1)
         {
              for (i = 0; i < natmspred; i++) { fclmpred[i].close(); }
         }
    }

    if (temflag == 1)
    {
         fclose(ftveg);
         fclose(fstxt);
         fclose(felev);
         if (atmsflag == 0)
         {
              fclose(fnirr);
              fclose(fpar);
         }
         fclose(ftair);
         fclose(fprec);
         if(telmnt[0].tem.ag.tlulcflag == 1)
         {
              fclose(flulc);
              if(telmnt[0].tem.ag.RAP0flag == 0) { fclose(fnpp); }
         }
         for (i = 0; i < ntempred; i++) { ftempred[i].close(); }
         if (stateflag == 1)
         {
              for (i = 0; i < MAXSTATE; i++) { fclose(fstate[i]); }
         }
    }

    flog1.close();

    MPI::Finalize();
    return 1;
};

/* *************************************************************
*********************END OF MAIN PROGRAM************************
************************************************************* */


/* **************************************************************
************************************************************** */

int askpred(char pvarname[MAXPRED][9], int spred)
{
  int i;
  int j;
  int k;
  int t;
  int cnt;
  int numpred;
  int length;

  int posspred = MAXPRED;
  if (atmsoutfg == 1 && temflag == 0) { posspred = NUMATMS+1; }

  cout << endl << endl << "           POSSIBLE OUTPUT VARIABLES:"
		<< endl << endl;
  for (i = 0; i < (posspred/5); i++)
  {
	for (t = 0; t < 5; t++)
	{
	  cout << pvarname[(5*i)+t] << " ";
	}
	cout << endl;
  }

  flog1 << endl << endl << "           POSSIBLE OUTPUT VARIABLES:"
		<< endl << endl;
  for (i = 0; i < (posspred/5); i++)
  {
     for (t = 0; t < 5; t++)
    {
      flog1 << pvarname[(5*i)+t] << ' ';
    }
    flog1 << endl;
  }

  cout << endl << endl << "How many variables are to be mapped (max "
	   << posspred << ") in output files?  ";
//  cin >> numpred;
  fpara >> numpred;
  cout << numpred << endl;

  flog1 << endl << endl << "How many variables are to be mapped (max "
		<< posspred << ") in output files?" << numpred << endl << endl;

  cout << "Please enter output variable: " << endl;
  flog1 << "Please enter output variable: " << endl;

  for (i = 0; i < numpred; i++)
  {
    cnt = i + 1;
    k = i+spred;
    cout << cnt << " ";
//	cin >> predmap[k];
	fpara >> predmap[k];
    length = strlen(predmap[k]);
    for (j = 0; j < length; j++) { predmap[k][j] = toupper(predmap[k][j]); }
	flog1 << cnt << " " << predmap[k] << endl;
  }

  return numpred;

};

/* **************************************************************
************************************************************** */
void read_dr(void)
{

	int i;
	int itype;
	end1 = 1;
	fatalerr = 0;
// Skip to the desired record in the GIS data sets

	if (elmnt.strtflag == 0)
	{
		for (i = 0; i < elmnt.numskip; i++)
		{
			if (atmsflag == 1)
			{
				end1 = telmnt[0].atmsgisin(fclds,cldflag,flonlat,telmnt[0].lonlatflag,
                                   numspin,spintime,RTIME);
				if (end1 == -1) { exit(0); }
			}

			if (temflag == 1)
			{
				end1 = telmnt[0].veggisin(flog1,fatalerr,atmsflag,temflag,ftveg);
				if (end1 == -1) { exit(0); }

				for (itype = 0; itype < telmnt[icount].maxtype; itype++)
				{
					if (temflag == 1 && fatalerr == 0)
					{
						end1 = telmnt[0].temgisin(flog1,fatalerr,atmsflag,
                                      itype,telmnt[icount].maxtype,
                                      kdinflg,stateflag,
	     	                      fstxt,felev,
                                      fnirr,fpar,
                                      ftair,fprec,
                                      flulc,fnpp,
                                      fkdin,fstate,
                                      numspin,spintime,RTIME);
						if (end1 == -1) { exit(0); }
					}
				}
			}
		}
	}

	icount = 0;
	while (icount < MAXGRID && end1 != -1)
	{


// Get spatially-explicit cloudiness data   *********************

		if (atmsflag == 1)
		{


			end1 = telmnt[icount].atmsgisin(fclds,cldflag,flonlat,
                                       telmnt[icount].lonlatflag,
                                       numspin,spintime,RTIME);
		}

// Get spatially-explicit vegetation data   *********************

		if (end1 != -1 && temflag == 1)
		{
			end1 = telmnt[icount].veggisin(flog1,fatalerr,atmsflag,temflag,ftveg);
		}

/* *************************************************************
	BEGIN VEGETATION MOSAIC LOOP
************************************************************* */

		if (end1 != -1 && temflag == 1)
		{
			for (itype = 0; itype < telmnt[icount].maxtype; itype++)
			{
          // Get spatially-explicit input data for TEM

				if (end1 != -1)
				{
					end1 = telmnt[icount].temgisin(flog1,fatalerr,
                                           atmsflag,itype,
                                           telmnt[icount].maxtype,
                                           kdinflg,stateflag,
											fstxt,felev,
                                           fnirr,fpar,
                                           ftair,fprec,
                                           flulc,fnpp,
                                           fkdin,fstate,
                                           numspin,spintime,RTIME);
				}
			}


		}

		if (end1 != -1) { ++icount; }

	}
	if(icount == 0)
	{
		cout <<"Read input data failed "<<endl;
		flog1 <<"Read input data failed "<<endl;
		exit(-1);
	}
};
//************************************************************** */

/* **************************************************************
************************************************************** */

void extrapolate(void)
{

    int i;
    int j;
    int k;
    int dm;
    int maxrec;
    int dyr;
    int itype;
//  double totamnt[NUMEQ+2];



//* *************************************************************



    end1 = 1;
    fatalerr = 0;

//* *************************************************************

// Skip to the desired record in the GIS data sets

    if (elmnt.strtflag == 0)
    {
         for (i = 0; i < elmnt.numskip; i++)
         {
              if (atmsflag == 1)
              {
                   end1 = telmnt[0].atmsgisin(fclds,cldflag,flonlat,telmnt[0].lonlatflag,
                                   numspin,spintime,RTIME);
                   if (end1 == -1) { exit(0); }
              }

              if (temflag == 1)
              {
                   end1 = telmnt[0].veggisin(flog1,fatalerr,atmsflag,temflag,ftveg);
                   if (end1 == -1) { exit(0); }

                   for (itype = 0; itype < telmnt[icount].maxtype; itype++)
                   {
	                if (temflag == 1 && fatalerr == 0)
                        {
	                     end1 = telmnt[0].temgisin(flog1,fatalerr,atmsflag,
                                      itype,telmnt[icount].maxtype,
                                      kdinflg,stateflag,
	     	                      fstxt,felev,
                                      fnirr,fpar,
                                      ftair,fprec,
                                      flulc,fnpp,
                                      fkdin,fstate,
                                      numspin,spintime,RTIME);
                             if (end1 == -1) { exit(0); }
                        }
                   }
              }
         }
    }

//************************************************************** */

    cout << endl;
    flog1 << endl << endl;

    telmnt[0].col = MISSING;
    telmnt[0].row = MISSING;
    telmnt[0].tem.totyr = -99;

    elmnt.show(flog1, telmnt[0].col, telmnt[0].row,
			 telmnt[0].tem.totyr, telmnt[0].tem.inittol);


    while (end1 != -1 && fatalerr == 0)
    {

         icount = 0;
	 while (icount < MAXGRID && end1 != -1)
	 {

//      cout <<endl<<"do loop"<<endl;
// Get spatially-explicit cloudiness data   *********************
    //  printf("read cld\n");
              if (atmsflag == 1)
	      {


		   end1 = telmnt[icount].atmsgisin(fclds,cldflag,flonlat,
                                        telmnt[icount].lonlatflag,
                                        numspin,spintime,RTIME);
	      }
  //   printf("read veg\n");
// Get spatially-explicit vegetation data   *********************

              if (end1 != -1 && temflag == 1)
              {
                   end1 = telmnt[icount].veggisin(flog1,fatalerr,atmsflag,temflag,ftveg);
              }
//	  cout <<endl<<"pause 2"<<endl<<"icount ="<< icount<<endl;
/* *************************************************************
		BEGIN VEGETATION MOSAIC LOOP
************************************************************* */

//      printf("read other data\n");
              if (end1 != -1 && temflag == 1)
              {
                   for (itype = 0; itype < telmnt[icount].maxtype; itype++)
		   {
          // Get spatially-explicit input data for TEM

	                if (end1 != -1)
                        {
                             end1 = telmnt[icount].temgisin(flog1,fatalerr,
                                           atmsflag,itype,
                                           telmnt[icount].maxtype,
                                           kdinflg,stateflag,
					   fstxt,felev,
                                           fnirr,fpar,
                                           ftair,fprec,
                                           flulc,fnpp,
                                           fkdin,fstate,
                                           numspin,spintime,RTIME);
                        }
                   }
//	  cout <<endl<<"pause 3"<<endl <<"icount ="<< icount<<endl;

	      }

	      if (end1 != -1) { ++icount; }
//	  cout <<endl<<"pause 4"<<endl<<"icount ="<< icount<<endl;

         }
         maxrec = icount;

  //   cout <<"run "<<maxrec<<endl;
         for (i = 0; i < maxrec; i++)
         {

      // Run TEM for "i" grid cells
//	  cout <<endl<<"begin run tem"<<endl;
              telmnt[i].runtem(flog1,predmap,
                       cldflag,atmsflag,atmsoutfg, natmspred, fclmpred,
                       temflag,kdinflg,ntempred,ftempred,
					   stateflag,equil,totsptime, RTIME,ispinout);
//	 cout <<endl<<"runtem"<<endl<<"pause 5"<<endl;

         }

/* *************************************************************
	  Write georeferenced data to output files
************************************************************* */

         for (i = 0; i < maxrec; i++)
         {
              for (dyr = 0; dyr < telmnt[i].outyr; dyr++)
              {
	           for (itype = 0; itype < telmnt[i].maxtype; itype++)
                   {
                        if (temflag == 1)
                        {
            // Write out spatially explicit kd parameters

	                     if (dyr == 0 && kdoutflg == 1)
                             {
	                          kdparam.outdel(fkdout, telmnt[i].col, telmnt[i].row,
                                    telmnt[i].tem.veg.temveg,
                                    telmnt[i].tem.veg.subtype[telmnt[i].mez][itype],
                                    telmnt[i].tem.microbe.kdsave[itype],
                                    telmnt[i].contnent);
	                     }
	                }
	           }
              }

              if (telmnt[i].tem.intflag > 0)
              {
                   if (elmnt.count < elmnt.grdcnt)
                   {
		        cout << "Integration terminated before attaining"
		             << " tolerance level" << endl;
	           }
		   flog1 << "Integration terminated before attaining"
		        << "s tolerance level" << endl;
              } 
              elmnt.show(flog1, telmnt[i].col, telmnt[i].row,
                 telmnt[i].tem.totyr, telmnt[i].tem.tol);
              if (elmnt.stopflag == 1) { end1 = -1; }
	 }

    }

		

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void startclm(void)
{

	int i;
	char ifilename[25];


	cldflag = 1;
	cout << "Do you have spatially explicit solar radiation data"
  	   << " or cloudiness data?:" << endl;
	cout << "Enter 0 for solar radiation data (W/ sq. m):" << endl;
	cout << "Enter 1 for cloudiness data (percent cloudiness): ";

	fpara >> cldflag;

	flog1 << "Do you have spatially explicit solar radiation"
        << " data or cloudiness data?:" << endl;
	flog1 << "Enter 0 for solar radiation data (W/ sq. m):" << endl;
	flog1 << "Enter 1 for cloudiness data (percent cloudiness): " << endl;
	flog1 << "cldflag = " << cldflag << endl;

	telmnt[0].clm.tcldsflag = 0;
	if (cldflag == 1)
	{
		cout << "Do you have transient cloudiness data?:" << endl;
		cout << "Enter 0 for no:" << endl;
		cout << "Enter 1 for yes: ";

		fpara >> telmnt[0].clm.tcldsflag;

		flog1 << endl << "Do you have transient cloudiness data?:" << endl;
		flog1 << "Enter 0 for no:" << endl;
		flog1 << "Enter 1 for yes: " << endl;
		flog1 << "telmnt[0].clm.tcldsflag = " << telmnt[0].clm.tcldsflag << endl << endl;

		cout << "Please enter the name of the file containing"
			<< " the mean monthly cloudiness data: " << endl;
		cout << "               (e.g., CLDS.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing the"
	             << " mean monthly cloudiness data: " << endl;
		flog1 << "               (e.g., CLDS.GIS) " << endl;
		flog1 << ifilename << endl << endl;
	}
	else
	{
		cout << "Do you have transient mean monthly solar radiation data?:" << endl;
		cout << "Enter 0 for no:" << endl;
		cout << "Enter 1 for yes: ";

		fpara >> telmnt[0].clm.tcldsflag;

		flog1 << endl << "Do you have transient mean monthly"
	      << " solar radiation data?:" << endl;
		flog1 << "Enter 0 for no:" << endl;
		flog1 << "Enter 1 for yes: " << endl;
		flog1 << "telmnt[0].clm.tcldsflag = " << telmnt[0].clm.tcldsflag
			<< endl << endl;

		cout << "Please enter the name of the file containing"
			<< " the mean monthly solar radiation data: " << endl;
		cout << "               (e.g., NIRR.GIS) " << endl;

		fpara >> ifilename;
		flog1 << "Please enter the name of the file containing"
		  << " the mean monthly solar radiation data: " << endl;
		flog1 << "               (e.g., NIRR.GIS) " << endl;
		flog1 << ifilename << endl << endl;
	}

	sprintf(fnmclds, "%s", ifilename);


	cout << "How do you locate your grid cells?" << endl;
	cout << "Enter 0 for column/row:" << endl;
	cout << "Enter 1 for longitude/latitude: ";

	fpara >> telmnt[0].lonlatflag;

	flog1 << "How do you locate your grid cells?" << endl;
	flog1 << "Enter 0 for column/row:" << endl;
	flog1 << "Enter 1 for longitude/latitude: " << telmnt[0].lonlatflag
  		<< endl << endl;

	if (telmnt[0].lonlatflag == 0)
	{
		cout << "Please enter the name of the file containing the latitude data: "
			<< endl;
		cout << "               (e.g., LAT.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the latitude data: " << endl;
		flog1 << "               (e.g., LAT.GIS) " << endl;
		flog1 << ifilename << endl << endl;

		sprintf(fnmlonlat, "%s", ifilename);
	}

	atmsoutfg = 0;
	cout << endl << "Do you wish output data from the irradiance model? " << endl;
	cout << "  Enter 0 for no" << endl;
	cout << "  Enter 1 for yes: ";

	fpara >> atmsoutfg;

	flog1 << endl << "Do you wish output data from the irradiance model? " << endl;
	flog1 << "  Enter 0 for no" << endl;
	flog1 << "  Enter 1 for yes: " << endl;
	flog1 << "atmsoutfg = " << atmsoutfg << endl << endl;

	if (atmsoutfg == 1)
	{
		totpred = natmspred = askpred(telmnt[0].clm.predstr, totpred);

		for (i = 0; i < natmspred; i++)
		{
			cout << endl << "Enter the name of the OUTPUT file to contain ";
			cout << predmap[i] << ":  ";
	
			fpara >> fnmclmpred[i];
	  
			flog1 << endl << "Enter the name of the OUTPUT file to contain ";
			flog1 << predmap[i] << ":  " << fnmclmpred[i] << endl;
			
		}
	}

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void starttem(void)
{

//read necessary parameters and name of the input files	
	int i;
	int numcmnt;
	char ifilename[40];
	char tempredfile[25];

	telmnt[0].tem.initrun(flog1,equil);
	telmnt[0].tem.ask(flog1);


// Get soil texture dependent parameters

	telmnt[0].tem.soil.getecd(flog1);


// Get vegetation type dependent parameters

	telmnt[0].tem.soil.getrootz(flog1);
	telmnt[0].tem.veg.getecd(flog1);//get kc ki kn1
	telmnt[0].tem.veg.getleafecd(flog1);
	telmnt[0].tem.microbe.getvegecd(flog1);//kn2


//Get calibration site specific parameters

	cout << "Please enter the number of community types with calibration data:";

	fpara >> numcmnt;

	flog1 << endl << endl << "Please enter the number of"
  	    << " community types with calibration data:";
	flog1 << numcmnt << endl;

	telmnt[0].tem.getsitecd(numcmnt, flog1);


// Get vegetation mosaic information

	telmnt[0].tem.veg.getvtype(flog1);


	cout << "Please enter the name of the file containing the soil texture data:";
	cout << endl;
	cout << "               (e.g., TEXTURE.GIS) " << endl;

	fpara >> ifilename;
  
	flog1 << "Please enter the name of the file containing"
  		<< " the soil texture data:";
	flog1 << endl;
	flog1 << "               (e.g., TEXTURE.GIS) " << endl;
	flog1 << ifilename << endl << endl;
	sprintf(fnmstxt, "%s", ifilename);

	cout << "Please enter the name of the file containing the"
  	   << " vegetation data: " << endl;
	cout << "               (e.g., TEMVEG.GIS) " << endl;

	fpara >> ifilename;

	flog1 << "Please enter the name of the file containing"
  		<< " the vegetation data: " << endl;
	flog1 << "               (e.g., TEMVEG.GIS) " << endl;
	flog1 << ifilename << endl << endl;
	sprintf(fnmtveg, "%s", ifilename);


	cout << "Please enter the name of the file containing"
  	   << " the elevation data: " << endl;
	cout << "               (e.g., ELEV.GIS) " << endl;

	fpara >> ifilename;

	flog1 << "Please enter the name of the file"
  		<< " containing the elevation data: " << endl;
	flog1 << "               (e.g., ELEV.GIS) " << endl;
	flog1 << ifilename << endl << endl;
	sprintf(fnmelev, "%s", ifilename);


	telmnt[0].tem.atms.ttairflag = 0;
	if (equil == 0)
	{
		cout << "Do you have transient air temperature data?:" << endl;
		cout << "Enter 0 for No:" << endl;
		cout << "Enter 1 for Yes: ";

		fpara >> telmnt[0].tem.atms.ttairflag;

		flog1 << endl << "Do you have transient air temperature data?:" << endl;
		flog1 << "Enter 0 for No:" << endl;
		flog1 << "Enter 1 for Yes: " << endl;
		flog1 << "telmnt[0].tem.atms.ttairflag = "
			<< telmnt[0].tem.atms.ttairflag << endl << endl;
	}

	cout << "Please enter the name of the file containing the"
  	   << " mean monthly air temperature data: " << endl;
	cout << "               (e.g., TAIR.GIS) " << endl;

	fpara >> ifilename;
	flog1 << "Please enter the name of the file containing"
  		<< " the mean monthly air temperature data: " << endl;
	flog1 << "               (e.g., TAIR.GIS) " << endl;
	flog1 << ifilename << endl << endl;
	sprintf(fnmtair, "%s", ifilename);


	telmnt[0].tem.atms.tprecflag = 0;
	if (equil == 0)
	{
		cout << "Do you have transient precipitation data?:" << endl;
		cout << "Enter 0 for No:" << endl;
		cout << "Enter 1 for Yes: ";

		fpara >> telmnt[0].tem.atms.tprecflag;

		flog1 << endl << "Do you have transient precipitation data?:" << endl;
		flog1 << "Enter 0 for No:" << endl;
		flog1 << "Enter 1 for Yes: " << endl;
		flog1 << "temnt[0].tem.atms.tprecflag = "
			<< telmnt[0].tem.atms.tprecflag << endl << endl;
	}

	cout << "Please enter the name of the file containing"
  	   << " the monthly precipitation data: " << endl;
	cout << "               (e.g., PREC.GIS) " << endl;

	fpara >> ifilename;
  
	flog1 << "Please enter the name of the file containing"
  	   << " the monthly precipitation data: " << endl;
	flog1 << "               (e.g., PREC.GIS) " << endl;
	flog1 << ifilename << endl << endl;
	sprintf(fnmprec, "%s", ifilename);


	if (atmsflag == 0)
	{
		cout << "Please enter the name of the file containing"
			<< " the net irradiance data: " << endl;
		cout << "               (e.g., NIRR.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the net irradiance data: " << endl;
		flog1 << "               (e.g., NIRR.GIS) " << endl;
		flog1 <<  ifilename << endl << endl;
		sprintf(fnmnirr,"%s",ifilename);


		cout << "Please enter the name of the file containing"
			<< " the mean monthly photosynthetically active radiation data: " << endl;
		cout << "               (e.g., PAR.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the mean monthly photosynthetically active radiation data: "
		  << endl;
		flog1 << "               (e.g., PAR.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmpar, "%s", ifilename);

	}

	cout << endl << endl << "Enter the initial concentration"
  	   << " of carbon dioxide in ppmv: ";

	fpara >> telmnt[0].tem.atms.initco2;

	flog1 << endl << endl << "Enter the initial concentration of"
		<< " carbon dioxide in ppmv: " << telmnt[0].tem.atms.initco2
		<< endl << endl;

	cout << endl << endl << "Enter the final equilibrium concentration"
  	   << " of carbon dioxide in ppmv: ";

	fpara >> telmnt[0].tem.atms.co2level;

	flog1 << endl << endl << "Enter the final equilibrium concentration"
		<< " of carbon dioxide in ppmv: " << telmnt[0].tem.atms.co2level
		<< endl << endl;

	telmnt[0].tem.atms.tco2flag = 0;
	if (equil == 0)
	{
		cout << "Do you have transient CO2 data?:" << endl;
		cout << "Enter 0 for No:" << endl;
		cout << "Enter 1 for Yes: ";

		fpara >> telmnt[0].tem.atms.tco2flag;

		flog1 << "Do you have transient CO2 data?:" << endl;
		flog1 << "Enter 0 for No:" << endl;
		flog1 << "Enter 1 for Yes: " << endl;
		flog1 << "telmnt[0].tem.atms.tco2flag = "
			<< telmnt[0].tem.atms.tco2flag << endl << endl;

		if (telmnt[0].tem.atms.tco2flag == 1)
		{
			telmnt[0].tem.atms.loadmpitCO2(flog1, totsptime, RTIME);
		}
	}

	cout << endl << endl << "Enter the factor for changing C:N"
  	   << " per ppmv of enhanced CO2:" << endl;
	cout << "                     (Enter 0.0 for no change): " << endl;

	fpara >> telmnt[0].tem.veg.dc2n;

	flog1 << endl << "Enter the factor for changing C:N per"
  		<< " ppmv of enhanced CO2:" << endl;
	flog1 << "                     (Enter 0.0 for no change): " << endl;
	flog1 << "telmnt[0].tem.veg.dc2n = " << telmnt[0].tem.veg.dc2n
  		<< endl << endl;

	telmnt[0].tem.ag.tlulcflag = 0;
	if (equil == 0)
	{
		cout << "Do you have transient land use data?:" << endl;
		cout << "Enter 0 for No:" << endl;
		cout << "Enter 1 for Yes: ";

		fpara >> telmnt[0].tem.ag.tlulcflag;

		flog1 << "Do you have transient land use data?:" << endl;
		flog1 << "Enter 0 for No:" << endl;
		flog1 << "Enter 1 for Yes: ";
		flog1 << "telmnt[0].tem.ag.tlulcflag = "
			<< telmnt[0].tem.ag.tlulcflag << endl << endl;;

		if (telmnt[0].tem.ag.tlulcflag == 1)
		{
			cout << "Please enter the name of the file containing"
		  	   << " the land use data: " << endl;
			cout << "               (e.g., LULC.GIS) " << endl;

			fpara >> ifilename;
	  
			flog1 << "Please enter the name of the file containing the"
	  			<< " land use data: " << endl;
			flog1 << "               (e.g., LULC.GIS) " << endl;
			flog1 << ifilename << endl << endl;
			sprintf(fnmlulc, "%s", ifilename);


			cout << "Will Relative Agricultural Production (RAP)"
	  			<< " always be 0.0?" << endl;
			cout << "Enter 0 for No:" << endl;
			cout << "Enter 1 for Yes: ";
	
			fpara >> telmnt[0].tem.ag.RAP0flag;

			flog1 << "Will Relative Agricultural Production (RAP)"
	  			<< " always be 0.0?" << endl;
			flog1 << "Enter 0 for No:" << endl;
			flog1 << "Enter 1 for Yes: " << endl;
			flog1 << "telmnt[0].tem.ag.RAP0flag = "
	  			<< telmnt[0].tem.ag.RAP0flag << endl;

			if(telmnt[0].tem.ag.RAP0flag == 0)
			{
				cout << "Please enter the name of the file containing"
					<< " the monthly potential NPP data: " << endl;
 				cout << "               (e.g., NPP.GIS) " << endl;
	
				fpara >> ifilename;

				flog1 << "Please enter the name of the file containing"
					<< " the monthly potential NPP data: " << endl;
 				flog1 << "               (e.g., NPP.GIS) " << endl;
				flog1 << ifilename << endl << endl;
				sprintf(fnmnpp, "%s", ifilename);

			}

			telmnt[0].tem.ag.getecd(flog1);
		}
	}

	kdinflg = 0;
	cout << "Do you want to use spatially explicit values of Kd?" << endl;
	cout << "Enter 0 for No:" << endl;
	cout << "Enter 1 for Yes: ";

	fpara >> kdinflg;

	flog1 << "Do you want to use spatially explicit values of Kd?" << endl;
	flog1 << "Enter 0 for No:" << endl;
	flog1 << "Enter 1 for Yes: " << endl;
	flog1 << "kdinflg = " << kdinflg << endl << endl;

	if (kdinflg == 1)
	{
		cout << "Please enter the name of the file containing"
			<< " the Kd values: " << endl;
		cout << "               (e.g., KDIN.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the Kd values: " << endl;
		flog1 << "               (e.g., KDIN.GIS) " << endl;
		flog1 << ifilename << endl << endl;

		sprintf(fnmkdin, "%s", ifilename);

	}

	kdoutflg = 0;
	cout << "Do you want to generate spatially explicit values of Kd?" << endl;
	cout << "Enter 0 for No:" << endl;
	cout << "Enter 1 for Yes: ";

	fpara >> kdoutflg;

	flog1 << endl << "Do you want to generate spatially"
  	    << " explicit values of Kd?" << endl;
	flog1 << "Enter 0 for No:" << endl;
	flog1 << "Enter 1 for Yes: " << endl;
	flog1 << "kdoutflg = " << kdoutflg << endl << endl;

	if (kdoutflg == 1)
	{
		cout << "Please enter the name of the file to contain the"
		  << " Kd values: " << endl;
		cout << "               (e.g., KDOUT.GIS) " << endl;

		fpara >> fnmkdout;

		flog1 << "Please enter the name of the file to contain the"
		  << " Kd values: " << endl;
		flog1 << "               (e.g., KDOUT.GIS) " << endl;
		flog1 << fnmkdout << endl << endl;
		
	}

	ntempred = askpred(telmnt[0].tem.predstr, totpred);
	totpred += ntempred;

	for (i = natmspred; i < totpred; i++)
	{
		cout << endl << "Enter the name of the OUTPUT file to contain "
			<< predmap[i] << ":  ";
	
		fpara >> fnmtempred[i - natmspred];
	
		flog1 << "Enter the name of the OUTPUT file to contain "
			<< predmap[i] << ":  " << fnmtempred[i - natmspred] << endl;
		
	}

	stateflag = 0;
	cout << endl << "Do you want to use spatially explicit"
  	   << " data for initial conditions? " << endl;
	cout << "  Enter 1 for YES:" << endl;
	cout << "  Enter 0 for NO:" << endl;
	cout << "  NOTE:  If YES, you will need spatially"
  	   << " explicit data for VEGC," << endl;
	cout << "         STRUCTN, SOLC, SOLN, AVLN, NSTORE,"
  	   << " AVAILH2O, RGRNDH2O," << endl;
	cout << "         SNOWPACK, SGRNDH2O, UNNORMLEAF, PET,"
       << " EET,(i.e. 13 data files)" << endl;

	fpara >> stateflag;

	flog1 << endl << "Do you want to use spatially"
  	   << " explicit data for intial conditions? " << endl;
	flog1 << "  Enter 1 for YES:" << endl;
	flog1 << "  Enter 0 for NO:" << endl;
	flog1 << "  NOTE:  If YES, you will need spatially"
		<< " explicit data for VEGC," << endl;
	flog1 << "         STRUCTN, SOLC, SOLN, AVLN, NSTORE,"
		<< " AVAILH2O, RGRNDH2O," << endl;
	flog1 << "         SNOWPACK, SGRNDH2O, UNNORMLEAF, PET,"
		<< " EET,(i.e. 13 data files)" << endl;
	flog1 << "stateflag = " << stateflag << endl << endl;

	if (stateflag == 1)
	{
		cout << "Please enter the name of the file containing"
			<< " the vegetation carbon data: " << endl;
		cout << "               (e.g., VEGC.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the vegetation carbon data: " << endl;
		flog1 << "               (e.g., VEGC.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[0], "%s", ifilename);

		cout << "Please enter the name of the file containing"
			<< " the vegetation structural nitrogen data: " << endl;
		cout << "               (e.g., STRN.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the vegetation structural nitrogen data: " << endl;
		flog1 << "               (e.g., STRN.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[1], "%s", ifilename);


		cout << "Please enter the name of the file containing"
			<< " the soil organic carbon data: " << endl;
		cout << "               (e.g., SOLC.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the soil organic carbon data: " << endl;
		flog1 << "               (e.g., SOLC.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[2], "%s", ifilename);


		cout << "Please enter the name of the file containing"
			<< " the soil organic nitrogen data: " << endl;
		cout << "               (e.g., SOLN.GIS) " << endl;

		fpara >> ifilename;
	
		flog1 << "Please enter the name of the file containing"
		  << " the soil organic nitrogen data: " << endl;
		flog1 << "               (e.g., SOLN.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[3], "%s", ifilename);


		cout << "Please enter the name of the file containing"
			<< " the soil available nitrogen data: " << endl;
		cout << "               (e.g., AVLN.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
			<< " the soil availablenitrogen data: " << endl;
		flog1 << "               (e.g., AVLN.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[4], "%s", ifilename);


		cout << "Please enter the name of the file containing"
			<< " the vegetation labile nitrogen data: " << endl;
		cout << "               (e.g., STON.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the vegetation labile nitrogen data: " << endl;
		flog1 << "               (e.g., STON.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[5], "%s", ifilename);


		cout << "Please enter the name of the file containing"
			<< " the available soil moisture data: " << endl;
		cout << "               (e.g., AVLW.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the available soil moisture data: " << endl;
		flog1 << "               (e.g., AVLW.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[6], "%s", ifilename);


		cout << "Please enter the name of the file containing"
			<< " the rain-derived groundwater data: " << endl;
		cout << "               (e.g., RGRW.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the rain-derived groundwater data: " << endl;
	        flog1 << "               (e.g., RGRW.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[7], "%s", ifilename);


		cout << "Please enter the name of the file containing"
			<< " the snowpack data: " << endl;
		cout << "               (e.g., SNWP.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the snowpack data: " << endl;
		flog1 << "               (e.g., SNWP.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[8],"%s", ifilename);


		cout << "Please enter the name of the file containing"
			<< " the snowpack-derived groundwater data: " << endl;
		cout << "               (e.g., SGRW.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the snowpack-derived groundwater data: " << endl;
		flog1 << "               (e.g., SGRW.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[9], "%s", ifilename);


		cout << "Please enter the name of the file containing"
			<< " the unnormalized leaf phenology data: " << endl;
		cout << "               (e.g., UNLF.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the unnormalized leaf phenology data: " << endl;
		flog1 << "               (e.g., UNLF.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[10], "%s",ifilename);


		cout << "Please enter the name of the file containing"
			<< " the potential evapotranspiration data: " << endl;
		cout << "               (e.g., PET.GIS) " << endl;

		fpara >> ifilename;

		flog1 << "Please enter the name of the file containing"
		  << " the potential evapotranspiration data: " << endl;
		flog1 << "               (e.g., PET.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[11], "%s", ifilename);

		cout << "Please enter the name of the file containing"
			<< " the estimated evapotranspiration data: " << endl;
		cout << "               (e.g., EET.GIS) " << endl;

		fpara >> ifilename;
	
		flog1 << "Please enter the name of the file containing"
			<< " the estimated evapotranspiration data: " << endl;
		flog1 << "               (e.g., EET.GIS) " << endl;
		flog1 << ifilename << endl << endl;
		sprintf(fnmstate[12], "%s", ifilename, "r");

	}

};



void cpflags(const int opt)
{
    int i;
//	cout <<"opt = "<<opt<<endl;
    if(opt == 1)
    {
//		cout <<"copy from global variable to struct"<<endl;
         tmflgs.stmflg = stmflg;   
         tmflgs.equil = equil;
         tmflgs.RTIME = RTIME;

         tmflgs.spinflag = spinflag;
         tmflgs.numspin = numspin;
         tmflgs.spintime = spintime;
         tmflgs.ispinout = ispinout; 
         tmflgs.totsptime = totsptime;
         tmflgs.transtime = transtime;



         tmflgs.atmsflag = atmsflag;
         tmflgs.atmsoutfg = atmsoutfg;
         tmflgs.temflag = temflag;
         tmflgs.stateflag = stateflag;

         tmflgs.cldflag = cldflag;
         tmflgs.natmspred = natmspred;
         tmflgs.ntempred = ntempred;
//		cout <<"ntempred = "<<ntempred<<" "<<tmflgs.ntempred<<endl;
         tmflgs.totpred = totpred;
         for(i = 0; i < MAXPRED; i ++)
              strcpy(tmflgs.predmap[i], predmap[i]);
			
         strcpy(tmflgs.fnmlonlat,fnmlonlat);
         strcpy(tmflgs.fnmstxt,fnmstxt);
         strcpy(tmflgs.fnmelev,fnmelev);
         strcpy(tmflgs.fnmtveg,fnmtveg);
         strcpy(tmflgs.fnmclds,fnmclds);
         strcpy(tmflgs.fnmtair,fnmtair);
         strcpy(tmflgs.fnmprec,fnmprec);
         strcpy(tmflgs.fnmnirr,fnmnirr);
         strcpy(tmflgs.fnmpar,fnmpar);
         strcpy(tmflgs.fnmlulc,fnmlulc);
         strcpy(tmflgs.fnmnpp,fnmnpp);
         if (stateflag == 1)
         {
              for (i = 0; i < MAXSTATE; i++)strcpy(tmflgs.fnmstate[i],fnmstate[i]);
         }

         for(i = 0; i < natmspred; i ++)
              strcpy(fnminout.fnmclmpred[i], fnmclmpred[i]);
         for(i = 0; i < ntempred; i ++)
              strcpy(fnminout.fnmtempred[i], fnmtempred[i]);
         strcpy(tmflgs.mdir, mdir);
         strcpy(fnminout.fnmkdout, fnmkdout);
         for(i = 0; i < MAXNPAR; i ++)
              tmflgs.lhflag[i] = lhflag[i];
         tmflgs.whichcmnt = whichcmnt;
         tmflgs.kdinflg = kdinflg;
         tmflgs.kdoutflg = kdoutflg;

         tmflgs.f1stnd = f1stnd;

         temstac.avlnflag = telmnt[0].tem.avlnflag;
         temstac.nfeed = telmnt[0].tem.nfeed;
         temstac.initbase = telmnt[0].tem.initbase;
         temstac.baseline = telmnt[0].tem.baseline;
         temstac.moistlim = telmnt[0].tem.moistlim;
         temstac.strteq = telmnt[0].tem.strteq;
         temstac.maxyears = telmnt[0].tem.maxyears;
         temstac.runsize = telmnt[0].tem.runsize;
         temstac.maxnrun = telmnt[0].tem.maxnrun;
         temstac.rheqflag = telmnt[0].tem.rheqflag;
         temstac.wtol = telmnt[0].tem.wtol;
         temstac.ctol = telmnt[0].tem.ctol;
         temstac.ntol = telmnt[0].tem.ntol;
         temstac.startyr = telmnt[0].tem.startyr;
         temstac.endyr = telmnt[0].tem.endyr;
         temstac.diffyr = telmnt[0].tem.diffyr;

         temstac.inittol = telmnt[0].tem.inittol;
         temstac.maxit = telmnt[0].tem.maxit;
         temstac.maxitmon = telmnt[0].tem.maxitmon;		
    }
    else
    {
         //copy from struct to global variables

         stmflg = tmflgs.stmflg;   
         equil = tmflgs.equil;		
         RTIME = tmflgs.RTIME;
//       cout <<"equil = "<<tmflgs.equil<<endl;
//       cout <<"stmflg = "<<tmflgs.stmflg<<endl;
//       cout <<"RTIME = "<<tmflgs.RTIME<<endl;

         spinflag = tmflgs.spinflag;
         numspin = tmflgs.numspin;
         spintime = tmflgs.spintime;
         ispinout = tmflgs.ispinout; 
         totsptime = tmflgs.totsptime;
         transtime = tmflgs.transtime;

         atmsflag = tmflgs.atmsflag;
         atmsoutfg = tmflgs.atmsoutfg;
         temflag = tmflgs.temflag;
         stateflag = tmflgs.stateflag;

         cldflag = tmflgs.cldflag;
         natmspred = tmflgs.natmspred;
         ntempred = tmflgs.ntempred;
         totpred = tmflgs.totpred;
         for(i = 0; i < MAXPRED; i ++)
              strcpy(predmap[i], tmflgs.predmap[i]);
         strcpy(fnmlonlat,tmflgs.fnmlonlat);
         strcpy(fnmstxt,tmflgs.fnmstxt);
         strcpy(fnmelev,tmflgs.fnmelev);
         strcpy(fnmtveg,tmflgs.fnmtveg);
         strcpy(fnmclds,tmflgs.fnmclds);
         strcpy(fnmtair,tmflgs.fnmtair);
         strcpy(fnmprec,tmflgs.fnmprec);
         strcpy(fnmnirr,tmflgs.fnmnirr);
         strcpy(fnmpar,tmflgs.fnmpar);
         strcpy(fnmlulc,tmflgs.fnmlulc);
         strcpy(fnmnpp,tmflgs.fnmnpp);

         if (stateflag == 1)
         {
              for (i = 0; i < MAXSTATE; i++)strcpy(fnmstate[i],tmflgs.fnmstate[i]);
         }

         for(i = 0; i < natmspred; i ++)
              strcpy(fnmclmpred[i], fnminout.fnmclmpred[i]);
         for(i = 0; i < ntempred; i ++)
              strcpy(fnmtempred[i], fnminout.fnmtempred[i]);

         strcpy(mdir, tmflgs.mdir);

         strcpy(fnmkdout, fnminout.fnmkdout);
         for(i = 0; i < MAXNPAR; i ++)
              lhflag[i] = tmflgs.lhflag[i];

         whichcmnt = tmflgs.whichcmnt;
         kdinflg = tmflgs.kdinflg;
         kdoutflg = tmflgs.kdoutflg;

         f1stnd = tmflgs.f1stnd;


         telmnt[0].tem.avlnflag = temstac.avlnflag;
         telmnt[0].tem.nfeed = temstac.nfeed ;
         telmnt[0].tem.initbase = temstac.initbase;
         telmnt[0].tem.baseline = temstac.baseline;
         telmnt[0].tem.moistlim = temstac.moistlim;
         telmnt[0].tem.strteq = temstac.strteq;
         telmnt[0].tem.maxyears = temstac.maxyears;
         telmnt[0].tem.runsize = temstac.runsize;
         telmnt[0].tem.maxnrun = temstac.maxnrun;
         telmnt[0].tem.rheqflag = temstac.rheqflag;
         telmnt[0].tem.wtol = temstac.wtol;
         telmnt[0].tem.ctol = temstac.ctol;
         telmnt[0].tem.ntol = temstac.ntol;
         telmnt[0].tem.startyr = temstac.startyr;
         telmnt[0].tem.endyr = temstac.endyr;
         telmnt[0].tem.diffyr = temstac.diffyr;

         telmnt[0].tem.inittol = temstac.inittol;
         telmnt[0].tem.maxit = temstac.maxit;
         telmnt[0].tem.maxitmon = temstac.maxitmon;
    }
}
