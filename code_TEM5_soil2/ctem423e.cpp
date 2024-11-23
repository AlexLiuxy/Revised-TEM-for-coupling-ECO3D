/* **************************************************************
*****************************************************************
CTEM423.CPP - Calibration version of the Terrestrial Ecosystem
			Model Version 4.2
*****************************************************************
************************************************************** */
//This file includes fixes for the reset function.  The model previously gave
//different results when run from the beginning compared to when run after using
//the HOME and END keys.

//This file also includes fixes for the Water Balance display screens.

//Both these changes were made by CFD on September 22 and 25, 2000.  To find the
//changes search for "CFD" and look for the dates of the comments.  


#include <stdio.h>
#include <conio.h>
#include <iostream.h>
#include <fstream.h>
#include <fcntl.h>
#include <iomanip.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <cstring.h>

const int CYCLE = 12;
const int MAXRTIME = 600;

const int MAXPRED = 85;
const int MAXCAL = 21;

const double MISSING = -99999.9;

//Modules representing climate and TEM

#include "potclm423e.cpp"      // Potsclm Class
#include "pctem423e.cpp"       // PCTEM Class

Potsclm clm;
PCTTEM tem;


enum ecalkey { C_CFALL, C_C2NB, C_CNEVN, C_CNMN,  C_CMAX,
               C_KRB,   C_KDC,  C_NMAX,  C_NFALL, C_MNUP,
               C_FPCMX, C_SLA,  C_CO2,   C_TOL,   C_LFMX,
               C_KLFC,  C_COV,  C_CNSOL, C_MNFIX, C_NLOS,
               C_MOPT,

               C_STXT,   C_AVLNFLG, C_NFEED, C_BLINE, C_MOISTL,
               C_CLM,    C_SEY0,    C_SEY1,  C_SEY2,  C_RESET,
               C_MANFLG, C_AGSTATE, C_CALWIND };

//added by CFD September 25, 2000
enum wcalkey {
				C_INT_TOL, 	C_NEW_TEXT, 	C_NEW_CLIM, 	C_SWY0,
				C_SWY1,			 C_SWY2, 		C_SWY3,    	C_SWY4,
				C_WRESET, 	C_WCALWIND
				 };
/* *************************************************************
		 Function Declarations
************************************************************* */

inline ecalkey& next(ecalkey& s)
{
  return s = (C_CALWIND == s) ? C_CFALL : ecalkey(s+1);
}

inline ecalkey& prev(ecalkey& s)
{
  return s = (C_CFALL == s) ? C_CALWIND : ecalkey(s-1);
}
//added by CFD September 25, 2000
inline wcalkey& next(wcalkey& current)
{
  	return current = (C_WCALWIND == current) ? C_INT_TOL : wcalkey(current + 1);
}
//added by CFD September 25, 2000
inline wcalkey& prev(wcalkey& current)
{
  	return current = (C_INT_TOL == current) ? C_WCALWIND : wcalkey(current - 1);
}
void calecdin(void);
void calclmin(void);
void displayECalibPar(ecalkey& dcal);
//added by CFD September 22, 2000
void displayWCalibPar (wcalkey &);
void ecalibrate(int k);
void setCalVar(void);
//void showlphn(void);
void showclm(void);
void updateECalibPar(ecalkey& dcal, const double& factor, int& t);
//added by CFD September 22, 2000
void updateWCalibPar (wcalkey & , double & factor, int & t);
void wcalibrate(int k);

// *************************************************************

int reset;
int equilsol;

int dyr;
int newclm = 0;
int newtext = 0;
double ecalvar[21];
//added by CFD September 22, 2000
double wcalvar[1];
ecalkey evar;
//added by CFD September 22, 2000
wcalkey wvar;
//removed by CFD
//int ivar;

int kdinflg = 0;
int stateflag = 0;
int manflag = 0;

double col;
double row;
double lon;
double lat;
double elev;

char calclm[80];
ofstream flog1;

/* *************************************************************
**********************START MAIN PROGRAM************************
************************************************************* */

int main()
{

  int i;
  int k;
  int key;

  reset = 1;
  clrscr();

/* Open log file
*  flog1.open("stem.log"); */

  cout << endl << endl << "Initialize conditions for TEM:" << endl << endl;

  tem.pcinitRun(flog1);

  calclmin();
  calecdin();

  for (i = 0; i < CYCLE; i++)
  {
    tem.atms.co2[i] = tem.atms.co2level;
  }

  tem.atms.initco2 = tem.atms.co2level;

  showclm();

  tem.soil.xtext(tem.veg.cmnt, tem.soil.pctsilt, tem.soil.pctclay);

  tem.ECDsetELMNTstate(tem.veg.cmnt, tem.soil.psiplusc);
  tem.setELMNTflux();
  tem.setPrevState(tem.prevy, tem.y);

  cout << "Begin to initialize ecd" << endl;

  tem.setELMNTecd(kdinflg,tem.veg.cmnt, tem.soil.psiplusc);
  tem.setELMNTevap(stateflag, tem.veg.cmnt, tem.atms.pet, tem.atms.tair);
  tem.pcdisplayStext();

  tem.mintflag = 0;
  tem.tol = tem.inittol;
  tem.retry = 0;
  tem.endgrid = 0;

  tem.firstcal = 1;
  tem.topwind = 0;

  clrscr();
  cout << "Do you want to start with the WBM screen? ";
  if (toupper(getch()) == 'Y')
  {
    tem.calwind = 1;
    cout << "Y" << endl;
  }
  else
  {
    tem.calwind = 2;
    cout << "N" << endl;
  }

  dyr = 0;
  //Removed by CFD September 25, 2000
  //ivar = 1;
  tem.endeq = 1;

  tem.retry = 0;
  tem.endgrid = 0;
  equilsol = 0;

  tem.ag.state = 0;
  tem.ag.prvstate = 0;

  while (dyr > -10)
  {
    tem.pcstepyr(dyr, tem.mintflag, tem.tol);

    // Update annual agricultural product pools and fluxes

    tem.ag.updateyr(dyr);

    ++dyr;

// Check to see if steady state conditions have been reached.

    if (equilsol == 1) { tem.endgrid = 1; }

    if (tem.equil == 1 && dyr >= tem.strteq && equilsol == 0)
    {
      if (tem.wtol >= fabs(tem.atms.yrrain + tem.atms.yrsnowfall \
                      - tem.atms.yreet - tem.soil.yrrrun - tem.soil.yrsrun) \
                      && tem.nfeed == 0 && tem.rheqflag == 0 \
                      && (tem.ctol >= fabs(tem.veg.yrnpp - tem.veg.yrltrc)) \
                      && (tem.ctol >= fabs(tem.veg.yrnpp - tem.veg.yrltrc)))
      {
 	equilsol = 1;
      }
      if (tem.wtol >= fabs(tem.atms.yrrain + tem.atms.yrsnowfall \
                      - tem.atms.yreet - tem.soil.yrrrun - tem.soil.yrsrun) \
                      && tem.nfeed == 0 && tem.rheqflag == 1 \
                      && (tem.ctol >= fabs(tem.yrnep))\
                      && (tem.ctol >= fabs(tem.veg.yrnpp - tem.veg.yrltrc)) \
                      && (tem.ctol >= fabs(tem.veg.yrltrc - tem.microbe.yrrh)))
      {
	equilsol = 1;
	  }
	  if (tem.wtol >= fabs(tem.atms.yrrain + tem.atms.yrsnowfall \
				   - tem.atms.yreet - tem.soil.yrrrun - tem.soil.yrsrun) \
				   && tem.nfeed == 1 && tem.rheqflag == 1 \
				   && (tem.ntol >= fabs(tem.soil.yrnin - tem.soil.yrnlost)) \
				   && (tem.ntol >= fabs(tem.veg.yrnup - tem.veg.yrltrn)) \
				   && (tem.ntol >= fabs(tem.veg.yrnup - tem.microbe.yrnmin)) \
				   && (tem.ntol >= fabs(tem.veg.yrltrn - tem.microbe.yrnmin)) \
				   && (tem.ctol >= fabs(tem.yrnep)) \
				   && (tem.ctol >= fabs(tem.veg.yrnpp - tem.veg.yrltrc)) \
				   && (tem.ctol >= fabs(tem.veg.yrltrc - tem.microbe.yrrh)))
      {
	equilsol = 1;
      }
    }

    if (tem.calwind == 1 )
    {
      tem.pcdisplayYrWBM();

      window(46,24,80,24);
      gotoxy(1,1);
      printf(" 4.2.3-15JUN00  YRNETW = %8.6lf  ",
            (tem.atms.yrrain+tem.atms.yrsnowfall-tem.atms.yreet
            -tem.soil.yrrrun-tem.soil.yrsrun));

    }
    else
    {
      if (tem.ag.state == 0) { tem.pcdisplayYrTEM(manflag); }
      else { tem.pcdisplayYrAgTEM(manflag); }

      tem.pcdisplayLimits();

      window(46,24,80,24);
      gotoxy(1,1);
      printf("4.2.3-15JUN00 YRNEP = %8.3lf     ",tem.yrnep);
    }

    if (tem.adapttol == 1 && (dyr+1) == tem.maxyears)
    {
      tem.endgrid = 1;
      tem.tolbomb = 1;
    }

    if (manflag == 1)
    {
      k = 79;
      if (tem.calwind == 1) { wcalibrate(k); }
      else { ecalibrate(k); }
    }

    if (kbhit() != 0)
    {
      key = getch();
      if (key == 0)
      {
	k = getch();
	if (k == 79)
        {
	  if (tem.calwind == 1) { wcalibrate(k); }
	  else { ecalibrate(k); }
	}
      }
    }

    if (tem.endgrid == 1)
    {
      k = 79;
      if (tem.calwind == 1) { wcalibrate(k); }
      else { ecalibrate(k); }
    }
  }

  return 0;

};

/* *************************************************************
*********************END OF MAIN PROGRAM************************
************************************************************* */


/* *************************************************************
************************************************************* */

void calclmin(void)
{

  int i;
  int j;
  ifstream fclm;
  char dummy[12];
  int dumveg;

  tem.soil.getecd(tem.soilfile);

  cout << endl << "Enter file with climatic data: ";
  cin >> calclm;
  fclm.open(calclm,ios::in);

  if (!fclm)
  {
    cerr << "\nCannot open " << calclm << " for data input" << endl;
    exit(-1);
  }

// Read in basic climate data

  fclm >> dummy >> col;
  fclm >> dummy >> row;
  fclm >> dummy >> tem.elev;
  fclm >> dummy >> dumveg;
  fclm >> dummy >> tem.soil.pctsand;
  fclm >> dummy >> tem.soil.pctsilt;
  fclm >> dummy >> tem.soil.pctclay;
  fclm >> dummy >> tem.soil.wsoil;
  fclm >> dummy >> tem.ag.RAP;

  fclm >> dummy >> dummy >> dummy >> dummy;

  for (i = 0; i < CYCLE; i++)
  {
    fclm >> j >> clm.clds[i] >> tem.atms.tair[i] >> tem.atms.prec[i];
  }

  fclm.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void calecdin(void)
{

  char ecd[80];

  tem.soil.getrootz(tem.rootfile);

  tem.veg.getecd(tem.vegfile);
  tem.veg.leafyrs = 10;
  tem.veg.getleafecd(tem.leaffile);

  tem.microbe.getvegecd(tem.mcrvfile);
  tem.ag.getecd(tem.agfile);

  cout << endl << "Enter name of file with parameter values: ";
  cin >> ecd;
  tem.getsitecd(0, ecd);

  tem.veg.adjc2n = 1.0 + (tem.veg.dc2n * (tem.atms.co2[11] - tem.atms.initco2));
  tem.veg.cneven[tem.veg.cmnt] = tem.veg.initcneven[tem.veg.cmnt]
                                 * tem.veg.adjc2n;

  tem.soil.xtext(tem.veg.cmnt, tem.soil.pctsilt, tem.soil.pctclay);

  tem.pcdisplayInitState(ecd,col,row);

  tem.pcdisplayVegECD();

  tem.pcdisplayOtherECD();

};

/***************************************************************
***************************************************************/


/* *************************************************************
************************************************************** */

void displayECalibPar(ecalkey& dcal)
{

  switch(dcal)
  {
	case C_CFALL:  printf("CFALL = %8.6lf   ",tem.veg.cfall[tem.veg.cmnt]);
				   break;
	case C_C2NB:   printf("PLNTCNB = %8.2lf ",tem.veg.c2nb[tem.veg.cmnt]);
				   break;
	case C_CNEVN:  printf("INCNEV = %8.2lf  ",tem.veg.initcneven[tem.veg.cmnt]);
				   break;
	case C_CNMN:   printf("CNMIN = %8.2lf   ",tem.veg.cnmin[tem.veg.cmnt]);
				   break;
	case C_CMAX:   printf("CMAX = %8.2lf    ",tem.veg.cmax);
				   break;
	case C_KRB:    printf("KRB = %8.6lf     ",tem.veg.krb[tem.veg.cmnt]);
				   break;
	case C_KDC:    printf("KDC = %8.6lf     ",tem.microbe.kdc);
				   break;
	case C_NMAX:   printf("NMAX = %8.4lf    ",tem.veg.nmax);
				   break;
	case C_NFALL:  printf("NFALL = %8.6lf   ",tem.veg.nfall[tem.veg.cmnt]);
				   break;
	case C_MNUP:   printf("MNUP = %8.4lf    ",tem.microbe.nup);
				   break;
	case C_FPCMX:  printf("FPCMX = %8.6lf   ",tem.veg.fpcmax[tem.veg.cmnt]);
				   break;
	case C_SLA:    printf("SLA = %8.4lf     ",tem.veg.sla[tem.veg.cmnt]);
				   break;
	case C_CO2:    printf("CO2 = %8.1lf     ",tem.atms.co2level);
				   break;
	case C_TOL:    printf("INT TOL = %8.6lf ",tem.tol);
				   break;
	case C_LFMX:   printf("LEAFMXC = %8.2lf ",tem.veg.leafmxc[tem.veg.cmnt]);
				   break;
	case C_KLFC:   printf("KLEAFC = %8.2lf   ",tem.veg.kleafc[tem.veg.cmnt]);
				   break;
	case C_COV:    printf("COV = %8.6lf     ",tem.veg.cov[tem.veg.cmnt]);
				   break;
	case C_CNSOL:  printf("CNSOIL = %8.2lf  ",tem.microbe.cnsoil[tem.veg.cmnt]);
				   break;
	case C_MNFIX:  printf("NFIX   = %8.1lf  ",tem.microbe.nfixpar[tem.veg.cmnt]);
				   break;
	case C_NLOS:   printf("NLOSS  = %8.1lf  ",tem.soil.nloss[tem.veg.cmnt]);
				   break;
	case C_MOPT:   printf("MOPT = %8.4lf    ",tem.microbe.moistopt[tem.veg.cmnt]);
				   break;
	case C_STXT:   printf("NEW TEXTURE? (%4.2lf)",tem.soil.psiplusc);
				   break;
	case C_AVLNFLG: if (tem.avlnflag == 1) { printf("AVAILN = ON"); }
				if (tem.avlnflag == 0) { printf("AVAILN = OFF"); }
				break;
	case C_NFEED:   if (tem.nfeed == 1) { printf("NFEED = ON");  }
				if (tem.nfeed == 0) { printf("NFEED = OFF"); }
				break;
	case C_BLINE:   if (tem.baseline == 1) { printf ("BLINE = ON"); }
				if (tem.baseline == 0) { printf ("BLINE = OFF"); }
				break;
	case C_MOISTL:  if (tem.moistlim == 1) { printf("MOISTLIM = ON");  }
				if (tem.moistlim == 0) { printf("MOISTLIM = OFF"); }
				break;
	case C_CLM:     printf("NEW CLIMATE?");
				break;
	case C_SEY0:    tem.displayOptionalEflx(tem.sey[0]);
				break;
	case C_SEY1:    tem.displayOptionalEflx(tem.sey[1]);
				break;
	case C_SEY2:    tem.displayOptionalEflx(tem.sey[2]);
				break;
	case C_RESET:   if (reset == 1) { cout << "reset = ON" << endl; }
				else { cout << "reset = OFF" << endl; }
				break;
	case C_MANFLG:  if (manflag == 1) { cout << "manual = ON" << endl; }
				else { cout << "manual = OFF" << endl; }
				break;
	case C_AGSTATE: if (tem.ag.state == 1) { cout << "agstate = ON" << endl; }
				else { cout << "agstate = OFF" <<endl; }
				break;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */
//added by CFD September 22, 2000
void displayWCalibPar (wcalkey & whichWcalib)
{
	switch(whichWcalib)
	{
		case C_INT_TOL:
		printf("INT TOL = %8.6lf ",tem.tol);
		break;
		case C_NEW_TEXT:
		printf("NEW TEXTURE? (%4.2lf)",tem.soil.psiplusc);
		break;
		case C_NEW_CLIM:
		printf("NEW CLIMATE?       ");
		break;
		case C_SWY0:
		tem.displayOptionalWflx(tem.swy[0]);
			break;
		case C_SWY1:
		tem.displayOptionalWflx(tem.swy[1]);
			break;
		case C_SWY2:
		tem.displayOptionalWflx(tem.swy[2]);
			break;
		case C_SWY3:
		tem.displayOptionalWflx(tem.swy[3]);
			break;
		case C_SWY4:
		tem.displayOptionalWflx(tem.swy[4]);
			break;
		case C_WRESET:
		if (reset == 1)
			cout << "reset = ON" << endl;
			else
			cout << "reset = OFF" << endl;
		 break;

	  default:
		break;
   }
}

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void ecalibrate(int t)
{

  double efactor;
  double invar;

  tem.endgrid = 0;

  window(43,23,80,23);
  printf("AN = %1d, NF = %1d, BL = %1d, MF = %1d, RS = %1d",tem.avlnflag,
		tem.nfeed,tem.baseline,tem.moistlim,reset);

  window(57,1,64,1);
  tem.displayOptionalEflx(tem.sey[0]);

  window(65,1,72,1);
  tem.displayOptionalEflx(tem.sey[1]);

  window(73,1,79,1);
  tem.displayOptionalEflx(tem.sey[2]);

  window(1,23,15,23);
  cout << "CALIBRATION:   " << endl;

  if (tem.firstcal == 1)
  {
	evar = C_CFALL;
	tem.firstcal = 0;
  }

  window(16,23,34,23);
  delline();
  displayECalibPar(evar);

  while (t == 79 || t == 72 || t == 75 || t == 77 || t == 80 || t == 82)
  {

	if (tem.tolbomb == 0 && tem.intbomb == 0)
	{
	  t = getch(); t = getch();
	}
	else
	{

	  // If simulation reaches RUNSIZE automatically adjust the
	  //   tolerance level (i.e. tol)
	  if (tem.tolbomb == 1)
	  {
		tem.retry = 1;
		tem.tolbomb = 0;
		tem.nattempt += 1;
		evar = C_TOL;
		if (tem.nattempt < tem.maxnrun)
		{
	  tem.tol /= 10.0;
	  t = 0;
		}
		else
		{
	  tem.tol = tem.inittol;
	  tem.nattempt = 0;
		}
	  }

	  // If simulation "black holes", allow user to adjust tolerance level
	  if (tem.intbomb == 1)
	  {
		tem.retry = 1;
		tem.intbomb = 0;
		evar = C_TOL;
	  }
	}

	if (tem.intbomb == 1) { evar = C_TOL; }

	// Set calvar[] with current TEM parameter values
	setCalVar();

	efactor = 1.0;
	window(35,23,42,23);
	delline();

	// Select or modify a specific calibration parameter
	//   represented by calvar[]

	switch (t)
	{
	  case 79: break;
	  case 72: efactor = 1.01; break;
	  case 80: efactor = 0.99; break;
	  case 77: next(evar); break;
	  case 75: prev(evar); break;
	  case 82: if (ecalkey(evar) < MAXCAL)
			   {
			 cin >> invar;
			 ecalvar[evar] = invar;
		   }
		   break;
	  default: t = 0; break;
	}

	window(35,23,42,23);
	delline();

	window(16,23,34,23);
	delline();

	// Update TEM calibration parameters with new calvar[] value
	updateECalibPar(evar, efactor, t);

	if (t == 0) break;
  }


  if (reset == 1 || tem.retry == 1)
  {
	tem.retry = 0;
	dyr = 0;
	tem.ECDsetELMNTstate(tem.veg.cmnt, tem.soil.psiplusc);
	tem.setELMNTflux();
	tem.setPrevState(tem.prevy, tem.y);
	tem.prevy[tem.I_UNRMLF] = 0.5;
	tem.prevy[tem.I_LEAF] = 0.5;
	//added by CFD to account for differences between running the model
	//before and after using the END and HOME keys - September 21, 2000
	tem.setELMNTecd(kdinflg,tem.veg.cmnt, tem.soil.psiplusc);
	tem.setELMNTevap(stateflag, tem.veg.cmnt, tem.atms.pet, tem.atms.tair);
  }


  if (newtext == 1 || newclm == 1)
  {
	window(1,1,80,24);
	gotoxy(1,1);
	clrscr();
	tem.setELMNTecd(kdinflg,tem.veg.cmnt,tem.soil.psiplusc);
	tem.pcdisplayStext();
	tem.topwind = 0;
	newclm = 0;
	newtext = 0;
  }


  if (equilsol == 1)
  {
	equilsol = 0;
	dyr = 0;
  }

};

/* *************************************************************
************************************************************** */


/***************************************************************
***************************************************************/

void setCalVar(void)
{

  ecalvar[C_CFALL]  =  tem.veg.cfall[tem.veg.cmnt];
  ecalvar[C_C2NB]  =  tem.veg.c2nb[tem.veg.cmnt];
  ecalvar[C_CNEVN]  =  tem.veg.initcneven[tem.veg.cmnt];
  ecalvar[C_CNMN]  =  tem.veg.cnmin[tem.veg.cmnt];
  ecalvar[C_CMAX]  =  tem.veg.cmax;
  ecalvar[C_KRB]  =  tem.veg.krb[tem.veg.cmnt];
  ecalvar[C_KDC]  =  tem.microbe.kdc;
  ecalvar[C_NMAX]  =  tem.veg.nmax;
  ecalvar[C_NFALL]  =  tem.veg.nfall[tem.veg.cmnt];
  ecalvar[C_MNUP]  =  tem.microbe.nup;
  ecalvar[C_FPCMX] =  tem.veg.fpcmax[tem.veg.cmnt];
  ecalvar[C_SLA] =  tem.veg.sla[tem.veg.cmnt];
  ecalvar[C_CO2] =  tem.atms.co2level;
  ecalvar[C_TOL] =  tem.tol;
  ecalvar[C_LFMX] =  tem.veg.leafmxc[tem.veg.cmnt];
  ecalvar[C_KLFC] =  tem.veg.kleafc[tem.veg.cmnt];
  ecalvar[C_COV] =  tem.veg.cov[tem.veg.cmnt];
  ecalvar[C_CNSOL] =  tem.microbe.cnsoil[tem.veg.cmnt];
  ecalvar[C_MNFIX] =  tem.microbe.nfixpar[tem.veg.cmnt];
  ecalvar[C_NLOS] =  tem.soil.nloss[tem.veg.cmnt];
  ecalvar[C_MOPT] =  tem.microbe.moistopt[tem.veg.cmnt];

}

/***************************************************************
***************************************************************/


/***************************************************************
***************************************************************/

void showclm()
{

  tem.pcdisplayClm(calclm, clm.clds);

  // Determine monthly gross irradiance (GIRR)

  lat = (double) row;

  clm.yrsumday = 0.0;
  for (int dm = 0; dm < CYCLE; dm++)
  {
	clm.girr[dm] = clm.xgirr(lat,dm,clm.yrsumday);

// Determine net irradiance (NIRR)

	tem.atms.nirr[dm] = clm.xnirr(clm.clds[dm],clm.girr[dm]);

// Determine photosynthetically active radiation (PAR)

	tem.atms.par[dm] = clm.xpar(clm.clds[dm],tem.atms.nirr[dm]);
  }

  tem.pcdisplayPAR(calclm, clm.girr);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void updateECalibPar(ecalkey& dcal, const double& factor, int& t)
{
  int i;

  switch (dcal)
  {
	case C_CFALL :  tem.veg.cfall[tem.veg.cmnt] = ecalvar[C_CFALL] * factor;
					printf("CFALL = %8.6lf   ",tem.veg.cfall[tem.veg.cmnt]);
				 break;
	case C_C2NB :   tem.veg.c2nb[tem.veg.cmnt] = ecalvar[C_C2NB] * factor;
					printf("PLNTCNB = %8.2lf ",tem.veg.c2nb[tem.veg.cmnt]);
				 break;
	case C_CNEVN :  tem.veg.initcneven[tem.veg.cmnt] = ecalvar[C_CNEVN] * factor;
				 printf("INCNEV = %8.2lf  ",tem.veg.initcneven[tem.veg.cmnt]);
				 break;
	case C_CNMN :   tem.veg.cnmin[tem.veg.cmnt] = ecalvar[C_CNMN] * factor;
				 printf("CNMIN = %8.2lf   ",tem.veg.cnmin[tem.veg.cmnt]);
				 break;
	case C_CMAX :   tem.veg.cmax = ecalvar[C_CMAX] * factor;
				 printf("CMAX = %8.2lf    ",tem.veg.cmax);
				 break;
	case C_KRB :    tem.veg.krb[tem.veg.cmnt] = ecalvar[C_KRB] * factor;
				 printf("KRB = %8.6lf     ",tem.veg.krb[tem.veg.cmnt]);
				 break;
	case C_KDC :    tem.microbe.kdc = ecalvar[C_KDC] * factor;
				 printf("KDC = %8.6lf     ",tem.microbe.kdc);
				 break;
	case C_NMAX :   tem.veg.nmax = ecalvar[C_NMAX] * factor;
				 printf("NMAX = %8.4lf    ",tem.veg.nmax);
				 break;
	case C_NFALL :  tem.veg.nfall[tem.veg.cmnt] = ecalvar[C_NFALL] * factor;
				 printf("NFALL = %8.6lf   ",tem.veg.nfall[tem.veg.cmnt]);
				 break;
	case C_MNUP :   tem.microbe.nup = ecalvar[C_MNUP] * factor;
				 printf("NUP = %8.4lf     ",tem.microbe.nup);
				 break;
	case C_FPCMX :  tem.veg.fpcmax[tem.veg.cmnt] = ecalvar[C_FPCMX] * factor;
				 printf("FPCMX = %8.6lf   ",tem.veg.fpcmax[tem.veg.cmnt]);
				 break;
	case C_SLA :    tem.veg.sla[tem.veg.cmnt] = ecalvar[C_SLA] * factor;
				 printf("SLA = %8.4lf    ",tem.veg.sla[tem.veg.cmnt]);
				 break;
	case C_CO2 :    tem.atms.co2level = ecalvar[C_CO2] * factor;
				 printf("CO2 = %8.1lf     ",tem.atms.co2level);
				 for (i = 0; i < CYCLE; i++)
					 {
				   tem.atms.co2[i] = tem.atms.co2level;
				 }
				 break;
	case C_TOL :    tem.tol = ecalvar[C_TOL] * factor;
				 printf("INT TOL = %8.6lf ",tem.tol);
				 break;
	case C_LFMX :   tem.veg.leafmxc[tem.veg.cmnt] = ecalvar[C_LFMX] * factor;
				 printf("LEAFMXC = %8.6lf  ",tem.veg.leafmxc[tem.veg.cmnt]);
				 break;
	case C_KLFC :   tem.veg.kleafc[tem.veg.cmnt] = ecalvar[C_KLFC] * factor;
				 printf("KLEAFC = %8.6lf   ",tem.veg.kleafc[tem.veg.cmnt]);
//             tem.microbe.decay = 0.26299 +
//						(1.14757*tem.microbe.procmntos[tem.veg.cmnt])
//                     - (0.42956*pow(tem.microbe.procmntos[tem.veg.cmnt],2.0));
				 break;
	case C_COV :     tem.veg.cov[tem.veg.cmnt] = ecalvar[C_COV] * factor;
				 printf("COV = %8.6lf    ",tem.veg.cov[tem.veg.cmnt]);
				 break;
	case C_CNSOL :   tem.microbe.cnsoil[tem.veg.cmnt]
					= ecalvar[C_CNSOL] * factor;
				 printf("CNSOIL = %8.2lf  ",tem.microbe.cnsoil[tem.veg.cmnt]);
				 break;
	case C_MNFIX :   tem.microbe.nfixpar[tem.veg.cmnt]
					= ecalvar[C_MNFIX] * factor;
				 printf("NFIX   = %8.1lf  ",tem.microbe.nfixpar[tem.veg.cmnt]);
				 break;
	case C_NLOS :    tem.soil.nloss[tem.veg.cmnt] = ecalvar[C_NLOS] * factor;
				 printf("NLOSS  = %8.1lf  ",tem.soil.nloss[tem.veg.cmnt]);
				 break;
	case C_MOPT :    tem.microbe.moistopt[tem.veg.cmnt]
					= ecalvar[C_MOPT] * factor;
				 printf("MOPT = %8.4lf   ", tem.microbe.moistopt[tem.veg.cmnt]);
				 break;
	case C_STXT :    if (factor != 1)
					 {
				   window(1,1,80,24);
				   gotoxy(1,1);
				   clrscr();
				   cout << "Enter proportion of silt plus clay: ";
				   cin >> tem.soil.psiplusc;
				   cout << "Enter percent clay: ";
				   cin >> tem.soil.pctclay;
				   tem.soil.pctsilt = (tem.soil.psiplusc * 100.0)
										  - tem.soil.pctclay;
				   tem.soil.pctsand = (1.0 - tem.soil.psiplusc)* 100.0;
				   tem.soil.xtext(tem.veg.cmnt,
									  tem.soil.pctsilt, tem.soil.pctclay);
				   newtext = 1;
				   reset = 1;
				   t = 0;
				 }
				 else
					 {
				   printf("NEW TEXTURE? (%4.2lf)",tem.soil.psiplusc);
				 }
				 break;
	case C_AVLNFLG : if (factor != 1) { tem.avlnflag += 1; }
				 if (tem.avlnflag > 1) { tem.avlnflag = 0; }
				 if (tem.avlnflag == 1) { cout << "AVAILN = ON"; }
				 if (tem.avlnflag == 0) { cout << "AVAILN = OFF"; }
				 break;
	case C_NFEED :   if (factor != 1) { tem.nfeed += 1; }
				 if (tem.nfeed > 1) { tem.nfeed = 0; }
				 if (tem.nfeed == 1) { cout << "NFEED = ON"; }
				 if (tem.nfeed == 0) { cout << "NFEED = OFF"; }
				 break;
	case C_BLINE :   if (factor != 1) { tem.baseline += 1; }
				 if (tem.baseline > 1) { tem.baseline = 0; }
				 if (tem.baseline == 1) { cout << "BLINE = ON"; }
				 if (tem.baseline == 0) { cout << "BLINE = OFF"; }
				 break;
	case C_MOISTL :  if (factor != 1) { tem.moistlim += 1; }
				 if(tem.moistlim > 1) { tem.moistlim = 0; }
				 if (tem.moistlim == 1) { cout << "MOISTLIM = ON"; }
				 if (tem.moistlim == 0) { cout << "MOISTLIM = OFF"; }
				 break;
	case C_CLM :     if (factor != 1)
					 {
				   window(1,1,80,24);
				   gotoxy(1,1);
				   clrscr();
				   calclmin();
				   showclm();
				   tem.soil.xtext(tem.veg.cmnt,
									  tem.soil.pctsilt, tem.soil.pctclay);
				   newclm = 1;
				   reset = 1;
				   t = 0;
				 }
				 else
					 {
				   cout << "NEW CLIMATE?       ";
				 }
				 break;
	case C_SEY0 :    if (factor > 1) { tem.next(tem.sey[0]); }
				 if (factor < 1) { tem.prev(tem.sey[0]); }
				 tem.displayOptionalEflx(tem.sey[0]);
					 break;
	case C_SEY1 :    if (factor > 1) { tem.next(tem.sey[1]); }
					 if (factor < 1) { tem.prev(tem.sey[1]); }
				 tem.displayOptionalEflx(tem.sey[1]);
				 break;
	case C_SEY2 :    if (factor > 1) { tem.next(tem.sey[2]); }
				 if (factor < 1) { tem.prev(tem.sey[2]); }
				 tem.displayOptionalEflx(tem.sey[2]);
				 break;
	case C_RESET :   if (factor != 1) { reset += 1; }
					 if (reset > 1) { reset = 0; }
				 if (reset == 1) { cout << "reset = ON"; }
				 else { cout << "reset = OFF"; }
				 break;
	case C_MANFLG :  if (factor != 1) { manflag += 1; }
				 if (manflag > 1) { manflag = 0; }
				 if (manflag == 1) { cout << "manual = ON"; }
				 else { cout << "manual = OFF"; }
				 break;
	case C_AGSTATE : if (factor != 1) { tem.ag.state += 1; }
				 if (tem.ag.state > 1) { tem.ag.state = 0; }
				 if (tem.ag.state == 1) { cout << "agstate = ON"; }
				 else { cout << "agstate = OFF"; }
				 break;
	case C_CALWIND : if (factor != 1) { tem.calwind += 1; }
				 if (tem.calwind > 2) { tem.calwind = 1; }
				 if (tem.calwind == 1)
					 {
				   tem.firstcal = 1;
				   cout << "screen = WBM" << endl;
				   tem.topwind = 0;
				 }
				 else
					 {
				   tem.firstcal = 0;
				   cout << "screen = TEM" << endl;
				   tem.topwind = 1;
				 }
				 break;
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */
//added by CFD September 22, 2000 - whole function
void updateWCalibPar (wcalkey & ivar, double & factor, int & t)
{
   switch (ivar)
   {
	case C_INT_TOL: //Integrator tolerance
			tem.tol = wcalvar[C_INT_TOL] * factor;
			//delline();
			printf("INT TOL = %8.6lf ",tem.tol);
		break;
	  case C_NEW_TEXT: //Psiplusc
		//if they want to change this - up or down arrow, doesn't matter
		if (factor != 1)
		{
			window(1,1,80,24);
			gotoxy(1,1);
			clrscr();
			cout << "Enter proportion of silt plus clay: ";
			cin >> tem.soil.psiplusc;
			cout << "Enter percent clay: ";
			cin >> tem.soil.pctclay;
			tem.soil.pctsilt = (tem.soil.psiplusc * 100.0) - tem.soil.pctclay;
			tem.soil.pctsand = (1.0 - tem.soil.psiplusc)* 100.0;
			tem.soil.xtext(tem.veg.cmnt, tem.soil.pctsilt, tem.soil.pctclay);
			newtext = 1;
			reset = 1;
			t = 0;
		}
		 //if the user doesn't want to change this
		else
		{
			delline();
			printf("NEW TEXTURE? (%4.2lf)",tem.soil.psiplusc);
		}
		break;
	  case C_NEW_CLIM:
		//if they want to change the climate data
		 //up or down arrow, doesn't matter
		if (factor != 1)
		{
			window(1,1,80,24);
			gotoxy(1,1);
			clrscr();
			calclmin();
			showclm();
			tem.soil.xtext(tem.veg.cmnt, tem.soil.pctsilt, tem.soil.pctclay);
			newclm = 1;
			reset = 1;
			t = 0;
		}
		else //if the user doesn't want to change the climate
		{
			delline();
			cout << "NEW CLIMATE?       ";
		}
		break;
	  case C_SWY0: //the 1st of 5 optional water balance model varibles that
				//can be displayed on a monthly scale (upper right corner)
		//up arrow
		if (factor > 1)
			tem.next(tem.swy[0]);
		 //down arrow
		if (factor < 1)
			tem.prev(tem.swy[0]);
		delline();
		tem.displayOptionalWflx(tem.swy[0]);
		break;
	  case C_SWY1: //the 2nd of 5 optional water balance model varibles that
				//can be displayed on a monthly scale (upper right corner)
		//up arrow
		if (factor > 1)
			tem.next(tem.swy[1]);
		 //down arrow
		if (factor < 1)
			tem.prev(tem.swy[1]);
		delline();
		tem.displayOptionalWflx(tem.swy[1]);
			break;
	case C_SWY2: //the 3rd of 5 optional water balance model varibles that
				//can be displayed on a monthly scale (upper right corner)
		//up arrow
		if (factor > 1)
			tem.next(tem.swy[2]);
		 //down arrow
		if (factor < 1)
			tem.prev(tem.swy[2]);
		delline();
		tem.displayOptionalWflx(tem.swy[2]);
		break;
	case C_SWY3: //the 4th of 5 optional water balance model varibles that
	//up arrow
		if (factor > 1)
			tem.next(tem.swy[3]);
		 //down arrow
		if (factor < 1)
			tem.prev(tem.swy[3]);
		delline();
		tem.displayOptionalWflx(tem.swy[3]);
		break;
	case C_SWY4: //the 5th of 5 optional water balance model varibles that
				//can be displayed on a monthly scale (upper right corner)
		//up arrow
		if (factor > 1)
			tem.next(tem.swy[4]);
		 //down arrow
		if (factor < 1)
			tem.prev(tem.swy[4]);
		delline();
		tem.displayOptionalWflx(tem.swy[4]);
		break;
	case C_WRESET: //the user wants to change if the model is reset or not
		//up or down arrow
		if (factor != 1)
			++reset;
		 //down arrow
		if (reset > 1)
			reset = 0;
		delline();
		 //display the choice
		if (reset == 1)
			cout << "reset = ON";
		else
			cout << "reset = OFF";
		break;
	  case C_WCALWIND: //Change the calibration screen
		// up or down arrow
		if (factor != 1)
			tem.calwind += 1;
		if (tem.calwind > 2)
			tem.calwind = 1;
		if (tem.calwind == 2)
		{
			tem.firstcal = 1;
			delline();
			cout << "screen = TEM" << endl;
			tem.topwind = 0;
		}
		else
		{
			tem.firstcal = 0;
			delline();
			cout << "screen = WBM" << endl;
			tem.topwind = 1;
		}
		break;
	default: cout<<"No ivar found\n";
   }
}


/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void wcalibrate(int t)
{
  double invar;
  //moved to global variable to match ecalibrate by CFD September 22, 2000
	//int calvar[1];
  double factor = 1.0;
  tem.endgrid = 0;

  window(43,23,80,23);
  printf("RS = %1d             WTOL = %8.6lf",reset,tem.wtol);

   //changed the next 5 pairs of lines - CFD September 25, 2000
	window(44,1,50,1);
	tem.displayOptionalWflx(tem.swy[0]);

	window(51,1,57,1);
	tem.displayOptionalWflx(tem.swy[1]);

	window(58,1,65,1);
	tem.displayOptionalWflx(tem.swy[2]);

	window(66,1,73,1);
	tem.displayOptionalWflx(tem.swy[3]);

	window(74,1,80,1);
	tem.displayOptionalWflx(tem.swy[4]);

  window(1,23,15,23);
  cout << "CALIBRATION:   " << endl;

  if (tem.firstcal == 1)
  {
	//changed September 25, 2000 - CFD
	wvar = C_INT_TOL;
	tem.firstcal = 0;
  }

  window(16,23,34,23);
  delline();
  //changed September 25, 2000 - CFD
   displayWCalibPar (wvar);

  while (t == 79 || t == 72 || t == 75 || t == 77 || t == 80 || t == 82)
  {

	if (tem.tolbomb == 0 && tem.intbomb == 0)
	{
	  t = getch(); t = getch();
	}
	else
	{
	  if (tem.tolbomb == 1)
	  {
	tem.retry = 1;
	tem.tolbomb = 0;
	tem.nattempt += 1;
	//changed September 25, 2000 - CFD
   wvar = C_INT_TOL;
	if (tem.nattempt < tem.maxnrun)
		{
	  tem.tol /= 10.0;
	  t = 0;
	}
	else
		{
	  tem.tol = tem.inittol;
	  tem.nattempt = 0;
	}
	  }
	  if (tem.intbomb == 1)
	  {
	tem.retry = 1;
	tem.intbomb = 0;
	//changed September 25, 2000 - CFD
   wvar = C_INT_TOL;
	  }
	}

	if (tem.intbomb == 1)
	{ 	//changed September 25, 2000 - CFD
		wvar = C_INT_TOL;
	}
	 //changed September 22, 2000 - CFD
	wcalvar[0] =  tem.tol;

	factor = 1.0;
	window(35,23,42,23);
	printf(" ");
	switch (t)
	{
	  case 79: break;
	  case 72: factor = 1.01; break;
	  case 80: factor = 0.99; break;
	  //changed September 25, 2000 - CFD
	  case 77: next(wvar); break;
		//changed September 25, 2000 - CFD
	  case 75: prev(wvar); break;
	  //changed September 25, 2000 - CFD
	  case 82: if (wcalkey(wvar) <= C_NEW_CLIM)
			   {
		 cin >> invar;
		 wcalvar[wvar] = invar;
		   }
		   break;
	  default: t = 0; break;
	}

	window(35,23,42,23);
	delline();

	window(16,23,34,23);
	delline();

	//added by CFD September 22, 2000
	//This function replaces a switch statement
	updateWCalibPar (wvar, factor, t);

	if (t == 0) break;
  }


  if (reset == 1 || tem.retry == 1)
  {
	tem.retry = 0;
	dyr = 0;
	tem.ECDsetELMNTstate(tem.veg.cmnt, tem.soil.psiplusc);
	tem.setELMNTflux();
	tem.setPrevState(tem.prevy, tem.y);
	tem.prevy[tem.I_UNRMLF] = 0.5;
	tem.prevy[tem.I_LEAF] = 0.5;
	//added by CFD September 22, 2000 - next two lines
	tem.setELMNTecd(kdinflg,tem.veg.cmnt, tem.soil.psiplusc);
	tem.setELMNTevap(stateflag, tem.veg.cmnt, tem.atms.pet, tem.atms.tair);
  }


  if (newtext == 1 || newclm == 1)
  {
	window(1,1,80,24);
	gotoxy(1,1);
	clrscr();
	tem.setELMNTecd(kdinflg,tem.veg.cmnt,tem.soil.psiplusc);
	tem.pcdisplayStext();
	tem.topwind = 0;
	newclm = 0;
	newtext = 0;
  }

  if (equilsol == 1)
  {
	equilsol = 0;
	dyr = 0;
  }

};


