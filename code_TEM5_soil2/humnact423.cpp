/* **************************************************************
*****************************************************************
HUMANACT42.CPP - describes human disturbances to natural ecosystems

Modifications:

19990821 - DWK changed char dummy[12] to char dummy[30] in
 void getecd(ofstream& rflog1)
19990821 - DWK changed char dummy[12] to char dummy[30] in void
getecd(char ecd[80])
20000102 - DWK added compiler directives

*****************************************************************
************************************************************** */

#if !defined(HUMNACT423_H)
  #include "humnact423.hpp"
#endif

Humanact::Humanact()
{

  int i;

  c2n = 54.29;
//  cfall = 0.20;
//  nfall = 0.20;

  state = 0;
  prvstate = 0;

  for (i = 0; i < 10; i++)
  {
    initPROD10[i].carbon = 0.0;
    initPROD10[i].nitrogen = 0.0;
  }
  for (i = 0; i < 100; i++)
  {
    initPROD100[i].carbon = 0.0;
    initPROD100[i].nitrogen = 0.0;
  }

};

/* **************************************************************
************************************************************** */



/* **************************************************************
************************************************************** */

void Humanact::conversion(int& ez, Tveg42& veg, Tsoil4& soil)
{

  int dm;

  formPROD10.carbon  = prod10par[ez] * veg.plant[CYCLE-1].carbon;
  formPROD10.nitrogen  = prod10par[ez] * (veg.strctrl[CYCLE-1].nitrogen
						+ veg.labile[CYCLE-1].nitrogen);
  formPROD100.carbon = prod100par[ez] * veg.plant[CYCLE-1].carbon;
  formPROD100.nitrogen = prod100par[ez] * (veg.strctrl[CYCLE-1].nitrogen
						+ veg.labile[CYCLE-1].nitrogen);

  PROD10.carbon += formPROD10.carbon;
  PROD10.nitrogen += formPROD10.nitrogen;
  PROD100.carbon += formPROD100.carbon;
  PROD100.nitrogen += formPROD100.nitrogen;
  for (dm = 0; dm < CYCLE; dm++)
  {
	slash[dm].carbon = slashpar[ez] * veg.plant[CYCLE-1].carbon
						/ (double) CYCLE;
	slash[dm].nitrogen = slashpar[ez]  * (veg.strctrl[CYCLE-1].nitrogen
						+ veg.labile[CYCLE-1].nitrogen) / (double) CYCLE;
	vconvrtflx[dm].carbon = (vconvert[ez] * veg.plant[CYCLE-1].carbon)
						/ (double) CYCLE;
	sconvrtflx[dm].carbon =  (sconvert[ez] * soil.org[CYCLE-1].carbon)
						/ (double) CYCLE;
    convrtflx[dm].carbon = vconvrtflx[dm].carbon + sconvrtflx[dm].carbon;
	vconvrtflx[dm].nitrogen = ((1.0 - nvretconv[ez]) * vconvert[ez]
						* (veg.strctrl[CYCLE-1].nitrogen
						+ veg.labile[CYCLE-1].nitrogen))  / (double) CYCLE;
	sconvrtflx[dm].nitrogen = ((1.0 - nsretconv[ez]) * sconvert[ez]
						* soil.org[CYCLE-1].nitrogen)  / (double) CYCLE;
	convrtflx[dm].nitrogen = vconvrtflx[dm].nitrogen + sconvrtflx[dm].nitrogen;
	nvretent[dm] = (nvretconv[ez] * vconvert[ez]
					* (veg.strctrl[CYCLE-1].nitrogen
					+ veg.labile[CYCLE-1].nitrogen)) / (double) CYCLE;
	nsretent[dm] = (nsretconv[ez] * sconvert[ez] * soil.org[CYCLE-1].nitrogen)
					 / (double) CYCLE;
	nretent[dm] = nvretent[dm] + nsretent[dm];

	yrvconvrtC += vconvrtflx[dm].carbon;
    yrvconvrtN += vconvrtflx[dm].nitrogen;
    yrnrent += nretent[dm];
    yrnvrent += nvretent[dm];

  }

};
 /* **************************************************************
************************************************************** */

void Humanact::getecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the data file (.ECD) with agricultural"
	   << " parameter values: " << endl;
//cin >> ecd;
  fpara >> ecd;

  rflog1 << "Enter name of the soil data file (.ECD) with agricultural"
  		 << " parameter values: " << ecd << endl;

  getecd(ecd);

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Humanact::getecd(char ecd[80])
{

  int NUMVAR = 13;
  int i;
  int dcmnt;
  char dummy[30];
  ifstream infile;

  long update;

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy; }
  for (dcmnt = 1; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> dummy >> dummy;
	infile >> slashpar[dcmnt] >> vconvert[dcmnt] >> prod10par[dcmnt]
			>> prod100par[dcmnt];
	infile >> sconvert[dcmnt] >> nvretconv[dcmnt] >> nsretconv[dcmnt]
		   >> vrespar[dcmnt];
    infile >> cfall[dcmnt] >> nfall[dcmnt];
    infile >> update;
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Humanact::resetPROD()
{

  PROD1.carbon = 0.0;
  PROD1.nitrogen = 0.0;
  PROD1decay.carbon = 0.0;
  PROD1decay.nitrogen = 0.0;

  formPROD10.carbon  = 0.0;
  formPROD10.nitrogen  = 0.0;
  PROD10decay.carbon  = 0.0;
  PROD10decay.nitrogen  = 0.0;
  PROD10.carbon = 0.0;
  PROD10.nitrogen = 0.0;

  formPROD100.carbon = 0.0;
  formPROD100.nitrogen = 0.0;
  PROD100decay.carbon = 0.0;
  PROD100decay.nitrogen = 0.0;
  PROD100.carbon = 0.0;
  PROD100.nitrogen = 0.0;

  TOTPROD.carbon = 0.0;
  TOTPROD.nitrogen = 0.0;

  formTOTPROD.carbon = 0.0;
  formTOTPROD.nitrogen = 0.0;

  TOTPRODdecay.carbon = 0.0;
  TOTPRODdecay.nitrogen = 0.0;

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Humanact::updateyr(const int& dyr)
{

  int i;
  int j;
  int k;
  double yrtrash;

  prvstate = state;
  PROD1.carbon = formPROD1.carbon;
  PROD1.nitrogen = formPROD1.nitrogen;
  PROD1decay.carbon = PROD1.carbon;
  PROD1decay.nitrogen = PROD1.nitrogen;

  j = dyr%10;
  initPROD10[j].carbon = formPROD10.carbon;
  initPROD10[j].nitrogen = formPROD10.nitrogen;
  PROD10decay.carbon  = 0.0;
  PROD10decay.nitrogen  = 0.0;
  for ( i = 0; i < 10; i++)
  {
    yrtrash = initPROD10[i].carbon * 0.10;
    PROD10decay.carbon += yrtrash;
    yrtrash = initPROD10[i].nitrogen * 0.10;
    PROD10decay.nitrogen += yrtrash;
  }
  PROD10.carbon -= PROD10decay.carbon;
  PROD10.nitrogen -= PROD10decay.nitrogen;

  k = dyr%100;
  initPROD100[k].carbon = formPROD100.carbon;
  initPROD100[k].nitrogen = formPROD100.nitrogen;
  PROD100decay.carbon = 0.0;
  PROD100decay.nitrogen = 0.0;
  for ( i = 0; i < 100; i++)
  {
    yrtrash = initPROD100[i].carbon * 0.01;
    PROD100decay.carbon += yrtrash;
    yrtrash = initPROD100[i].nitrogen * 0.01;
    PROD100decay.nitrogen += yrtrash;
  }

  PROD100.carbon -= PROD100decay.carbon;
  PROD100.nitrogen -= PROD100decay.nitrogen;

  TOTPROD.carbon = PROD1.carbon + PROD10.carbon + PROD100.carbon;
  TOTPROD.nitrogen = PROD1.nitrogen + PROD10.nitrogen + PROD100.nitrogen;

  formTOTPROD.carbon = formPROD1.carbon + formPROD10.carbon + formPROD100.carbon;
  formTOTPROD.nitrogen = formPROD1.nitrogen + formPROD10.nitrogen
			+ formPROD100.nitrogen;

  TOTPRODdecay.carbon = PROD1decay.carbon + PROD10decay.carbon
			+ PROD100decay.carbon;
  TOTPRODdecay.nitrogen = PROD1decay.nitrogen + PROD10decay.nitrogen
  			+ PROD100decay.nitrogen;

};
