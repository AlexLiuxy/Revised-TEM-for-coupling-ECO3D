/* **************************************************************
*****************************************************************
TSOIL423.CPP - object describing general characteristics of soil
            - modified by DWK on 20000102
*****************************************************************
************************************************************** */

#if !defined(TSOIL423_H)
  #include "tsoil423.hpp"
#endif

/* **************************************************************
************************************************************** */

Tsoil4::Tsoil4(void)
{

  text  = -99;
  wsoil = -99;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tsoil4::getecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the soil (.ECD) data file with parameter values: ";
  cout << endl;
//cin >> ecd;
  fpara >> ecd;

  rflog1 << "Enter name of the soil (.ECD) data file with parameter values: ";
  rflog1 << ecd << endl;

  getecd(ecd);

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Tsoil4::getecd (char ecd[80])
{

  char dummy[12];
  ifstream infile;

  long update;

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  infile >> dummy >> dummy >> dummy;
  infile >> dummy >> pctpora >> update;
  infile >> dummy >> pctporb >> update;
  infile >> dummy >> fldcapa >> update;
  infile >> dummy >> fldcapb >> update;
  infile >> dummy >> wiltpta >> update;
  infile >> dummy >> wiltptb >> update;

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::getrootz(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the data file containing the rooting depths:";
  cout << endl;
  cout << "               (e.g., ROOTZVEG.ECD)" << endl;
//cin >> ecd;
  fpara >> ecd;

  rflog1 << "Enter name of the data file containing the rooting depths:";
  rflog1 << endl;
  rflog1 << "               (e.g., ROOTZVEG.ECD)" << endl;
  rflog1 << ecd << endl;

  getrootz(ecd);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::getrootz(char ecd[80])
{

  const int NUMVAR = 7;
  char dummy[NUMVAR][10];
  ifstream infile;

  int i;
  int dcmnt;
  int  rootveg[MAXCMNT];
  long update[MAXCMNT];
  char vegname[MAXCMNT][31];

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> rootveg[dcmnt] >> vegname[dcmnt];
    infile >> rootza[dcmnt] >> rootzb[dcmnt] >> rootzc[dcmnt];
    infile >> minrootz[dcmnt] >> update[dcmnt];
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::lake(double& tair,double& prec,double& rain,double& snowfall,
       		  double& pet, double& eet, int& dm)
{

  rgrndh2o[dm] = 0.0;
  sperc[dm] = 0.0;
  snowpack[dm] = 0.0;
  sgrndh2o[dm] = 0.0;
  moist[dm] = 0.0;

  if (tair >= -1.0)
  {
   rain = prec;
    snowfall = 0.0;
  }
  else
  {
    rain = 0.0;
    snowfall = prec;
  }

  eet = pet;
  h2oyld[dm] = prec - pet;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Tsoil4::percol(double& rain, double& snowinf, double& eet,
                    double& avlh2o, const int& dm)
{

  double extra;
  double recharge;
  sperc[dm] = 0.0;
  rperc[dm] = 0.0;

  recharge = rain + snowinf;
  if (recharge <= 0.0) { recharge = 0.001; }
  if ((avlh2o + rain + snowinf - eet) > awcapmm)
  {
    extra = rain + snowinf + avlh2o - awcapmm - eet;
    sperc[dm] = snowinf * extra / recharge;
    rperc [dm] = rain * extra / recharge;
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

double Tsoil4::rrunoff(const double& rgrndh2o, const double& rperc)
{

  double rrunof;

  rrunof = 0.5 * (rgrndh2o + rperc);

  return rrunof;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void Tsoil4::showecd(void)
{

  cout << endl << "                   SOIL CHARACTERISTICS OF SITE";
  cout << endl << endl;
  printf("PSAND    = %5.2lf      PSILT = %5.2lf      PCLAY = %5.2lf\n",
         pctsand, pctsilt, pctclay);

  printf("POROSITY = %5.2lf   PCFLDCAP = %5.2lf   PCWILTPT = %5.2lf\n",
         pctpor, pcfldcap, pcwiltpt);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

double Tsoil4::snowmelt(const double& elev, const double& tair,
                        const double& prevtair, const double& snowpack)
{

  double snowflux = 0.0;

  if (tair >= -1.0)
  {
    if (elev <= 500.0) { snowflux = snowpack;}
    else
    {
      if (prevtair < -1) { snowflux = 0.5 * snowpack; }
      else { snowflux = snowpack; }
    }
  }

  return snowflux;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

double Tsoil4::srunoff(const double& elev, const double& tair,
                       const double& prevtair, const double& prev2tair,
                       const double& sgrndh2o, const double& sperc)
{

  double srunof = 0.0;

  if (tair >= -1.0)
  {
    if (prevtair < -1.0) { srunof = 0.1 * (sgrndh2o + sperc); }
    else
    {
      if (prev2tair < -1)
      {
	if (elev <= 500.0) { srunof = 0.5 * (sgrndh2o + sperc); }
	else { srunof = 0.25 * (sgrndh2o + sperc); }
      }
      else { srunof = 0.5 * (sgrndh2o + sperc); }
    }
  }

  return srunof;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void Tsoil4::xtext(int& cmnt, double& pctsilt, double& pctclay)
{

  totpor = fldcap = wiltpt = MISSING;
  awcapmm =  MISSING;

  psiplusc = (pctsilt + pctclay) * 0.01;
  if (psiplusc < 0.01) { psiplusc = 0.01; }

  rootz = (rootza[cmnt] * pow(psiplusc, 2.0)) + (rootzb[cmnt] * psiplusc)
          + rootzc[cmnt];
  if (rootz < minrootz[cmnt]) { rootz = minrootz[cmnt]; }

  pctpor = (pctpora * psiplusc) + pctporb;
  pcfldcap = (fldcapa * psiplusc) + fldcapb;
  pcwiltpt = (wiltpta * psiplusc) + wiltptb;

  totpor  = rootz * pctpor * 10.0;
  fldcap  = rootz * pcfldcap * 10.0;
  wiltpt  = rootz * pcwiltpt * 10.0;

  awcapmm = fldcap - wiltpt;

};

