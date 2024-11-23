/* **************************************************************
*****************************************************************
TVEG423.CPP - object describing characteristics of vegetation used
	      in the Terrestrial Ecosystem Model (TEM)
            - Bug fix added to hpp file by DWK on 19991028
            - modified from TVEG42.CPP by DWK on 20000102
20000614 - DWK changes innpp to innpp[CYCLE] and ingpp to ingpp[CYCLE] in
           tveg423e.hpp

2001048 -- Q. Z. modify GPP with freezing index
*****************************************************************
************************************************************** */

#if !defined(TVEG423_H)
  #include "tveg423e.hpp"
#endif

/* *************************************************************
************************************************************* */

Tveg42::Tveg42() : Biome()
{

  leafyrs = 10;
  topt = -999.9;

};

/* **************************************************************
************************* Public Functions **********************
************************************************************** */


/* *************************************************************
************************************************************* */

double Tveg42::deltaleaf(const int& dcmnt, double& eet,
                         double& prveetmx, double& prvleaf)
{

  double normeet;
  double unnormleaf;

  if (prveetmx <= 0.0) { prveetmx = 1.0; }
  normeet = eet / prveetmx;
  unnormleaf = (aleaf[dcmnt] * normeet) + (bleaf[dcmnt] * prvleaf)
               + cleaf[dcmnt];
  if (unnormleaf < (0.5 * minleaf[dcmnt]))
  {
    unnormleaf = 0.5 * minleaf[dcmnt];
  }

  return unnormleaf;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tveg42::getecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the file with the vegetation "
	   << "parameter values (.ECD):";
  cout << endl;
//cin >> ecd;
  fpara >> ecd;

  rflog1 << "Enter name of the file with the vegetation parameter "
		 << "values (.ECD):";
  rflog1 << ecd << endl << endl;

  getecd(ecd);

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tveg42::getecd(char ecd[80])
{

  const int NUMVAR = 21;
  char dummy[NUMVAR][10];
  ifstream infile;
  int i;
  int dcmnt;

  int vegid[MAXCMNT];
  char vegname[MAXCMNT][31];
  long update[MAXCMNT];

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << endl << "Cannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> vegid[dcmnt] >> vegname[dcmnt];
    infile >> kc[dcmnt];
    infile >> ki[dcmnt];
    infile >> gva[dcmnt];

    infile >> tmin[dcmnt];
    infile >> toptmin[dcmnt];
    infile >> toptmax[dcmnt];
    infile >> tmax[dcmnt];

    infile >> raq10a0[dcmnt];
    infile >> raq10a1[dcmnt];
    infile >> raq10a2[dcmnt];
    infile >> raq10a3[dcmnt];

    infile >> kn1[dcmnt];
    infile >> labncon[dcmnt];

    infile >> leafmxc[dcmnt];
    infile >> kleafc[dcmnt];
    infile >> sla[dcmnt];
    infile >> cov[dcmnt];
    infile >> fpcmax[dcmnt];

    infile >> update[dcmnt];
  }
  infile.close();

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tveg42::getleafecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the file with leaf parameter values (.ECD):" << endl;
//cin >> ecd;
  fpara >> ecd;

  rflog1 << "Enter name of the file with leaf parameter values (.ECD):";
  rflog1 << ecd << endl << endl;

  getleafecd(ecd);

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

void Tveg42::getleafecd(char ecd[80])
{

  const int NUMVAR = 7;
  char dummy[NUMVAR][8];
  ifstream infile;
  int i;
  int dcmnt;

  int lfvegid[MAXCMNT];
  char lfvegname[MAXCMNT][31];
  int update[MAXCMNT];

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> lfvegid[dcmnt] >> lfvegname[dcmnt];
    infile >> minleaf[dcmnt];
    infile >> aleaf[dcmnt];
    infile >> bleaf[dcmnt];
    infile >> cleaf[dcmnt];
    infile >> update[dcmnt];
  }

  infile.close();

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

// modified for thawing
double Tveg42::gppxclm(int& dcmnt, double& co2, double& par,
                       double& temp, double& gv,
                       double& leaf, double& foliage,double& thawpercent)

//double Tveg42::gppxclm(int& dcmnt, double& co2, double& par,
//                       double& temp, double& gv,
//                       double& leaf, double& foliage)
{

  double gpp;


/* **************************************************************
   gpp:    gpp as influenced by carbon dioxide (co2), moisture
	   (gv), phenology (leaf), photosynthetically active
	   radiation (par), and air temperature (temp)
************************************************************** */

  gpp  = co2 * gv;
  gpp *= cmax * foliage / (kc[dcmnt] + gpp);
  gpp *= leaf * par / (ki[dcmnt] + par);
  gpp *= temp;

  gpp *= thawpercent; // addition by Q.Z.

  return gpp;

};

/* **************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void Tveg42::leafinit(ofstream& rflog1)
{

  int lfswtch;

  cout << "Enter number of years for model run:  " << endl;
  cin >> leafyrs;

  cout << "Do you have a file containing the phenology parameters?:";
  cout << endl;
  cout << "Enter 0 for no:" << endl;
  cout << "Enter 1 for yes:" << endl;
//cin >> lfswtch;
  fpara >> lfswtch;

  rflog1 << "Enter number of years for model run:  " << leafyrs << endl;
  rflog1 << "Do you have a file containing the phenology parameters?:";
  rflog1 << endl;
  rflog1 << "Enter 0 for no:" << endl;
  rflog1 << "Enter 1 for yes:" << endl;
  rflog1 << lfswtch << endl;

  if (lfswtch == 0 )
  {
    cout << "Enter regression coefficient for relative EET (i.e. 'a'):  ";
//  cin >> aleaf[0];
	fpara >> aleaf[0];

    cout << "Enter regression coefficient for LAI(t-1) (i.e. 'b'):  ";
//  cin >> bleaf[0];
	fpara >> bleaf[0];

	cout << "Enter regression intercept (i.e. 'c'):  ";
//  cin >> cleaf[0];
	fpara >> cleaf[0];

    rflog1 << "Enter regression coefficient for relative EET (i.e. 'a'): ";
    rflog1 <<  aleaf[0] << endl;
    rflog1 << "Enter regression coefficient for LAI(t-1) (i.e. 'b'):  ";
    rflog1 << bleaf[0] << endl;
    rflog1 << "Enter regression intercept (i.e. 'c'):  ";
    rflog1 << cleaf[0] << endl;

    minleaf[0] = 2.0;
    while (minleaf[0] >= 1.0)
    {
      cout << "Enter minimum LAI (must be less than 1.0):  ";
      cin >> minleaf[0];

      rflog1 << "Enter minimum LAI (must be less than 1.0):  ";
      rflog1 << minleaf[0];
      rflog1 << endl << endl;
    }
  }
  else { getleafecd(rflog1); }

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

double Tveg42::nupxclm(int& dcmnt, double& soilh2o, double& availn,
                       double& respq10, double& ksoil, double& foliage)
{

/* **************************************************************
   nuptake:  uptake of nitrogen by plants as influenced by
	     available nitrogen concentration (availn), moisture
	     (ksoil), and air temperature (respq10)
************************************************************** */

  double nuptake;

   nuptake  = (availn * ksoil) / soilh2o;


//   if (nuptake > 0.00002) nuptake = 0.000018;  // debug by Q. zhuang, 2/Nov/2001
//   if (nuptake < 0.00001) nuptake = 0.000018;  // debug by Q. zhuang, 2/Nov/2001

//   printf(" %f",nuptake);
   
   nuptake *= nmax * foliage / (kn1[dcmnt] + nuptake);
   nuptake *= respq10;

  return nuptake;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

double Tveg42::rmxclm(int& dcmnt, double& vegc, double& respq10)
{

/* **************************************************************
   rm: plant maintenance respiration as influenced by plant
       biomass (vegc) and air temperature (respq10)
************************************************************** */

  double rmaint;

  kr = exp((kra[dcmnt]*vegc) + krb[dcmnt]);
  rmaint  = kr * vegc;
  rmaint *= respq10;

  return rmaint;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tveg42::showecd(int& dcmnt)
{

  cout << endl << "             VEGETATION PARAMETERS INFLUENCED BY CLIMATE";
  cout << endl << endl;
  printf("     KI = %6.2lf   KC = %6.2lf   KN1 = %6.4lf   GVA = %8.4lf\n",
         ki[dcmnt], kc[dcmnt], kn1[dcmnt], gva[dcmnt]);

  printf("   TMIN = %5.1lf    TOPTMIN = %5.1lf   TOPTMAX = %5.1lf   \
		TMAX = %5.1lf\n",
         tmin[dcmnt], toptmin[dcmnt], toptmax[dcmnt], tmax[dcmnt]);

  printf(" RAQ10A0 = %7.5lf  RAQ10A1 = %7.5lf  RAQ10A2 = %7.5lf  \
  		RAQ10A3 = %7.5lf\n",
         raq10a0[dcmnt], raq10a1[dcmnt], raq10a2[dcmnt], raq10a3[dcmnt]);

};

/* **************************************************************
************************************************************** */

void Tveg42::showleaf(int& dcmnt)
{

  cout << endl << "         PARAMETERS FOR THE LEAF PHENOLOGY MODEL";
  cout << endl << endl;
  printf("     ALEAF = %7.5lf     BLEAFC = %7.5lf         CLEAF = %8.5lf\n",
         aleaf[dcmnt], bleaf[dcmnt], cleaf[dcmnt]);
  printf("   MINLEAF = %4.2lf       MAXLEAF = %7.4lf      UNLEAF12 = %7.4lf\n",
         minleaf[dcmnt], prvleafmx[dcmnt], unleaf12[dcmnt]);

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tveg42::updateC2N(const int& dcmnt, double& yreet, double& yrpet,
                       double& currentco2, double& initco2)
{
  if (yrpet != 0.0)
  {
    c2n = c2nb[dcmnt] + c2na[dcmnt]*(yreet/yrpet);
  }
  else { c2n = c2nb[dcmnt]; }

  if (c2n < c2nmin[dcmnt]) { c2n = c2nmin[dcmnt]; }
  adjc2n = 1.0 + (dc2n * (currentco2 - initco2));
  c2n *= adjc2n;
  cneven[dcmnt] = initcneven[dcmnt] * adjc2n;

};

/* **************************************************************
************************* Private Functions *********************
************************************************************** */

double Tveg42::rq10(int& dcmnt, double& tair)
{

  double raq10;

/* **************************************************************
 rq10: effect of temperature on plant respiration
************************************************************** */

  raq10 = raq10a0[dcmnt] + (raq10a1[dcmnt]*tair)
          + (raq10a2[dcmnt]*pow(tair,2.0))
          + (raq10a3[dcmnt]*pow(tair,3.0));

  return pow(raq10,tair/10.0);

};

