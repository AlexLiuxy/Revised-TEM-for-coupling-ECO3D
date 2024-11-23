/* **************************************************************
*****************************************************************
TMCRB423.CPP - object describing characteristics of soil microbes
               used in the Terrestrial Ecosystem Model (TEM)
*****************************************************************
************************************************************** */

#if !defined(TMCRB423_H)
  #include "tmcrb423.hpp"
#endif

/* **************************************************************
************************* Public Functions **********************
************************************************************** */

/* **************************************************************
************************************************************** */

void Tmicrobe4::getvegecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the file with microbe parameter values (.ECD)";
  cout << endl;
  cout << "dependent upon vegetation : " << endl;
//cin >> ecd;
  fpara >> ecd;

  rflog1 << "Enter name of the file with microbe parameter values (.ECD)";
  rflog1 << endl;
  rflog1 << "dependent upon vegetation: " << ecd << endl << endl;

  getvegecd(ecd);

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe4::getvegecd(char ecd[80])
{

  const int NUMVAR = 8;
  char dummy[NUMVAR][12];
  ifstream infile;
  int i;
  int dcmnt;

  int vegid[MAXCMNT];
  char vegname[MAXCMNT][31];
  long update[MAXCMNT];

  infile.open(ecd, ios::in);

  if (!infile)
  {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 1; dcmnt < MAXCMNT; dcmnt++)
  {
     infile >> vegid[dcmnt] >> vegname[dcmnt];
     infile >> rhq10[dcmnt];
     infile >> kn2[dcmnt];
     infile >> moistmin[dcmnt] >> moistopt[dcmnt] >> moistmax[dcmnt];
     infile >> update[dcmnt];
  }

  infile.close();

};

/* **************************************************************
************************************************************** */



/* **************************************************************
************************************************************** */

double Tmicrobe4::nminxclm(const int& dcmnt, const int& dm, double& soilh2o,
                           double& soilorgc, double& soilorgn, double& availn,
                           double& decay, double& rh, double& ksoil)
{

  double nmin;
  double tcnsoil;

  tcnsoil = cnsoil[dcmnt];

  nuptake[dm] = 0.0;
  nmin = 0.0;
  if (soilorgc > 0.0 && soilorgn > 0.0)
  {
    nuptake[dm]  = (availn * ksoil) / soilh2o;

//   if (nuptake[dm] > 0.00002) nuptake[dm] = 0.000018;
// debug by Q. zhuang, 2/Nov/2001
//   if (nuptake[dm] < 0.00001) nuptake[dm] = 0.000018;
// debug by Q. zhuang, 2/Nov/2001

    nuptake[dm] /= (kn2[dcmnt] + nuptake[dm]);
    nmin   = ((soilorgn / soilorgc) - (nup * nuptake[dm] * decay)) * rh;

    if (nmin >= 0.0) { nmin *= (soilorgn/soilorgc) * tcnsoil; }
    else { nmin *= (soilorgc/soilorgn) / tcnsoil; }
  }

  return nmin;

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

double Tmicrobe4::rhxclm(double& soilorgc, double& dq10, double& moist)
{

  return kd * soilorgc * moist * dq10;

//  if (soilorgc > 0.001) { return kd * soilorgc * moist * dq10; }
//  else { return 0.0; }

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Tmicrobe4::showecd(const int& dcmnt)
{

  cout << endl << "   MICROBIAL PARAMETERS INFLUENCED BY CLIMATE";
  cout << endl << endl;
  printf("          KN2 = %6.4lf\n", kn2[dcmnt]);
  printf("        RHQ10 = %6.2lf\n", rhq10[dcmnt]);
  printf("     MOISTMIN = %8.6lf\n", moistmin[dcmnt]);
  printf("     MOISTOPT = %8.6lf\n", moistopt[dcmnt]);
  printf("     MOISTMAX = %8.6lf\n", moistmax[dcmnt]);
  printf("       CNSOIL = %5.2lf\n", cnsoil[dcmnt]);

};

/* **************************************************************
************************************************************** */

/* **************************************************************
************************************************************** */

double Tmicrobe4::yrkd(int nfeed, double& yrltrc, double& yrltrn, int dcmnt)
{

  double yrkd;

  if (yrltrn <= 0.0) { yrltrn = 0.001; }
  if (yrltrc < 0.0)
  {
    cout << "YRLTRC is < 0.0 in microbe.yrkd()" << endl;
    exit(-1);
  }
  if (nfeed == 0) { yrkd = kdc; }
  else
  {
    yrkd = kdc * pow((yrltrc/yrltrn),-0.784) / pow(lcclnc[dcmnt],-0.784);
  }

  return yrkd;
};

