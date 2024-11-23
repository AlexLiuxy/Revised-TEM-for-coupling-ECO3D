/* **************************************************************
*****************************************************************
POTCLM42.CPP - describes physical characteristics of the atmosphere
             - modified by DWK on 20000102
20000616 - DWK adds the enum variable sradkey to declaration in
           potsclm423e.hpp
*****************************************************************
************************************************************** */

#if !defined(POTCLM423_H)
  #include "potclm423e.hpp"
#endif

/* **************************************************************
************************************************************** */

Potsclm::Potsclm() : Atmosphere()
{

// Initialize predstr array

  strcpy(predstr[0],"GIRR");
  strcpy(predstr[1],"NIRR");
  strcpy(predstr[2],"PAR");
  strcpy(predstr[3],"CLDINESS");

};

/* **************************************************************
		     Public Functions
************************************************************** */

double Potsclm::mkclds(const double& girr, const double& nirr)
{

  double clouds;

  if (nirr >= (0.71 * girr)) { return clouds = 0.0; }
  else
  {
    clouds = 1.0 - (((nirr/girr) - 0.23)/0.48);
    clouds *= 100.0;
  }
  if (clouds > 100.0) { clouds = 100.0; }

  return clouds;

};


/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

double Potsclm::xgirr(const double& lat, const int& dm, double& sumday)
{

  const double pi = 3.141592654;                // Greek "pi"
  const double sp = 1368.0 * 3600.0 / 41860.0;  // solar constant

  double lambda;
  double sumd;
  double sig;
  double eta;
  double sinbeta;
  double sb;
  double sotd;
  int day;
  int hour;
  double gross;

  ndays[0] = 31;
  ndays[1] = 28;
  ndays[2] = 31;
  ndays[3] = 30;
  ndays[4] = 31;
  ndays[5] = 30;
  ndays[6] = 31;
  ndays[7] = 31;
  ndays[8] = 30;
  ndays[9] = 31;
  ndays[10] = 30;
  ndays[11] = 31;

  lambda = lat * pi / 180.0;

  gross = 0.0;
  for (day = 0; day < ndays[dm]; day++)
  {
    ++sumday;
    sumd = 0;
    sig = -23.4856*cos(2 * pi * (sumday + 10.0)/365.25);
    sig *= pi / 180.0;

    for (hour = 0; hour < 24; hour++)
    {
      eta = (double) ((hour+1) - 12) * pi / 12.0;
      sinbeta = sin(lambda)*sin(sig) + cos(lambda)*cos(sig)*cos(eta);
      sotd = 1 - (0.016729 * cos(0.9856 * (sumday - 4.0) * pi / 180.0));
      sb = sp * sinbeta / pow(sotd,2.0);
      if (sb >= 0.0) { sumd += sb; }
    }

    gross += sumd;
  }

  gross /= (double) ndays[dm];

  return gross;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Potsclm::xnirr(const double& clds, const double& girr)
{

  double nirr;

  if (clds >= 0.0)
  {
    nirr = girr * (0.251 + (0.509*(1.0 - clds/100.0)));
  }
  else { nirr = MISSING; }

  return nirr;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double Potsclm::xpar(const double& clds, const double& nirr)
{

  double par;

  if (clds >= 0.0)
  {
      par = nirr * ((0.2 * clds / 100.0) + 0.45);
  }
  else { par = MISSING; }

  return par;

};


