/* ****************************************************************
ATMS423.CPP - object describes physical characteristics of the
	      atmosphere

           - modified from atms41d.cpp by DWK 20000102
*****************************************************************
************************************************************** */

#if !defined(ATMS423_H)
 #include "atms423.hpp"
#endif

Atmosphere::Atmosphere()
{

  co2level = 0;

// Initialize number of days per month for each month of year

  daze[0] = daze[2] = daze[4] = daze[6] = daze[7] = 31.0;
  daze[9] = daze[11] = 31.0;
  daze[3] = daze[5] = daze[8] = daze[10] = 30.0;
  daze[1] = 28.0;


};

/* **************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

double Atmosphere::petjh(const double& nirr, const double& tair, const int& dm)
{

  double f;
  double rt;
  double pet;

  f = ((9.0/5.0) * tair) + 32.0;
  rt = nirr * 0.016742;
  pet = ((0.014*f) - 0.37) * rt * daze[dm];

  if (pet < 0.0) { pet = 0.0; }

  return pet;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Atmosphere::precsplt(double& prec, double& tair, double& rain, double& snowfall)
{


/* *************************************************************
	Willmott's assumptions on snow/rain split:
************************************************************** */

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

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

double Atmosphere::xeet(const double& rain, const double& snowinf, const double& pet, const double& avlh2o, const double& awcapmm, const int& dm)
{

  const double edpar = 5.0;
  double gm;
  double ep;
  double def;
  double prob;
  double rbar;
  double dsm;
  double aet;

  if ((rain+snowinf) >= pet)
  {
    aet = pet;
  }
  else
  {
    gm = (1.0 - exp(-edpar * avlh2o/awcapmm)) / (1.0 - exp(-edpar));
    ep = pet / daze[dm];
    def = ep + awcapmm - avlh2o;
    prob = 1.0 - exp(-0.005*(rain + snowinf));
    if (prob != 0.0) { rbar = (rain + snowinf) / (daze[dm] * prob); }
    else { rbar = 0.0; }

    if (rbar != 0.0)
    {
      dsm = rbar*prob*(gm + ((1.0-gm) * exp(-ep/rbar)) - exp(-def/rbar)) - (ep*gm);
    }
    else {
      dsm = -ep*gm;
    }

    dsm *= daze[dm];

    aet = rain + snowinf - dsm;
    if (aet > pet) { aet = pet; }

  }

  return aet;
};


