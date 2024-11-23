/* **************************************************************
*****************************************************************
ATMS423.HPP - object describes physical characteristics of the
	      atmosphere

           - modified from atms41d.hpp by DWK 20000102


2001048 - Q. Z. Soil thermal model
*****************************************************************
************************************************************** */

// Class uses global constants CYCLE and MAXRTIME

#if !defined(ATMS423_H)
#define ATMS423_H
class Atmosphere
{

  public:

     Atmosphere();

/* **************************************************************
		 Public Functions
************************************************************** */


     void precsplt(double& prec, double& tair, double& rain, double& snowfall);
     double petjh(const double& nirr, const double& tair, const int& dm);
     double xeet(const double& rain, const double& snowinf, const double& pet, const double& avlh2o, const double& awcapmm, const int& dm);



/* **************************************************************
		 Public Variables
************************************************************** */

// Gaseous components of air

     double initco2;               // initial CO2 concentration (ppmv)
     double co2level;              // constant CO2 concentration (ppmv)
     double co2[CYCLE];            // variable CO2 concentration (ppmv)
     double tco2[MAXRTIME][CYCLE]; // transient CO2 concentration (ppmv)          
     long co2year[MAXRTIME];       // year of CO2 data

 //  Light Variables

     double girr[CYCLE];        // Gross Irradiance (cal/(sq. cm * day))
     double nirr[CYCLE];        // Net Irradiance   (cal/(sq. cm * day))
     double nirr1[CYCLE];        // Net Irradiance   (cal/(sq. cm * day))
     double par[CYCLE];         // Photosynthetically Active Radiation
     double par1[CYCLE];         // Photosynthetically Active Radiation
				      // (cal/(sq.cm * day))
//   Evapotranspiration Variables

     double pet[CYCLE];         // Potential Evapotranspiration (mm)     
     double eet[CYCLE];         // Estimated Actual Evapotranspiration (mm)


//   Climatic Variables

     double tair[CYCLE];        // Surface Air Temperature (degrees C)
     double ttair[MAXRTIME][CYCLE];
     long tairyear[MAXRTIME];       // year of air temperature data

     double msoil[CYCLE];        // Surface Soil Moisture (degrees C)  XL add for soil moisture
     double tmsoil[MAXRTIME][CYCLE];
     long msoilyear[MAXRTIME];       // year of soil moisture data

     double frontd[CYCLE],thawbe[CYCLE],thawend[CYCLE]; // front frozen depth
     double tsoil[CYCLE];       // Top 20cm soil temperature (degrees C)
     double dst5[CYCLE],dst10[CYCLE],dst20[CYCLE],dst50[CYCLE],dst100[CYCLE],dst200[CYCLE]; // for different depth soil temperature

     double clds[CYCLE];        // Cloudiness (%)
     double tclds[MAXRTIME][CYCLE];
     long cldsyear[MAXRTIME];       // year of cloudiness data

 //    double sun[CYCLE];         // percent sunshine duration

     double  prec[CYCLE];       // Total Precipitation (mm)
     double  tprec[MAXRTIME][CYCLE];
     long precyear[MAXRTIME];       // year of precipitation data

     double  snowfall[CYCLE];   // Snow (mm)
     double  rain[CYCLE];       // Rainfall (mm)

     double prevtair;           // Previous Month's Air Temperature
				//   (degrees C)
     double prev2tair;          // Previous 2 Month's Air Temperature
				//   (degrees C)

     double prvpetmx;           // Maximum PET of previous year
     double prveetmx;           // Maximum EET of previous year

     double avetair;            // Mean annual air temperature (degrees C)
     double mxtair;             // Maximum monthly air temperature (degrees C)
     double mxttair[MAXRTIME];
     double yrprec;             // Annual sum of total precipitation (mm)
     double yrtprec[MAXRTIME];
     double yrrain;             // Annual sum of rainfall (mm)
     double yrpet;              // Annual sum of potential evapotranspiration (mm)
     double yreet;              // Annual sum of estimated actual evapotranspiration (mm)
     double yrsnowfall;         // Annual sum of snow

     double yrfrontd,yrthawbegin,yrthawend; // for soil thermal model
     double yrtsoil; // for soil temperature
     double yrdst5; // for soil temperature
     double yrdst10; // for soil temperature
     double yrdst20; // for soil temperature
     double yrdst50; // for soil temperature
     double yrdst100; // for soil temperature
     double yrdst200; // for soil temperature

     int tco2flag;
     int ttairflag;
     int tprecflag;
     int tcldsflag;
     int tmsoilflag; //XL add for soil moisture


 // private:      // changed private into public by Q. Z. for soil thermal model

/* **************************************************************
		      Private Variables
************************************************************** */

     double daze[CYCLE];          // number of days per month

};

#endif
