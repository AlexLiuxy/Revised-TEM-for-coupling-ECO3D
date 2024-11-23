/* **************************************************************
*****************************************************************
HATMS423.HPP - uses the physical characteristics of the atmosphere
	       as described by TEM, but uses transient CO2, specified
	       every 6 months, from an ASCII file

               20000102 - DWK added compiler directives
*****************************************************************
************************************************************** */

#if !defined(ATMS423_H)
  #include "atms423.cpp"    // MPIatms inherits Atmosphere class
#endif

#if !defined(TCO2DAT423_H)
  #include "tco2dat423.cpp" // MPIatms uses CO2data class
#endif

// MPIatms also uses the global constants CYCLE and MAXRTIME

#if !defined(HATMS423_H)
#define HATMS423_H

class MPIatms : public Atmosphere {

  public:

     MPIatms();

/* **************************************************************
		 Public Functions
************************************************************** */

     void loadfosfuel(ofstream& rflog1, const int& RTIME);
     void loadmpitCO2(ofstream& rflog1, const int& totsptime, const int& RTIME);
     void loadmpitCO2(char ifilename[80], const int& totsptime, const int& RTIME);


/* **************************************************************
		 Public Variables
************************************************************** */

     double ffuel[MAXRTIME];
     double maxffuel;
     int ffuelyear[MAXRTIME];

};

#endif
