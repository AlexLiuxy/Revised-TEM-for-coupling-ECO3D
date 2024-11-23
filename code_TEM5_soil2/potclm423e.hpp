/* **************************************************************
*****************************************************************
POTCLM423.HPP - describes physical characteristics of the atmosphere
             - modified by DWK on 20000102
20000616 - DWK adds the enum variable sradkey to declaration in
           potsclm423e.hpp
*****************************************************************
************************************************************** */

// Global constants

const int NUMATMS = 4;

#if !defined(ATMS423_H)
  #include "atms423.cpp"    // Potsclm inherits Atmosphere class
#endif

// Potsclm also uses the global constants CYCLE, MISSING and NUMATMS

#if !defined(POTCLM423_H)
#define POTCLM423_H

class Potsclm : public Atmosphere {

  public:

     Potsclm();

     enum sradkey { I_GIRR, I_NIRR, I_PAR, I_CLDS };

/* **************************************************************
		 Public Functions
************************************************************** */
	  
     double mkclds(const double& girr, const double& nirr);
     double xgirr(const double& lat, const int& dm, double& sumday);
     double xnirr(const double& clds, const double& girr);
     double xpar(const double& clds, const double& nirr);

/* **************************************************************
		 Public Variables
************************************************************** */
	  
     char predstr[NUMATMS][9];
     double yrsumday;
     

  private:

/* **************************************************************
		 Private Variables
************************************************************** */

     int ndays[CYCLE];                 // number of days per month
};

#endif
