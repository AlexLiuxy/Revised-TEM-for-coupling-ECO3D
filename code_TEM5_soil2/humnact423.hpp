/* **************************************************************
*****************************************************************
HUMNACT423.HPP - describes human disturbances to natural ecosystems

                20000102 - DWK added compiler directives
*****************************************************************
************************************************************** */

#if !defined(BIOMS423_H)
  #include "bioms423.hpp"  // Humnact uses Biomass class
#endif

// Humnact also uses the global constants CYCLE, NUMMSAC, NUMVEG
//   and MAXRTIME

#if !defined(HUMNACT423_H)
#define HUMNACT423_H

class Humanact {

  public:

  Humanact();

/* **************************************************************
		 Public Functions
************************************************************** */

     void conversion(int& ez, Tveg42& veg, Tsoil4& soil);
     void getecd (ofstream& rflog1);
     void getecd (char ecd[80]);
     void resetPROD();
     void updateyr(const int& dyr);



  /* **************************************************************
		 Public Variables
************************************************************** */

    Biomass PROD1;

    Biomass initPROD10[10];
    Biomass PROD10;

    Biomass initPROD100[100];
    Biomass PROD100;

    Biomass TOTPROD;
    double totgC[CYCLE];

    Biomass convrtflx[CYCLE];
    double yrconvrtC;
    double yrconvrtN;

    Biomass vconvrtflx[CYCLE];
    double yrvconvrtC;
    double yrvconvrtN;

    Biomass sconvrtflx[CYCLE];
    double yrsconvrtC;
    double yrsconvrtN;



    double totnpp[CYCLE];
    Biomass formPROD1;
    Biomass formPROD10;
    Biomass formPROD100;
    Biomass formTOTPROD;

    Biomass PROD1decay;
    Biomass PROD10decay;
    Biomass PROD100decay;
    Biomass TOTPRODdecay;

    Biomass slash[CYCLE];
    double yrslashC;
    double yrslashN;
    double nretent[CYCLE];
    double nvretent[CYCLE];
    double nsretent[CYCLE];
    double yrnrent;
    double yrnvrent;
    double yrnsrent;
    double potnpp[CYCLE];
    double natsoil;
    double kd;
    double tpotnpp[NUMMSAC][MAXRTIME][CYCLE];
    Biomass npp[CYCLE];
    double yrnppC;
    double yrnppN;
    Biomass ltrfal[CYCLE];
    double yrltrc;
    double yrltrn;
    double fertn[CYCLE];
    double yrfertn;

    double slashpar[NUMVEG];
    double vconvert[NUMVEG];
    double prod10par[NUMVEG];
    double prod100par[NUMVEG];
    double sconvert[NUMVEG];
    double nvretconv[NUMVEG];
    double nsretconv[NUMVEG];
    double vrespar[NUMVEG];
    double cfall[NUMVEG];
    double nfall[NUMVEG];

    double RAP;
    double tRAP[MAXRTIME];
    double c2n;

    int state;
    int tstate[MAXRTIME];
    int distflag;
    int tlulcflag;
    int lulcyear[MAXRTIME];

    // Is RAP always going to be zero?
    // (if so, we don't need potential NPP data)
    int RAP0flag;

    int prvstate;
//    char predstr[35][9];
};

#endif

