/* *************************************************************
****************************************************************
PCTEM423.HPP - Class adds DOS-specific interface to the core
                 Terrestrial Ecosystem Model Version 4.2 to allow
     	         calibration of the model on a PC

Modifications:

2000103 - DWK added compiler directives
20000615 - DWK changes y[2] to y[I_SOLC] and y[3] to y[I_SOLN]
           at bottom of pcstepyr() in pctem423e.cpp
***************************************************************
************************************************************** */

// Global constants

const int WSY = 5;
const int ESY = 3;

#if !defined(TTEM423_H)
  #include "ttem423e1.cpp"      // PCTEM inherits TTEM class
#endif

#if !defined(PCTEM423_H)
#define PCTEM423_H

class PCTTEM : public TTEM {

  public:

    PCTTEM();

    enum seykey { NOEKEY,     GET_GPP,   GET_GPR,    GET_LTRC, GET_RH,
                  GET_LTRN,   GET_NMIN,  GET_NLST,   GET_NINP, GET_NEP,
                  GET_INGPP,  GET_INNPP, GET_INNUP,  GET_VEGC, GET_STRN,
                  GET_AGNPPC, GET_NMBL,  GET_NRSRB,  GET_NPP,  GET_VNUP,
                  GET_VSUP,   GET_VLUP,  GET_POTNPP, GET_L2SN, GET_FPC,
                  GET_LAI,    GET_LEAF,  GET_AGNPPN, GET_AGFRTN };

    enum swykey { NOWKEY,     GET_RAIN,  GET_RPERC, GET_RRUN, GET_SNWFALL,
                  GET_SNWINF, GET_SPERC, GET_SRUN,  GET_PET,  GET_EET,
                  GET_SH2O,   GET_PCTP,  GET_VSM,   GET_WYLD };

/* **************************************************************
			Public Functions
************************************************************** */

    inline seykey& next(seykey& s) { return s = (GET_AGFRTN== s) ? GET_GPP : seykey(s+1); }
    inline seykey& prev(seykey& s) { return s = (GET_GPP == s) ? GET_AGFRTN : seykey(s-1); }

    inline swykey& next(swykey& s) { return s = (GET_WYLD== s) ? GET_RAIN : swykey(s+1); }
    inline swykey& prev(swykey& s) { return s = (GET_RAIN == s) ? GET_WYLD : swykey(s-1); }

    void displayOptionalEflx(seykey& s);
    void displayOptionalWflx(swykey& s);
    void displayState(const int& dyr, const int& dm, char* Climit, double pstate[]);

    double getOptionalEflx(const int& dm, const int& optflx, double pstate[]);
    double getOptionalWflx(const int& optflx, double pstate[]);

    int pcadapt(const int& numeq, double pstate[], double& ptol, const int& pdm);

    void pcdisplayClm(char* calclm, double clouds[CYCLE]);
    void pcdisplayInitState(char* ecd, const double& col, const double& row);
    void pcdisplayLimits();
    void pcdisplayMonth(const int& dyr, const int& dm, double pstate[]);
    void pcdisplayOtherECD();
    void pcdisplayPAR(char* calclm, double girr[CYCLE]);
    void pcdisplayStext();
    void pcdisplayVegECD();
    void pcdisplayYrAgTEM(const int& manflag);
    void pcdisplayYrTEM(const int& manflag);
    void pcdisplayYrWBM();

    void pcinitRun(ofstream& rflog1);
    int pcstepyr(const int& dyr, int& intflag, double& tol);


/* **************************************************************
			 Public Variables
************************************************************** */


    seykey sey[ESY];
    swykey swy[WSY];

    int endgrid;
    int adapttol;
    int intbomb;
    int tolbomb;
    int mintflag;
    int topwind;
    int calwind;
    int firstcal;

// Names of files containing parameter values determined from the
//   literature

    char soilfile[40];
    char rootfile[40];
    char vegfile[40];
    char leaffile[40];
    char mcrvfile[40];
    char agfile[40];

    ofstream flog2;

  private:

/* **************************************************************
			Private Functions
************************************************************** */

    void pcdisplayDT(const double& tottime, const double& deltat);
    void pcdisplayODEerr(const int& test, double pstate[]);

    double outflux1;
    double outflux2;
    double outflux3;
    double outflux4;
    double outflux5;
    double outflux6;
    double outflux7;
    double outflux8;

    double eqsolc;
};

#endif
