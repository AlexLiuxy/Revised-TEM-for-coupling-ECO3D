/* **************************************************************
*****************************************************************
TELM423.HPP - Runs TEM for a single grid cell

Modifications:

20000202 - DWK changed initstat[MAXSTATE] to initstat[MAXSTATE+1];
20000219 - DWK changed ez to cmnt in temgisqc()
20000614 - DWK adds innpp[] and innpp[] to temwritepred() to allow these
           variables to be written out to spatially explicit files.  NOTE:
           THE EFFECT OF THESE CHANGES HAVE NOT BEEN EVALUATED YET.
20000614 - DWK resorts function names alphabetically to allow function
           descriptions to be found easier
20000616 - DWK changes the function of atmswritepred() in telm423e.cpp to write
           out solar radiation output variables directly to spatially explicit
           data sets instead of a 2-dimensional array formerly known as
           atmspred[][].  NOTE: THIS FUNCTION HAS NOT BEEN EVALUATED YET -
           PLACEHOLDER.
20000616 - DWK adds atmswritemiss() to write out solar radiation output variables
           for grid cells that have missing data.  NOTE: THIS FUNCTION HAS NOT
           BEEN EVALUATED YET - PLACEHOLDER.
20000616 - DWK changes the function of temwritepred() in telm423e.cpp to write
           out TEM output variables directly to spatially explicit data sets
           instead of a 3-dimensional array formerly known as tempred[][][].
20000616 - DWK adds temwritemiss() to write out TEM output for grid cells that
           either have missing data or have environmental conditions that will
           not allow biological organisms to survive (i.e., ling-term MAXTAIR <
           -1.0 degree Celsius or long-term mean annual precipitation equal to
           0.0).
20000629 - DWK changed dyr to outyr in call to temwritepred() in the "equilibrium
           portion" of runtem() in telm423e.cpp
20000718 - DWK modified runtem() in telm423e.cpp so that grid cells with missing
           data would still report output to spatially explicit files (i.e.,
           bug fix)
20000719 - DWK adds the global constant TQCZEROFLAG to flag when environmental
           conditions of a grid cell cannot support life.  In addition, the
           appropriate conditional statements related to tqc in
           runtem() were modified in telm423e.cpp to use this global constant.
20000719 - DWK modified the assignment of tem.veg.cmnt in runtem() to equal
           zero if tem.veg.temveg was less than zero or greater than (NUMVEG+1)
           (i.e., grid cell is covered with water) in telm423e.cpp
20000719 - DWK added a conditional statement to the transient portion of
           runtem() that if "tqc = ACCEPT" and "qc = REJECT", then
           "tqc = REJECT" in telm423e.cpp
20000719 - DWK fixed bug in temwritemiss() by switching the positions of
           ntempred and natmspred in the function call in telm423e.cpp.
20000719 - DWK assigns "ACCEPT" to tqc at the beginning of runtem() in
           telm423e.cpp
20000719 - DWK modifies the assignment of tqc from transqc() so that this
           assignment occurs only if tqc is currently equal ACCEPT in runtem()
           (i.e., it did not get changed because of unsuitable long-term
           environmental conditions in the grid cell in telm423e.cpp
20070306 J. Tang added ispinout switch to runtem           

*****************************************************************
************************************************************** */

//Modules representing climate and TEM

#if !defined(POTCLM423_H)
  #include "potclm423e.cpp"       // TEMelmnt uses the Potsclm class
#endif

#if !defined(TTEM423E1_H)
  #include "ttem423e1.cpp"      // TEMelmnt uses the TTEM class
#endif



//Modules describing the interface with spatially explicit data sets

#if !defined(TCLMDAT423_H)
  #include "tclmdat423.cpp"  //TEMelmnt uses the Clmdata class
#endif

#if !defined(TVEGDAT423_H)
  #include "tvegdat423.cpp"  //TEMelmnt uses the Vegdata class
#endif

#if !defined(TSOLDAT423_H)
  #include "tsoldat423.cpp"  //TEMelmnt uses the Soildata class
#endif

#if !defined(TELVDAT423_H)
  #include "telvdat423.cpp"  //TEMelmnt uses the Elevdata class
#endif

#if !defined(LULCDAT423_H)
  #include "lulcdat423.cpp"  //TEMelmnt uses the Lulcdata class
#endif

#if !defined(TTEMDAT423_H)
  #include "ttemdat423.cpp"  //TEMelmnt uses the Temdata class
#endif

#if !defined(TKDDAT423_H)
  #include "tkddat423.cpp"   //TEMelmnt uses the KDdata class
#endif

#if !defined(LATDAT423_H)
  #include "latdat423.cpp"   //TEMelmnt uses the Latdata class
#endif

#ifndef TELM423E1_H
#define TELM423E1_H

// Following global constant added by DWK on 20000719 to flag
// when long-term environmental conditions of a grid cell cannot
// support life
const int TQCZEROFLAG = 31;

class TEMelmnt {

  public:

     TEMelmnt();

/* *************************************************************
		 Public Function Declarations
************************************************************* */

     int atmsgisin(FILE* fclds, const int& cldflag, FILE* flonlat,
                   const int& lonlatflag, const int& numspin,
                   const int& spintime, const int& RTIME);
     void atmswritemiss(ofstream fout[NUMATMS], char predname[MAXPRED][9],
                        const int& dyr, const int& natmspred,
                        const double value);
     void atmswritepred(ofstream fout[NUMATMS], char predname[MAXPRED][9],
                        const int& dyr, const int& natmspred);
     void GISsetELMNTstate(const int& ez, Temdata initstat[NUMMSAC][MAXSTATE+1]);
     void runtem(ofstream& rflog1, char predmap[MAXPRED][9],
                 const int& cldflag, const int& atmsflag,
                 const int& atmsoutfg, int& natmspred,
                 ofstream fsradout[NUMATMS],
                 const int& temflag, const int& kdinflg,
                 int& ntempred, ofstream ftemout[MAXPRED],
                 const int& stateflag,
                 const int& equil, const int& totsptime,
                 const int& RTIME, const int& ispinout);
     int temgisin(ofstream& flog1, int& ftlerr, const int& atmsflag,
		  const int& itype, const int& maxtype, const int& kdinflg,
                  const int& stateflag,
		  FILE* fstxt, FILE*felev, FILE* fnirr, FILE* fpar,
		  FILE* ftair, FILE* fprec, FILE* fmsoil, FILE* flulc, FILE* fnpp,
                  FILE* fkdin, FILE* fstate[MAXSTATE+1],
		  const int& numspin, const int& spintime, const int& RTIME); //XL add for msoil
        void temwritemiss(ofstream fout[MAXPRED], char predname[MAXPRED][9],
                            const int& dyr, const int& itype,
                            const int& natmspred, const int& ntempred,
                            const double value);
     void temwritepred(ofstream fout[MAXPRED], char predname[MAXPRED][9],
                       const int& dyr, const int& itype,
                       const int& ntempred, const int& natmspred);
     int veggisin(ofstream& flog1, int& ftlerr, const int& atmsflag,
		  const int& temflag, FILE* ftveg);

     float col;
     float row;
     int lonlatflag;
     double lon;
     double lat;
     int carea;
     char contnent[9];

     int mez;
     int maxtype;

     int outyr;
     int itype;

     long atmstotyr[MAXRTIME];
     long ttotyr[MAXRTIME][NUMMSAC];

     Potsclm clm;
     TTEM tem;


/* *************************************************************
		 Private Function Declarations
************************************************************* */

  private:

     int coregerr(ofstream& rflog1, char varname1[9], const float& col1,
		  const float& row1, char varname2[9], const float& col2,
		  const float& row2);
     int loadteclds(FILE* fclds, Clmdata clds, const int& numspin,
                    const int& spintime, const int& RTIME);
     int loadtelulc(FILE* flulc, Lulcdata lulc, const int& numspin,
                    const int& spintime, const int& RTIME, int& ftlerr,
		    ofstream& flog1);
     int loadtenpp(FILE* fnpp, Temdata npp, const int& maxtype,
                   const int& RTIME, int& ftlerr, ofstream& flog1);
     int loadteprec(FILE* fprec, Clmdata prec, const int& numspin,
                    const int& spintime, const int& RTIME, int& ftlerr,
		    ofstream& flog1);
     int loadtemsoil(FILE* fmsoil, Clmdata msoil, const int& numspin,
                    const int& spintime, const int& RTIME, int& ftlerr,
		    ofstream& flog1);   //XL add for msoil
     int loadtetair(FILE* ftair, Clmdata tair, const int& numspin,
                    const int& spintime, const int& RTIME, int& ftlerr,
		    ofstream& flog1);
     int temgisqc(const int& stateflag, const double& pctsilt,
                  const double& pctclay, const int& cmnt, const double& elev,
                  double nirr[CYCLE], double par[CYCLE],
		  double tair[CYCLE], double& mxtair,
                  double prec[CYCLE],double msoil[CYCLE], double& yrprec,
		  Temdata initstat[NUMMSAC][MAXSTATE+1]);//XL add for msoil
     int transqc(int& maxyears, long& totyr, Biomass plant[CYCLE]);


// *************************************************************

     Temdata initstat[NUMMSAC][MAXSTATE+1];

     int qc;
     int lowtair;
     int noprec;

};

#endif

