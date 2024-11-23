/* **************************************************************
*****************************************************************
TELM423.CPP - Runs TEM for a single grid cell

Modifications:

20000214 - DWK changes initstat[itype][0].ave to
           initstat[itype][tem.I_VEGC].ave in runtem()
20000214 - DWK changes tem.ez to tem.veg.cmnt int runtem()
20000214 - DWK added lowtair if condition to temgisin()
20000219 - DWK set mez = 0 when mez < 0 or mez > NUMVEG in veggisin()
20000219 - DWK changed ez to cmnt in temgisqc()
20000614 - DWK adds innpp[] and innpp[] to temwritepred() to allow these
           variables to be written out to spatially explicit files.  NOTE:
           THE EFFECT OF THESE CHANGES HAS NOT BEEN EVLAUATED YET.
20000614 - DWK resorts function names alphabetically to allow function
           descriptions to be found easier
20000616 - DWK changes the function of atmswritepred() to write
           out solar radiation output variables directly to spatially explicit
           data sets instead of a 2-dimensional array formerly known as
           atmspred[][].  NOTE: THIS FUNCTION HAS NOT BEEN EVALUATED YET -
           PLACEHOLDER.
20000616 - DWK adds atmswritemiss() to write out solar radiation output variables
           for grid cells that have missing data.  NOTE: THIS FUNCTION HAS NOT
           BEEN EVALUATED YET - PLACEHOLDER.
20000616 - DWK changes the function of temwritepred() to write out TEM output
           variables directly to spatially explicit data sets instead of a 3-
           dimensional array formerly known as tempred[][][].
20000616 - DWK adds temwritemiss() to write out TEM output for grid cells that
           either have missing data or have environmental conditions that will
           not allow biological organisms to survive (i.e., ling-term MAXTAIR <
           -1.0 degree Celsius or long-term mean annual precipitation equal to
           0.0).
           TEM output files instead of
20000629 - DWK changed dyr to outyr in call to temwritepred in the "equilibrium
           portion" of runtem()
20000718 - DWK modified runtem() so that grid cells with missing data would
           still report output to spatially explicit files (i.e., bug fix)
20000719 - DWK added the global constant TQCZEROFLAG to telm423e.hpp to flag
           when environmental conditions of a grid cell cannot support life.  In
           addition, the appropriate conditional statements related to tqc in
           runtem() and transqc() were modified to use this global constant.
20000719 - DWK modified the assignment of tem.veg.cmnt in runtem() to equal
           zero if tem.veg.temveg was less than zero or greater than (NUMVEG+1)
           (i.e., grid cell is covered with water)
20000719 - DWK added a conditional statement to the transient portion of
           runtem() that if "tqc = ACCEPT" and "qc = REJECT", then
           "tqc = REJECT"
20000719 - DWK fixed bug in temwritemiss() by switching the positions of
           ntempred and natmspred in the function call.
20000719 - DWK assigns "ACCEPT" to tqc at the beginning of runtem()
20000719 - DWK modifies the assignment of tqc from transqc() so that this
           assignment occurs only if tqc is currently equal ACCEPT in runtem()
           (i.e., it did not get changed because of unsuitable long-term
           environmental conditions in the grid cell.
20010315 - JSC adjust in equilibrium and transient run parts for the
				use of transient solar radiation input
            	LOOK for JSC

20010418 - Q. Z. Soil thermal model

20020202 Kick added the following in the temwritepred ()
         //strcmp(predname[k],tem.predstr[tem.I_FPC]) == 0
      
*****************************************************************
************************************************************** */

#if !defined(TELM423E1_H)
  #include "telm423e1.hpp"
#endif

/* ************************************************************ */

TEMelmnt::TEMelmnt()
{

  col = MISSING;
  row = MISSING;
  carea = -999;
  lowtair = 0;
  noprec = 0;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

int TEMelmnt::atmsgisin(FILE* fclds, const int& cldflag,
                        FILE* flonlat, const int& lonlatflag,
                        const int& numspin, const int& spintime,
                        const int& RTIME)
{

  int i;
  int gisend;
  Clmdata clds;
  Latdata lonlat;

  if (clm.tcldsflag == 1)
  {
    gisend = loadteclds(fclds, clds, numspin, spintime, RTIME);
    if (gisend == -1) { return gisend; }
  }
  else
  {
    gisend = clds.getdel(fclds);
    if (gisend == -1)
    {
      col = MISSING;
      row = MISSING;
      return gisend;
    }
    col = clds.col;
    row = clds.row;
    carea = clds.carea;
    for (i = 0; i < CYCLE; i++)
    {
      if (cldflag == 1) { clm.clds[i] = clds.mon[i]; }
      else { clm.nirr[i] = clds.mon[i]; }
    }
    strcpy(contnent,clds.contnent);
  }

  if (lonlatflag == 1) { lat = (double) row; }
  else
  {
    gisend = lonlat.getdel(flonlat);
    if (gisend != -1) { lat = lonlat.lat; }
    else
    {
      col = MISSING;
      row = MISSING;
      return gisend;
    }
  }

  return gisend;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt::atmswritemiss(ofstream fout[NUMATMS], char predname[MAXPRED][9],
                             const int& dyr, const int& natmspred,
                             const double value)
{

  int i;
  int dm;
  Clmdata atmspred;

  for (i = 0; i < natmspred; i++)
  {
    for (dm = 0; dm < CYCLE; dm++)
    {
      atmspred.mon[dm] = value;
    }

	 atmspred.outdel(fout[i], col, row, predname[i], carea, atmstotyr[dyr],
                    atmspred.mon, contnent);
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TEMelmnt::atmswritepred(ofstream fout[NUMATMS], char predname[MAXPRED][9],
                             const int& dyr, const int& natmspred)
{

  int i;
  int dm;
  Clmdata atmspred;

  for (i = 0; i < natmspred; i++)
  {
    if (strcmp(predname[i],clm.predstr[clm.I_PAR]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++) { atmspred.mon[dm] = clm.par[dm]; }
    }
    else if (strcmp(predname[i],clm.predstr[clm.I_NIRR]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++) { atmspred.mon[dm] = clm.nirr[dm]; }
   }
    else if (strcmp(predname[i],clm.predstr[clm.I_GIRR]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++) { atmspred.mon[dm] = clm.girr[dm]; }
    }
    else if (strcmp(predname[i],clm.predstr[clm.I_CLDS]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++) { atmspred.mon[dm] = clm.clds[dm]; }
    }

	 atmspred.outdel(fout[i], col, row, predname[i], carea, atmstotyr[dyr],
                    atmspred.mon, contnent);
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int TEMelmnt::coregerr(ofstream& rflog1, char varname1[9],
                       const float& col1, const float& row1,
                       char varname2[9], const float& col2, const float& row2)
{

  int fatalerr = 0;

  if (col1 != col2 || row1 != row2)
  {
    fatalerr = 1;

    cout << "ERROR:  " << varname1 << " data and ";
    cout << varname2 << "data are not coregistered." << endl;
    cout << "COL = " << col1 << " and ROW = " << row1 << " in ";
    cout << varname1 << " data" << endl;
    cout << "COL = " << col2 << " and ROW = " << row2;
    cout << " in " << varname2 << " data" << endl;

    rflog1 << "ERROR:  " << varname1 << " data and ";
    rflog1 << varname2 << "data are not coregistered." << endl;
    rflog1 << "COL = " << col1 << " and ROW = " << row1 << " in ";
    rflog1 << varname1 << " data" << endl;
    rflog1 << "COL = " << col2 << " and ROW = " << row2 << " in ";
    rflog1 << varname2 << " data" << endl;
  }

  return fatalerr;

};

/* **************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt::GISsetELMNTstate(const int& dcmnt,
                                Temdata initstat[NUMMSAC][MAXSTATE+1])
{

/* **************************************************************
  Function initializes TEM state variables from spatially
explicit data sets
************************************************************** */

  int dm;

  // Initialize ODE carbon, nitrogen, and water state variables
  //   for element

  // Carbon
  tem.y[tem.I_VEGC] = initstat[itype][tem.I_VEGC].mon[CYCLE-1];
  tem.y[tem.I_SOLC] = initstat[itype][tem.I_SOLC].mon[CYCLE-1];

  // Nitrogen
  tem.y[tem.I_STRN] = initstat[itype][tem.I_STRN].mon[CYCLE-1];
  tem.y[tem.I_STON] = initstat[itype][tem.I_STON].mon[CYCLE-1];
  tem.y[tem.I_STON] *= 0.001;
  tem.y[tem.I_SOLN] = initstat[itype][tem.I_SOLN].mon[CYCLE-1];
  tem.y[tem.I_AVLN] = initstat[itype][tem.I_AVLN].mon[CYCLE-1];
  tem.y[tem.I_AVLN] *= 0.001;

  // Water
  tem.y[tem.I_AVLW] = initstat[itype][tem.I_AVLW].mon[CYCLE-1];
  tem.y[tem.I_RGRW] = initstat[itype][tem.I_RGRW].mon[CYCLE-1];
  tem.y[tem.I_SNWPCK] = initstat[itype][tem.I_SNWPCK].mon[CYCLE-1];
  tem.y[tem.I_SGRW] = initstat[itype][tem.I_SGRW].mon[CYCLE-1];
  tem.y[tem.I_SM] = tem.y[tem.I_AVLW] + tem.soil.wiltpt;
  tem.y[tem.I_PCTP] = 100.0 * tem.y[tem.I_SM] / tem.soil.totpor;
  tem.y[tem.I_VSM] = tem.y[tem.I_SM] / (tem.soil.rootz * 1000.0);
  if (tem.y[tem.I_VSM] <= 0.0) {
    tem.y[tem.I_VSM] = 0.001;
  }

  // Phenology
  tem.y[tem.I_UNRMLF] = initstat[itype][MAXSTATE].mon[CYCLE-1];
  tem.y[tem.I_UNRMLF] *= 0.01;
  tem.y[tem.I_LEAF] = tem.y[tem.I_UNRMLF] / tem.veg.prvleafmx[dcmnt];


  // Initialize monthly carbon, nitrogen, and water pools for element
  for (dm = 0; dm < CYCLE; dm++)
  {
    // Carbon pools
    tem.veg.plant[dm].carbon = tem.y[tem.I_VEGC];
    tem.soil.org[dm].carbon = tem.y[tem.I_SOLC];
    tem.totalc[dm] = 0.0;

    // Nitrogen pools
    tem.veg.strctrl[dm].nitrogen = tem.y[tem.I_STRN];
    tem.veg.labile[dm].nitrogen = tem.y[tem.I_STON];
    tem.veg.plant[dm].nitrogen = 0.0;
    tem.soil.org[dm].nitrogen = tem.y[tem.I_SOLN];
    tem.soil.availn[dm] = tem.y[tem.I_AVLN];

    // Water pools
    tem.soil.avlh2o[dm] = tem.y[tem.I_AVLW];
    tem.soil.rgrndh2o[dm] = tem.y[tem.I_RGRW];
    tem.soil.snowpack[dm] = tem.y[tem.I_SNWPCK];
    tem.soil.sgrndh2o[dm] = tem.y[tem.I_SGRW];
    tem.soil.moist[dm] = tem.y[tem.I_SM];
    tem.soil.pctp[dm] = tem.y[tem.I_PCTP];
    tem.soil.vsm[dm] = tem.y[tem.I_VSM];

    // Phenology
    tem.veg.unnormleaf[dm] = tem.y[tem.I_UNRMLF];
    tem.veg.leaf[dm] = tem.y[tem.I_LEAF];

  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::loadteclds(FILE* fclds, Clmdata clds,
                         const int& numspin, const int& spintime,
                         const int& RTIME)
{

  int i;
  int j;
  int k;
  int gisend = 1;
  int dm;
  int dyr;
  int totspin;

// Get CLDINESS from transient data set

  totspin = numspin * spintime;
  for (i = totspin; i < RTIME; i++)
  {
    k = i - totspin;
    gisend = clds.getdel(fclds);
    if (gisend == -1)
    {
      col = MISSING;
      row = MISSING;
      return gisend;
    }

// Determine CLDINESS for equilibrium conditions

    if ( k == 0)
    {
      col = clds.col;
      row = clds.row;
      carea = clds.carea;
      strcpy(contnent,clds.contnent);
      clm.cldsyear[0] = clds.year - totspin;
      for (j = 0; j < CYCLE; j++)
      {
        clm.tclds[0][j] = clds.mon[j];
      }
    }

// Determine CLDINESS for transient conditions

    else
    {
      for (j = 0; j < CYCLE; j++)
      {
        clm.tclds[i][j] = clds.mon[j];
      }
      clm.cldsyear[i] = clds.year;
    }
  }

// Determine CLDINESS during spin up

  for (i = 0; i < numspin; i++)
  {
    dyr = totspin + 1;
    for (j = 0; j < spintime; j++)
    {
      k = (i * spintime) + j + 1;
      clm.cldsyear[k] = clm.cldsyear[0] + k;
      for (dm = 0; dm < CYCLE; dm++)
      {
        clm.tclds[k][dm] = clm.tclds[dyr][dm];
      }
      ++dyr;
    }
  }

  return gisend;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::loadtelulc(FILE* flulc, Lulcdata lulc,
                         const int& numspin, const int& spintime,
                         const int& RTIME, int& ftlerr, ofstream& flog1)
{

  int i;
  int j;
  int k;
  int gisend = 1;
  int totspin;

// Get LULC from transient data set

  totspin = numspin*spintime;
  for (i = totspin; i < RTIME; i++)
  {
    k = i - totspin;
    gisend = lulc.getdel(flulc);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "LULC", lulc.col, lulc.row);

// Determine LULC for equilibrium conditions

    if ( k == 0)
    {
      tem.ag.lulcyear[0] = lulc.year - totspin;
      tem.ag.tstate[0] = lulc.agstate;
      tem.ag.tRAP[0] = lulc.RAP;
    }

// Determine LULC for transient conditions

    else
    {
      tem.ag.lulcyear[i] = lulc.year;
      tem.ag.tstate[i] = lulc.agstate;
      tem.ag.tRAP[i] = lulc.RAP;
    }
  }

// Determine LULC during spin up

  for (i = 0; i < numspin; i++)
  {
    for (j = 0; j < spintime; j++)
    {
      k = (i * spintime) + j + 1;
      tem.ag.lulcyear[k] = tem.ag.lulcyear[0] + k;
      tem.ag.tstate[k] = tem.ag.tstate[0];
      tem.ag.tRAP[k] = tem.ag.tRAP[0];
    }
  }

 return gisend;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::loadtenpp(FILE* fnpp, Temdata npp, const int& maxtype,
                        const int& RTIME, int& ftlerr, ofstream& flog1)
{

//Note: numspin and spintime removed from function call by D. Kicklighter 990724

  int dv;
  int gisend = 1;
  int dm;
  int dyr;
//  int totspin;

// Get POTNPP from transient data set

  for (dyr = 0; dyr < RTIME; dyr++)
  {
    for (dv = 0; dv < maxtype; dv++)
    {
      gisend = npp.getdel(fnpp);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "POTNPP", npp.col, npp.row);
      for (dm = 0; dm < CYCLE; dm++)
      {
        tem.ag.tpotnpp[dv][dyr][dm] = npp.mon[dm];
      }
    }
  }

  return gisend;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::loadteprec(FILE* fprec, Clmdata prec,
                         const int& numspin, const int& spintime,
                         const int& RTIME, int& ftlerr, ofstream& flog1)
{

  int i;
  int j;
  int k;
  int gisend = 1;
  int dm;
  int dyr;
  int totspin;


// Get PREC from transient data set

  noprec = 0;
  totspin = numspin * spintime;
  for (i = totspin; i < RTIME; i++)
  {
    k = i - totspin;
    gisend = prec.getdel(fprec);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "PREC", prec.col, prec.row);

// Determine PREC for equilibrium conditions

    if ( k == 0)
    {
      for (j = 0; j < CYCLE; j++) { tem.atms.tprec[0][j] = prec.mon[j]; }
      tem.atms.precyear[0] = prec.year - totspin;
      tem.atms.yrtprec[0] = prec.total;
      if (prec.total <= 0.0) { noprec = 1; }
    }

// Determine PREC for transient conditions

    else
    {
      for (j = 0; j < CYCLE; j++) { tem.atms.tprec[i][j] = prec.mon[j]; }
      tem.atms.precyear[i] = prec.year;
      tem.atms.yrtprec[i] = prec.total;
    }
  }

// Determine PREC during spin up

  for (i = 0; i < numspin; i++)
  {
    dyr = totspin + 1;
    for (j = 0; j < spintime; j++)
    {
      k = (i * spintime) + j + 1;
      tem.atms.precyear[k] = tem.atms.precyear[0] + k;
      tem.atms.yrtprec[k] = tem.atms.yrtprec[dyr];
      for (dm = 0; dm < CYCLE; dm++)
      {
        tem.atms.tprec[k][dm] = tem.atms.tprec[dyr][dm];
      }
      ++dyr;
    }
  }

  return gisend;

};

/* *************************************************************
************************************************************** */
//XL add for msoil

/* *************************************************************
************************************************************** */

int TEMelmnt::loadtemsoil(FILE* fmsoil, Clmdata msoil,
                         const int& numspin, const int& spintime,
                         const int& RTIME, int& ftlerr, ofstream& flog1)
{

  int i;
  int j;
  int k;
  int gisend = 1;
  int dm;
  int dyr;
  int totspin;


// Get MSOL from transient data set


  totspin = numspin * spintime;
  for (i = totspin; i < RTIME; i++)
  {
    k = i - totspin;
    gisend = msoil.getdel(fmsoil);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "MSOL", msoil.col, msoil.row);

// Determine MSOL for equilibrium conditions

    if ( k == 0)
    {
      for (j = 0; j < CYCLE; j++) { tem.atms.tmsoil[0][j] = msoil.mon[j]; }
      // tem.atms.precyear[0] = msoil.year - totspin;
      // tem.atms.yrtprec[0] = msoil.total;
      // tem.atms.tairyear[0] = tair.year - totspin;
      // tem.atms.mxttair[0] = tair.max;
      // if (tair.max < -1.0) { lowtair = 1; }
      // for (j = 0; j < CYCLE; j++)
      // {
      //   tem.atms.ttair[0][j] = tair.mon[j];
      // }
    }

// Determine MSOL for transient conditions

    else
    {
      for (j = 0; j < CYCLE; j++) { tem.atms.tmsoil[i][j] = msoil.mon[j]; }
      tem.atms.msoilyear[i] = msoil.year;
      // tem.atms.yrtprec[i] = prec.total;
    }
  }

// Determine PREC during spin up

  for (i = 0; i < numspin; i++)
  {
    dyr = totspin + 1;
    for (j = 0; j < spintime; j++)
    {
      k = (i * spintime) + j + 1;
      tem.atms.msoilyear[k] = tem.atms.msoilyear[0] + k;
      // tem.atms.yrtmsoil[k] = tem.atms.yrtmsoil[dyr];
      for (dm = 0; dm < CYCLE; dm++)
      {
        tem.atms.tmsoil[k][dm] = tem.atms.tmsoil[dyr][dm];
      }
      ++dyr;
    }
  }

  return gisend;

};

/* *************************************************************
************************************************************** */

/* *************************************************************
************************************************************** */

int TEMelmnt::loadtetair(FILE* ftair, Clmdata tair,
                         const int& numspin, const int& spintime,
                         const int& RTIME, int& ftlerr, ofstream& flog1)
{

  int i;
  int j;
  int k;
  int gisend = 1;
  int dm;
  int dyr;
  int totspin;

// Get TAIR from transient data set

  lowtair = 0;
  totspin = numspin * spintime;
  for (i = totspin; i < RTIME; i++)
  {
    k = i - totspin;
    gisend = tair.getdel(ftair);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "TAIR", tair.col, tair.row);

// Determine TAIR for equilibrium conditions

    if ( k == 0)
    {
      tem.atms.tairyear[0] = tair.year - totspin;
      tem.atms.mxttair[0] = tair.max;
      if (tair.max < -1.0) { lowtair = 1; }
      for (j = 0; j < CYCLE; j++)
      {
        tem.atms.ttair[0][j] = tair.mon[j];
      }
    }

// Determine TAIR for transient conditions

    else
    {
      tem.atms.tairyear[i] = tair.year;
      tem.atms.mxttair[i] = tair.max;
      for (j = 0; j < CYCLE; j++) { tem.atms.ttair[i][j] = tair.mon[j]; }
    }
  }
// Determine TAIR during spin up

  for (i = 0; i < numspin; i++)
  {
    dyr = totspin + 1;
    for (j = 0; j < spintime; j++)
    {
      k = (i * spintime) + j + 1;
      tem.atms.tairyear[k] = tem.atms.tairyear[0] + k;
      tem.atms.mxttair[k] = tem.atms.mxttair[dyr];
      for (dm = 0; dm < CYCLE; dm++)
      {
        tem.atms.ttair[k][dm] = tem.atms.ttair[dyr][dm];
      }
      ++dyr;
    }
  }

  return gisend;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TEMelmnt::runtem(ofstream& rflog1, char predmap[MAXPRED][9],
                      const int& cldflag, const int& atmsflag,
                      const int& atmsoutfg, int& natmspred,
                      ofstream fsradout[NUMATMS],
                      const int& temflag, const int& kdinflg,
                      int& ntempred, ofstream ftemout[MAXPRED],
                      const int& stateflag,
                      const int& equil, const int& totsptime,
                      const int& RTIME, const int& ispinout)
{

//  int i;
//  int j;
//  int k;
//  int dt;
  int dm;
  int dyr = 0;
  int wrtyr;
//  int maxtime;

  // Assignment of "ACCEPT" to tqc added by DWK on 20000719
  int tqc = ACCEPT;
  qc = ACCEPT;
  tem.totyr = 0;


  if (atmsflag == 1)
  {

// Check cloudiness input for valid data

    for (dm = 0; dm < CYCLE; dm++)
    {
      if (cldflag == 0 && clm.nirr[dm] <= -99.0) { qc = 1; }
      if (cldflag == 1 && clm.tcldsflag == 0 && clm.clds[dm] <= -99.0)
      {
        qc = 2;
      }
      if (cldflag == 1 && clm.tcldsflag == 1 && clm.tclds[0][dm] <= -99.0)
      {
        qc = 3;
      }
    }

    if (qc == ACCEPT)
    {

// Calculate GIRR from POTSCLM model

      clm.yrsumday = 0.0;
      for (dm = 0; dm < CYCLE; dm++)
      {
        clm.girr[dm] = clm.xgirr(lat,dm,clm.yrsumday);
      }
    }
    else { rflog1 << "cldqc = " << qc << endl; }
  }


//  cout <<endl<<"OK 1"<<endl;
//  system("pause");
/* *************************************************************
		BEGIN VEGETATION MOSAIC LOOP
************************************************************* */

  for (itype = 0; itype < maxtype; itype++)
  {
	qc = ACCEPT;
//	cout <<"OK 2"<<endl<<"maxtype="<<maxtype<<endl;
//	system("pause");

    //Following "if" statement added by DWK on 20000719 to modify
    // the assignment of tem.veg.cmnt

    if (tem.veg.temveg > 0 && tem.veg.temveg < (NUMVEG+1))
    {

//     printf(" %d", tem.veg.temveg);

    tem.veg.cmnt = tem.veg.subtype[mez][itype];

    }
	else { tem.veg.cmnt = 0; }

//	cout <<"OK 3"<<endl;


// Calculate NIRR, CLDINESS and PAR for Equilibrium Conditions with POTSCLM

    if (atmsflag == 1)
	{
      for (dm = 0; dm < CYCLE; dm++)
      {
 //	     if (clm.tcldsflag == 1) { clm.clds[dm] = clm.tclds[0][dm]; }
 //	     if (cldflag == 1)

 //       {
 //	       clm.nirr[dm] = clm.xnirr(clm.clds[dm], clm.girr[dm]);
 //       }
      if (clm.tcldsflag == 1)
       if (cldflag == 1) { clm.clds[dm] = clm.tclds[0][dm]; }
         else { clm.nirr[dm] = clm.tclds[0][dm]; }
//      else { clm.nirr[dm] = clm.tclds[0][dm] * 2.387; }
//     }
     if (cldflag == 1) {
      clm.nirr[dm] = clm.xnirr(clm.clds[dm], clm.girr[dm]);
      }
     //  JSC 3/15/01 added above from xtem41e for trans radiation
		else
        {
	       clm.clds[dm] = clm.mkclds(clm.girr[dm], clm.nirr[dm]);
        }
	     clm.par[dm]  = clm.xpar(clm.clds[dm], clm.nirr[dm]);
	  }
//	cout <<"OK 4"<<endl;
//	system("pause");

// Output POTSCLM steady state results to files

      if (atmsoutfg == 1 && itype == 0)
      {
	     dyr = 0;
	     atmswritepred(fsradout, predmap, dyr, natmspred);
	  }

	}
//	cout <<"OK 5"<<endl;
//	system("pause");


// Initialize TEM output variables to missing value

    if (temflag == 1)
    {

      tem.microbe.kdsave[itype] = -9.9999;
      if (kdinflg == 1) { tem.microbe.kdc = tem.microbe.kdin[itype]; }

// "Hand-off" POTSCLM equilibrium output and equilibrium input data to TEM

      for (dm = 0; dm < CYCLE; dm++)
      {

	     if (atmsflag == 1)
        {
	       tem.atms.nirr[dm] = clm.nirr[dm];
	       tem.atms.par[dm] = clm.par[dm];
	     }
	     tem.atms.co2[dm] = tem.atms.co2level;

        // Update TAIR with current year data from transient data files
	     if (tem.atms.ttairflag == 1)
        {
          tem.atms.tair[dm] = tem.atms.ttair[0][dm];
        }

//XL add for msoil
        // Update TAIR with current year data from transient data files
	     if (tem.atms.tmsoilflag == 1)
        {
          tem.atms.msoil[dm] = tem.atms.tmsoil[0][dm];
        }


        // Update PREC with current year data from transient data files
	     if (tem.atms.tprecflag == 1)
        {
          tem.atms.prec[dm] = tem.atms.tprec[0][dm];
        }

        // Update LULC with current year data from transient data files
   	  if (tem.ag.tlulcflag == 1)
        {
          // if rap is 0 (RAP0flag == 1), then we don't need tpotnpp
	       if(tem.ag.RAP0flag == 0)
          {
  	         tem.ag.potnpp[dm] = tem.ag.tpotnpp[itype][0][dm];
	       }
          else { tem.ag.potnpp[dm] = 0.0; }
	     }
	  }
//	  cout <<"OK 6"<<endl;
//	  system("pause");

      if (tem.atms.ttairflag == 1)
      {
        tem.atms.mxtair = tem.atms.mxttair[0];
        if (tem.atms.mxtair < -1.0) { lowtair = 1; }
      }
      if (tem.atms.tprecflag == 1)
      {
        tem.atms.yrprec = tem.atms.yrtprec[0];
        if (tem.atms.yrprec <= 0.0) { noprec = 1; }
      }

// Check TEM input for valid data

      qc = temgisqc(stateflag,tem.soil.pctsilt, tem.soil.pctclay,
                    tem.veg.cmnt, tem.elev, tem.atms.nirr, tem.atms.par,
                    tem.atms.tair, tem.atms.mxtair, tem.atms.prec,tem.atms.msoil,
					tem.atms.yrprec, initstat);
//	  cout <<"OK 7"<<endl;
//	  system("pause");


// Check TEM parameters for specific vegetation types

      // If "qc == REJECT", write code of underlying reason to log file
      if (qc != ACCEPT) { rflog1 << "temgisqc = " << qc << endl; }
	  else { qc = tem.ecdqc(tem.veg.cmnt); }

      // "qc == ACCEPT" condition moved ahead of "qc != ACCEPT".  In addition
      //  "qc != ACCEPT" condition modified to write out grid cells with
      // missing values (or 0.0) to spatially explicit files
      if (qc == ACCEPT)
      {

/* *******************************************************************
Start the Terrestrial Ecosystem Model (TEM) for Equilibrium Conditions
******************************************************************* */

        tem.soil.xtext(tem.veg.cmnt,tem.soil.pctsilt,tem.soil.pctclay);


// Initialize TEM parameters to grid cell's vegetation type, soil texture,
//   PET and AET

	     tem.setELMNTecd(kdinflg,tem.veg.cmnt, tem.soil.psiplusc);
        tem.setELMNTevap(stateflag, tem.veg.cmnt, tem.atms.pet, tem.atms.tair);


// "While" loop to allow adaptive integrator tolerance (i.e. tem.tol) to be
// reduced if chaotic behavior occurs

        tem.ag.state = 0;
        tem.ag.prvstate = 0;

	     for (dyr = 0; dyr < RTIME; dyr++) { tem.qualcon[dyr][itype] = 0; }

	     tem.nattempt = 0;
	     tem.tol = tem.inittol;
		 tem.baseline = tem.initbase;
//		 cout <<"OK 8"<<endl;
//		 system("pause");
	     while (tem.nattempt < tem.maxnrun)
		{

// Initialize standing stocks of carbon and nitrogen from spatially-explicit
// data sets (i.e. stateflag = 1) or calibration ("ECD") data

	       if (stateflag == 1) { GISsetELMNTstate(tem.veg.cmnt, initstat); }
	       else { tem.ECDsetELMNTstate(tem.veg.cmnt, tem.soil.psiplusc); }
	       tem.setELMNTflux();
//		   cout <<"OK 9"<<endl;
//         system("pause");
// Run TEM until steady state conditions occur (equilibrium)

		   tem.nattempt = tem.equilibrium(itype,tem.tol);
//   		   cout <<"OK 10"<<endl;
//		   system("pause");
		   if (tem.nattempt < tem.maxnrun) { tem.tol /= 10.0; }

		}
//		cout <<"OK 11"<<endl;
//		system("pause");
	     outyr = 0;
	     tem.qualcon[outyr][itype] += tem.nattempt;
	     if (tem.nattempt == tem.maxnrun && tem.totyr == tem.maxyears)
        {
	       tem.qualcon[outyr][itype] += 20;
	     }
      } // end of "qc == ACCEPT"



/* *************************************************************
Output TEM steady state results to files.  Assign the baseline equilibrium
conditions to the first year of the output data if a transient analysis is being
conducted.  Otherwise, report the number of years it took to reach an
equilibrium state for local environmental conditions
************************************************************* */

      // (The following code was added 20000718 by DWK and incorporates
      // components from both the former "qc == ACCEPT" and "qc == REJECT"
      // conditions if the eariler code).
	   if (equil == 0)
      {
        ttotyr[outyr][itype] = tem.startyr - totsptime - 1;
      }
	   else
      {
        if (qc == ACCEPT) { ttotyr[outyr][itype] = tem.totyr; }
        else {

/* ***********************************************************************
 If long-term environmental conditions of the grid cell cannot support the
 terrestrial biota (as indicated by input data), assume no living or dead
 organic matter exists in the grid cell and equilibrium conditions were
 achieved in a single year (i.e., tem.totyr = 1).  Otherwise, assign a missing
 value to flag the grid cell to provide information about potential problems
 with the input data sets.
 *********************************************************************** */

         // initstat[itype][0].ave changed to initstat[itype][tem.I_VEGC].ave by
        //    DWK on 20000214
        // tem.ez changed to tem.veg.cmnt by DWK on 20000214
          if ((tem.veg.cmnt > 0) && ((noprec == 1) || (lowtair == 1)
             || (stateflag == 1 && initstat[itype][tem.I_VEGC].ave < 0.1
             && initstat[itype][tem.I_VEGC].ave > MISSING)))
          {
	         ttotyr[outyr][itype] = 1;
          }
          else { ttotyr[outyr][itype] = -999; }
        }
      }



      // Write TEM output to spatially explicit data files
      if (qc == ACCEPT) // i.e., Good data
      {
        // DWK changed dyr to outyr in temwritepred on 20000629
        if(outyr > totsptime || ispinout == 1 || tem.totyr == (tem.startyr -1))

        {temwritepred(ftemout,predmap,outyr,itype,ntempred,natmspred);}
      }
      else // Something's not right!
      {

         // If "qc == REJECT", write code of underlying reason to log file
        rflog1 << "temecdqc = " << qc << endl;

        // Assign zero to all TEM output variables if long-term environmental
        // conditions cannot support life.

        if ((tem.veg.cmnt > 0) && ((noprec == 1) || (lowtair == 1)
           || (stateflag == 1 && initstat[itype][tem.I_VEGC].ave < 0.1
           && initstat[itype][tem.I_VEGC].ave > MISSING)))
        {

          // tqc is a flag used to assign zero to all TEM variables below for
          // both the initialization and the transient simulations of TEM
          tqc = TQCZEROFLAG;
			if(outyr >= totsptime || ispinout == 1 || tem.totyr == (tem.startyr - 1))
          temwritemiss(ftemout,predmap,outyr,itype,ntempred,natmspred,ZERO);
        }
        else
        {
          // Skip TEM if input data and/or parameters are missing
			if(outyr >= totsptime || ispinout == 1 || tem.totyr == (tem.startyr - 1))
          temwritemiss(ftemout,predmap,outyr,itype,ntempred,natmspred,MISSING);
        }
      } // end of "qc = REJECT"

 	   if (equil == 0)
      {

/* *******************************************************************
                     Begin Transient Conditions
******************************************************************* */

	     tem.microbe.kdsave[itype] = tem.microbe.kd;
	     tem.totyr = tem.startyr - totsptime;

        // Following "if" conditional statement added by DWK on 20000719
        if (tqc == ACCEPT)
        {
	       tqc = transqc(tem.maxyears,tem.totyr, tem.veg.plant);
        }

        dyr = 1;
        outyr = 1;
        tem.totyr = 0;
        tem.baseline = 0;
        tem.wrtyr = -99;

        while (tem.totyr < tem.endyr)
        {

          //Following "if" statement added by DWK on 20000719
          if (tqc == ACCEPT && qc != ACCEPT) { tqc = qc; }

          // Following "if" statement moved within year loop by DWK on
          // 20000718

//llc test
//  printf("=2=%i,  %i, \n",tqc,ACCEPT);
  	       if (tqc == ACCEPT)
          {

// Calculate NIRR, CLDINESS and PAR for Spin up and Transient Conditions
//   with POTSCLM  -- also, hand off NIRR and PAR to TEM

	         if (atmsflag == 1 && clm.tcldsflag == 1)
            {
              tem.totyr = clm.cldsyear[dyr];
		        for (dm = 0; dm < CYCLE; dm++)
             {
 //		          clm.clds[dm] = clm.tclds[dyr][dm];
 //		          if (cldflag == 1)
 //               {
 //		            clm.nirr[dm] = clm.xnirr(clm.clds[dm], clm.girr[dm]);
             if (cldflag == 1) { clm.clds[dm] = clm.tclds[dyr][dm]; }
              else { clm.nirr[dm] = clm.tclds[dyr][dm]; }
 //           else { clm.nirr[dm] = clm.tclds[dyr][dm] * 2.387; }

  				 if (cldflag == 1) {
             clm.nirr[dm] = clm.xnirr(clm.clds[dm], clm.girr[dm]);
 // JSC added 3/15/01 for transient radiation input
 		          }
  		          else
                {
 		            clm.clds[dm] = clm.mkclds(clm.girr[dm], clm.nirr[dm]);
 		          }
 		          tem.atms.par[dm]  = clm.xpar(clm.clds[dm], clm.nirr[dm]);
 		          tem.atms.nirr[dm] = clm.nirr[dm];

		        }

// Output POTSCLM transient results to files

		        if (atmsoutfg == 1 && itype == 0)
              {
 	             atmswritepred(fsradout, predmap, dyr, natmspred);
		        }

	         }



// Run the Terrestrial Ecosystem Model (TEM) under transient conditions

 	         wrtyr = tem.transient(dyr,itype, tem.tol,RTIME);

// Output TEM transient results for specified years to files




	         if ((wrtyr%tem.diffyr) == 0)
            {
		        ttotyr[outyr][itype] = tem.totyr;
			if(wrtyr >=0 || ispinout ==1 || tem.totyr == (tem.startyr - 1))        
              temwritepred(ftemout,predmap,outyr,itype,ntempred,natmspred);
		        ++outyr;
	        }
	         ++dyr;
  	       } // End of tqc = ACCEPT "if" statement
	       else
          {
            if (outyr == 1) { rflog1 << "tqc = " << tqc << endl; }
            tem.totyr = tem.startyr - totsptime - 1
                        + (outyr * tem.diffyr);

            ttotyr[outyr][itype] =  tem.totyr;

            // Assign zero to all TEM output variables if initial long-term
            // environmental conditions were unable to support life.
            if (tqc == TQCZEROFLAG)
            {
				if(outyr >= totsptime || ispinout == 1 || tem.totyr == (tem.startyr - 1))
			  temwritemiss(ftemout,predmap,outyr,itype,ntempred,
			  	natmspred,ZERO);
            }
            else
            {
              // Skip TEM if input data and/or parameters are missing
				if(outyr >= totsptime || ispinout == 1 || tem.totyr == (tem.startyr - 1))
			  temwritemiss(ftemout,predmap,outyr,itype,ntempred,
			  natmspred,MISSING);
            }
            ++outyr;

	       } // End of tqc = REJECT
        } // End of While totyr < endyr
	   } // End of Transient else
    } // End of the Terrestrial Ecosystem Model (TEM)
  } // End of Vegetation Mosaic Loop
//llc test

};

/* *************************************************************
************************************************************* */


/* **************************************************************
************************************************************** */

int TEMelmnt::temgisin(ofstream& flog1, int& ftlerr, const int& atmsflag,
		       const int& itype, const int& maxtype,
                       const int& kdinflg, const int& stateflag,
                       FILE* fstxt, FILE*felev, FILE* fnirr, FILE* fpar,
		       FILE* ftair, FILE* fprec, FILE* fmsoil, FILE* flulc, FILE* fnpp,
                       FILE* fkdin, FILE* fstate[MAXESTAT+7],
		       const int& numspin, const int& spintime,
                       const int& RTIME) //Xl add for soil moisture
{

  int i;
  int gisend = 1;

  Soildata fao;
  Elevdata elv;
  Clmdata nirr;
  Clmdata par;
  Clmdata tair;
  Clmdata prec;
  Clmdata msoil;
  Lulcdata lulc;
  Temdata npp;
  Temdata leaf;
  KDdata kddat;

  if (itype == 0)
  {

    noprec = 0;
    lowtair = 0;

    gisend = fao.getdel(fstxt);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "TEXTURE", fao.col, fao.row);
    carea = fao.carea;
    tem.soil.pctsilt = fao.pctsilt;
    tem.soil.pctclay = fao.pctclay;
    tem.soil.wsoil = fao.wsoil;

    gisend = elv.getdel(felev);
    if (gisend == -1) { return gisend; }
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "ELEV", elv.col, elv.row);
    tem.elev = elv.elev;

    if (atmsflag == 0)
    {
      gisend = nirr.getdel(fnirr);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "NIRR", nirr.col, nirr.row);
      for (i = 0; i < CYCLE; i++) { tem.atms.nirr[i] = nirr.mon[i]; }

      gisend = par.getdel(fpar);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "PAR", par.col, par.row);
      for (i = 0; i < CYCLE; i++) { tem.atms.par[i] = par.mon[i]; }
    }

    if (tem.atms.ttairflag == 1)
    {
      gisend = loadtetair(ftair, tair, numspin, spintime, RTIME, ftlerr, flog1);
      if (gisend == -1) { return gisend; }
      for (i = 0; i < CYCLE; i++) { tem.atms.tair[i] = tem.atms.ttair[0][i]; }
    }
    else
    {
      gisend = tair.getdel(ftair);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "TAIR", tair.col, tair.row);
      tem.atms.mxtair = tair.max;

      // lowtair condition added by DWK on 20000214
      if (tem.atms.mxtair < -1.0) { lowtair = 1; }
      for (i = 0; i < CYCLE; i++) { tem.atms.tair[i] = tair.mon[i]; }
    }

    if (tem.atms.tprecflag == 1)
    {
      gisend = loadteprec(fprec, prec, numspin, spintime, RTIME, ftlerr, flog1);
      if (gisend == -1) { return gisend; }
      for (i = 0; i < CYCLE; i++) { tem.atms.prec[i] = tem.atms.tprec[0][i]; }
    }
    else
    {
      gisend = prec.getdel(fprec);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "PREC", prec.col, prec.row);
      tem.atms.yrprec = prec.total;
      if (prec.total <= 0.0) { noprec = 1; }
      for (i = 0; i < CYCLE; i++) { tem.atms.prec[i] = prec.mon[i]; }
    }

//XL add for msoil
    if (tem.atms.tmsoilflag == 1)
    {
      gisend = loadtemsoil(fmsoil, msoil, numspin, spintime, RTIME, ftlerr, flog1);
      if (gisend == -1) { return gisend; }
      for (i = 0; i < CYCLE; i++) { tem.atms.msoil[i] = tem.atms.tmsoil[0][i]; }
    }
    else
    {
      gisend = msoil.getdel(fmsoil);
      if (gisend == -1) { return gisend; }
      ftlerr = coregerr(flog1,"TEMVEG", col, row, "MSOL", msoil.col, msoil.row);
      // tem.atms.mxmsoil = msoil.max;

      // // lowtair condition added by DWK on 20000214
      // if (tem.atms.mxtair < -1.0) { lowtair = 1; }
      // for (i = 0; i < CYCLE; i++) { tem.atms.tair[i] = tair.mon[i]; }
    }

    if (tem.ag.tlulcflag == 1)
    {
      gisend = loadtelulc(flulc, lulc, numspin, spintime, RTIME, ftlerr, flog1);
      if (gisend == -1) { return gisend; }
      tem.ag.state = tem.ag.tstate[0];
      tem.ag.RAP = tem.ag.tRAP[0];

//Note: numspin and spintime removed from function call by D. Kicklighter 990724

      if(tem.ag.RAP0flag == 0)  // only get potnpp if RAP > 0
      {
	gisend = loadtenpp(fnpp, npp, maxtype, RTIME, ftlerr, flog1);
      }
      if (gisend == -1) { return gisend; }
    }
  }


  if (kdinflg == 1)
  {
    kddat.getdel(fkdin);
    ftlerr = coregerr(flog1,"TEMVEG", col, row, "KD", kddat.col, kddat.row);
    tem.microbe.kdin[itype] = kddat.kd;
  }

/* Read in initialization data */

  if (stateflag == 1)
  {
    for (i = 0; i < MAXSTATE; i++)
    {
      gisend = initstat[itype][i].getdel(fstate[i]);
      if (gisend == -1) { return gisend; }
	  ftlerr = coregerr(flog1,"TEMVEG", col, row, "Initial",
	   initstat[itype][i].col, initstat[itype][i].row);
      if (ftlerr != 0)
      {
	cout  << "Initial data set " << (i+1) << endl;
	flog1 << "Initial data set " << (i+1) << endl;
      }
    }

	tem.veg.prvleafmx[initstat[itype][MAXESTAT+4].subtveg-1]
		= initstat[itype][MAXESTAT+4].max * 0.01;
    tem.atms.yrpet = initstat[itype][MAXESTAT+5].total;
    tem.atms.prvpetmx = initstat[itype][MAXESTAT+5].max;
    tem.atms.yreet = initstat[itype][MAXESTAT+6].total;
    tem.atms.prveetmx = initstat[itype][MAXESTAT+6].max;
    for (i = 0; i < CYCLE; i++)
    {
      tem.atms.pet[i] = initstat[itype][MAXESTAT+5].mon[i];
      tem.veg.unnormleaf[i] = initstat[itype][MAXESTAT+4].mon[i] * 0.01;
    }
  }

//  for (i = 0; i < RTIME; i++) {
//    flog1 << tem.atms.tairyear[i];
//    for (int dm = 0; dm < CYCLE; dm++) {
//      flog1 << " " << setprecision(1) << tem.atms.ttair[i][dm];
//    }
//    flog1 << endl;
//  }

  return gisend;

};

/* *************************************************************
************************************************************** */


/* * *************************************************************
************************************************************** */

int TEMelmnt::temgisqc(const int& stateflag, const double& pctsilt,
			   const double& pctclay,
			   const int& cmnt, const double& elev, double nirr[CYCLE],
			   double par[CYCLE],double tair[CYCLE], double& mxtair,
			   double prec[CYCLE],double msoil[CYCLE], double& yrprec,
		       Temdata initstat[NUMMSAC][MAXSTATE+1])
{

  int i;
  int qc;

  qc = ACCEPT;

  if (pctsilt < 0.0) { return qc = 1; }
  if (pctclay < 0.0) { return qc = 2; }
// Following condition modified by DWK on 20000219
//  if (ez < 0 || ez >= NUMVEG) { return qc = 3; }
  if (cmnt < 1 || cmnt > NUMVEG) { return qc = 3; }
  if (elev <= -999.0) { return qc = 4;}
  if (mxtair < -1.0) { return qc = 5; }
  if (yrprec <= 0.0) { return qc = 6; }

  for (i = 0; i < CYCLE; i++)
  {
    if (nirr[i] <= -1.0) { return qc = 7; }
    if (par[i] <= -1.0) { return qc = 8; }
    if (tair[i] <= -99.0) { return qc = 9; }
    if (prec[i] <= -1.0) { return qc = 10; }
    if (msoil[i] <= -1.0) { return qc = 11; }//XL add for msoil
  }


  if (stateflag == 1)
  {
    for (i = 0; i < (MAXESTAT+5); i++)
    {
      if (initstat[0][i].mon[11] <= MISSING) { return qc = (11 + i); }
    }
    if (initstat[0][0].ave < 0.1) { return qc = 17 + MAXESTAT; }
  }

  return qc;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt::temwritemiss(ofstream fout[MAXPRED], char predname[MAXPRED][9],
                            const int& dyr, const int& itype,
                            const int& ntempred, const int& natmspred,
                            const double value)
{

  int i;
  int k;
  int dm;
  Temdata tempred;

  for (i = 0; i < ntempred; i++)
  {
    k = i + natmspred;

    for (dm = 0; dm < CYCLE; dm++)
    {
   	tempred.mon[dm] = value;
    }

    // Write output data to files

    tempred.outdel(fout[i], col, row, predname[k],
                   tem.veg.temveg, tem.veg.subtype[mez][itype],
                   (100.0 * tem.soil.psiplusc),
                   tem.qualcon[dyr][itype],
                   carea,
                   ttotyr[dyr][itype],
                   tempred.mon,
                   contnent);

  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TEMelmnt::temwritepred(ofstream fout[MAXPRED], char predname[MAXPRED][9],
                            const int& dyr, const int& itype,
                            const int& ntempred, const int& natmspred)
{

  int i;
  int k;
  int dm;
  Temdata tempred;


  for (i = 0; i < ntempred; i++)
  {
    k = i + natmspred;

    if (strcmp(predname[k],tem.predstr[tem.I_NPP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
    	  tempred.mon[dm] = tem.veg.npp[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VEGC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.plant[dm].carbon;  // VEGC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SOLC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.org[dm].carbon;  // SOLC

      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_STRN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.strctrl[dm].nitrogen;  // VSTRUCTN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_STON]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
  	     tempred.mon[dm] = tem.veg.labile[dm].nitrogen * 1000.0;  // VSTOREN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SOLN]) == 0)
    {
    for (dm = 0; dm < CYCLE; dm++)
      {
     	  tempred.mon[dm] = tem.soil.org[dm].nitrogen;  // SOILORGN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AVLN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.availn[dm] * 1000.0;  // AVAILN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_NMIN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.microbe.netnmin[dm] * 1000.0;  // NETNMIN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_NEP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.nep[dm];  // NEP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_GPP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.gpp[dm];  // GPP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_INGPP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.ingpp[dm];  // INGPP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_INNPP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.innpp[dm];  // INNPP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RH]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.microbe.rh[dm];  // RH
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_NLST]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.nlost[dm] * 1000.0;  // NLOST
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_NINP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.ninput[dm] * 1000.0;  // NINPUT
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RVMNT]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.rm[dm];  // RVMAINT
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RVGRW]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.rg[dm];  // RVGRWTH
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_GPR]) == 0)
    {

      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.gpr[dm];  // GPR
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_LTRC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.ltrfal[dm].carbon;  // LTRC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VNUP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.nuptake[dm] * 1000.0;  // VEGNUP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_LTRN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.ltrfal[dm].nitrogen * 1000.0;  // LTRN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_MNUP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.microbe.nuptake[dm] * 1000.0;  // MICRONUP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VNMBL]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.nmobil[dm] * 1000.0;  // VNMOBIL
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VNRSRB]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.nresorb[dm] * 1000.0;  // VNRESORB
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VSUP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.suptake[dm] * 1000.0;  // VEGSUP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VLUP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.luptake[dm] * 1000.0;  // VEGLUP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_UNRMLF]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.unnormleaf[dm] * 100.0;  // UNRMLEAF
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_LEAF]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.veg.leaf[dm] * 100.0;      // LEAF
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_CNVRTC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.convrtflx[dm].carbon;   // CONVERTC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_CNVRTN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.convrtflx[dm].nitrogen * 1000.0; // CONVERTN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRDF10C]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
		 tempred.mon[dm] = tem.ag.formPROD10.carbon / (double) CYCLE;
		 // PRDF10C
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRDF10N]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.formPROD10.nitrogen / (double) CYCLE) *
                          1000.0; // PRDF10N
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRDF100C]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
		 tempred.mon[dm] = tem.ag.formPROD100.carbon / (double) CYCLE;
		 // PRDF100C
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRDF100N]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.formPROD100.nitrogen / (double) CYCLE) *
                          1000.0; // PRDF100N
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PROD10C]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD10.carbon; // PROD10C
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PROD10N]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD10.nitrogen * 1000.0; // PROD10N
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PROD100C]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD100.carbon; // PROD100C
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PROD100N]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD100.nitrogen * 1000.0; // PROD100N
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRD10FC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
		 tempred.mon[dm] = tem.ag.PROD10decay.carbon / (double) CYCLE;
		 // PRD10FC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRD10FN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.PROD10decay.nitrogen / (double) CYCLE) *
                          1000.0; // PRD10FN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRD100FC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
		 tempred.mon[dm] = tem.ag.PROD100decay.carbon / (double) CYCLE;
		 // PRD100FC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PRD100FN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.PROD100decay.nitrogen / (double) CYCLE) *
                          1000.0; // PRD100FN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGNPPC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.npp[dm].carbon; // AGNPPC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGNPPN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.npp[dm].nitrogen * 1000.0; // AGNPPN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGFPRDC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
		 tempred.mon[dm] = tem.ag.formPROD1.carbon / (double) CYCLE;
		 // AGFPRODC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGFPRDN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.formPROD1.nitrogen / (double) CYCLE) *
                          1000.0; // AGFPRODN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGFRTN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.fertn[dm] * 1000.0; // AGFERTN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGLTRC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.ltrfal[dm].carbon; // AGLTRFC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGLTRN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.ltrfal[dm].nitrogen * 1000.0; // AGLTRFN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SLASHC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.slash[dm].carbon; // SLASHC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SLASHN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.slash[dm].nitrogen; // SLASHN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGPRDC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD1.carbon; // AGPRODC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGPRDN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.PROD1.nitrogen * 1000.0; // AGPRODN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGPRDFC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
		 tempred.mon[dm] = tem.ag.PROD1decay.carbon / (double) CYCLE;
		 // AGPRODFC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AGPRDFN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.PROD1decay.nitrogen / (double) CYCLE) *
                          1000.0; // AGPRODFN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTFPRDC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
		 tempred.mon[dm] = tem.ag.formTOTPROD.carbon / (double) CYCLE;
		 // TOTFPRDC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTFPRDN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.formTOTPROD.nitrogen / (double) CYCLE) *
                          1000.0; // TOTFPRDN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTPRDC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.TOTPROD.carbon; // TOTPRODC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTPRDN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.TOTPROD.nitrogen * 1000.0; // TOTPRODN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTPRDFC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
		 tempred.mon[dm] = tem.ag.TOTPRODdecay.carbon / (double) CYCLE;
		 // TOTPRDFC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTPRDFN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = (tem.ag.TOTPRODdecay.nitrogen / (double) CYCLE) *
                          1000.0; // TOTPRDFN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTEC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.totalc[dm]; // TOTEC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTNPP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.totnpp[dm]; // TOTNPP
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTGC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
		 tempred.mon[dm] = tem.totalc[dm] + tem.ag.TOTPROD.carbon;
		 // TOTGC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_CFLX]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.cflux[dm]; // CFLUX
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SCNVRTC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.sconvrtflx[dm].carbon; // SCONVRTC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SCNVRTN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.sconvrtflx[dm].nitrogen; // SCONVRTN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_NSRTNT]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.ag.nsretent[dm]; // NSRETENT
      }
    }
     else if (strcmp(predname[k],tem.predstr[tem.I_VSM]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.vsm[dm] * 100.0;      // VSM
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_PET]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.atms.pet[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_EET]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.atms.eet[dm];              // EET
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RAIN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.atms.rain[dm];             // RAIN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SNWINF]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.snowinf[dm];          // SNOWINF
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_WYLD]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.h2oyld[dm];          // H2OYIELD
      }
    }

       // for soil temperature
    else if (strcmp(predname[k],tem.predstr[tem.I_TSOIL]) == 0)
     {
      for (dm = 0; dm < CYCLE; dm++)
         {
	        tempred.mon[dm] = tem.atms.tsoil[dm];          // Tsoil
         }
      }
    else if (strcmp(predname[k],tem.predstr[tem.I_DST5]) == 0)
     {
      for (dm = 0; dm < CYCLE; dm++)
         {
	        tempred.mon[dm] = tem.atms.dst5[dm];          // Tsoil
         }
      }
   else if (strcmp(predname[k],tem.predstr[tem.I_DST10]) == 0)
     {
      for (dm = 0; dm < CYCLE; dm++)
         {
	        tempred.mon[dm] = tem.atms.dst10[dm];          // Tsoil
         }
      }
   else if (strcmp(predname[k],tem.predstr[tem.I_DST20]) == 0)
     {
      for (dm = 0; dm < CYCLE; dm++)
         {
	        tempred.mon[dm] = tem.atms.dst20[dm];          // Tsoil
         }
      }
   else if (strcmp(predname[k],tem.predstr[tem.I_DST50]) == 0)
     {
      for (dm = 0; dm < CYCLE; dm++)
         {
	        tempred.mon[dm] = tem.atms.dst50[dm];          // Tsoil
         }
      }
   else if (strcmp(predname[k],tem.predstr[tem.I_DST100]) == 0)
     {
      for (dm = 0; dm < CYCLE; dm++)
         {
	        tempred.mon[dm] = tem.atms.dst100[dm];          // Tsoil
         }
      }
   else if (strcmp(predname[k],tem.predstr[tem.I_DST200]) == 0)
     {
      for (dm = 0; dm < CYCLE; dm++)
         {
	        tempred.mon[dm] = tem.atms.dst200[dm];          // Tsoil
         }
      }
   else if (strcmp(predname[k],tem.predstr[tem.I_FRONTD]) == 0)
     {
      for (dm = 0; dm < CYCLE; dm++)
         {
	        tempred.mon[dm] = tem.atms.frontd[dm];          // Tsoil
         }
      }
   else if (strcmp(predname[k],tem.predstr[tem.I_THAWBE]) == 0)
     {
      for (dm = 0; dm < CYCLE; dm++)
         {
	        tempred.mon[dm] = tem.atms.thawbe[dm];          // Tsoil
         }
      }
  else if (strcmp(predname[k],tem.predstr[tem.I_THAWEND]) == 0)
       {
      for (dm = 0; dm < CYCLE; dm++)
         {
	        tempred.mon[dm] = tem.atms.thawend[dm];          // Tsoil
         }
      }
   // end of soil temperature


    else if (strcmp(predname[k],tem.predstr[tem.I_PCTP]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.pctp[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SM]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.moist[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_AVLW]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.avlh2o[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RRUN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.rrun[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SNWFAL]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.atms.snowfall[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SNWPCK]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.snowpack[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SRUN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.srun[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RGRW]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.rgrndh2o[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SGRW]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.sgrndh2o[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_RPERC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.rperc[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_SPERC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.soil.sperc[dm];
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_TOTC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
	     tempred.mon[dm] = tem.totalc[dm];  // TOTALC
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_VEGN]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
 	     tempred.mon[dm] = tem.veg.plant[dm].nitrogen;  // VEGN
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_LAI]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
 	     tempred.mon[dm] = tem.veg.lai[dm];  // LAI
      }
    }
    else if (strcmp(predname[k],tem.predstr[tem.I_FPC]) == 0)
    {
      for (dm = 0; dm < CYCLE; dm++)
      {
 	     tempred.mon[dm] = tem.veg.fpc[dm] * 100.0;  // FPC
      }
    }

    // Write output data to files

    if (strcmp(predname[k],tem.predstr[tem.I_VSM]) == 0
        || strcmp(predname[k],tem.predstr[tem.I_PCTP]) == 0
        || strcmp(predname[k],tem.predstr[tem.I_FPC]) == 0
        || strcmp(predname[k],tem.predstr[tem.I_LEAF]) == 0)
    {
      tempred.poutdel(fout[i], col, row, predname[k],
                     tem.veg.temveg, tem.veg.cmnt,
                     (100.0 * tem.soil.psiplusc),
                     tem.qualcon[dyr][itype],
                     carea,
                     ttotyr[dyr][itype],
                     tempred.mon,
                     contnent);
    }
    else
    {
      tempred.outdel(fout[i], col, row, predname[k],
                     tem.veg.temveg, tem.veg.cmnt,
                     (100.0 * tem.soil.psiplusc),
                     tem.qualcon[dyr][itype],
                     carea,
                     ttotyr[dyr][itype],
                     tempred.mon,
                     contnent);
    }
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::transqc(int& maxyears, long& totyr, Biomass plant[CYCLE])
{

  int i;
  int qc;
  double sumcarbon = 0.0;
  qc = ACCEPT;

  if (totyr < 0 || totyr >= (long) maxyears) { return qc = 30; }
  for (i = 0; i < CYCLE; i++) { sumcarbon += plant[i].carbon; }
  // "qc = 31" changed to "qc = TQCZEROFLAG" by DWK on 20000719
  if (sumcarbon <= 0.0) { return qc = TQCZEROFLAG; }

  return qc;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

int TEMelmnt::veggisin(ofstream& flog1, int& ftlerr, const int& atmsflag,
		       const int& temflag, FILE* ftveg)
{

  int gisend;

  Vegdata tveg;

  gisend = tveg.getdel(ftveg);
  if (gisend == -1) { return gisend; }

  if (atmsflag == 0)
  {
    col = tveg.col;
    row = tveg.row;
  }
  else
  {
    ftlerr = coregerr(flog1,"CLOUDS", col, row, "TEMVEG", tveg.col, tveg.row);
  }

  strcpy(contnent, tveg.contnent);
  mez = tveg.temveg - 1;
  if (temflag == 1)
  {
    maxtype = tem.veg.numtype[mez];
    tem.veg.temveg = tveg.temveg;
  }
  if (mez < 0 || mez >= NUMVEG)
  {
    mez = 0;  // Added by DWK on 20000219
    maxtype = 1;
  }

  return gisend;

};

