/* *************************************************************
****************************************************************
PCTEM423.CPP - Class adds DOS-specific interface to the core
                 Terrestrial Ecosystem Model Version 4.2 to allow
                 calibration of the model on a PC

Modifications:

20000103 - DWK added compiler directives
20000615 - DWK changes y[2] to y[I_SOLC] and y[3] to y[I_SOLN]
           at bottom of pcstepyr()
20010610 - QZ added STM

***************************************************************
************************************************************** */

#if !defined(PCTEM423_H)
  #include "pctem423e.hpp"
#endif

/* *********************************************************** */

PCTTEM::PCTTEM() : TTEM()
{

  strcpy(soilfile,"tsoil423.ecd");
  strcpy(rootfile,"trotz423.ecd");
  strcpy(vegfile,"tveg423.ecd");
  strcpy(leaffile,"tleaf423.ecd");
  strcpy(mcrvfile,"tmcrv423.ecd");
  strcpy(agfile,"ag423.ecd");

  //added for soil thermal model
  strcpy(snowfile,"qtsp44a2.ecd");
  strcpy(soillfile,"qtla44a2.ecd");
  strcpy(soiltfile,"qtst44a2.ecd");

  DisplayFlag = true;
  rheqflag = 0;
  sey[0] = GET_LAI;
  sey[1] = GET_NPP;
  sey[2] = GET_VNUP;

  swy[0] = GET_RAIN;
  swy[1] = GET_SNWINF;
  swy[2] = GET_VSM;
  swy[3] = GET_PET;
  swy[4] = GET_EET;

//  flog2.open("kick.log", ios::out);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::displayOptionalEflx(seykey& s)
{
  switch(s)
  {
    case GET_GPP:  cout << "   GPP  "; break;
    case GET_GPR:  cout << "    RA  "; break;
    case GET_LTRC:  cout << "   LTRC "; break;
    case GET_RH:  cout << "    RH  "; break;
    case GET_LTRN:  cout << "   LTRN "; break;
    case GET_NMIN:  cout << "   NMIN "; break;
    case GET_NLST:  cout << "  NLOST "; break;
    case GET_NINP:  cout << " NINPUT "; break;
    case GET_NEP:  cout << "   NEP  "; break;
    case GET_INGPP: cout << "  INGPP "; break;
    case GET_INNPP: cout << "  INNPP "; break;
    case GET_INNUP: cout << "  INNUP "; break;
    case GET_VEGC: cout << "   VEGC "; break;
    case GET_STRN: cout << " STRUCN "; break;
    case GET_AGNPPC: cout << " AGNPPC "; break;
    case GET_NMBL: cout << " NMOBIL "; break;
    case GET_NRSRB: cout << " NRTRAN "; break;
    case GET_NPP: cout << "   NPP  "; break;
    case GET_VNUP: cout << " UPTAKE "; break;
    case GET_VSUP: cout << " SUPTAK "; break;
    case GET_VLUP: cout << " LUPTAK "; break;
    case GET_POTNPP: cout << " POTNPP "; break;
    case GET_L2SN: cout << "   LCON "; break;
    case GET_FPC: cout << "    FPC "; break;
    case GET_LAI: cout << "    LAI "; break;
    case GET_LEAF: cout << "   LEAF "; break;
    case GET_AGNPPN: cout << " AGNPPN "; break;
    case GET_AGFRTN: cout << "AGFERTN "; break;

    case GET_TSOIL: cout << " TSOIL "; break; // added for soil thermal model
    case GET_FRONTD: cout << " FRONTD "; break; // added for soil thermal model

  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::displayOptionalWflx(swykey& s)
{
  switch (s)
  {
    case GET_RAIN:  cout << "   RAIN"; break;
    case GET_RPERC:  cout << "  RPERC"; break;
    case GET_RRUN:  cout << "   RRUN"; break;
    case GET_SNWFALL:  cout << " SNWFAL"; break;
    case GET_SNWINF:  cout << " SNWINF"; break;
    case GET_SPERC:  cout << "  SPERC"; break;
    case GET_SRUN:  cout << "   SRUN"; break;
    case GET_PET:  cout << "    PET"; break;
    case GET_EET:  cout << "    EET"; break;
    case GET_SH2O: cout << " SMOIST"; break;
    case GET_PCTP: cout << "   PCTP"; break;
    case GET_VSM: cout << "    VSM"; break;
    case GET_WYLD: cout << "   WYLD"; break;

                    //added for soil thermal model
  	 case GET_TOP20ST: cout << "TOP20ST"; break;
    case GET_DST5: cout << "   DST5"; break;
 	 case GET_DST10: cout << "  DST10"; break;
 	 case GET_DST20: cout << "  DST20"; break;
 	 case GET_DST50: cout << "  DST50"; break;
 	 case GET_DST100: cout << " DST100"; break;
 	 case GET_DST200: cout << " DST200"; break;
    case GET_FFRONTD: cout << " FFRONTD"; break;
    case GET_TTHAWB: cout << " TTHAWB"; break;
    case GET_TTHAWE: cout << " TTHAWE"; break;

  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::displayState(const int& dyr, const int& dm,
                          char* Climit, double pstate[])
{
  printf("%3d-%2d %8.2lf %7.2lf %9.2lf %7.2lf",
		  dyr,(dm+1),
		  pstate[I_VEGC],
		  pstate[I_STRN],
		  pstate[I_SOLC],
		  pstate[I_SOLN]);
  printf(" %6.3lf %6.3lf %7.3lf %6.2lf %s %7.3lf",
		  pstate[I_AVLN],
		  pstate[I_STON],
		  outflux1,
		  outflux2,
		  Climit,
          outflux3);

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

double PCTTEM::getOptionalEflx(const int& dm, const int& optflx,
                               double pstate[])
{

  double outflux;

  switch (optflx)
  {
    case GET_GPP:    outflux = pstate[I_GPP]; break;
    case GET_GPR:    outflux = pstate[I_GPR]; break;
    case GET_LTRC:   outflux = pstate[I_LTRC]; break;
    case GET_RH:     outflux = pstate[I_RH]; break;
    case GET_LTRN:   outflux = pstate[I_LTRN]; break;
    case GET_NMIN:   outflux = pstate[I_NMIN]; break;
    case GET_NLST:   outflux = pstate[I_NLST]; break;
    case GET_NINP:   outflux = pstate[I_NINP]; break;
    case GET_NEP:    outflux = pstate[I_NEP]; break;
    case GET_INGPP:  outflux = pstate[I_INGPP]; break;
    case GET_INNPP:  outflux = pstate[I_INNPP]; break;
    case GET_INNUP:  outflux = pstate[I_INNUP]; break;
    case GET_VEGC:   outflux = pstate[I_VEGC]; break;
    case GET_STRN:   outflux = pstate[I_STRN]; break;
    case GET_AGNPPC: outflux = pstate[I_AGNPPC]; break;
    case GET_NMBL:   outflux = pstate[I_VNMBL]; break;
    case GET_NRSRB:  outflux = pstate[I_VNRSRB]; break;
    case GET_NPP:    outflux = pstate[I_NPP]; break;
    case GET_VNUP:   outflux = pstate[I_VNUP]; break;
    case GET_VSUP:   outflux = pstate[I_VSUP]; break;
    case GET_VLUP:   outflux = pstate[I_VLUP]; break;
    case GET_POTNPP: outflux = ag.potnpp[dm]; break;
    case GET_L2SN:   if (pstate[I_STRN] != 0.0)
                     {
                       outflux = pstate[I_STON]/pstate[I_STRN];
                     }
                     else { outflux = MISSING; } break;
    case GET_FPC:    outflux = veg.fpc[dm]; break;
    case GET_LAI:    outflux = veg.lai[dm]; break;
    case GET_LEAF:   outflux = pstate[I_LEAF]; break;
    case GET_AGNPPN: outflux = pstate[I_AGNPPN]; break;
    case GET_AGFRTN: outflux = pstate[I_AGFRTN]; break;

        // added for soil thermal model
    case GET_TSOIL: outflux = pstate[I_TSOIL]; break;
    case GET_FRONTD: outflux = pstate[I_FRONTD]; break;

    default:         outflux = MISSING;
  }

  return outflux;

};

/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

double PCTTEM::getOptionalWflx(const int& optflx, double pstate[])
{

  double outflux;

  switch (optflx)
  {
    case GET_RAIN:    outflux = pstate[I_RAIN]; break;
    case GET_RPERC:   outflux = pstate[I_RPERC]; break;
    case GET_RRUN:    outflux = pstate[I_RRUN]; break;
    case GET_SNWFALL: outflux = pstate[I_SNWFAL]; break;
    case GET_SNWINF:  outflux = pstate[I_SNWINF]; break;
    case GET_SPERC:   outflux = pstate[I_SPERC]; break;
    case GET_SRUN:    outflux = pstate[I_SRUN]; break;
    case GET_PET:     outflux = pstate[I_PET]; break;
    case GET_EET:     outflux = pstate[I_EET]; break;
    case GET_SH2O:    outflux = pstate[I_SM]; break;
    case GET_PCTP:    outflux = pstate[I_PCTP]; break;
    case GET_VSM:     outflux = 100.0*pstate[I_VSM]; break;
    case GET_WYLD:    outflux = pstate[I_WYLD]; break;

    // addition for soil thermal model
  	 case GET_TOP20ST: outflux = pstate[I_TSOIL]; break;
    case GET_DST5:    outflux = pstate[I_DST5]; break;
 	 case GET_DST10:   outflux = pstate[I_DST10]; break;
 	 case GET_DST20:   outflux = pstate[I_DST20]; break;
 	 case GET_DST50:   outflux = pstate[I_DST50]; break;
 	 case GET_DST100:  outflux = pstate[I_DST100]; break;
 	 case GET_DST200:  outflux = pstate[I_DST200]; break;
    case GET_FFRONTD: outflux = pstate[I_FRONTD]; break;
    case GET_TTHAWB:  outflux = pstate[I_THAWBE]; break;
    case GET_TTHAWE:  outflux = pstate[I_THAWEND]; break;

    default:          outflux = MISSING;
  }

  return outflux;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int PCTTEM::pcadapt(const int& numeq, double pstate[],
                    double& ptol, const int& pdm)
{

  int i;
  double ipart;
  double fpart;
  double time = 0.0;
  double dt = 1.0;
  int mflag = 0;
  long nintmon = 0;
  double oldstate[NUMEQ];

  blackhol = 0;
  while (time != 1.0)
  {
    test = REJECT;
    if (syint == 1)
    {
      while (test != ACCEPT)
      {
	rkf(numeq,pstate,dt,pdm);
	test = boundcon(dum4,error,ptol);
	if (test != ACCEPT)
        {
          // Display ODE errors to DOS screen

          pcdisplayODEerr(test, pstate);
        }

        if (dt <= pow(0.5,maxit))
        {
	  test = ACCEPT;
	  mflag = 1;
          if (nintmon == 0)
          {
            for(i = 0; i < numeq;i++) { oldstate[i] = pstate[i]; }
          }
	  ++nintmon;
        }

        if (test == ACCEPT)
        {
          for(i = 0; i < numeq;i++) { pstate[i] = dum4[i]; }
          time += dt;

          // Display time updates to the DOS screen

          pcdisplayDT(time, dt);

          fpart = modf((0.01 + (time/(2.0*dt))),&ipart);
          if ( fpart < 0.1 && dt < 1.0) { dt *= 2.0; }
        }
        else { dt *= 0.500; }

        if (nintmon == maxitmon)
        {
          time = 1.0;
          blackhol = 1;
          for(i = 0; i < numeq;i++) { pstate[i] = oldstate[i]; }
        }
      } /* end test while */
    }   /* end rkf integrator */
  }     /* end time while */

  return mflag;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayClm(char* calclm, double clouds[CYCLE])
{
  int dm;

  double yrclds = 0;
  double aveclds;
  double yrprec = 0;
  double yrtair = 0;

  for (dm = 0; dm < CYCLE; dm++)
  {
    yrclds += clouds[dm];
    yrtair += atms.tair[dm];
    yrprec += atms.prec[dm];
  }

  aveclds = yrclds / 12.0;
  atms.avetair = yrtair / 12.0;

  window(1,1,80,25);
  gotoxy(1,1);
  clrscr();
  cout << endl << "                ECD Climate File = "
  << calclm << endl << endl;
  cout << "           MONTH      CO2       CLDS    TAIR    PRECIP"
  	<< endl << endl;
  for (dm = 0; dm < CYCLE; dm++)
  {
    printf("            %3d     %6.1lf    %6.1lf  %6.1lf   %6.1lf\n",
          (dm+1),atms.co2[dm],clouds[dm],atms.tair[dm],atms.prec[dm]);
  }
  printf("\n        AVERAGE               %6.1f  %6.1lf\n",
        aveclds, atms.avetair);
  printf("                                ANNUAL PREC = %8.2lf mm/yr\n\n",
        yrprec);

  cout << endl << "Press any key to continue . . ." << endl;
  while (getch() == '\0');

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayDT(const double& tottime, const double& deltat)
{

  window(1,15,39,15);
  gotoxy(1,1);
  printf("TIME = %10.8lf   DT = %10.8lf    ", tottime, deltat);

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayInitState(char* ecd, const double& col, const double& row)
{

  window(1,1,80,25);
  gotoxy(1,1);
  clrscr();
  cout << "                      ECD File = " << ecd << endl << endl;
  printf("  COL = %6.1lf          ROW = %6.1lf          ELEV = %6.0lf\n\n",
        col,row,elev);
  cout << veg.cmnt_name;
  printf("    CMNT = %d        WETSOIL = %d\n\n", veg.cmnt,soil.wsoil);

  soil.showecd();



  cout << endl;
  cout << endl;
  cout << "                PARAMETERS FOR STATE VARIABLES" << endl << endl;
  cout << "           PLANTC  STRUCTN     SOILC    SOILN   AVAILN  LABILN";
  cout << endl;
  printf("YICA    %9.2lf %8.2lf %9.2lf %8.2lf %8.3lf %7.3lf\n",
        vegca[veg.cmnt],strna[veg.cmnt],
        solca[veg.cmnt],solna[veg.cmnt],
        avlna[veg.cmnt],stona[veg.cmnt]);
  printf("YICB    %9.2lf %8.2lf %9.2lf %8.2lf %8.3lf %7.3lf\n",
        vegcb[veg.cmnt],strnb[veg.cmnt],
        solcb[veg.cmnt],solnb[veg.cmnt],
        avlnb[veg.cmnt],stonb[veg.cmnt]);

  cout << endl << endl <<  "Press any key to continue . . ." << endl;
  while (getch() == '\0');

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayLimits()
{
  window(1,24,45,24);
  gotoxy(1,1);
  delline ();
  if (veg.yrinnpp > veg.yrnpp && veg.yrinnup == veg.yrnup)
  {
    cout << "NPP LIMITED BY N UPTAKE, C UPTAKE RESTRICTED ";
  }
  if (veg.yrinnpp == veg.yrnpp && veg.yrinnup == veg.yrnup)
  {
    cout << "NPP NOT LIMITED BY EITHER C OR N UPTAKE      ";
  }
  if (veg.yrinnpp == veg.yrnpp && veg.yrinnup > veg.yrnup)
  {
    cout << "NPP LIMITED BY C UPTAKE, N UPTAKE RESTRICTED ";
  }
  if (veg.yrinnpp > veg.yrnpp && veg.yrinnup > veg.yrnup)
  {
    cout << "NPP LIMITED BY BOTH C AND N UPTAKE           ";
  }
  if (veg.yrinnpp < veg.yrnpp || veg.yrinnup < veg.yrnup)
  {
    cout << "ERROR-NPP OR N UPTAKE > INITIAL VALUES       ";
  }

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayMonth(const int& dyr, const int& dm, double pstate[])
{

  char limitC[2];

  if (topwind == 0)
  {
     window(1,1,80,24);
     clrscr();
     if (calwind == 1)
     {
        window(1,1,43,1);
        cout << "    TIME   AVAILW   RGRNDW  SNOWPCK SRGRNDW";
    }
    else
    {
       window(1,1,49,1);
       cout << "   TIME   VEG C    STR N   SOIL C   SOIL N  AVALN LABILN";
    }

    topwind = 1;
  }

  // Assign values of optional variables to outflux? for later screen display

  outflux1 = getOptionalEflx(dm, sey[0], pstate);

  outflux2 = getOptionalEflx(dm, sey[1], pstate);

  outflux3 = getOptionalEflx(dm, sey[2], pstate);

  outflux4 = getOptionalWflx(swy[0], pstate);

  outflux5 = getOptionalWflx(swy[1], pstate);

  outflux6 = getOptionalWflx(swy[2], pstate);

  outflux7 = getOptionalWflx(swy[3], pstate);

  outflux8 = getOptionalWflx(swy[4], pstate);

  window(1,2,80,13);
  gotoxy(1,1);
  delline();
  gotoxy(1,12);

  // Display monthly values for selected C and N pools and fluxes

  if (calwind == 1)
  {
	printf("%4d-%2d %9.2lf %8.2lf %8.2lf %7.2lf",
			dyr,(dm+1),
			pstate[I_AVLW],
			pstate[I_RGRW],
			pstate[I_SNWPCK],
			pstate[I_SGRW]);
	printf(" %6.1lf %6.1lf %6.1lf %6.1lf %7.1lf",
			outflux4,
			outflux5,
            outflux6,
            outflux7,
            outflux8);
  }
  else
  {

    // Productivity is nitrogen limited

    if ((pstate[I_INNPP] > pstate[I_NPP])
       && (pstate[I_INNUP] == pstate[I_VNUP]))
    {
      strcpy(limitC, "N");
    }

    // Productivity is either carbon or nitrogen limited

    if ((pstate[I_INNPP] == pstate[I_NPP])
       && (pstate[I_INNUP] == pstate[I_VNUP]))
    {
      strcpy(limitC, "E");
    }

    // Productivity is carbon limited (climate)

    if ((pstate[I_INNPP] == pstate[I_NPP])
       && (pstate[I_INNUP] > pstate[I_VNUP]))
    {
      strcpy(limitC, "C");
    }

    // Productivity is limited by both carbon and nitrogen

    if ((pstate[I_INNPP] > pstate[I_NPP])
       && (pstate[I_INNUP] > pstate[I_VNUP]))
    {
      strcpy(limitC, "B");
    }

    // Unknown limits on productivity

    if (pstate[I_INNPP] < pstate[I_NPP]
       || pstate[I_INNUP] < pstate[I_VNUP])
//    else
    {
      strcpy(limitC, "?");
    }

    displayState(dyr, dm, limitC, pstate);
  }

  window(1,14,80,14);
  gotoxy(1,1);
  delline();
//  printf("                                                                         ");

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayODEerr(const int& test, double pstate[])
{
  int errtest;

  window(40,15,80,15);
  gotoxy(1,1);
  if (test != ACCEPT)
  {
    errtest = test - 1;
    switch(errtest)
    {
      case   I_VEGC: printf("   VEGC = %8.2lf  Error = %11.8lf  ",
                           pstate[I_VEGC],error[I_VEGC]);
                     break;
      case   I_STRN: printf("   STRN = %8.2lf  Error = %11.8lf  ",
                           pstate[I_STRN],error[I_STRN]);
                     break;
      case   I_SOLC: printf("  SOILC = %8.2lf  Error = %11.8lf  ",
                           pstate[I_SOLC],error[I_SOLC]);
                     break;
      case   I_SOLN: printf("  SOILN = %8.2lf  Error = %11.8lf  ",
                           pstate[I_SOLN],error[I_SOLN]);
                     break;
      case   I_AVLN: printf(" AVAILN = %8.5lf  Error = %11.8lf  ",
                           pstate[I_AVLN],error[I_AVLN]);
                     break;
      case   I_STON: printf("LABILEN = %8.3lf  Error = %11.8lf  ",
                           pstate[I_STON],error[I_STON]);
                     break;
      case    I_GPP: printf("    GPP = %8.1lf  Error = %11.8lf  ",
                           pstate[I_GPP],error[I_GPP]);
                     break;
      case    I_NPP: printf("    NPP = %8.1lf  Error = %11.8lf  ",
                           pstate[I_NPP],error[I_NPP]);
                     break;
      case   I_VNUP: printf(" UPTAKE = %8.3lf  Error = %11.8lf  ",
                           pstate[I_VNUP],error[I_VNUP]);
                     break;
      case    I_GPR: printf("     RA = %8.1lf  Error = %11.8lf  ",
                           pstate[I_GPR],error[I_GPR]);
                     break;
      case   I_LTRC: printf("   LTRC = %8.1lf  Error = %11.8lf  ",
                           pstate[I_LTRC],error[I_LTRC]);
                     break;
      case     I_RH: printf("     RH = %8.1lf  Error = %11.8lf  ",
                           pstate[I_RH],error[I_RH]);
                     break;
      case   I_LTRN: printf("   LTRN = %8.3lf  Error = %11.8lf  ",
                           pstate[I_LTRN],error[I_LTRN]);
                     break;
      case   I_NMIN: printf("   NMIN = %8.3lf  Error = %11.8lf  ",
                           pstate[I_NMIN],error[I_NMIN]);
                     break;
      case   I_NLST: printf("  NLOST = %8.3lf  Error = %11.8lf  ",
                           pstate[I_NLST],error[I_NLST]);
                     break;
      case   I_NINP: printf(" NINPUT = %8.3lf  Error = %11.8lf  ",
                           pstate[I_NINP],error[I_NINP]);
                     break;
      case  I_RVMNT: printf("     RM = %8.1lf  Error = %11.8lf  ",
                           pstate[I_RVMNT],error[I_RVMNT]);
                     break;
      case  I_RVGRW: printf("     RG = %8.1lf  Error = %11.8lf  ",
                           pstate[I_RVGRW],error[I_RVGRW]);
                     break;
      case   I_MNUP: printf("MCRONUP = %8.3lf  Error = %11.8lf  ",
                           pstate[I_MNUP],error[I_MNUP]);
                     break;
      case    I_NEP: printf("    NEP = %8.3lf  Error = %11.8lf  ",
                           pstate[I_NEP],error[I_NEP]);
                     break;
      case  I_INGPP: printf("INITGPP = %8.1lf  Error = %11.8lf  ",
                           pstate[I_INGPP],error[I_INGPP]);
                     break;
      case  I_INNPP: printf("INITNPP = %8.1lf  Error = %11.8lf  ",
                           pstate[I_INNPP],error[I_INNPP]);
                     break;
      case  I_INNUP: printf("INUPTAK = %8.3lf  Error = %11.8lf  ",
                           pstate[I_INNUP],error[I_INNUP]);
                     break;
      case  I_VNMBL: printf(" NMOBIL = %8.3lf  Error = %11.8lf  ",
                           pstate[I_VNMBL],error[I_VNMBL]);
                     break;
      case   I_VSUP: printf("SUPTAKE = %8.3lf  Error = %11.8lf  ",
                           pstate[I_VSUP],error[I_VSUP]);
                     break;
      case   I_VLUP: printf("LUPTAKE = %8.3lf  Error = %11.8lf  ",
                           pstate[I_VLUP],error[I_VLUP]);
                     break;
      case I_UNRMLF: printf(" UNLEAF = %8.3lf  Error = %11.8lf  ",
                           pstate[I_UNRMLF],error[I_UNRMLF]);
                     break;
      case   I_LEAF: printf("   LEAF = %8.3lf  Error = %11.8lf  ",
                           pstate[I_LEAF],error[I_LEAF]);
                     break;


      case   I_AVLW: printf(" AVAILW = %8.2lf  Error = %11.8lf  ",
                           pstate[I_AVLW],error[I_AVLW]);
                     break;
      case   I_RGRW: printf(" RGRNDW = %8.2lf  Error = %11.8lf  ",
                           pstate[I_RGRW],error[I_RGRW]);
                     break;
      case I_SNWPCK: printf(" SNWPCK = %8.2lf  Error = %11.8lf  ",
                           pstate[I_SNWPCK],error[I_SNWPCK]);
                     break;
      case   I_SGRW: printf(" SGRNDW = %8.2lf  Error = %11.8lf  ",
                           pstate[I_SGRW],error[I_SGRW]);
                     break;
      case     I_SM: printf(" SMOIST = %8.5lf  Error = %11.8lf  ",
                           pstate[I_SM],error[I_SM]);
                     break;
      case   I_PCTP: printf("   PCTP = %8.3lf  Error = %11.8lf  ",
                           pstate[I_PCTP],error[I_PCTP]);
                     break;
      case    I_VSM: printf("    VSM = %8.3lf  Error = %11.8lf  ",
                           pstate[I_VSM],error[I_VSM]);
                     break;
      case   I_RAIN: printf("   RAIN = %8.5lf  Error = %11.8lf  ",
                           pstate[I_RAIN],error[I_RAIN]);
                     break;
      case  I_RPERC: printf("  RPERC = %8.3lf  Error = %11.8lf  ",
                           pstate[I_RPERC],error[I_RPERC]);
                     break;
      case   I_RRUN: printf("   RRUN = %8.1lf  Error = %11.8lf  ",
                           pstate[I_RRUN],error[I_RRUN]);
                     break;
      case I_SNWFAL: printf("SNOWFAL = %8.1lf  Error = %11.8lf  ",
                           pstate[I_SNWFAL],error[I_SNWFAL]);
                     break;
      case I_SNWINF: printf("SNOWINF = %8.3lf  Error = %11.8lf  ",
                           pstate[I_SNWINF],error[I_SNWINF]);
                     break;
      case  I_SPERC: printf("  SPERC = %8.1lf  Error = %11.8lf  ",
                           pstate[I_SPERC],error[I_SPERC]);
                     break;
      case   I_SRUN: printf("   SRUN = %8.1lf  Error = %11.8lf  ",
                           pstate[I_SRUN],error[I_SRUN]);
                     break;
      case    I_PET: printf("    PET = %8.1lf  Error = %11.8lf  ",
                           pstate[I_PET],error[I_PET]);
                     break;
      case    I_EET: printf("    EET = %8.3lf  Error = %11.8lf  ",
                           pstate[I_EET],error[I_EET]);
                     break;
      case   I_WYLD: printf(" H2OYLD = %8.3lf  Error = %11.8lf  ",
                           pstate[I_WYLD],error[I_WYLD]);
                     break;
    }
  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayOtherECD()
{
  window(1,1,80,25);
  gotoxy(1,1);
  clrscr();

  veg.showecd(veg.cmnt);

  veg.showleaf(veg.cmnt);

  microbe.showecd(veg.cmnt);
  cout << endl << "\nPress any key to continue . . ." << endl;
  while (getch() == '\0');

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayPAR(char* calclm, double girr[CYCLE])
{

  int dm;
  float pday;
  double yrpar = 0;

  window(1,1,80,25);
  gotoxy(1,1);
  clrscr();
  cout << "                ECD Climate File = " << calclm << endl << endl;
  cout << "                     Irradiance Data" << endl << endl;
  cout << "           MONTH     GIRR     NIRR     PAR" << endl << endl;
  for(dm = 0; dm < CYCLE; dm++)
  {
    printf("            %3d     %6.1lf  %6.1lf   %6.1lf\n",
		  (dm+1),girr[dm],atms.nirr[dm],atms.par[dm]);

    switch(dm)
    {
      case 0:  pday = 31.0; break;
      case 1:  pday = 28.25; break;
      case 2:  pday = 31.0; break;
      case 3:  pday = 30.0; break;
      case 4:  pday = 31.0; break;
      case 5:  pday = 30.0; break;
      case 6:  pday = 31.0; break;
      case 7:  pday = 31.0; break;
      case 8:  pday = 30.0; break;
      case 9:  pday = 31.0; break;
      case 10: pday = 30.0; break;
      case 11: pday = 31.0; break;
    }
    yrpar += (atms.par[dm] * pday);
  }

  yrpar /= 1000.0;
  printf("\n                       ANNUAL PAR = %8.2lf kcal/cm2/yr\n\n", yrpar);

  cout << endl << "Press any key to continue . . ." << endl;
  while (getch() == '\0');

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayStext()
{

  eqsolc = (-0.827*atms.avetair) + (0.0224*pow(atms.avetair,2.0))
           + (0.127*(atms.yrrain+soil.yrsnowinf)/10.0);
  eqsolc += (-0.000938*pow((atms.yrrain+soil.yrsnowinf)/10.0,2.0))
            + (0.000899*soil.pctsilt*(atms.yrrain+soil.yrsnowinf)/10.0);
  eqsolc += (0.0006*soil.pctclay*(atms.yrrain+soil.yrsnowinf)/10.0) + 4.09;
  eqsolc *= 1000.0;

  window(1,1,80,25);
  gotoxy(1,1);
  clrscr();
  cout << "         INITIAL CONDITIONS AND CALCULATED PARAMETERS";
  cout << endl << endl;

  printf("PSIPLUSC = %6.4lf     ROOT DEPTH = %8.4lf\n\n",
        soil.psiplusc, soil.rootz);

  cout << "       PLANTC  STRUCTN     SOILC    SOILN   AVAILN  STOREN" << endl;
  printf("    %9.2lf %8.2lf %9.2lf %8.2lf %8.3lf %7.3lf\n\n\n",\
	  y[I_VEGC],y[I_STRN],y[I_SOLC],y[I_SOLN],y[I_AVLN],y[I_STON]);

  printf("CMAX =  %8.2lf           KDC =  %8.6lf\n", veg.cmax, microbe.kdc);
  printf("NMAX = %9.6lf            NUP = %9.6lf\n", veg.nmax, microbe.nup);


  cout << endl << endl <<  "Press any key to continue . . ." << endl;
  while (getch() == '\0');

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayVegECD()
{
  window(1,1,80,25);
  gotoxy(1,1);
  clrscr();
  cout << "              VEGETATION PARAMETERS FOR THE CARBON CYCLE";
  cout << endl << endl;

  printf("CMAXCUT = %4.2lf     C1A = %7.2lf    C1B = %7.2lf  ",
		veg.cmaxcut[veg.cmnt],veg.cmax1a[veg.cmnt],veg.cmax1b[veg.cmnt]);
  printf("C2A = %7.2lf  C2B = %7.2lf\n",
		veg.cmax2a[veg.cmnt],veg.cmax2b[veg.cmnt]);

  printf("  CFALL = %8.6lf     KRA = %10.7lf  KRB = %10.7lf\n",
        veg.cfall[veg.cmnt],veg.kra[veg.cmnt],veg.krb[veg.cmnt]);

  printf("    KDA = %8.6lf     KDB =  %8.6lf  LCCLNC = %8.6lf PFTOS = %6.4lf\n",
        microbe.kda[veg.cmnt],microbe.kdb[veg.cmnt],
        microbe.lcclnc[veg.cmnt],microbe.propftos[veg.cmnt]);

  cout << endl << "             VEGETATION PARAMETERS FOR THE NITROGEN CYCLE";
  cout << endl << endl;

  printf("NMAXCUT = %4.2lf      N1A = %8.4lf  N1B = %8.4lf  ",
		veg.nmaxcut[veg.cmnt], veg.nmax1a[veg.cmnt],veg.nmax1b[veg.cmnt]);
  printf("N2A = %8.4lf  N2B = %8.4lf\n",
		veg.nmax2a[veg.cmnt],veg.nmax2b[veg.cmnt]);

  printf("  NFALL = %8.6lf     NUPA = %8.4lf     NUPB = %8.4lf\n",
        veg.nfall[veg.cmnt],microbe.nupa[veg.cmnt],microbe.nupb[veg.cmnt]);

  printf("  NLOSS = %8.6lf  NFIXPAR = %8.6lf  LABNCON = %8.6lf\n",
        soil.nloss[veg.cmnt],microbe.nfixpar[veg.cmnt],veg.labncon[veg.cmnt]);

  cout << endl;
  cout << "             VEGETATION PARAMETERS FOR NITROGEN FEEDBACK";
  cout << endl << endl;
  printf("INITCNEVEN = %6.2lf       CNMIN = %6.2lf\n",
		veg.initcneven[veg.cmnt],veg.cnmin[veg.cmnt]);
  printf("PLANTCNA = %8.4lf  PLANTCNB = %8.4lf PLANTCNMIN = %8.4lf\n",
  		veg.c2na[veg.cmnt],veg.c2nb[veg.cmnt],veg.c2nmin[veg.cmnt]);

  cout << endl << "Press any key to continue . . ." << endl;
  while (getch() == '\0');

}

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayYrAgTEM(const int& manflag)
{
  window(40,16,80,16);
  gotoxy(1,1);
  printf("PROD10C = %8.2lf PROD10N = %7.3lf      ",
        ag.PROD10.carbon,ag.PROD10.nitrogen);

  window(40,17,80,17);
  gotoxy(1,1);
  printf("PRD100C = %8.2lf PRD100N = %7.3lf      ",
        ag.PROD100.carbon,ag.PROD100.nitrogen);

  window(40,18,80,18);
  gotoxy(1,1);
  printf("YRRH10C = %6.1lf   YRRH10N = %7.3lf     ",
        ag.PROD10decay.carbon,ag.PROD10decay.nitrogen);

  window(1,16,39,16);
  gotoxy(1,1);
  printf("AGFPRDC = %8.2lf  YRAVLN =  %8.3lf ",
        ag.formPROD1.carbon+ag.formPROD10.carbon+ag.formPROD100.carbon,
        soil.yravln);

  window(1,17,39,17);
  gotoxy(1,1);
  printf("MF = %1d     AGS = %1d   YRSOCN =  %7.3lf ",
        manflag,ag.state,soil.yrc2n);

  window(1,18,19,18);
  gotoxy(1,1);
  printf("YRSOLC = %9.2lf   YRSOLN =  %7.3lf ",soil.yrorgc,soil.yrorgn);

  window(40,19,73,19);
  gotoxy(1,1);
  printf("YRH100C = %6.1lf   YRH100N = %7.3lf      ",
        ag.PROD100decay.carbon,ag.PROD100decay.nitrogen);

  window(1,19,39,19);
  gotoxy(1,1);
  printf("YRAGPRDC = %7.2lf   YRAGPN =  %7.3lf ",
        ag.formPROD1.carbon,ag.formPROD1.nitrogen);

  window(40,20,73,20);
  gotoxy(1,1);
  printf(" YAGNPP = %6.1lf  YRAGNUP = %8.3lf     ",ag.yrnppC,ag.yrnppN);

  window(1,20,39,20);
  gotoxy(1,1);
  printf("YRNRETN  = %8.3lf  YRCONVN = %7.3lf ",ag.yrnrent,ag.yrconvrtN);

  window(40,21,73,21);
  gotoxy(1,1);
  printf(" YAGLTR = %6.1lf YRAGLTRN = %8.3lf     ",ag.yrltrc,ag.yrltrn);

  window(40,22,73,22);
  gotoxy(1,1);
  printf("   YRRH = %6.1lf   YRNMIN = %8.3lf     ",
        microbe.yrrh,microbe.yrnmin);

  window(1,21,39,21);
  gotoxy(1,1);
  printf("YRAGFERT = %8.3lf  YRNLOST = %7.3lf ",ag.yrfertn,soil.yrnlost);

  window(1,22,39,22);
  gotoxy(1,1);
  printf("YRNIN  =   %8.3lf YRENLOST = %7.3lf ",
        soil.yrnin,(soil.yrnlost+ag.formPROD1.nitrogen+ag.yrconvrtN));

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayYrTEM(const int& manflag)
{
  window(40,16,80,16);
  gotoxy(1,1);
  printf("  YRSUP = %8.3lf  YRLUP = %8.3lf     ",veg.yrsup,veg.yrlup);

  window(40,17,80,17);
  gotoxy(1,1);
  printf("YRINGPP = %6.1lf YRNMOBIL = %8.3lf     ",veg.yringpp,veg.yrnmobil);

  window(40,18,80,18);
  gotoxy(1,1);
  printf("  YRGPP = %6.1lf YRNRSORB = %8.3lf     ",veg.yrgpp,veg.yrnrsorb);

  window(1,16,39,16);
  gotoxy(1,1);
  printf("YRVEGC = %9.2lf    YRVEGN =  %6.2lf ",veg.yrcarbon,veg.yrnitrogen);

  window(1,17,39,17);
  gotoxy(1,1);
  printf("MF = %1d     AGS = %1d    YRSTRN =  %6.2lf ",
        manflag,ag.state,veg.yrstructn);

  window(1,18,19,18);
  gotoxy(1,1);
  printf("YRSOLC = %9.2lf    YRLABN =  %7.3lf",soil.yrorgc,veg.yrstoren);

  window(40,19,73,19);
  gotoxy(1,1);
  printf("YRINNPP = %6.1lf   YRINUP = %8.3lf     ",veg.yrinnpp,veg.yrinnup);

  window(1,19,39,19);
  gotoxy(1,1);
  printf("YRVGCN = %9.2lf    YRSOLN = %7.2lf ",veg.yrc2n,soil.yrorgn);

  window(40,20,73,20);
  gotoxy(1,1);
  printf("  YRNPP = %6.1lf    YRNUP = %8.3lf     ",veg.yrnpp,veg.yrnup);

  window(1,20,39,20);
  gotoxy(1,1);
  printf("YRSOCN = %9.2lf    YRAVLN = %8.3lf",soil.yrc2n,soil.yravln);

  window(40,21,73,21);
  gotoxy(1,1);
  printf("  YRLTR = %6.1lf   YRLTRN = %8.3lf     ",veg.yrltrc,veg.yrltrn);

  window(40,22,73,22);
  gotoxy(1,1);
  printf("   YRRH = %6.1lf   YRNMIN = %8.3lf     ",
        microbe.yrrh,microbe.yrnmin);

  window(1,21,39,21);
  gotoxy(1,1);
  printf("YDECAY =   %8.3lf   FOLIAG = %8.3lf",microbe.decay,veg.foliage);

  window(1,22,39,22);
  gotoxy(1,1);
  printf("YRNIN  =   %8.3lf  YRNLOST = %8.3lf",soil.yrnin,soil.yrnlost);

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcdisplayYrWBM()
{

  window(42,16,80,16);
  gotoxy(1,1);
  printf(" YRPREC = %8.1lf  YRSFALL = %8.1lf  ",
        atms.yrrain+atms.yrsnowfall,atms.yrsnowfall);

  window(42,17,80,17);
  gotoxy(1,1);
  printf(" YRRAIN = %8.1lf YRSNWINF = %8.1lf  ",atms.yrrain,soil.yrsnowinf);

  window(42,18,80,18);
  gotoxy(1,1);
  printf("YRRPERC = %8.1lf  YRSPERC = %8.1lf  ",soil.yrrperc,soil.yrsperc);

  window(1,16,41,16);
  gotoxy(1,1);
  printf(" YRAVLW  = %8.1lf YRSNWPCK =  %8.1lf ",
        soil.yravlh2o,soil.yrsnowpack);

  window(1,17,41,17);
  gotoxy(1,1);
  printf("YRRGRNDW = %8.1lf YRSGRNDW =  %8.1lf ",
        soil.yrrgrndh2o,soil.yrsgrndh2o);

  window(1,18,41,18);
  gotoxy(1,1);
  printf("YRSMOIST = %8.1lf   YRPCTP =  %8.1lf",soil.yrsmoist,soil.yrpctp);

  window(42,19,73,19);
  gotoxy(1,1);
  printf(" YRRRUN = %8.1lf   YRSRUN = %8.1lf",soil.yrrrun,soil.yrsrun);

  window(1,19,41,19);
  gotoxy(1,1);
  printf("   YRVSM = %9.2lf                      ",soil.yrvsm*100.0);

  window(42,20,73,20);
  gotoxy(1,1);
  printf(" YRWYLD =   %6.1lf                    ",soil.yrh2oyld);

  window(1,20,41,20);
  gotoxy(1,1);
  printf("   YRPET = %9.2lf   YREET =  %8.1lf",atms.yrpet,atms.yreet);

  window(42,21,73,21);
  gotoxy(1,1);
  printf("YRULEAF = %8.1lf   YRLEAF = %8.1lf",
        veg.yrunleaf*100.0,veg.yrleaf*100.0);

  window(42,22,73,22);
  gotoxy(1,1);
  printf("                                       ");

  window(1,21,41,21);
  gotoxy(1,1);
  printf(" AWCAPMM = %8.1lf  ROOTZMM =  %8.1lf ",
        soil.awcapmm,soil.rootz*1000.0);

  window(1,22,41,22);
  gotoxy(1,1);
  printf("WILTPTMM = %8.1lf FLDCAPMM =  %8.1lf    ",
        soil.wiltpt,soil.fldcap);

  window(1,24,45,24);
  gotoxy(1,1);
  delline ();

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void PCTTEM::pcinitRun(ofstream& rflog1)
{

/* **************************************************************
	     Run Model with Nitrogen Limitation?
************************************************************** */

  cout << "Do you want to start calibration allowing"
  	   << " available N to fluctuate?  ";
  if (toupper(getch()) == 'Y')
  {
    avlnflag = 1;
    cout << "Y" << endl;
  }
  else 
  {
    avlnflag = 0;
    cout << "N" << endl;
  }

  baseline = 0;
  cout << "Do you want to start calibration with N feedback on GPP?  ";
  if (toupper(getch()) == 'Y')
  {
    nfeed = 1;
    cout << "Y" << endl;
    cout << "Do you want to solve for baseline soil nitrogen?  ";
    if (toupper(getch()) == 'Y')
    {
       baseline = 1;
       cout << "Y" << endl;
    }
    else { cout << "N" << endl; }
  }
  else 
  {
     nfeed = 0;
     cout << "N" << endl;
  }

/* **************************************************************
		 Run Model with Moisture Limitation?
************************************************************** */

  cout << "Do you want to start calibration with moisture limitation?  ";
  if (toupper(getch()) == 'Y') 
  {
    moistlim = 1;
    cout << "Y" << endl;
  }
  else 
  {  
     moistlim = 0;
     cout << "N" << endl;
  }

/* **************************************************************
		 Run Model Until Steady State Conditions?
************************************************************** */

  equil = 0;
  adapttol = 0;
  intbomb = 0;
  tolbomb = 0;
  cout << "Do you want the model to stop at steady state conditions?  ";
  if (toupper(getch()) == 'Y')
  {
    equil = 1;
    cout << "Y" << endl;
	cout << endl << "How many years do you want to wait before"
		 << " checking equilibrium conditions?  ";
//  cin >> strteq;
	fpara >> strteq;
    strteq *= 12;

    if (nfeed == 0) 
    {
       cout << "Do you want decomposition to come into equilibrium?  ";
       if (toupper(getch()) == 'Y') 
       {
	  rheqflag = 1;
	  cout << "Y" << endl;
       }
       else { cout << "N" << endl; }
    }

    cout << "Enter the absolute tolerance for the water cycle?  ";
//  cin >> wtol;
	fpara >> wtol;

    cout << "Enter the absolute tolerance for the carbon cycle?  ";
//  cin >> ctol;
	fpara >> ctol;

    if (nfeed == 1) 
    {
      rheqflag = 1;
      cout << "Enter the absolute tolerance for the nitrogen cycle?  ";
//    cin >> ntol;
      fpara >> ntol;
    }

	cout << "Do you want to have the integrator tolerance"
		 << " adjusted automatically?  ";
    if (toupper(getch()) == 'Y') 
    {
      cout << "Y" << endl;
      adapttol = 1;
      ask(rflog1);
	  cout << "Enter the maximum number of years for the model to run:  ";
//	  cin >> maxyears;
	  fpara >> maxyears;
	  cout << "Enter the maximum number of attempts to reach a solution:  ";
//    cin >> maxnrun;
	  fpara >> maxnrun;
    }
    else { cout << "N" << endl; }
  }
  else { cout << "N" << endl << endl; }

  if (adapttol == 0) { ask(rflog1); }

  cout << endl << "Enter the level of carbon dioxide in ppmv: ";
//cin >> atms.co2level;
  fpara >> atms.co2level;
  cout << endl << endl;

  cout << endl << "Enter the factor for changing C:N per ppmv"
	   << " of enhanced CO2: " << endl;
  cout << "               (Enter 0.0 for no change)" << endl;
//cin >> veg.dc2n;
  fpara >> veg.dc2n;
  cout << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

int PCTTEM::pcstepyr(const int& dyr, int& intflag, double& tol)
{

  int dm;
  int mintflag;

  if (dyr == 0) { microbe.kd = microbe.kdc; }
  else
  {
    if (ag.state == 0 && ag.prvstate == 0)
    {
      microbe.kd = microbe.yrkd(nfeed, veg.yrltrc, veg.yrltrn, veg.cmnt);
      ag.kd = microbe.kd;
      ag.natsoil = soil.org[CYCLE-1].carbon;
    }
    else
    {
      if (soil.org[CYCLE-1].carbon < ag.natsoil)
      {
        microbe.kd = ag.kd * soil.org[CYCLE-1].carbon / ag.natsoil;
      }
      else { microbe.kd = ag.kd; }
    }
  }

// Initialize annual rates

  resetYrFlux();

  if (ag.state == 1 && ag.prvstate == 0)
  {
    y[I_VEGC] = ag.vrespar[veg.cmnt] * veg.plant[CYCLE-1].carbon;
    y[I_STRN] = ag.vrespar[veg.cmnt] * veg.strctrl[CYCLE-1].nitrogen;
    y[I_STON] = ag.vrespar[veg.cmnt] * veg.labile[CYCLE-1].nitrogen;

    // Create PROD10 and PROD100 from conversion to agriculture

    ag.conversion(veg.cmnt, veg, soil);

  }

  // Revert to natural vegetation after cropland abandonment

  if (ag.state == 0 && ag.prvstate == 1)
  {
    veg.plant[CYCLE-1].carbon = y[I_VEGC];
    veg.strctrl[CYCLE-1].nitrogen = y[I_STRN];
    veg.labile[CYCLE-1].nitrogen = y[I_STON];
  }

  // Begin seasonal cycle

   for (dm = 0; dm < CYCLE; dm++)
  {

    // Get environmental conditions for month "dm"

	 atms.precsplt(atms.prec[dm], atms.tair[dm], atms.rain[dm],
	 		 atms.snowfall[dm]);

    switch (dm) {
      case 0:  atms.prevtair = atms.tair[CYCLE-1];
	       atms.prev2tair = atms.tair[CYCLE-2];
          break;
      case 1:  atms.prevtair = atms.tair[0];
	       atms.prev2tair = atms.tair[CYCLE-1];
	       break;
      default: atms.prevtair = atms.tair[dm-1];
	       atms.prev2tair = atms.tair[dm-2];
	       break;
      }
	soil.snowinf[dm] = soil.snowmelt(elev, atms.tair[dm], atms.prevtair,
		 y[I_SNWPCK]);
    soil.snowpack[dm]= atms.snowfall[dm] - soil.snowinf[dm];
    if (soil.snowpack[dm] < 0.0) { soil.snowpack[dm] = 0.0; }

   if (stmflg == 1) {
     switch (dm) {
	 case 0: sthermal.airt19= (atms.tair[CYCLE-1]+atms.tair[0])/2.0;
	 // satisfy soil thermal model
             sthermal.airt29= atms.tair[0];
             sthermal.airt39= (atms.tair[0]+atms.tair[1])/2.0;
             break;
	 case 11:sthermal.airt19= (atms.tair[10]+atms.tair[11])/2.0;
	  // satisfy soil thermal model
             sthermal.airt29= atms.tair[11];
             sthermal.airt39= (atms.tair[11]+atms.tair[0])/2.0;
            break;
	 default:sthermal.airt19= (atms.tair[dm-1]+atms.tair[dm])/2.0;
	 // satisfy soil thermal model
             sthermal.airt29= atms.tair[dm];
             sthermal.airt39= (atms.tair[dm]+atms.tair[dm+1])/2.0;
            break;
         }

    switch (dm) {
	 case 0: sthermal.hsnow19=(soil.snowpack[CYCLE-1]+soil.snowpack[0])
	 		/2.0/1000.0; // satisfy soil thermal model
             sthermal.hsnow29= soil.snowpack[0]/1000.0;
             sthermal.hsnow39= (soil.snowpack[1]+soil.snowpack[0])/2.0/1000.0;
             break;
	 case 11: sthermal.hsnow19=(soil.snowpack[11]+soil.snowpack[10])/2.0/1000.0;
	  // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[11]/1000.0;
              sthermal.hsnow39= (soil.snowpack[11]+soil.snowpack[0])/2.0/1000.0;
            break;
	 default: sthermal.hsnow19=(soil.snowpack[dm-1]+soil.snowpack[dm])
	 			/2.0/1000.0; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[dm]/1000.0;
			  sthermal.hsnow39= (soil.snowpack[dm]+soil.snowpack[dm+1])
			  /2.0/1000.0;
            break;
         }


  sthermal.calcu_cdsnow = 0.2;
   // dummy value, it will be changed in the STM use prescribe value

  sthermal.water2 = soil.pctp[dm]/100.0;
   // assume moss and organic have the same water content, 29/march/2001

   if (fswitch==0) {
		if (sthermal.soiltemp_(&sthermal.water2, &sthermal.calcu_cdsnow,
		 &sthermal.airt19, &sthermal.airt29, &sthermal.airt39,
		  &sthermal.hsnow19, &sthermal.hsnow29, &sthermal.hsnow39,
		  sthermal.weight9, &sthermal.smass9, &sthermal.is9, &sthermal.ism19,
		   sthermal.x9, sthermal.dx9, sthermal.xfb9, sthermal.xfa9,
		sthermal.water9, sthermal.t9, &sthermal.tsoil,&sthermal.frontd,
		&sthermal.thawbegin,&sthermal.thawend,sthermal.diffsoilt,veg.cmnt) !=0)
      { printf("bad tem"); getch();};

   atms.tsoil[dm]=sthermal.tsoil;

   ttsoil[dm]=sthermal.tsoil;
   ddst5[dm]=sthermal.diffsoilt[0];
   ddst10[dm]=sthermal.diffsoilt[1];
   ddst20[dm]=sthermal.diffsoilt[2];
   ddst50[dm]=sthermal.diffsoilt[3];
   ddst100[dm]=sthermal.diffsoilt[4];
   ddst200[dm]=sthermal.diffsoilt[5];
   ffrontd[dm]=sthermal.frontd;

   atms.frontd[dm]=sthermal.frontd;
   atms.dst5[dm]=sthermal.diffsoilt[0];
   atms.dst10[dm]=sthermal.diffsoilt[1];
   atms.dst20[dm]=sthermal.diffsoilt[2];
   atms.dst50[dm]=sthermal.diffsoilt[3];
   atms.dst100[dm]=sthermal.diffsoilt[4];
   atms.dst200[dm]=sthermal.diffsoilt[5];

   atms.frontd[dm]=ffrontd[dm];
   atms.thawbe[dm]=tthawb[dm];
   atms.thawend[dm]=tthawe[dm];

    }
    else {
     atms.tsoil[dm]=ttsoil[dm];
     atms.frontd[dm]=ffrontd[dm];
//     atms.thawbegin[dm]=tthawb[dm];
//     atms.thawend[dm]=tthawe[dm];
     atms.dst5[dm]=ddst5[dm] ;
     atms.dst10[dm]=ddst10[dm];
     atms.dst20[dm]=ddst20[dm];
     atms.dst50[dm]=ddst50[dm];
     atms.dst100[dm]=ddst100[dm];
     atms.dst200[dm]=ddst200[dm];
    }
  } // stmflg ===1


   // * calculate the thawing-proportion in early month
   // and later fall new algorithm, 28/Nov/2001
 if (atms.dst5[dm] ==0.0) atms.dst5[dm] =0.001;
 if (atms.dst5[dm-1] ==0.0) atms.dst5[dm-1] =0.002;
 if (atms.dst5[dm+1] ==0.0) atms.dst5[dm+1] =0.003;

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] > 0.0)
		&& (atms.dst5[dm+1] > 0.0))
	 veg.thawpercent[dm] = 1- (atms.dst5[dm]
		/ (atms.dst5[dm] - atms.dst5[dm-1]));

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] < 0.0)
	&& (atms.dst5[dm+1] < 0.0))
	 veg.thawpercent[dm] = 0.0;

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] < 0.0)
	&& (atms.dst5[dm+1] > 0.0))
	 veg.thawpercent[dm] = 1- (atms.dst5[dm+1]
	 / (atms.dst5[dm+1] - atms.dst5[dm]));

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] > 0.0)
	&& (atms.dst5[dm+1] < 0.0))
     veg.thawpercent[dm] = atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm+1]);

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] > 0.0)
	&& (atms.dst5[dm+1] > 0.0))
     veg.thawpercent[dm] = 1.0;

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] < 0.0)
	&& (atms.dst5[dm+1] < 0.0))
     veg.thawpercent[dm] = atms.dst5[dm-1] / (atms.dst5[dm-1] - atms.dst5[dm]);

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] > 0.0)
	&& (atms.dst5[dm+1] < 0.0))
	veg.thawpercent[dm] = ((1- (atms.dst5[dm]
		/ (atms.dst5[dm] - atms.dst5[dm-1])))
        + (atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm+1])))/2.0;

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] < 0.0)
	&& (atms.dst5[dm+1] > 0.0))
   veg.thawpercent[dm] = ((atms.dst5[dm-1] / (atms.dst5[dm-1] - atms.dst5[dm]))
       + ( 1- (atms.dst5[dm+1] / (atms.dst5[dm+1] - atms.dst5[dm]))))/2.0;

 if (veg.thawpercent[dm] > 1.0) veg.thawpercent[dm] = 1.0;
 if (veg.thawpercent[dm] < 0.0) veg.thawpercent[dm] = 0.0;

/*
 if (atms.dst5[dm] > 0.0)
     {
      if (atms.dst5[dm-1] > 0.0)
         {
          if (atms.dst5[dm+1] > 0.0)
            {
             veg.thawpercent[dm] = 1.0; // summer months
			 // To overcome wired JUNE phenomena, Q. Z. Sept. 03
			 /2001, delete the following dm ==5 sentences
			 if (dm == 5) {
			   if (atms.frontd[5] < 1.5) veg.thawpercent[5]
			   = abs(atms.frontd[5] / 1.5);

//              veg.thawpercent[4] = abs(atms.dst5[3]
				/ (atms.dst5[4] - atms.dst5[3])); // may
//              veg.thawpercent[dm] = 1.1 * veg.thawpercent[4];
			 // Assume Jun thawing has not reached the rooting depth yet
			 }
			}
		  else //atms.dst5[dm+1] < 0.0)
			{
		 //   veg.thawpercent[dm] = atms.dst5[dm]
		 /(atms.dst5[dm] - atms.dst5[dm+1]); // sept.
		   veg.thawpercent[dm] = 1.0;
			}
		 }
		 else // atms.dst5[dm-1]<0.0
		 {
		   if (atms.dst5[dm+1] > 0.0)
		   {
//            veg.thawpercent[dm] = abs(atms.dst5[dm-1]
			/ (atms.dst5[dm] - atms.dst5[dm-1])); // may
			veg.thawpercent[dm] =1.0;

		   }
		   else // //atms.dst5[dm+1] < 0.0)
			{
			 veg.thawpercent[dm] =0.0;  // rare to happen
			}
		 }
	 }
  else // (atms.dst5[dm] < 0.0)
   {
	 if (atms.dst5[dm-1] > 0.0)
		 {
		  if (atms.dst5[dm+1] > 0.0)
			{
			 veg.thawpercent[dm] = 1.0; // rare happen
			}
		  else
			{
				if (atms.dst5[dm-1] - atms.dst5[dm] == 0)
				{
				   veg.thawpercent[dm] = 1.0;
				}
				else
				{
				   veg.thawpercent[dm] = atms.dst5[dm-1]
				   /(atms.dst5[dm-1] - atms.dst5[dm]); // october
				}
//           veg.thawpercent[dm] = atms.dst5[dm-1]
			/(atms.dst5[dm-1] - atms.dst5[dm]); // October
//              veg.thawpercent[dm] = 1.0;
			 }
		 }
		 else // atms.dst5[dm-1]<0.0
		 {
		   if (atms.dst5[dm+1] > 0.0)
		   {
				if (atms.dst5[dm] - atms.dst5[dm-1] == 0)
				{
				   veg.thawpercent[dm] = 0.01;
				}
				else
				{
				   veg.thawpercent[dm] = abs(atms.dst5[dm]
				   / (atms.dst5[dm] - atms.dst5[dm-1])); // april
				}
		   }
		   else // //atms.dst5[dm+1] < 0.0)
			{
			 veg.thawpercent[dm] =0.0;  // rare to happen
			}
		 }
	  }

  */
	getenviron(dm);

	// Call adaptive Runge-Kutta integrator extended to write to a DOS window

	mintflag = pcadapt(NUMEQ,y,tol,dm);

	if (mintflag == 1) { intflag = 1; }

	if (adapttol == 1 && blackhol == 1 && endgrid == 0)
	{
//      endgrid = 1;
	  intbomb = 1;
	}

	// Check mass balance
	massbal(y,prevy);

	// Save current months pools and fluxes to check next month's
	//   mass balance

	setPrevState(prevy, y);

	// Update carbon, nitrogen and water pools and
	// fluxes from integrator results

	setMonth(dm, y);

//    cout << endl << " INNPP = " << veg.innpp[dm];
//    cout << " NPP = " << veg.npp[dm];
//    exit(-1);

	// Display monthly results to a DOS window
	if (DisplayFlag)
		pcdisplayMonth(dyr, dm, y);


//    cout << endl << " UNRMLEAF = " << veg.unnormleaf[dm];
//    cout << " LEAF = " << veg.leaf[dm];
//    exit(-1);

	resetODEflux(y);
  }

  //  Determine maximum eet, maximum pet, and new gpp optimum temperature

  for (dm = 0; dm < CYCLE; dm++)
  {
	if (dm == 0)
	{
	  atms.prveetmx = atms.eet[0];
	  atms.prvpetmx = atms.pet[0];
	  veg.topt   = atms.tair[0];
	  veg.prvleafmx[veg.cmnt] = veg.unnormleaf[0];
	}
	else
	{
	  if (atms.eet[dm] > atms.prveetmx) { atms.prveetmx = atms.eet[dm]; }
	  if (atms.pet[dm] > atms.prvpetmx) { atms.prvpetmx = atms.pet[dm]; }
	  if (veg.aleaf[veg.cmnt] == 0.0 && veg.bleaf[veg.cmnt] == 0.0
		 && veg.cleaf[veg.cmnt] == 1.0)
	  {
		if (atms.tair[dm] > veg.topt) { veg.topt = atms.tair[dm]; }
	  }
	  else
	  {
		if (veg.unnormleaf[dm] > veg.prvleafmx[veg.cmnt])
		{
	  veg.prvleafmx[veg.cmnt] = veg.unnormleaf[dm];
	  veg.topt = atms.tair[dm];
	}
	  }
	}
  }

  if (veg.topt > veg.toptmax[veg.cmnt]) { veg.topt = veg.toptmax[veg.cmnt]; }
  if (veg.topt < veg.toptmin[veg.cmnt]) { veg.topt = veg.toptmin[veg.cmnt]; }


// Determine vegetation C/N parameter as a function of vegetation type,
// annual PET, and annual EET

  veg.updateC2N(veg.cmnt, atms.yreet, atms.yrpet, atms.co2[11], atms.initco2);


  soil.yravlh2o /= 12.0;
  soil.yrrgrndh2o /= 12.0;
  soil.yrsnowpack /= 12.0;
  soil.yrsgrndh2o /= 12.0;
  soil.yrsmoist /= 12.0;
  soil.yrpctp /= 12.0;
  soil.yrvsm /= 12.0;

  veg.yrcarbon  /= 12.0;
  veg.yrnitrogen  /= 12.0;
  veg.yrstructn /= 12.0;

  if (veg.yrstructn != 0.0)
  {
	veg.yrc2n  = veg.yrcarbon / veg.yrstructn;
  }
  else { veg.yrc2n = MISSING; }               // specific to pcstepyr?

  soil.yrorgc /= 12.0;
  soil.yrorgn /= 12.0;

  if (soil.yrorgn != 0.0)
  {
	soil.yrc2n = soil.yrorgc / soil.yrorgn;
  }
  else { soil.yrc2n = MISSING; }             // specific to pcstepyr?

  soil.yravln  /= 12.0;
  veg.yrstoren /= 12.0;
  veg.yrunleaf /= 12.0;
  veg.yrleaf /= 12.0;
  veg.yrfpc /= 12.0;


  // y[2] changed to y[I_SOLC] by DWK on 20000615 and
  // y[3] changed to y[I_SOLN] by DWK on 20000615
  if (baseline == 1)
  {
	soil.yrnin = 0.0;
	soil.yrnlost = 0.0;
	if (y[I_SOLC]/microbe.cnsoil[veg.cmnt] >= y[I_SOLN])
	{
	  soil.yrnin = (y[I_SOLC]/microbe.cnsoil[veg.cmnt]) - y[I_SOLN];
	}
	else
	{
	  soil.yrnlost = y[I_SOLN] - (y[I_SOLC]/microbe.cnsoil[veg.cmnt]);
	}
	y[I_SOLN] = y[I_SOLC]/microbe.cnsoil[veg.cmnt];
  }

  return intflag;

};

