/* *************************************************************
****************************************************************
TTEM423.CPP - Terrestrial Ecosystem Model Version 4.2
****************************************************************

Modifications:

19991028 - DWK add bug fixes
20000107 - DWK adds compiler directive
20000107 - DWK renames functions to better describe purpose
20000130 - DWK changes prevleafmx to initleafmx in fourth getsiteecd()
20000130 - DWK changes y[2] to y[I_SOLC] and y[3] to y[I_SOLN]
           at bottom of stepyr()
20000201 - DWK changes pstate[1] to pstate[I_STRN] in delta()
20000207 - DWK added a minimum VEGC and SOLC (0.00001 gC m-2) condition
           to massbal(); Commented out by DWK on 20000210
20000614 - DWK veg.ingpp changed to veg.ingpp[] in delta(), setELMNTflux(), and
           setmonth()
20000614 - DWK veg.inpp changed to veg.innpp[] in delta(), setELMNTflux(), and
           setmonth()
20000614 - DWK resorts function names alphabetically to allow function
           descriptions to be found easier
20000616 - CFD removed Rh and Nmin from boundcon(), added RvMaint to boundcon()
20000616 - CFD restructured massbal()
20000630 - DWK changes missing values from -99.99 to -999.99 in ecdqc()
20010418 - Q. Z. Soil thermal model
20020202 - Kick and Qianlai, add leaf to calculation of LAI in delta()
20020202 - Q. Z change the include from ttem423e.cpp to ttem423e1.cpp

****************************************************************
************************************************************** */

#ifndef TTEM423E1_H
  #include "ttem423e1.hpp"
#endif

/* *********************************************************** */

TTEM::TTEM() : Odeint4()
{

  nfert = 1.00000;
  tol = inittol;
  totyr = -99;

// Identify potential output variables from TEM

// Ecosystem carbon pools ***************************************

  strcpy(predstr[I_VEGC],"VEGC");       // vegetation carbon
  strcpy(predstr[I_SOLC],"SOILORGC");   // soil organic carbon


  strcpy(predstr[I_TOTC],"TOTALC");     // total carbon


// Ecosystem nitrogen pools *************************************

  // total nitrogen stored in vegetation
  strcpy(predstr[I_VEGN],"VEGN");

  // vegetation structural nitrogen
  strcpy(predstr[I_STRN],"VSTRUCTN");

  // vegetation labile nitrogen
  strcpy(predstr[I_STON],"VSTOREN");

  strcpy(predstr[I_SOLN],"SOILORGN");  // soil organic nitrogen
  strcpy(predstr[I_AVLN],"AVAILN");    // soil available nitrogen


// Carbon and nitrogen pools associated with human products *****

  // carbon in agricultural products
  strcpy(predstr[I_AGPRDC],"AGPRODC");

  // nitrogen in agricultural products
  strcpy(predstr[I_AGPRDN],"AGPRODN");

  // carbon pool of products that decompose in 10 years
  strcpy(predstr[I_PROD10C],"PROD10C");

  // nitrogen pool of products that decompose in 10 years
  strcpy(predstr[I_PROD10N],"PROD10N");

  // carbon pool of products that decompose in 100 years
  strcpy(predstr[I_PROD100C],"PROD100C");

  // nitrogen pool of products that decompose in 100 years
  strcpy(predstr[I_PROD100N],"PROD100N");

  // carbon in all product pools
  strcpy(predstr[I_TOTPRDC],"TOTPRODC");

  // nitrogen in all product pools
  strcpy(predstr[I_TOTPRDN],"TOTPRODN");

  // total carbon pool found in ecosystem excluding products
  strcpy(predstr[I_TOTEC],"TOTEC");

  // total carbon pool found in ecosystem including products
  strcpy(predstr[I_TOTGC],"TOTGC");


// Ecosystem water pools ****************************************

  // available soil moisture
  strcpy(predstr[I_AVLW],"AVAILH2O");

  // groundwater pool resulting from rainfall
  strcpy(predstr[I_RGRW],"RGRNDH2O");

  strcpy(predstr[I_SNWPCK],"SNOWPACK");  // snowpack

  // groundwater pool resulting from snow melt
  strcpy(predstr[I_SGRW],"SGRNDH2O");

  strcpy(predstr[I_SM],"SOILH2O");       // soil moisture

   // added on 15/NOV/98 for soil temperature
   strcpy(predstr[I_TSOIL],"TSOIL");
   strcpy(predstr[I_DST5],"DST5");
   strcpy(predstr[I_DST10],"DST10");
   strcpy(predstr[I_DST20],"DST20");
   strcpy(predstr[I_DST50],"DST50");
   strcpy(predstr[I_DST100],"DST100");
   strcpy(predstr[I_DST200],"DST200");
   strcpy(predstr[I_FRONTD],"FRONTD");
   strcpy(predstr[I_THAWBE],"THAWBE");
   strcpy(predstr[I_THAWEND],"THAWEND");
// end of ...



  // soil moisture expressed as percent total porosity
  strcpy(predstr[I_PCTP],"PCTP");

  strcpy(predstr[I_VSM],"VSM");      // volumetric soil moisture


// Carbon fluxes for natural ecosystems *************************

  // GPP not limited by nutrient availability
  strcpy(predstr[I_INGPP],"VEGINGPP");

  strcpy(predstr[I_GPP],"GPP");      // gross primary production

  // NPP not limited by nutrient availability
  strcpy(predstr[I_INNPP],"VEGINNPP");

  strcpy(predstr[I_NPP],"NPP");      // net primary production
  strcpy(predstr[I_GPR],"GPR");      // gross plant respiration

  // vegetation maintenance respiration
  strcpy(predstr[I_RVMNT],"RVMAINT");

  // vegetation growth respiration
  strcpy(predstr[I_RVGRW],"RVGRWTH");

  strcpy(predstr[I_LTRC],"LTRC");   // litterfall carbon
  strcpy(predstr[I_RH],"RH");       // heterotrophic respiration
  strcpy(predstr[I_NEP],"NEP");     // net ecosystem production


// Nitrogen fluxes for natural ecosystems ***********************

  // total nitrogen inputs into ecosystem
  strcpy(predstr[I_NINP],"NINPUT");

  // VEGNUP not limited by carbon availability
  strcpy(predstr[I_INNUP],"VEGINNUP");

  // nitrogen uptake by vegetation
  strcpy(predstr[I_VNUP],"VEGNUP");

  // vegetation nitrogen uptake for structural components
  strcpy(predstr[I_VSUP],"VEGSUP");

  // vegetation nitrogen uptake for labile components
  strcpy(predstr[I_VLUP],"VEGLUP");

  // nitrogen mobilization by vegetation
  strcpy(predstr[I_VNMBL],"VNMOBIL");

  strcpy(predstr[I_VNRSRB],"VNRESORB"); // nitrogen resorption by vegetation
  strcpy(predstr[I_LTRN],"LTRN");       // litterfall nitrogen
  strcpy(predstr[I_MNUP],"MICRONUP");   // nitrogen uptake by microbes
  strcpy(predstr[I_NMIN],"NETNMIN");    // net nitrogen mineralization
  strcpy(predstr[I_NLST],"NLOST");      // total nitrogen losses from ecosystem


// Water fluxes *************************************************

  strcpy(predstr[I_RAIN],"RAIN");        // rainfall

  // percolation of rainwater through soil profile
  strcpy(predstr[I_RPERC],"RPERC");

  strcpy(predstr[I_RRUN],"RRUN");        // runoff of rainwater
  strcpy(predstr[I_SNWFAL],"SNOWFALL");  // snowfall

  // infiltration into the soil of water from snowmelt
  strcpy(predstr[I_SNWINF],"SNOWINF");

  // percolation of snowmelt through soil profile
  strcpy(predstr[I_SPERC],"SPERC");

  strcpy(predstr[I_SRUN],"SRUN");        // runoff of snowmelt
  strcpy(predstr[I_PET],"PET");          // potential evapotranspiration
  strcpy(predstr[I_EET],"EET");          // estimated evapotranspiration
  strcpy(predstr[I_WYLD],"H2OYIELD");    // water yield


// Phenology variables for natural ecosystems *******************

  // un-normalized relative phenology variable
  strcpy(predstr[I_UNRMLF],"UNRMLEAF");

  // normalized relative phenology variable (0 - 1.0)
  strcpy(predstr[I_LEAF],"LEAF");

  strcpy(predstr[I_FPC],"FPC");          // foliar projected cover
  strcpy(predstr[I_LAI],"LAI");          // leaf area index


// Carbon and nitrogen fluxes associated with agricultural
// conversion ***************************************************

  strcpy(predstr[I_CNVRTC],"CONVERTC"); // carbon loss during conversion
  strcpy(predstr[I_CNVRTN],"CONVERTN"); // nitrogen loss during conversion
  strcpy(predstr[I_SCNVRTC],"SCONVRTC");
  strcpy(predstr[I_SCNVRTN],"SCONVRTN");
  strcpy(predstr[I_NVRTNT],"NVRETENT");
  strcpy(predstr[I_NSRTNT],"NSRETENT");
  strcpy(predstr[I_NRETNT],"NRETENT");

  // carbon associated with slash left after conversion
  strcpy(predstr[I_SLASHC],"SLASHC");

  // nitrogen associated with slash left after conversion
  strcpy(predstr[I_SLASHN],"SLASHN");

  // carbon loss to formation of products that decompose in 10 years
  strcpy(predstr[I_PRDF10C],"PRDF10C");

  // nitrogen loss to formation of products that decompose in 10 years
  strcpy(predstr[I_PRDF10N],"PRDF10N");

  // carbon loss to formation of products that decompose in 100 years
  strcpy(predstr[I_PRDF100C],"PRDF100C");

  // nitrogen loss to formation of products that decompose in 100 years
  strcpy(predstr[I_PRDF100N],"PRDF100N");


// Carbon and nitrogen fluxes associated with agriculture *******

  // agricultral net primary production
  strcpy(predstr[I_AGNPPC],"AGNPPC");

  // nitrogen uptake by crops associated with AGNPPC
  strcpy(predstr[I_AGNPPN],"AGNPPN");

  // carbon loss to formation of agricultural products
  strcpy(predstr[I_AGFPRDC],"AGFPRODC");

  // nitrogen loss to formation of agricultural products
  strcpy(predstr[I_AGPRDN],"AGFPRODN");

  strcpy(predstr[I_AGFRTN],"AGFERTN"); // nitrogen fertilization
  strcpy(predstr[I_AGLTRC],"AGLTRFC"); // litterfall carbon from crops
  strcpy(predstr[I_AGLTRN],"AGLTRFN"); // litterfall nitrogen from crops


// Carbon and nitrogen fluxes associated with human products ****

  // carbon loss to the formation of all products

  strcpy(predstr[I_TOTFPRDC],"TOTFPRDC");

  // nitrogen loss to the formation of all products

  strcpy(predstr[I_TOTFPRDN],"TOTFPRDN");

  // carbon loss to resulting from decomposition of agricultural
  //   products

  strcpy(predstr[I_AGPRDFC],"AGPRODFC");

  // nitrogen loss resulting from decomposition of agricultural
  //   products

  strcpy(predstr[I_AGPRDFN],"AGPRODFN");

  // carbon loss resulting from decomposition of PROD10C

  strcpy(predstr[I_PRD10FC],"PRD10FC");

  // nitrogen loss resulting from decomposition of PROD10N

  strcpy(predstr[I_PRD10FN],"PRD10FN");

  // carbon loss resulting from decomposition of PROD100C

  strcpy(predstr[I_PRD100FC],"PRD100FC");

  // nitrogen loss resulting from decomposition of PROD100N

  strcpy(predstr[I_PRD100FN],"PRD100FN");

  // carbon loss resulting from decomposition of all products

  strcpy(predstr[I_TOTPRDFC],"TOTPRDFC");

  // nitrogen loss resulting from decomposition of all products

  strcpy(predstr[I_TOTPRDFN],"TOTPRDFN");


// Integrated carbon fluxes *************************************

  // total net primary production (NPP+AGNPPC)

  strcpy(predstr[I_TOTNPP],"TOTNPP");

  // carbon flux from ecosystem (NEP+CONVERTC)

  strcpy(predstr[I_CFLX],"CFLUX");

};

/* **************************************************************
************************* Public Functions **********************
************************************************************** */


/* **************************************************************
************************************************************** */
//Rh and Nmin removed by CFD 20000616

int TTEM::boundcon(double ptstate[], double err[], double& ptol)
{

  int test = ACCEPT;

// Check carbon and nitrogen state variables
  if (err[I_VEGC] > fabs(ptol * ptstate[I_VEGC]))
  {
    return test = temkey(I_VEGC)+1;
  }
  if (nfeed == 1 && err[I_STRN] > fabs(ptol * ptstate[I_STRN]))
  {
    return test = temkey(I_STRN)+1;
  }
  if (err[I_SOLC] > fabs(ptol * ptstate[I_SOLC]))
  {
    return test = temkey(I_SOLC)+1;
  }
  if (nfeed == 1 && err[I_SOLN] > fabs(ptol * ptstate[I_SOLN]))
  {
    return test = temkey(I_SOLN)+1;
  }
  if (nfeed == 1 && err[I_AVLN] > fabs(ptol * ptstate[I_AVLN]))
  {
    return test = temkey(I_AVLN)+1;
  }
  if (err[I_GPP] > fabs(ptol * ptstate[I_GPP]))
  {
    return test = temkey(I_GPP)+1;
  }
  if (err[I_NPP] > fabs(ptol * ptstate[I_NPP]))
  {
    return test = temkey(I_NPP)+1;
  }
  if (nfeed == 1 && err[I_VNUP] > fabs(ptol * ptstate[I_VNUP]))
  {
    return test = temkey(I_VNUP)+1;
  }
  if (nfeed == 1 && err[I_VSUP] > fabs(ptol * ptstate[I_VSUP]))
  {
    return test = temkey(I_VSUP)+1;
  }
  if (nfeed == 1 && err[I_STON] > fabs(ptol * ptstate[I_STON]))
  {
    return test = temkey(I_STON)+1;
  }
  if (nfeed == 1 && err[I_VNMBL] > fabs(ptol * ptstate[I_VNMBL]))
  {
    return test = temkey(I_VNMBL)+1;
  }
  //20000616 - CFD added RvMaint
  if (nfeed == 1 && err[I_RVMNT] > fabs(ptol * ptstate[I_RVMNT]))
  {
    return test = temkey(I_RVMNT)+1;
  }
  //20000616 - CFD removed Rh and Nmin

  // Check water state variables

  if (err[I_AVLW]  > fabs(ptol * ptstate[I_AVLW]))
  {
    return test = temkey(I_AVLW)+1;
  }
  if (err[I_RGRW]  > fabs(ptol * ptstate[I_RGRW]))
  {
    return test = temkey(I_RGRW)+1;
  }
  if (err[I_SNWPCK]  > fabs(ptol * ptstate[I_SNWPCK]))
  {
    return test = temkey(I_SNWPCK)+1;
  }
  if (err[I_SGRW]  > fabs(ptol * ptstate[I_SGRW]))
  {
    return test = temkey(I_SGRW)+1;
  }
  if (err[I_RPERC]  > fabs((ptol) * ptstate[I_RPERC]))
  {
    return test = temkey(I_RPERC)+1;
  }
  if (err[I_EET]  > fabs((ptol) * ptstate[I_EET]))
  {
    return test = temkey(I_EET)+1;
  }


  return test;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM::delta(const int& dm, double pstate[], double pdstate[])
{

  nfert = 1.0000000;

  atms.eet[dm] = atms.xeet(atms.rain[dm], soil.snowinf[dm], atms.pet[dm],
                           pstate[I_AVLW], soil.awcapmm, dm);

  soil.percol(atms.rain[dm], soil.snowinf[dm], atms.eet[dm],
              pstate[I_AVLW],dm);

  if ((pstate[I_AVLW] + soil.snowinf[dm] + atms.rain[dm] - atms.eet[dm]
     - soil.rperc[dm] - soil.sperc[dm]) < 0.0)
  {
    atms.eet[dm] = pstate[I_AVLW] + soil.snowinf[dm] + atms.rain[dm]
                   - soil.rperc[dm] - soil.sperc[dm];
  }


  soil.rrun[dm] = soil.rrunoff(pstate[I_RGRW], soil.rperc[dm]);
  soil.srun[dm] = soil.srunoff(elev, atms.tair[dm], atms.prevtair,
                               atms.prev2tair, pstate[I_SGRW], soil.sperc[dm]);

  soil.h2oyld[dm] = soil.rrun[dm] + soil.srun[dm];


  if (moistlim == 0)
  {
    veg.unnormleaf[dm] = veg.deltaleaf(veg.cmnt, atms.pet[dm],
                                       atms.prvpetmx, prevy[I_UNRMLF]);
  }
  else
  {
    veg.unnormleaf[dm] = veg.deltaleaf(veg.cmnt, atms.eet[dm],
                                       atms.prveetmx, prevy[I_UNRMLF]);
  }
  if (veg.prvleafmx[veg.cmnt] <= 0.0) { veg.leaf[dm] = 0.0; }
  else { veg.leaf[dm] = veg.unnormleaf[dm]/veg.prvleafmx[veg.cmnt]; }
  if (veg.leaf[dm] < veg.minleaf[veg.cmnt])
  {
    veg.leaf[dm] = veg.minleaf[veg.cmnt];
  }
  if (veg.leaf[dm] > 1.0) { veg.leaf[dm] = 1.0; }

  veg.alleaf = veg.leafmxc[veg.cmnt]/(1.0 + veg.kleafc[veg.cmnt]
               * exp(veg.cov[veg.cmnt]*pstate[I_VEGC]));
  veg.lai[dm] = veg.sla[veg.cmnt] * veg.alleaf * veg.leaf[dm];
   // modified 02/02/2002 Kick and Qianlai
  veg.foliage = veg.alleaf / veg.leafmxc[veg.cmnt];
  veg.fpc[dm] = 1.0 - exp(-0.5*veg.lai[dm]);
  // modified by Kick and Qialai 02/02/2002

  deltaxclm(veg.cmnt, soil.pcfldcap, dm);

  if (ag.state == 0)
  {

      // added veg.thawpercent for thaw-frozen mechanism by qianlai 28/08/2000
    veg.ingpp[dm] = veg.gppxclm(veg.cmnt, atms.co2[dm], atms.par[dm],
							temp, gv, veg.leaf[dm], veg.foliage,
							veg.thawpercent[dm] );
//    veg.ingpp[dm] = veg.gppxclm(veg.cmnt, atms.co2[dm], atms.par[dm],
//                            temp, gv, veg.leaf[dm], veg.foliage);
    if (veg.ingpp[dm] < 0.0) { veg.ingpp[dm] = 0.0; }

    veg.inuptake = veg.nupxclm(veg.cmnt, pstate[I_SM], pstate[I_AVLN],
                               respq10, ksoil, veg.foliage);
  }
  else
  {
    veg.ingpp[dm] = 0.0;
    veg.inuptake = 0.0;
  }
  microbe.rh[dm] = microbe.rhxclm(pstate[I_SOLC], dq10, rhmoist);
  if (microbe.rh[dm] < 0.0) { microbe.rh[dm] = 0.0; }

  soil.ninput[dm] = 0.0;

  microbe.netnmin[dm] = microbe.nminxclm(veg.cmnt, dm, pstate[I_SM],
                                         pstate[I_SOLC], pstate[I_SOLN],
                                         pstate[I_AVLN], microbe.decay,
                                         microbe.rh[dm], ksoil);

  if (ag.state == 0)
  {
    veg.ltrfal[dm].carbon = veg.cfall[veg.cmnt] * pstate[I_VEGC];
    if (veg.ltrfal[dm].carbon < 0.0) { veg.ltrfal[dm].carbon = 0.0; }

    veg.ltrfal[dm].nitrogen = veg.nfall[veg.cmnt] * pstate[I_STRN];
    if (veg.ltrfal[dm].nitrogen < 0.0) { veg.ltrfal[dm].nitrogen = 0.0; }

    veg.rm[dm] = veg.rmxclm(veg.cmnt,pstate[I_VEGC],respq10);
    if (veg.rm[dm] < 0.0) { veg.rm[dm] = 0.0; }

    veg.innpp[dm] = veg.ingpp[dm] - veg.rm[dm];

    veg.rg[dm] = 0;
    if (veg.innpp[dm] > 0.0)
    {
      veg.rg[dm]  = 0.2 * veg.innpp[dm];
      veg.innpp[dm] *= 0.8;
    }

    if (veg.inuptake > pstate[I_AVLN] + (nfert*soil.ninput[dm])
    + microbe.netnmin[dm])
    {
      veg.inuptake = pstate[I_AVLN] + (nfert*soil.ninput[dm])
                     + microbe.netnmin[dm];
    }

    if (veg.inuptake < 0.0) { veg.inuptake = 0.0; }

    veg.nuptake[dm] = veg.inuptake;
    veg.suptake[dm] = veg.nuptake[dm];
    veg.luptake[dm] = 0.0;
    veg.gpp[dm] = veg.ingpp[dm];
    veg.npp[dm] = veg.innpp[dm];
    veg.nmobil[dm] = 0.0;
    veg.nresorb[dm] = 0.0;

// Nitrogen feedback of GPP (nfeed == 1)

    if (nfeed == 1)
    {
      if (veg.inuptake == 0.0) { veg.inuptake = 0.000001; }
      veg.inprodcn = veg.innpp[dm] / (veg.inuptake + pstate[I_STON]);

      if (veg.ltrfal[dm].nitrogen <= veg.ltrfal[dm].carbon/veg.cneven[veg.cmnt])
      {
	veg.nresorb[dm] = veg.ltrfal[dm].carbon/veg.cneven[veg.cmnt]
                              - veg.ltrfal[dm].nitrogen;
      }
      else
      {
	veg.ltrfal[dm].nitrogen = veg.ltrfal[dm].carbon/veg.cneven[veg.cmnt];
	veg.nresorb[dm] = 0.0;
      }
      if (pstate[I_VEGC] > 0.0)
      {
        veg.nresorb[dm] *= (pstate[I_STRN]/pstate[I_VEGC]) * veg.c2n;
      }

      if (veg.inprodcn > veg.cneven[veg.cmnt])
      {
	veg.npp[dm] = veg.cneven[veg.cmnt] * (veg.nuptake[dm] + pstate[I_STON]);
	if (veg.npp[dm] < 0.0) { veg.npp[dm] = 0.0; }
	veg.rg[dm] = 0.25 * veg.npp[dm];
  	veg.gpp[dm] = veg.npp[dm] + veg.rg[dm] + veg.rm[dm];
	if (veg.gpp[dm] < 0.0) { veg.gpp[dm] = 0.0; }

	veg.nmobil[dm] = pstate[I_STON];
      }

      if (veg.inprodcn <= veg.cneven[veg.cmnt])
      {
        veg.nuptake[dm] = veg.inuptake * (veg.inprodcn - veg.cnmin[veg.cmnt])
                          * (veg.inprodcn - 2*veg.cneven[veg.cmnt]
                          + veg.cnmin[veg.cmnt]);
	veg.nuptake[dm] /= ((veg.inprodcn - veg.cnmin[veg.cmnt])
                           * (veg.inprodcn - 2*veg.cneven[veg.cmnt]
                           + veg.cnmin[veg.cmnt])) - pow(veg.inprodcn
                           - veg.cneven[veg.cmnt],2.0);
	if (veg.nuptake[dm] < 0.0) { veg.nuptake[dm] = 0.0; }
	if (pstate[I_STON] >= veg.npp[dm]/veg.cneven[veg.cmnt])
        {
	   veg.nmobil[dm] = veg.npp[dm]/veg.cneven[veg.cmnt];
	   if (veg.nmobil[dm] < 0.0 && pstate[I_VEGC] > 0.0)
           {
	     veg.nmobil[dm] *= (pstate[I_STRN]/pstate[I_VEGC]) * veg.c2n;
	   }
	   veg.suptake[dm] = 0.0;
	}
	else
        {
	  veg.nmobil[dm] = pstate[I_STON];
	  veg.suptake[dm] = (veg.npp[dm]/veg.cneven[veg.cmnt])
                             - veg.nmobil[dm];
	  if (veg.suptake[dm] < 0.0) { veg.suptake[dm] = 0.0; }
	  if (veg.suptake[dm] > veg.nuptake[dm])
          {
            veg.suptake[dm] = veg.nuptake[dm];
          }
	}

    // ptate[1] changed to pstate[I_STRN] by DWK on 20000201
	if ((pstate[I_STON] + veg.nuptake[dm] - veg.suptake[dm]
           + veg.nresorb[dm] - veg.nmobil[dm]) < (veg.labncon[veg.cmnt]
           * (pstate[I_STRN] + veg.suptake[dm] - veg.ltrfal[dm].nitrogen
           - veg.nresorb[dm] + veg.nmobil[dm])))
        {
	  veg.luptake[dm] = veg.nuptake[dm] - veg.suptake[dm];
	}
	else
        {
	   veg.luptake[dm] = (veg.labncon[veg.cmnt] * (pstate[I_STRN]
                             + veg.suptake[dm] - veg.ltrfal[dm].nitrogen
                             - veg.nresorb[dm] + veg.nmobil[dm]))
                             - (pstate[I_STON] + veg.nresorb[dm]
                             - veg.nmobil[dm]);
	   if (veg.luptake[dm] < 0.0) { veg.luptake[dm] = 0.0; }
	   veg.nuptake[dm] = veg.suptake[dm] + veg.luptake[dm];
	}
      }
    }

    veg.gpr[dm] = veg.rm[dm] + veg.rg[dm];
    nep[dm] = veg.npp[dm] - microbe.rh[dm];
    cflux[dm] = nep[dm];

    ag.npp[dm].carbon = 0.0;
    ag.npp[dm].nitrogen = 0.0;
    ag.fertn[dm] = 0.0;
    ag.ltrfal[dm].carbon = 0.0;
    ag.ltrfal[dm].nitrogen = 0.0;
    ag.slash[dm].carbon = 0.0;
    ag.slash[dm].nitrogen = 0.0;
    ag.sconvrtflx[dm].carbon = 0.0;
    ag.sconvrtflx[dm].nitrogen = 0.0;
    ag.nsretent[dm] = 0.0;
  }
  else
  {
    veg.ltrfal[dm].carbon = 0.0;
    veg.ltrfal[dm].nitrogen = 0.0;
    veg.gpr[dm] = 0.0;
    veg.rm[dm] = 0.0;
    veg.rg[dm] = 0.0;
    veg.innpp[dm] = 0.0;
    veg.nuptake[dm] = 0.0;
    veg.suptake[dm] = 0.0;
    veg.luptake[dm] = 0.0;
    veg.gpp[dm] = 0.0;
    veg.npp[dm] = 0.0;
    veg.nmobil[dm] = 0.0;
    veg.nresorb[dm] = 0.0;

    ag.npp[dm].carbon = ag.RAP * ag.potnpp[dm];
    if (ag.npp[dm].carbon < 0.0) { ag.npp[dm].carbon = 0.0; }
    ag.npp[dm].nitrogen = ag.npp[dm].carbon / ag.c2n;
    ag.fertn[dm] = 0.0;
    ag.ltrfal[dm].carbon = ag.npp[dm].carbon * ag.cfall[veg.cmnt];
    ag.ltrfal[dm].nitrogen = ag.npp[dm].nitrogen * ag.nfall[veg.cmnt];
    nep[dm] = ag.npp[dm].carbon - microbe.rh[dm];
    cflux[dm] = nep[dm] - ag.convrtflx[dm].carbon;
  }

// Changes in available nitrogen in soil

  if (avlnflag == 1)
  {
    if (ag.state == 0)
    {
      soil.nlost[dm] = pstate[I_AVLN] / pstate[I_SM];
      soil.nlost[dm] *= ((atms.rain[dm] + soil.snowinf[dm] - atms.eet[dm])
                        + (soil.rootz * 1000.0)) / (soil.rootz * 1000.0);

//     if  (soil.nlost[dm] > 1000.0 ) soil.nlost[dm] = 0.002771;
// denugged by Q. Z. 04/Nov/2001
//     if  (soil.nlost[dm] < 0.000001 ) soil.nlost[dm] = 0.002771;
// denugged by Q. Z. 04/Nov/2001

      soil.nlost[dm] *= soil.nloss[veg.cmnt];

      if (soil.nlost[dm] > pstate[I_AVLN] - veg.nuptake[dm]
         + microbe.netnmin[dm] + soil.ninput[dm])
      {
        soil.nlost[dm] = pstate[I_AVLN] - veg.nuptake[dm] + microbe.netnmin[dm]
                         + soil.ninput[dm];
      }
      if (soil.nlost[dm] < 0.0)
      {
        soil.nlost[dm]  = 0.0;
        microbe.netnmin[dm] = soil.nlost[dm] + veg.nuptake[dm] - soil.ninput[dm]
                              - pstate[I_AVLN];
      }
    }
    else
    {
      soil.ninput[dm] = ag.nretent[dm];
      soil.nlost[dm] = 0.0;

      // Condition added to limit addition of agfertn to only
      //   periods of crop growth
      if (ag.npp[dm].nitrogen > 0.0)
      {
        if (pstate[I_AVLN] + soil.ninput[dm] + microbe.netnmin[dm]
            - soil.nlost[dm] < ag.npp[dm].nitrogen)
        {
   	   ag.fertn[dm] = ag.npp[dm].nitrogen + soil.nlost[dm] - pstate[I_AVLN]
           - soil.ninput[dm] - microbe.netnmin[dm];
	   if (ag.fertn[dm] < 0.0)
           {
	      ag.fertn[dm] = 0.0;
	      microbe.netnmin[dm] = soil.nlost[dm] + ag.npp[dm].nitrogen
                                    - soil.ninput[dm] - pstate[I_AVLN];
	   }
	   soil.ninput[dm] += ag.fertn[dm];
        }
      }
    }
  }
  else
  {
    soil.nlost[dm] = soil.ninput[dm] - veg.nuptake[dm] - ag.npp[dm].nitrogen
                     + microbe.netnmin[dm];
  }

// Describe monthly changes to carbon pools and fluxes for ODE state variables
// (i.e., pdstate)

  // Carbon pools in natural ecosystems

  pdstate[I_VEGC] = veg.gpp[dm] - veg.gpr[dm] - veg.ltrfal[dm].carbon;
  pdstate[I_SOLC] = veg.ltrfal[dm].carbon + ag.ltrfal[dm].carbon
                    + ag.slash[dm].carbon - ag.sconvrtflx[dm].carbon
                    - microbe.rh[dm];
//llc test

//  printf("=1=%f, =2=%f, =3=%f, =4=%f, =4=%f, =5=%f,\n",veg.ltrfal[dm].carbon,ag.ltrfal[dm].carbon,ag.slash[dm].carbon,ag.sconvrtflx[dm].carbon,microbe.rh[dm],pdstate[I_SOLC] );

  // Nitrogen pools in natural ecosystems

  pdstate[I_STRN] = veg.suptake[dm] - veg.ltrfal[dm].nitrogen - veg.nresorb[dm]
                    + veg.nmobil[dm];
  pdstate[I_SOLN] = veg.ltrfal[dm].nitrogen + ag.ltrfal[dm].nitrogen
                    + ag.slash[dm].nitrogen - ag.sconvrtflx[dm].nitrogen
                    - ag.nsretent[dm] - microbe.netnmin[dm];
  pdstate[I_AVLN] = soil.ninput[dm] - soil.nlost[dm] + microbe.netnmin[dm]
                    - veg.nuptake[dm] - ag.npp[dm].nitrogen;
  pdstate[I_STON] = veg.luptake[dm] + veg.nresorb[dm] - veg.nmobil[dm];

  // Human product pools

  pdstate[I_AGPRDC] = 0.0;
  pdstate[I_AGPRDN] = 0.0;
  pdstate[I_PROD10C] = 0.0;
  pdstate[I_PROD10N] = 0.0;
  pdstate[I_PROD100C] = 0.0;
  pdstate[I_PROD100N] = 0.0;
  pdstate[I_TOTPRDC] = 0.0;
  pdstate[I_TOTPRDN] = 0.0;
  pdstate[I_TOTEC] = 0.0;
  pdstate[I_TOTGC] = 0.0;

  // Water pools

  pdstate[I_AVLW] = soil.snowinf[dm] + atms.rain[dm] - atms.eet[dm]
                    - soil.rperc[dm] - soil.sperc[dm];
  pdstate[I_RGRW] = soil.rperc[dm] - soil.rrun[dm];
  pdstate[I_SNWPCK] = atms.snowfall[dm] - soil.snowinf[dm];
  pdstate[I_SGRW] = soil.sperc[dm] - soil.srun[dm];
  pdstate[I_SM] = pdstate[I_AVLW];
  pdstate[I_PCTP] = 100.0 * pdstate[I_SM]/soil.totpor;
  pdstate[I_VSM] = pdstate[I_SM]/(soil.rootz*1000.0);
  if (pstate[I_VSM]+pdstate[I_VSM] <= 0.0) {
    pdstate[I_VSM] = 0.001 - pstate[I_VSM];
  }

  // Carbon fluxes in natural ecosystems

  pdstate[I_INGPP] = veg.ingpp[dm];
  pdstate[I_GPP] = veg.gpp[dm];
  pdstate[I_INNPP] = veg.innpp[dm];
  pdstate[I_NPP] = veg.npp[dm];
  pdstate[I_GPR] = veg.gpr[dm];
  pdstate[I_RVMNT] = veg.rm[dm];
  pdstate[I_RVGRW] = veg.rg[dm];
  pdstate[I_LTRC] = veg.ltrfal[dm].carbon;
  pdstate[I_RH] = microbe.rh[dm];
  pdstate[I_NEP] = nep[dm];

  // Nitrogen fluxes in natural ecosystems

  pdstate[I_NINP] = soil.ninput[dm];
  pdstate[I_INNUP] = veg.inuptake;
  pdstate[I_VNUP] = veg.nuptake[dm];
  pdstate[I_VSUP] = veg.suptake[dm];
  pdstate[I_VLUP] = veg.luptake[dm];
  pdstate[I_VNMBL] = veg.nmobil[dm];
  pdstate[I_VNRSRB] = veg.nresorb[dm];
  pdstate[I_LTRN] = veg.ltrfal[dm].nitrogen;
  pdstate[I_MNUP] = microbe.nuptake[dm];
  pdstate[I_NMIN] = microbe.netnmin[dm];
  pdstate[I_NLST] = soil.nlost[dm];

  // Water fluxes

  pdstate[I_RAIN] = atms.rain[dm];
  pdstate[I_RPERC] = soil.rperc[dm];
  pdstate[I_RRUN] = soil.rrun[dm];
  pdstate[I_SNWFAL] = atms.snowfall[dm];
  pdstate[I_SNWINF] = soil.snowinf[dm];
  pdstate[I_SPERC] = soil.sperc[dm];
  pdstate[I_SRUN] = soil.srun[dm];
  pdstate[I_PET] = atms.pet[dm];
  pdstate[I_EET] = atms.eet[dm];
  pdstate[I_WYLD] = soil.rrun[dm] + soil.srun[dm];

       // for soil temperature
  pdstate[I_TSOIL] = atms.tsoil[dm];
  pdstate[I_DST5] = atms.dst5[dm];
  pdstate[I_DST10] = atms.dst10[dm];
  pdstate[I_DST20] = atms.dst20[dm];
  pdstate[I_DST50] = atms.dst50[dm];
  pdstate[I_DST100] = atms.dst100[dm];
  pdstate[I_DST200] = atms.dst200[dm];
  pdstate[I_FRONTD] = atms.frontd[dm];
  pdstate[I_THAWBE] = atms.thawbe[dm];
  pdstate[I_THAWEND] = atms.thawend[dm];
// end of ...


  // Phenology

  pdstate[I_UNRMLF] = veg.unnormleaf[dm];
  pdstate[I_LEAF] = veg.leaf[dm];

  // Carbon and nitrogen fluxes from agricultural conversion

  pdstate[I_CNVRTC] = ag.convrtflx[dm].carbon;
  pdstate[I_CNVRTN] = ag.convrtflx[dm].nitrogen;
  pdstate[I_SCNVRTC] = ag.sconvrtflx[dm].carbon;
  pdstate[I_SCNVRTN] = ag.sconvrtflx[dm].nitrogen;
  pdstate[I_NVRTNT] = ag.nvretent[dm];
  pdstate[I_NSRTNT] = ag.nsretent[dm];
  pdstate[I_NRETNT] = ag.nretent[dm];
  pdstate[I_SLASHC] = ag.slash[dm].carbon;
  pdstate[I_SLASHN] = ag.slash[dm].nitrogen;
  pdstate[I_PRDF10C] = ag.formPROD10.carbon / (double) CYCLE;
  pdstate[I_PRDF10N] = ag.formPROD10.nitrogen / (double) CYCLE;
  pdstate[I_PRDF100C] = ag.formPROD100.carbon / (double) CYCLE;
  pdstate[I_PRDF100N] = ag.formPROD100.nitrogen / (double) CYCLE;

  // Carbon and nitrogen fluxes in agricultural ecosystems

  pdstate[I_AGNPPC] = ag.npp[dm].carbon;
  pdstate[I_AGNPPN] = ag.npp[dm].nitrogen;
  pdstate[I_AGFPRDC] = ag.npp[dm].carbon - ag.ltrfal[dm].carbon;
  pdstate[I_AGFPRDN] = ag.npp[dm].nitrogen - ag.ltrfal[dm].nitrogen;
  pdstate[I_AGLTRC] = ag.ltrfal[dm].carbon;
  pdstate[I_AGLTRN] = ag.ltrfal[dm].nitrogen;
  pdstate[I_AGFRTN] = ag.fertn[dm];

  // Carbon and nitrogen fluxes from human product pools

  pdstate[I_TOTFPRDC] = ag.formTOTPROD.carbon / (double) CYCLE;
  pdstate[I_TOTFPRDN] = ag.formTOTPROD.nitrogen / (double) CYCLE;
  pdstate[I_AGPRDFC] = ag.PROD1decay.carbon / (double) CYCLE;
  pdstate[I_AGPRDFN] = ag.PROD1decay.nitrogen / (double) CYCLE;
  pdstate[I_PRD10FC] = ag.PROD10decay.carbon / (double) CYCLE;
  pdstate[I_PRD10FN] = ag.PROD10decay.nitrogen / (double) CYCLE;
  pdstate[I_PRD100FC] = ag.PROD100decay.carbon / (double) CYCLE;
  pdstate[I_PRD100FN] = ag.PROD100decay.nitrogen / (double) CYCLE;
  pdstate[I_TOTPRDFC] = ag.TOTPRODdecay.carbon / (double) CYCLE;
  pdstate[I_TOTPRDFN] = ag.TOTPRODdecay.nitrogen / (double) CYCLE;

  // Integrated carbon fluxes

  pdstate[I_TOTNPP] = veg.npp[dm] + ag.npp[dm].carbon;
  pdstate[I_CFLX] = cflux[dm];

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM::deltaxclm(const int& dcmnt, const double& pcfldcap, const int& dm)
{

  double vfc;


  vfc = pcfldcap * 0.01;


/* gv: effect of moisture on primary productivity */

  if (moistlim == 0)
  {
    gv = 1.0;
  }
  else
  {
    if (atms.eet[dm]/atms.pet[dm] <= 0.1)
    {
      gv = (-10.0 * pow((atms.eet[dm]/atms.pet[dm]),2.0))
           + (2.9 * (atms.eet[dm]/atms.pet[dm]));
      if (gv < 0.0) { gv = 0.0; }
    }
    else
    {
      gv = 0.1 + (0.9 * atms.eet[dm] / atms.pet[dm]);
    }
  }


/* ksoil: effect of soil moisture on nitrogen uptake by plants
			 and microbes */

  if (moistlim == 0) { ksoil = pow(vfc,3.0); }
  else { ksoil = pow(y[I_VSM],3.0); }


/* rhmoist: effect of moisture on decomposition */

  if (moistlim == 0)
  {
    rhmoist = (vfc - microbe.moistmin[dcmnt])
              * (vfc - microbe.moistmax[dcmnt]);
    rhmoist /= rhmoist - pow((vfc - microbe.moistopt[dcmnt]),2.0);
  }
  else
  {
    rhmoist = (y[I_VSM] - microbe.moistmin[dcmnt])
              * (y[I_VSM] - microbe.moistmax[dcmnt]);
    rhmoist /= rhmoist - pow((y[I_VSM] - microbe.moistopt[dcmnt]),2.0);
  }
  if (rhmoist < 0.0) { rhmoist = 0.0; }

};

/* *************************************************************
************************************************************* */



/* *************************************************************
************************************************************* */

int TTEM::ecdqc(const int& dcmnt)
{

  int qc = ACCEPT;

  if (vegca[dcmnt] <= -999.99) { return qc = REJECT; }
  if (vegcb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (strna[dcmnt] <= -999.99) { return qc = REJECT; }
  if (strnb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solca[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solcb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solna[dcmnt] <= -999.99) { return qc = REJECT; }
  if (solnb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (avlna[dcmnt] <= -999.99) { return qc = REJECT; }
  if (avlnb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (stona[dcmnt] <= -999.99) { return qc = REJECT; }
  if (stonb[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.unleaf12[dcmnt] <= -9.99) { return qc = REJECT; }
  if (veg.prvleafmx[dcmnt] <= -9.99) { return qc = REJECT; }
  if (veg.cmaxcut[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.cmax1a[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cmax1b[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cmax2a[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cmax2b[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cfall[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.kra[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.krb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.kda[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.kdb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.lcclnc[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.propftos[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.nmaxcut[dcmnt] <= -99.99) { return qc = REJECT; }
  // DWK changed missing values from -99.99 to -999.99 for "nmax"s
  if (veg.nmax1a[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.nmax1b[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.nmax2a[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.nmax2b[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.nfall[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.nupa[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.nupb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (soil.nloss[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.nfixpar[dcmnt] <= -99.9) { return qc = REJECT; }
  if (veg.cneven[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.cnmin[dcmnt] <= -999.99) { return qc = REJECT; }
  if (veg.c2na[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.c2nb[dcmnt] <= -99.99) { return qc = REJECT; }
  if (veg.c2nmin[dcmnt] <= -99.99) { return qc = REJECT; }
  if (microbe.cnsoil[dcmnt] <= -999.99) { return qc = REJECT; }

  return qc;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TTEM::ECDsetELMNTstate(const int& dcmnt, const double& psiplusc)
{

  int dyr;
  int dm;

  for (dm = 0; dm < CYCLE; dm++)
  {
    y[I_AVLW] = soil.avlh2o[dm] = soil.awcapmm;
    y[I_RGRW] = soil.rgrndh2o[dm] = 0.0;
    y[I_SNWPCK] = soil.snowpack[dm] = 0.0;
    y[I_SGRW] = soil.sgrndh2o[dm] = 0.0;
    y[I_SM] = soil.moist[dm] = soil.awcapmm + soil.wiltpt;
    y[I_PCTP] = soil.pctp[dm] = 100.0 * soil.moist[dm] / soil.totpor;
    soil.vsm[dm] = soil.moist[dm] / (soil.rootz * 1000.0);
    if (soil.vsm[dm] <= 0.0)
    {
      soil.vsm[dm] = 0.001;
    }
    y[I_VSM] = soil.vsm[dm];

    veg.plant[dm].carbon = vegca[dcmnt] * psiplusc + vegcb[dcmnt];
    y[I_VEGC] = veg.plant[dm].carbon;

    veg.strctrl[dm].nitrogen = strna[dcmnt] * psiplusc + strnb[dcmnt];
    y[I_STRN] = veg.strctrl[dm].nitrogen;

    soil.org[dm].carbon = solca[dcmnt] * psiplusc + solcb[dcmnt];
    y[I_SOLC] = soil.org[dm].carbon;

//llc test

//  printf("=1=%f, =2=%f, =3=%f, =4=%f,\n",y[I_SOLC],solca[dcmnt],psiplusc,solcb[dcmnt] );

    soil.org[dm].nitrogen = solna[dcmnt] * psiplusc + solnb[dcmnt];
    y[I_SOLN] = soil.org[dm].nitrogen;

    soil.availn[dm] = avlna[dcmnt] * psiplusc + avlnb[dcmnt];
    y[I_AVLN] = soil.availn[dm];

    veg.labile[dm].nitrogen = stona[dcmnt] * psiplusc + stonb[dcmnt];
    y[I_STON] = veg.labile[dm].nitrogen;

    y[I_UNRMLF] = veg.unnormleaf[dm] = veg.unleaf12[dcmnt];

    veg.plant[dm].nitrogen = 0.0;

    veg.alleaf = veg.leafmxc[dcmnt]/(1.0 + veg.kleafc[dcmnt]
                 * exp(veg.cov[dcmnt]*y[I_VEGC]));
    veg.lai[dm] = veg.sla[dcmnt] * veg.alleaf;
    veg.foliage = veg.alleaf / veg.leafmxc[dcmnt];
    veg.fpc[dm] = 1.0;

    y[I_PROD10C] = ag.PROD10.carbon = 0.0;
    y[I_PROD10N] = ag.PROD10.nitrogen = 0.0;
    y[I_PROD100C] = ag.PROD100.carbon = 0.0;
    y[I_PROD100N] = ag.PROD100.nitrogen = 0.0;
    y[I_AGPRDC] = ag.PROD1.carbon = 0.0;
    y[I_AGPRDN] = ag.PROD1.nitrogen = 0.0;
    y[I_TOTPRDC] = ag.TOTPROD.carbon = 0.0;
    y[I_TOTPRDN] = ag.TOTPROD.nitrogen = 0.0;
    y[I_TOTEC] = totalc[dm] = veg.plant[dm].carbon + soil.org[dm].carbon;
    y[I_TOTGC] = totalc[dm];
  }

  for (dyr = 0; dyr < 10; dyr++)
  {
    ag.initPROD10[dyr].carbon = 0.0;
    ag.initPROD10[dyr].nitrogen = 0.0;
  }

  for (dyr = 0; dyr < 100; dyr++)
  {
    ag.initPROD100[dyr].carbon = 0.0;
    ag.initPROD100[dyr].nitrogen = 0.0;
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

int TTEM::equilibrium(const int& itype, double& tol)
{

  int dyr = 0;

  setPrevState(prevy,y);
//  cout <<"run ttem:: equilibrium"<<endl;
//  system("pause");

  totyr = 0;
  endeq = 0;
  intflag = 0;
  initFlag = 0;
  while ((dyr < runsize) && (endeq < 2))
  {
//    cout <<"stepyr()"<<endl;
	endeq = stepyr(dyr,itype,intflag,tol);
    ++dyr;
    ++totyr;

 // Reset product fluxes and pools to zero
//	cout <<"ag.resetPROD()"<<endl;
	ag.resetPROD();
//	cout <<"run ag.resetPROD()"<<endl;
//	system("pause");

// Check to see if steady state conditions have been reached.

    if (dyr >= strteq && endeq == 0)
    {
      if (nfeed == 0 && rheqflag == 0
         && (wtol >= fabs(atms.yrrain + atms.yrsnowfall - atms.yreet
         - soil.yrh2oyld))
         && ((ctol >= fabs(veg.yrnpp - veg.yrltrc))
         || ( 0.001 >= fabs(veg.yrnpp - veg.yrltrc))))
      {
	     endeq = 1;
      }
      if (nfeed == 0 && rheqflag == 1 && (wtol >= fabs(atms.yrrain
         + atms.yrsnowfall - atms.yreet - soil.yrh2oyld))
         && (ctol >= fabs(yrnep))
         && (ctol >= fabs(veg.yrnpp - veg.yrltrc))
         && (ctol >= fabs(veg.yrltrc - microbe.yrrh)))
      {
	     endeq = 1;
      }
      if (nfeed == 1 && rheqflag == 1 && (wtol >= fabs(atms.yrrain
         + atms.yrsnowfall - atms.yreet - soil.yrh2oyld))
         && (ntol >= fabs(soil.yrnin - soil.yrnlost))
         && (ntol >= fabs(veg.yrnup - veg.yrltrn))
         && (ntol >= fabs(veg.yrnup - microbe.yrnmin))
         && (ntol >= fabs(veg.yrltrn - microbe.yrnmin))
         && (ctol >= fabs(yrnep)) && (ctol >= fabs(veg.yrnpp - veg.yrltrc))
         && (ctol >= fabs(veg.yrltrc - microbe.yrrh)))
      {
        endeq = 1;
      }
    }
  }
//  cout <<"OK 12"<<endl;
//  system("pause");

  if (endeq == 2)
  {
    nattempt = maxnrun;
    initFlag = 1;
  }

  if (dyr >= runsize && endeq < 2) { ++nattempt; }

  return nattempt;

};
/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM::getenviron(const int& dm)
{

  // Determine monthly potential evapotranspiration

  atms.pet[dm] = atms.petjh(atms.nirr[dm], atms.tair[dm], dm);
  if (atms.pet[dm] <= 0.0) { atms.pet[dm] = 0.001; }

  // Determine if monthly precipitation occurs as rain or snow

  atms.precsplt(atms.prec[dm],atms.tair[dm],atms.rain[dm],atms.snowfall[dm]);

  // Determine previous two month's air temperatures for following
  //   snowmelt calculations

  switch (dm)
  {
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

  // Determine contribution of snowmelt to soil moisture

  soil.snowinf[dm] = soil.snowmelt(elev, atms.tair[dm],
                                   atms.prevtair, y[I_SNWPCK]);

  monthxclm(veg.cmnt, veg.topt, dm);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(ofstream& rflog1)
{

  char ecd[80];

  cout << "Enter name of the site (.ECD) data file with the parameter values:";
  cout << endl;
// cin >> ecd;
   fpara >> ecd;

  rflog1 << "Enter name of the site (.ECD)"
  		 << " data file with the parameter values:";
  rflog1 << ecd << endl << endl;

  getsitecd(ecd);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(const int& numcmnt, ofstream& rflog1)
{

  int dv;
  char ecd[80];

  cout << "Enter name of the site (.ECD) data"
  	   << " file with the parameter values  for cmnt" << endl;
  rflog1 << "Enter name of the site (.ECD) data file with"
  		<< " the parameter values cmnt" << endl;
  for (dv = 0; dv < numcmnt; dv++)
  {
    cout << (dv+1) << ": ";
//  cin >> ecd;
	fpara >> ecd;

    rflog1 << (dv+1) << ": " << ecd << endl;

    getsitecd(dv, ecd);
  }
  rflog1 << endl;

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(char ecd[80])
{

  const int NUMVAR = 52;
  char dummy[NUMVAR][10];
  ifstream infile;
  int i;
  int dcmnt;

  int sitetveg[MAXCMNT];
  int sitewsoil[MAXCMNT];
  long update[MAXCMNT];
  char sitevegtype[MAXCMNT][31];
  char sitename[MAXCMNT][17];
  char sitetext[MAXCMNT][14];
  float sitecol[MAXCMNT];
  float siterow[MAXCMNT];

  infile.open(ecd, ios::in);

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (dcmnt = 0; dcmnt < MAXCMNT; dcmnt++)
  {
    infile >> sitetveg[dcmnt] >> sitevegtype[dcmnt];
    infile >> sitecol[dcmnt] >> siterow[dcmnt];
    infile >> sitename[dcmnt] >> sitetext[dcmnt] >> sitewsoil[dcmnt];
    infile >> vegca[dcmnt] >> vegcb[dcmnt];
    infile >> strna[dcmnt] >> strnb[dcmnt];
    infile >> solca[dcmnt] >> solcb[dcmnt];
    infile >> solna[dcmnt] >> solnb[dcmnt];
    infile >> avlna[dcmnt] >> avlnb[dcmnt];
    infile >> stona[dcmnt] >> stonb[dcmnt];
    infile >> veg.unleaf12[dcmnt];
    infile >> veg.initleafmx[dcmnt];      // Changed by DWK on 19991028
    infile >> veg.cmaxcut[dcmnt];
    infile >> veg.cmax1a[dcmnt] >> veg.cmax1b[dcmnt];
    infile >> veg.cmax2a[dcmnt] >> veg.cmax2b[dcmnt];
    infile >> veg.cfall[dcmnt];
    infile >> veg.kra[dcmnt] >> veg.krb[dcmnt];
    infile >> microbe.kda[dcmnt] >> microbe.kdb[dcmnt];
    infile >> microbe.lcclnc[dcmnt] >> microbe.propftos[dcmnt];
    infile >> veg.nmaxcut[dcmnt];
    infile >> veg.nmax1a[dcmnt] >> veg.nmax1b[dcmnt];
    infile >> veg.nmax2a[dcmnt] >> veg.nmax2b[dcmnt];
    infile >> veg.nfall[dcmnt];
    infile >> microbe.nupa[dcmnt] >> microbe.nupb[dcmnt];
    infile >> soil.nloss[dcmnt];
    infile >> microbe.nfixpar[dcmnt];
    infile >> veg.initcneven[dcmnt] >> veg.cnmin[dcmnt];
    infile >> veg.c2na[dcmnt] >> veg.c2nb[dcmnt] >> veg.c2nmin[dcmnt];
    infile >> microbe.cnsoil[dcmnt];
    infile >> update[dcmnt];

//    veg.initcneven[i] = veg.cneven[i];
//    veg.adjc2n = 1.0 + (veg.dc2n * (atms.co2[11] - atms.initco2));
//    veg.cneven[i] = veg.initcneven[i] * veg.adjc2n;
  }

  infile.close();

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::getsitecd(const int& dv, char ecd[80])
{
  // Function added by DWK on 20000102

  char dummy[12];
  // string changed to char[80] by DWK on 20000210
  char sitename[80];
  float sitecol;
  float siterow;
  long updated;

  fecd[dv].open(ecd,ios::in);

  if (!fecd[dv])
  {
    cerr << endl << "Cannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  fecd[dv] >> dummy >> veg.cmnt;
  fecd[dv] >> dummy >> veg.cmnt_name;
  fecd[dv] >> dummy >> sitename;
  fecd[dv] >> dummy >> sitecol;
  fecd[dv] >> dummy >> siterow;
  fecd[dv] >> dummy >> updated;

  fecd[dv] >> dummy >> vegca[veg.cmnt];
  fecd[dv] >> dummy >> vegcb[veg.cmnt];
  fecd[dv] >> dummy >> strna[veg.cmnt];
  fecd[dv] >> dummy >> strnb[veg.cmnt];
  fecd[dv] >> dummy >> solca[veg.cmnt];
  fecd[dv] >> dummy >> solcb[veg.cmnt];
  fecd[dv] >> dummy >> solna[veg.cmnt];
  fecd[dv] >> dummy >> solnb[veg.cmnt];
  fecd[dv] >> dummy >> avlna[veg.cmnt];
  fecd[dv] >> dummy >> avlnb[veg.cmnt];
  fecd[dv] >> dummy >> stona[veg.cmnt];
  fecd[dv] >> dummy >> stonb[veg.cmnt];

  fecd[dv] >> dummy >> veg.unleaf12[veg.cmnt];
  // veg.prevleafmx changed to veg.initleafmx by DWK on 20000130
  fecd[dv] >> dummy >> veg.initleafmx[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmaxcut[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax1a[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax1b[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax2a[veg.cmnt];
  fecd[dv] >> dummy >> veg.cmax2b[veg.cmnt];
  fecd[dv] >> dummy >> veg.cfall[veg.cmnt];
  fecd[dv] >> dummy >> veg.kra[veg.cmnt];
  fecd[dv] >> dummy >> veg.krb[veg.cmnt];
  fecd[dv] >> dummy >> microbe.kda[veg.cmnt];
  fecd[dv] >> dummy >> microbe.kdb[veg.cmnt];
  fecd[dv] >> dummy >> microbe.lcclnc[veg.cmnt];
  fecd[dv] >> dummy >> microbe.propftos[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmaxcut[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax1a[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax1b[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax2a[veg.cmnt];
  fecd[dv] >> dummy >> veg.nmax2b[veg.cmnt];
  fecd[dv] >> dummy >> veg.nfall[veg.cmnt];
  fecd[dv] >> dummy >> microbe.nupa[veg.cmnt];
  fecd[dv] >> dummy >> microbe.nupb[veg.cmnt];
  fecd[dv] >> dummy >> soil.nloss[veg.cmnt];
  fecd[dv] >> dummy >> microbe.nfixpar[veg.cmnt];
  fecd[dv] >> dummy >> veg.initcneven[veg.cmnt];
  fecd[dv] >> dummy >> veg.cnmin[veg.cmnt];
  fecd[dv] >> dummy >> veg.c2na[veg.cmnt];
  fecd[dv] >> dummy >> veg.c2nb[veg.cmnt];
  fecd[dv] >> dummy >> veg.c2nmin[veg.cmnt];
  fecd[dv] >> dummy >> microbe.cnsoil[veg.cmnt];

  fecd[dv].close();

};

/* *************************************************************
************************************************************** */
void TTEM::initrun(Temstac& temstac)
{
	avlnflag=temstac.avlnflag;
	nfeed=temstac.nfeed;
	initbase=temstac.initbase;
	baseline=temstac.baseline;
	
	moistlim=temstac.moistlim;
	strteq=temstac.strteq;
	maxyears=temstac.maxyears;
	runsize=temstac.runsize;
	maxnrun=temstac.maxnrun;
	rheqflag=temstac.rheqflag;
		
	
	wtol=temstac.wtol;
	ctol=temstac.ctol;
 	ntol=temstac.ntol;
 	startyr=temstac.startyr;
 	endyr=temstac.endyr;
 	diffyr=temstac.diffyr;
    	

	
};


/* **************************************************************
************************************************************** */

void TTEM::initrun(ofstream& rflog1, const int& equil,Temstac& temstac)
{

  avlnflag = nfeed = rheqflag = 0;

/* **************************************************************
		  Run Model with Nitrogen Limitation?
************************************************************** */

  cout << endl << "Do you want to allow available N to fluctuate?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes: ";
//cin >> avlnflag;
  fpara >> avlnflag;
	temstac.avlnflag=avlnflag;
  rflog1 << endl << "Do you want to allow available N to fluctuate?" << endl;
  rflog1 << "  Enter 0 for No" << endl;
  rflog1 << "  Enter 1 for Yes: " << endl;
  rflog1 << "avlnflag = " << avlnflag << endl << endl;

  cout << endl << "Do you want nitrogen feedback on GPP?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes: ";
//cin >> nfeed;
  fpara >> nfeed;
	temstac.nfeed=nfeed;
  rflog1 << endl << "Do you want nitrogen feedback on GPP?" << endl;
  rflog1 << "  Enter 0 for No" << endl;
  rflog1 << "  Enter 1 for Yes: " << endl;
  rflog1 << "nfeed = " << nfeed << endl << endl;

  baseline = initbase = 0;
  if (nfeed == 1)
  {
    cout << endl << "Do you want to solve for baseline soil nitrogen?" << endl;
    cout << "  Enter 0 for No" << endl;
    cout << "  Enter 1 for Yes: ";
//	cin >> initbase;
    fpara >> initbase;
    baseline = initbase;

	rflog1 << endl << "Do you want to solve for"
		   << " baseline soil nitrogen?" << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "baseline = " << baseline << endl << endl;
  }
  temstac.initbase=initbase;
  temstac.baseline=baseline;

/* **************************************************************
			 Run Model with Moisture Limitation?
************************************************************** */

  moistlim = 0;
  cout << endl << "Do you want to run the model with"
  	   << " moisture limitation?" << endl;
  cout << "  Enter 0 for No" << endl;
  cout << "  Enter 1 for Yes: ";
//cin >> moistlim;
  fpara >> moistlim;
	temstac.moistlim=moistlim;

  rflog1 << endl << "Do you want to run the model"
  		 << " with moisture limitation?" << endl;
  rflog1 << "  Enter 0 for No" << endl;
  rflog1 << "  Enter 1 for Yes: " << endl;
  rflog1 << "moistlim = " << moistlim << endl << endl;


/* ***************************************************************
	       Details for Steady State Conditions
************************************************************** */


  maxyears = 0;
  maxnrun = 0;

  cout << endl << "How many years do you want to wait"
  		<< " before checking equilibrium conditions? ";
//cin >> strteq;
  fpara >> strteq;
	temstac.strteq=strteq;
  rflog1 << endl;
  rflog1 << "How many years do you want to wait before"
  		 << " checking equilibrium conditions? ";
  rflog1 << endl;
  rflog1 << "strteq = " << strteq << endl << endl;

  cout << endl << "Enter the maximum number of years for the model to run: ";
//cin >> maxyears;
  fpara >> maxyears;
	temstac.maxyears=maxyears;

  rflog1 << endl << "Enter the maximum number of years for the model to run: ";
  rflog1 << endl;
  rflog1 << "maxyears = " << maxyears << endl << endl;

  runsize = maxyears;
  temstac.runsize=runsize;


  cout << endl << "Enter the maximum number of attempts to reach a solution: ";
//cin >> maxnrun;
  fpara >> maxnrun;
	temstac.maxnrun=maxnrun;

  rflog1 << endl;
  rflog1 << "Enter the maximum number of attempts to reach a solution: ";
  rflog1 << endl;
  rflog1 << "maxnrun = " << maxnrun << endl << endl;

  if (nfeed == 0)
  {
    cout << endl << "Do you want decomposition to come into equilibrium? ";
    cout << "  Enter 0 for No" << endl;
    cout << "  Enter 1 for Yes: ";
//  cin >> rheqflag;
	fpara >> rheqflag;

    rflog1 << endl;
    rflog1 << "Do you want decomposition to come into equilibrium? " << endl;
    rflog1 << "  Enter 0 for No" << endl;
    rflog1 << "  Enter 1 for Yes: " << endl;
    rflog1 << "rheqflag = " << rheqflag << endl << endl;
  }
	temstac.rheqflag=rheqflag;

  wtol = 1000.0;
  cout << endl;
  cout << "What absolute tolerance do you want to"
  	   << " use for checking equilibrium";
  cout << endl;
  cout << "of the water cycle? ";
//cin >> wtol;
  fpara >> wtol;
	temstac.wtol=wtol;

  rflog1 << endl;
  rflog1 << "What absolute tolerance do you want to"
         << " use for checking equilibrium";
  rflog1 << endl;
  rflog1 << "of the water cycle? wtol = " << wtol << endl << endl;

  ctol = 1000.0;
  cout << endl;
  cout << "What absolute tolerance do you want to"
  	   << " use for checking equilibrium";
  cout << endl;
  cout << "of the carbon cycle? ";
//cin >> ctol;
  fpara >> ctol;
	temstac.ctol=ctol;

  rflog1 << endl;
  rflog1 << "What absolute tolerance do you want to use"
  		 << " for checking equilibrium";
  rflog1 << endl;
  rflog1 << "of the carbon cycle?" << endl;
  rflog1 << "ctol = " << ctol << endl << endl;

  ntol = 1000.0;
  if (nfeed == 1)
  {
    rheqflag = 1;
	temstac.rheqflag=rheqflag;
    cout << endl;
	cout << "What absolute tolerance do you want to"
		 << " use for checking equilibrium";
    cout << endl;
    cout << "of the nitrogen cycle? ";
//  cin >> ntol;
	fpara >> ntol;

    rflog1 << endl;
	rflog1 << "What absolute tolerance do you want to"
		   << " use for checking equilibrium";
    rflog1 << endl;
    rflog1 << "of the nitrogen cycle?" << endl;
    rflog1 << "ntol = " << ntol << endl << endl;
  }
	temstac.ntol=ntol;
  if (equil == 0)
  {

    cout << endl << endl;
    cout << "What year do you want to start collecting output data? ";
//  cin >> startyr;
	fpara >> startyr;

    rflog1 << endl << endl;
    rflog1 << "What year do you want to start collecting output data? ";
    rflog1 << "startyr = " << startyr << endl;

    cout << endl << endl;
    cout << "What year do you want to stop collecting output data? ";
//  cin >> endyr;
	fpara >> endyr;

    rflog1 << endl << endl;
    rflog1 << "What year do you want to stop collecting output data? ";
    rflog1 << "endyr = " << endyr << endl;

	cout << "How often (x years) should data be"
	     << " collected after the initial year? ";
//  cin >> diffyr;
	fpara >> diffyr;

	rflog1 << "How often (x years) should data be collected"
		   << " after the initial year? ";
    rflog1 << "diffyr = " << diffyr << endl;

  }
	temstac.startyr=startyr;
	temstac.endyr=endyr;
	temstac.diffyr=diffyr;	  

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM::massbal(double y[NUMEQ], double prevy[NUMEQ])
{


  if (y[I_SM] != y[I_AVLW] + soil.wiltpt)
  {
    y[I_SM] = y[I_AVLW] + soil.wiltpt;
  }

  if (y[I_PCTP] != 100.0 * y[I_SM] / soil.totpor)
  {
    y[I_PCTP] = 100.0 * y[I_SM] / soil.totpor;
  }

  if (y[I_VSM] != y[I_SM] / (soil.rootz * 1000.0))
  {
    y[I_VSM] = y[I_SM] / (soil.rootz * 1000.0);
    if (y[I_VSM] <= 0.0) { y[I_VSM] = 0.001; }
  }

  if ((y[I_SNWPCK] - prevy[I_SNWPCK]) != (y[I_SNWFAL] - y[I_SNWINF]))
  {
    y[I_SNWINF] = y[I_SNWFAL] - y[I_SNWPCK] + prevy[I_SNWPCK];
  }

  if ((y[I_AVLW] - prevy[I_AVLW]) != (y[I_SNWINF] + y[I_RAIN]
       - y[I_RPERC] - y[I_EET] - y[I_SPERC]))
  {
    y[I_SPERC] = y[I_SNWINF] + y[I_RAIN] - y[I_RPERC] - y[I_EET]
                 - y[I_AVLW] + prevy[I_AVLW];
  }

  if ((y[I_RGRW] - prevy[I_RGRW]) != (y[I_RPERC] - y[I_RRUN]))
  {
    y[I_RRUN] = y[I_RPERC] - y[I_RGRW] + prevy[I_RGRW];
  }

  if ((y[I_SGRW] - prevy[I_SGRW]) != (y[I_SPERC] - y[I_SRUN]))
  {
    y[I_SRUN] = y[I_SPERC] - y[I_SGRW] + prevy[I_SGRW];
  }

  if (y[I_WYLD] != y[I_RRUN] + y[I_SRUN])
  {
    y[I_WYLD] = y[I_RRUN] + y[I_SRUN];
  }
/************************* Carbon Cycle Balances **************************/
  if (y[I_INNPP] < y[I_NPP]) { y[I_INNPP] = y[I_NPP]; }

  if (y[I_INGPP] < y[I_GPP]) { y[I_INGPP] = y[I_GPP]; }

  if (y[I_GPR] != y[I_GPP] - y[I_NPP])
  {
    y[I_GPR] = y[I_GPP] - y[I_NPP];
  }

  if (y[I_GPR] != y[I_RVMNT] + y[I_RVGRW])
  {
    y[I_RVGRW] = y[I_GPR] - y[I_RVMNT];
  }

  if (ag.state == 0)
  {
    // Minimum VEGC added by DWK on 20000207
    //  if (y[I_VEGC] < 0.00001) { y[I_VEGC] = 0.00001; }
    if (y[I_VEGC] - prevy[I_VEGC] != y[I_NPP] - y[I_LTRC])
    {
      y[I_LTRC] = y[I_NPP] - y[I_VEGC] + prevy[I_VEGC];
    }
  }
  else
  {
    if (y[I_AGNPPC] != y[I_AGLTRC] + y[I_AGFPRDC])
    {
      y[I_AGLTRC] = y[I_AGNPPC] - y[I_AGFPRDC];
    }
  }

  // Minimum SOLC added by DWK on 20000207
  //  if (y[I_SOLC] < 0.00001) { y[I_SOLC] = 0.00001; }
  if (y[I_SOLC] - prevy[I_SOLC] != y[I_LTRC] + y[I_AGLTRC]
      + y[I_SLASHC] - y[I_SCNVRTC] - y[I_RH])
  {
    y[I_RH] = y[I_LTRC] + y[I_AGLTRC] + y[I_SLASHC] - y[I_SCNVRTC]
    - y[I_SOLC] + prevy[I_SOLC];
  }

  if (y[I_NEP] != y[I_NPP] + y[I_AGNPPC] - y[I_RH])
  {
    y[I_NEP] = y[I_NPP] + y[I_AGNPPC] - y[I_RH];
  }

  if (y[I_CFLX] != y[I_NEP] - y[I_CNVRTC])
  {
    y[I_CFLX] = y[I_NEP] - y[I_CNVRTC];
  }

  /*********************Nitrogen Cycle Balances**********************/

  if (y[I_VNUP] < 0.0 ) { y[I_VNUP] = 0.0; }

  //if (y[I_VNUP] > y[I_INNUP]) { y[I_VNUP] = y[I_INNUP]; }
  if (y[I_INNUP] < y[I_VNUP]) { y[I_INNUP] = y[I_VNUP]; }

  if (y[I_VSUP] < 0.0) { y[I_VSUP] = 0.0; }

  if (y[I_VSUP] > y[I_VNUP]) { y[I_VSUP] = y[I_VNUP]; }

  if (y[I_VLUP] != y[I_VNUP] - y[I_VSUP])
  {
    y[I_VLUP] = y[I_VNUP] - y[I_VSUP];
  }

  if (ag.state == 0)
  {
    if (y[I_STON] - prevy[I_STON] != y[I_VLUP] + y[I_VNRSRB] - y[I_VNMBL])
    {
      y[I_VNRSRB] = y[I_STON] - prevy[I_STON] + y[I_VNMBL] - y[I_VLUP];
    }

    if (y[I_STRN] - prevy[I_STRN] != y[I_VSUP] - y[I_LTRN]
        - y[I_VNRSRB] + y[I_VNMBL])
    {
      y[I_LTRN] = y[I_VSUP] - y[I_STRN] + prevy[I_STRN]
                  - y[I_VNRSRB] + y[I_VNMBL];
    }
  }
  else
  {
    if (y[I_AGNPPN] != y[I_AGLTRN] + y[I_AGFPRDN])
    {
      y[I_AGLTRN] = y[I_AGNPPN] - y[I_AGFPRDN];
    }
  }

  if (y[I_SOLN] - prevy[I_SOLN] != y[I_LTRN] + y[I_AGLTRN] + y[I_SLASHN]
      - y[I_NMIN] - y[I_SCNVRTN] - y[I_NSRTNT])
  {
    y[I_NMIN] = y[I_LTRN] + y[I_AGLTRN] + y[I_SLASHN] - y[I_SCNVRTN]
                - y[I_NSRTNT] - y[I_SOLN] + prevy[I_SOLN];
  }


  if (y[I_NLST] < 0.0) { y[I_NLST] = 0.0; }

  if (ag.state == 0)
  {
    if (y[I_NINP] < 0.0) { y[I_NINP] = 0.0; }

    if (y[I_AVLN] - prevy[I_AVLN] != y[I_NINP] - y[I_NLST] + y[I_NMIN]
        - y[I_VNUP])
    {
      if (y[I_NINP] + y[I_NMIN] - y[I_VNUP] - y[I_AVLN] + prevy[I_AVLN] > 0.0)
      {
        y[I_NLST] =  y[I_NINP] + y[I_NMIN] - y[I_VNUP] - y[I_AVLN]
                     + prevy[I_AVLN];
      }
      else
      {
        y[I_NINP] = y[I_NLST] - y[I_NMIN] + y[I_VNUP] + y[I_AVLN]
                    - prevy[I_AVLN];
      }
    }
  }

  if(ag.state == 1)
  {

// Following condition added to correct for negative AGFERTN
// D. Kicklighter 19990721

    if (y[I_NRETNT] != y[I_NVRTNT] + y[I_NSRTNT])
    {
      y[I_NRETNT] = y[I_NVRTNT] + y[I_NSRTNT];
    }

    if (y[I_AGFRTN] < 0.0) { y[I_AGFRTN] = 0.0; }

    if (y[I_NINP] != y[I_AGFRTN] + y[I_NRETNT])
    {
      y[I_NINP] = y[I_AGFRTN] + y[I_NRETNT];
    }

    if (y[I_NINP] < 0.0) { y[I_NINP] = 0.0; }

    if (y[I_AVLN] - prevy[I_AVLN] != y[I_NINP] - y[I_NLST] + y[I_NMIN]
       - y[I_AGNPPN])
    {
      y[I_NLST] =  prevy[I_AVLN] - y[I_AVLN] + y[I_NINP] + y[I_NMIN]
                   - y[I_AGNPPN];
    }
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::monthxclm(const int& dcmnt, const double& tgppopt, const int& dm)
{

double raq10;


/* temp: effect of temperature on primary productivity */

  if (atms.tair[dm] <= veg.tmin[dcmnt] || atms.tair[dm] >= veg.tmax[dcmnt])
  {
    temp = 0.0;
  }
  else
  {
    if (atms.tair[dm] >= tgppopt && atms.tair[dm] <= veg.toptmax[dcmnt])
    {
      temp = 1.0;
    }
    else
    {
      if (atms.tair[dm] > veg.tmin[dcmnt] && atms.tair[dm] < tgppopt)
      {
	temp = (atms.tair[dm]-veg.tmin[dcmnt])*(atms.tair[dm]-veg.tmax[dcmnt])
               /((atms.tair[dm]-veg.tmin[dcmnt])*(atms.tair[dm]-veg.tmax[dcmnt])
               - pow((atms.tair[dm]-tgppopt),2.0));
      }
      else
      {
	temp = (atms.tair[dm]-veg.tmin[dcmnt])*(atms.tair[dm]-veg.tmax[dcmnt])
               /((atms.tair[dm]-veg.tmin[dcmnt])*(atms.tair[dm]-veg.tmax[dcmnt])
               - pow((atms.tair[dm]-veg.toptmax[dcmnt]),2.0));
      }
    }
  }


/* respq10: effect of temperature on plant respiration  */

  raq10 = veg.raq10a0[dcmnt];
  //raq10 = veg.raq10a0[dcmnt] + (veg.raq10a1[dcmnt]*atms.tair[dm])
  //        + (veg.raq10a2[dcmnt]*pow(atms.tair[dm],2.0))
  //       + (veg.raq10a3[dcmnt]*pow(atms.tair[dm],3.0));
  respq10 = pow(raq10,atms.tair[dm]/10.0);


/* dq10: effect of temperature on decomposition */

// 19990821-previous year snowpack (soil.snowpack[dm]) changed
//          to previous month snowpack (y[NUMEEQ+2]) by Jim Long



 if (tmflgs.stmflg==1) {

	dq10 = pow(microbe.rhq10[dcmnt],atms.tsoil[dm]/10.0);
	// Replace air Temperature with Soil Temperature
    }
  else
   {

//  if (y[I_SNWPCK] > 0.0) { dq10 = 1.0; }
//  else
 // {
     dq10 = pow(microbe.rhq10[dcmnt],atms.tair[dm]/10.0);
   }

//  if (y[I_SNWPCK] > 0.0) { dq10 = 1.0; }
//  else
//  {
//    dq10 = pow(microbe.rhq10[dcmnt],atms.tair[dm]/10.0);
//  }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

void TTEM::resetODEflux(double y[])
{
  int i;

  for (i = MAXSTATE; i < NUMEQ; i++) { y[i] = 0.0; }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM::resetYrFlux(void)
{

  int dm;

  // Annual carbon storage
  veg.yrcarbon = 0.0;
  soil.yrorgc = 0.0;

  // Annual nitrogen storage
  veg.yrnitrogen = 0.0;
  veg.yrstructn = 0.0;
  veg.yrc2n = 0.0;
  veg.yrstoren = 0.0;
  soil.yrorgn = 0.0;
  soil.yrc2n = 0.0;
  soil.yravln = 0.0;

  // Annual carbon & nitrogen storage in agricultural ecosystems
  ag.PROD1.carbon = 0.0;
  ag.PROD1.nitrogen = 0.0;

  // Annual water storage
  soil.yravlh2o = 0.0;
  soil.yrrgrndh2o = 0.0;
  soil.yrsnowpack = 0.0;
  soil.yrsgrndh2o = 0.0;
  soil.yrsmoist = 0.0;
  soil.yrpctp = 0.0;
  soil.yrvsm = 0.0;

  // Annual carbon fluxes
  veg.yringpp = 0.0;
  veg.yrgpp = 0.0;
  veg.yrinnpp = 0.0;
  veg.yrnpp = 0.0;
  veg.yrltrc = 0.0;
  microbe.yrrh = 0.0;
  yrnep = 0.0;

  // Annual nitrogen fluxes
  soil.yrnin = 0.0;
  veg.yrinnup = 0.0;
  veg.yrnup = 0.0;
  veg.yrsup = 0.0;
  veg.yrlup = 0.0;
  veg.yrnmobil = 0.0;
  veg.yrnrsorb = 0.0;
  veg.yrltrn = 0.0;
  microbe.yrnmin = 0.0;
  soil.yrnlost = 0.0;

  // Annual water fluxes
  atms.yrrain = 0.0;
  soil.yrrperc = 0.0;
  soil.yrrrun = 0.0;
  atms.yrsnowfall = 0.0;
  soil.yrsnowinf = 0.0;
  soil.yrsperc = 0.0;
  soil.yrsrun = 0.0;
  atms.yrpet = 0.0;
  atms.yreet = 0.0;
  soil.yrh2oyld = 0.0;

     // for soil temperature
  atms.yrfrontd =0.0;
  atms.yrthawbegin =0.0;
  atms.yrthawend =0.0;
  atms.yrtsoil = 0.0;
  atms.yrdst5 = 0.0;
  atms.yrdst10 = 0.0;
  atms.yrdst20 = 0.0;
  atms.yrdst50 = 0.0;
  atms.yrdst100 = 0.0;
  atms.yrdst200 = 0.0;


  // Phenology
  veg.yrunleaf = 0.0;
  veg.yrleaf = 0.0;
  veg.yrfpc = 0.0;

  // Annual carbon and nitrogen fluxes from agricultural
  // conversion
  ag.yrconvrtC = 0.0;
  ag.yrvconvrtC = 0.0;
  ag.yrsconvrtC = 0.0;
  ag.yrconvrtN = 0.0;
  ag.yrvconvrtN = 0.0;
  ag.yrsconvrtN = 0.0;
  ag.yrnrent = 0.0;
  ag.yrnvrent = 0.0;
  ag.yrnsrent = 0.0;
  ag.yrslashC = 0.0;
  ag.yrslashN = 0.0;
  ag.formPROD10.carbon  = 0.0;
  ag.formPROD10.nitrogen  = 0.0;
  ag.formPROD100.carbon = 0.0;
  ag.formPROD100.nitrogen = 0.0;

  // Annual carbon and nitrogen fluxes from agriculture
  ag.yrnppC = 0.0;
  ag.yrnppN = 0.0;
  ag.formPROD1.carbon = 0.0;
  ag.formPROD1.nitrogen = 0.0;
  ag.yrltrc = 0.0;
  ag.yrltrn = 0.0;
  ag.yrfertn = 0.0;

  // Annual carbon and nitrogen fluxes from human products
  ag.formTOTPROD.carbon = 0.0;
  ag.formTOTPROD.nitrogen = 0.0;

  // Monthly carbon and nitrogen fluxes
  for (dm = 0; dm < CYCLE; dm++)
  {
    soil.ninput[dm] = 0.0;
    soil.nlost[dm] = 0.0;

    ag.slash[dm].carbon = 0.0;
    ag.slash[dm].nitrogen = 0.0;
    ag.vconvrtflx[dm].carbon = 0.0;
    ag.sconvrtflx[dm].carbon = 0.0;
    ag.convrtflx[dm].carbon = 0.0;
    ag.vconvrtflx[dm].nitrogen = 0.0;
    ag.sconvrtflx[dm].nitrogen = 0.0;
    ag.convrtflx[dm].nitrogen = 0.0;
    ag.nvretent[dm] = 0.0;
    ag.nsretent[dm] = 0.0;
    ag.nretent[dm] = 0.0;
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM::setELMNTecd(const int& kdinflg, const int& dcmnt,
                       const double& psiplusc)
{


// Initialize TEM parameters dependent upon a grid cell's soil texture

  if (psiplusc <= veg.cmaxcut[dcmnt])
  {
    veg.cmax = (veg.cmax1a[dcmnt] * psiplusc) + veg.cmax1b[dcmnt];
  }
  else
  {
    veg.cmax = (veg.cmax2a[dcmnt] * psiplusc) + veg.cmax2b[dcmnt];
  }

  if (kdinflg == 0)
  {
    microbe.kdc = (microbe.kda[dcmnt] / psiplusc) + microbe.kdb[dcmnt];
  }

  if (psiplusc <= veg.nmaxcut[dcmnt])
  {
    veg.nmax = (veg.nmax1a[dcmnt] * psiplusc) + veg.nmax1b[dcmnt];
  }
  else
  {
    veg.nmax = (veg.nmax2a[dcmnt] * psiplusc) + veg.nmax2b[dcmnt];
  }

  microbe.nup = (microbe.nupa[dcmnt] / psiplusc) + microbe.nupb[dcmnt];

// Set initial maximum relative leaf area

  veg.prvleafmx[dcmnt] = veg.initleafmx[dcmnt];  // Changed by DWK on 19991028


// Determine the "decay" parameter

  microbe.decay = 0.26299 + (1.14757*microbe.propftos[dcmnt])
                  - (0.42956*pow(microbe.propftos[dcmnt],2.0));

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM::setELMNTevap(const int& stateflg, const int& dcmnt,
                        double pet[CYCLE], double tair[CYCLE])
{

  int i;

// Determine initial values for atms.prvpetmx, atms.prveetmx,
//   and veg.topt

  if (stateflg == 1)
  {
    for (i = 0; i < CYCLE; i++)
    {
      if (pet[i] >= atms.prvpetmx)
      {
	     veg.topt = tair[i];
      }
      if (veg.aleaf[dcmnt] == 0.0 && veg.bleaf[dcmnt] == 0.0
         && veg.cleaf[dcmnt] == 1.0)
      {
	     if (tair[i] > veg.topt) { veg.topt = tair[i]; }
      }
      else
      {
	     if (veg.unnormleaf[i] >= veg.prvleafmx[dcmnt])
        {
	       veg.topt = tair[i];
	     }
      }
    }

    if (veg.topt > veg.toptmax[dcmnt]) { veg.topt = veg.toptmax[dcmnt]; }
    if (veg.topt < veg.toptmin[dcmnt]) { veg.topt = veg.toptmin[dcmnt]; }
  }

  else
  {
    for (i = 0; i < CYCLE; i++)
    {
      atms.pet[i] = atms.petjh(atms.nirr[i], tair[i], i);
      if (i == 0)
      {
	atms.prvpetmx = atms.prveetmx = atms.pet[0];
	veg.topt  = tair[0];
      }
      else
     {
	if (atms.pet[i] > atms.prvpetmx)
        {
	  atms.prvpetmx = atms.prveetmx = atms.pet[i];
	  veg.topt = tair[i];
	}
      }
    }

    atms.yrpet = 1.0;
    atms.yreet = 1.0;
  }


// Determine vegetation C/N parameter as a function of vegetation type,
// annual PET, and annual EET (annual EET initially equal to yrpet)

  veg.updateC2N(dcmnt,atms.yreet,atms.yrpet,atms.co2[11],atms.initco2);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void TTEM::setELMNTflux(void)
{

  // Initialize carbon, nitrogen and water fluxes including
  // ODE state variables (i.e., y[])
  int dm;

  for (dm = 0; dm < CYCLE; dm++)
  {
    // Initialize carbon fluxes in natural ecosystems to zero

    y[I_INGPP] = veg.ingpp[dm] = 0.0;
    y[I_GPP] = veg.gpp[dm] = 0.0;
    y[I_INNPP] = veg.innpp[dm] = 0.0;
    y[I_NPP] = veg.npp[dm] = 0.0;
    y[I_GPR] = veg.gpr[dm] = 0.0;
    y[I_RVMNT] = veg.rm[dm] = 0.0;
    y[I_RVGRW] = veg.rg[dm] = 0.0;
    y[I_LTRC] = veg.ltrfal[dm].carbon = 0.0;
    y[I_RH] = microbe.rh[dm] = 0.0;
    y[I_NEP] = nep[dm] = 0.0;

    // Initialize nitrogen fluxes in natural ecosystems to zero

    y[I_NINP] = soil.ninput[dm] = 0.0;
    y[I_INNUP] = veg.inuptake = 0.0;
    y[I_VNUP] = veg.nuptake[dm] = 0.0;
    y[I_VSUP] = veg.suptake[dm] = 0.0;
    y[I_VLUP] = veg.luptake[dm] = 0.0;
    y[I_VNMBL] = veg.nmobil[dm] = 0.0;
    y[I_VNRSRB] = veg.nresorb[dm] = 0.0;
    y[I_LTRN] = veg.ltrfal[dm].nitrogen = 0.0;
    y[I_MNUP] = microbe.nuptake[dm] = 0.0;
    y[I_NMIN] = microbe.netnmin[dm] = 0.0;
    y[I_NLST] = soil.nlost[dm] = 0.0;

    // Initialize water fluxes to zero

    y[I_RAIN] = 0.0;
    y[I_RPERC] = soil.rperc[dm] = 0.0;
    y[I_RRUN] = soil.rrun[dm] = 0.0;
    y[I_SNWFAL] = 0.0;
    y[I_SNWINF] = soil.snowinf[dm] = 0.0;
    y[I_SPERC] = soil.sperc[dm] = 0.0;
    y[I_SRUN] = soil.srun[dm] = 0.0;
    y[I_PET] = 0.0;
    y[I_EET] = atms.eet[dm] = 0.0;
    y[I_WYLD] = soil.h2oyld[dm] = 0.0;

// for soil temperature
  y[I_TSOIL] = atms.tsoil[dm]=0.0;
  y[I_DST5] = atms.dst5[dm]=0.0;
  y[I_DST10] = atms.dst10[dm]=0.0;
  y[I_DST20] = atms.dst20[dm]=0.0;
  y[I_DST50] = atms.dst50[dm]=0.0;
  y[I_DST100] = atms.dst100[dm]=0.0;
  y[I_DST200] = atms.dst200[dm]=0.0;
  y[I_FRONTD] = atms.frontd[dm]=0.0;
  y[I_THAWBE] = atms.thawbe[dm]=0.0;
  y[I_THAWEND] = atms.thawend[dm]=0.0;
// end of ...


    // Initialize carbon and nitrogen fluxes during conversion
    //  to zero

    y[I_CNVRTC] = ag.convrtflx[dm].carbon = 0.0;
    y[I_CNVRTN] = ag.convrtflx[dm].nitrogen = 0.0;
    y[I_SCNVRTC] = ag.sconvrtflx[dm].carbon = 0.0;
    y[I_SCNVRTN] = ag.sconvrtflx[dm].nitrogen = 0.0;
    y[I_NVRTNT] = ag.nvretent[dm] = 0.0;
    y[I_NSRTNT] = ag.nsretent[dm] = 0.0;
    y[I_NRETNT] = ag.nretent[dm] = 0.0;
    y[I_SLASHC] = ag.slash[dm].carbon = 0.0;
    y[I_SLASHN] = ag.slash[dm].nitrogen = 0.0;
    y[I_PRDF10C] = ag.formPROD10.carbon = 0.0;
    y[I_PRDF10N] = ag.formPROD10.nitrogen = 0.0;
    y[I_PRDF100C] = ag.formPROD100.carbon = 0.0;
    y[I_PRDF100N] = ag.formPROD100.nitrogen = 0.0;

    // Initialize carbon and nitrogen in agricultural ecosystems
    //   to zero

    y[I_AGNPPC] = ag.npp[dm].carbon = 0.0;
    y[I_AGNPPN] = ag.npp[dm].nitrogen = 0.0;
    y[I_AGFPRDC] = ag.formPROD1.carbon = 0.0;
    y[I_AGFPRDN] = ag.formPROD1.nitrogen = 0.0;
    y[I_AGLTRC] = ag.ltrfal[dm].carbon = 0.0;
    y[I_AGLTRN] = ag.ltrfal[dm].nitrogen = 0.0;
    y[I_AGFRTN] = ag.fertn[dm] = 0.0;

    // Initialize carbon and nitrogen from human products to zero

    y[I_TOTFPRDC] = ag.formTOTPROD.carbon = 0.0;
    y[I_TOTFPRDN] = ag.formTOTPROD.nitrogen = 0.0;
    y[I_AGPRDFC] = ag.PROD1decay.carbon = 0.0;
    y[I_AGPRDFN] = ag.PROD1decay.nitrogen = 0.0;
    y[I_PRD10FC] = ag.PROD10decay.carbon = 0.0;
    y[I_PRD10FN] = ag.PROD10decay.nitrogen = 0.0;
    y[I_PRD100FC] = ag.PROD100decay.carbon = 0.0;
    y[I_PRD100FN] = ag.PROD100decay.nitrogen = 0.0;
    y[I_TOTPRDFC] = ag.TOTPRODdecay.carbon = 0.0;
    y[I_TOTPRDFN] = ag.TOTPRODdecay.nitrogen = 0.0;

    // Initialize integrated carbon fluxes to zero

    y[I_TOTNPP] = ag.totnpp[dm] = 0.0;
    y[I_CFLX] = cflux[dm] = 0.0;

  }
};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************* */

void TTEM::setMonth(int& dm, double y[])
{

  // Carbon pools
  veg.plant[dm].carbon = y[I_VEGC];
  soil.org[dm].carbon = y[I_SOLC];
  totalc[dm] = veg.plant[dm].carbon + soil.org[dm].carbon;

  // Nitrogen pools
  veg.strctrl[dm].nitrogen = y[I_STRN];
  veg.labile[dm].nitrogen = y[I_STON];
  veg.plant[dm].nitrogen = veg.strctrl[dm].nitrogen + veg.labile[dm].nitrogen;
  soil.org[dm].nitrogen = y[I_SOLN];
  soil.availn[dm] = y[I_AVLN];

  // Water pools
  soil.avlh2o[dm] = y[I_AVLW];
  soil.rgrndh2o[dm] = y[I_RGRW];
  soil.snowpack[dm] = y[I_SNWPCK];
  soil.sgrndh2o[dm] = y[I_SGRW];
  soil.moist[dm] = y[I_SM];
  soil.pctp[dm] = y[I_PCTP];
  soil.vsm[dm] = y[I_VSM];

  // Monthly carbon fluxes in natural ecosystems
  veg.ingpp[dm] = y[I_INGPP];
  veg.gpp[dm] = y[I_GPP];
  veg.innpp[dm] = y[I_INNPP];
  veg.npp[dm] = y[I_NPP];
  veg.gpr[dm] = y[I_GPR];
  veg.rm[dm] = y[I_RVMNT];
  veg.rg[dm] = y[I_RVGRW];
  veg.ltrfal[dm].carbon = y[I_LTRC];
  microbe.rh[dm] = y[I_RH];
  nep[dm] = y[I_NEP];

  // Monthly nitrogen fluxes in natural ecosystems
  soil.ninput[dm] = y[I_NINP];
  veg.inuptake = y[I_INNUP];
  veg.nuptake[dm] = y[I_VNUP];
  veg.suptake[dm] = y[I_VSUP];
  veg.luptake[dm] = y[I_VLUP];
  veg.nmobil[dm] = y[I_VNMBL];
  veg.nresorb[dm] = y[I_VNRSRB];
  veg.ltrfal[dm].nitrogen = y[I_LTRN];
  microbe.nuptake[dm] = y[I_MNUP];
  microbe.netnmin[dm] = y[I_NMIN];
  soil.nlost[dm] = y[I_NLST];

  // Monthly water fluxes
  atms.rain[dm] = y[I_RAIN];
  soil.rperc[dm] = y[I_RPERC];
  soil.rrun[dm] = y[I_RRUN];
  atms.snowfall[dm] = y[I_SNWFAL];
  soil.snowinf[dm] = y[I_SNWINF];
  soil.sperc[dm] = y[I_SPERC];
  soil.srun[dm] = y[I_SRUN];
  atms.pet[dm] = y[I_PET];
  atms.eet[dm] = y[I_EET];
  soil.h2oyld[dm] = y[I_WYLD];

  // added for soil thermal model

  atms.tsoil[dm]=  y[I_TSOIL];
  atms.dst5[dm]=   y[I_DST5];
  atms.dst10[dm]=  y[I_DST10];
  atms.dst20[dm]=  y[I_DST20];
  atms.dst50[dm]=  y[I_DST50];
  atms.dst100[dm]= y[I_DST100];
  atms.dst200[dm]= y[I_DST200];
  atms.frontd[dm]= y[I_FRONTD];
  atms.thawbe[dm]= y[I_THAWBE];
  atms.thawend[dm]= y[I_THAWEND];

// end of adding


  // Monthly phenology
  veg.unnormleaf[dm] = y[I_UNRMLF];
  veg.leaf[dm] = y[I_LEAF];

  // Monthly carbon and nitrogen fluxes associated with
  //  agricultural conversion
  ag.convrtflx[dm].carbon = y[I_CNVRTC];
  ag.convrtflx[dm].nitrogen = y[I_CNVRTN];
  ag.sconvrtflx[dm].carbon = y[I_SCNVRTC];
  ag.sconvrtflx[dm].nitrogen = y[I_SCNVRTN];
  ag.nvretent[dm] = y[I_NVRTNT];
  ag.nsretent[dm] = y[I_NSRTNT];
  ag.nretent[dm] = y[I_NRETNT];

  // Monthly carbon and nitrogen fluxes from agricultural
  //   ecosystems
  ag.npp[dm].carbon = y[I_AGNPPC];
  ag.npp[dm].nitrogen = y[I_AGNPPN];
  ag.ltrfal[dm].carbon = y[I_AGLTRC];
  ag.ltrfal[dm].nitrogen = y[I_AGLTRN];
  ag.fertn[dm] = y[I_AGFRTN];

  // Monthly integrated carbon fluxes
  ag.totnpp[dm] = y[I_TOTNPP];
  cflux[dm] = y[I_CFLX];

  // Update sum of annual carbon storage
  veg.yrcarbon  += y[I_VEGC];
  soil.yrorgc += y[I_SOLC];

  // Update sum of annual nitrogen storage
  veg.yrnitrogen  += y[I_STRN] + y[I_STON];
  veg.yrstructn += y[I_STRN];
  soil.yrorgn += y[I_SOLN];
  soil.yravln  += y[I_AVLN];
  veg.yrstoren += y[I_STON];

  // Update sum of annual water storage
  soil.yravlh2o += y[I_AVLW];
  soil.yrrgrndh2o += y[I_RGRW];
  soil.yrsnowpack += y[I_SNWPCK];
  soil.yrsgrndh2o += y[I_SGRW];
  soil.yrsmoist += y[I_SM];
  soil.yrpctp += y[I_PCTP];
  soil.yrvsm += y[I_VSM];

  // Update sum of annual carbon fluxes in natural ecosystems
  veg.yringpp += y[I_INGPP];
  veg.yrgpp   += y[I_GPP];
  veg.yrinnpp += y[I_INNPP];
  veg.yrnpp   += y[I_NPP];
  veg.yrltrc  += y[I_LTRC];
  microbe.yrrh    += y[I_RH];
  yrnep   += y[I_NEP];

 // Update sum of annual nitrogen fluxes in natural ecosystems
  soil.yrnin   += y[I_NINP];
  veg.yrinnup += y[I_INNUP];
  veg.yrnup   += y[I_VNUP];
  veg.yrsup    += y[I_VSUP];
  veg.yrlup    += y[I_VLUP];
  veg.yrnmobil += y[I_VNMBL];
  veg.yrnrsorb += y[I_VNRSRB];
  veg.yrltrn  += y[I_LTRN];
  microbe.yrnmin  += y[I_NMIN];
  soil.yrnlost += y[I_NLST];

   // Update sum of annual water fluxes
  atms.yrrain += y[I_RAIN];
  soil.yrrperc += y[I_RPERC];
  soil.yrrrun += y[I_RRUN];
  atms.yrsnowfall += y[I_SNWFAL];
  soil.yrsnowinf += y[I_SNWINF];
  soil.yrsperc += y[I_SPERC];
  soil.yrsrun += y[I_SRUN];
  atms.yrpet += y[I_PET];
  atms.yreet += y[I_EET];
  soil.yrh2oyld += y[I_WYLD];

  // Update sum of annual phenology in natural ecosystems
  veg.yrunleaf += y[I_UNRMLF];
  veg.yrleaf += y[I_LEAF];
  veg.yrfpc += veg.fpc[dm];

 // Update sum of annual carbon and nitrogen fluxes from
 //   agricultural conversion
  ag.yrconvrtC += y[I_CNVRTC];
  ag.yrconvrtN += y[I_CNVRTN];
  ag.yrsconvrtC += y[I_SCNVRTC];
  ag.yrsconvrtN += y[I_SCNVRTN];
  ag.yrslashC += y[I_SLASHC];
  ag.yrslashN += y[I_SLASHN];

 // Update sum of annual carbon and nitrogen fluxes from
 //   agricultural ecosystems
  ag.yrnppC += y[I_AGNPPC];
  ag.yrnppN += y[I_AGNPPN];
  ag.formPROD1.carbon += y[I_AGFPRDC];
  ag.formPROD1.nitrogen += y[I_AGFPRDN];
  ag.yrltrc += y[I_AGLTRC];
  ag.yrltrn += y[I_AGLTRN];
  ag.yrfertn += y[I_AGFRTN];

   // Update sum of annual integrated carbon fluxes
  yrcflux += y[I_CFLX];

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void TTEM::setPrevState(double prevState[],double currentState[])
{
  for (int i = 0; i < NUMEQ; i++) { prevState[i] = currentState[i]; }

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************** */

int TTEM::stepyr(const int& dyr, const int& itype, int& intflag, double& tol)
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

  // Reset annual fluxes to zero
//  cout <<"resetYrFlux()"<<endl;
  resetYrFlux();
//  cout <<"run resetYrFlux()"<<endl;

  // Convert natural vegetation to agriculture

  if (ag.state == 1 && ag.prvstate == 0)
  {

    y[I_VEGC] = ag.vrespar[veg.cmnt] * veg.plant[CYCLE-1].carbon;
    y[I_STRN] = ag.vrespar[veg.cmnt] * veg.strctrl[CYCLE-1].nitrogen;
    y[I_STON] = ag.vrespar[veg.cmnt] * veg.labile[CYCLE-1].nitrogen;

    // Create PROD10 and PROD100 from conversion to agriculture

    ag.conversion(veg.cmnt, veg, soil);

  }
//  cout <<"ag.conversion()"<<endl;

  // Revert to natural vegetation after cropland abandonment

  if (ag.state == 0 && ag.prvstate == 1)
  {
    veg.plant[CYCLE-1].carbon = y[I_VEGC];
    veg.strctrl[CYCLE-1].nitrogen = y[I_STRN];
    veg.labile[CYCLE-1].nitrogen = y[I_STON];
  }

  for (dm = 0; dm < CYCLE; dm++)
  {
    if (initFlag == 1)
    {
      if (atms.tco2flag != 0) { atms.co2[dm] = atms.tco2[dyr][dm]; }
      if (atms.ttairflag != 0) { atms.tair[dm] = atms.ttair[dyr][dm]; }
      if (atms.tprecflag != 0) { atms.prec[dm] = atms.tprec[dyr][dm]; }
      if (atms.tmsoilflag != 0) { atms.msoil[dm] = atms.tmsoil[dyr][dm]; }//XL ad for msoil
      if (ag.tlulcflag != 0)
      {
        if (ag.RAP0flag == 0)
        {
	  ag.potnpp[dm] = ag.tpotnpp[itype][dyr][dm];
        }
        else { ag.potnpp[dm] = 0.0; }
      }
    }

    // Get environmental conditions for month "dm"

     // calculate snowpack for the specific year

//    atms.pet[dm] = atms.petjh(atms.nirr[dm], atms.tair[dm], dm);
//    if (atms.pet[dm] <= 0.0) { atms.pet[dm] = 0.001; }

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
	soil.snowinf[dm] = soil.snowmelt(elev, atms.tair[dm],
						 atms.prevtair, y[I_SNWPCK]);
    soil.snowpack[dm]= atms.snowfall[dm] - soil.snowinf[dm];
    if (soil.snowpack[dm] < 0.0) { soil.snowpack[dm] = 0.0; }

   if (tmflgs.stmflg == 1) {
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
             sthermal.soilm= atms.msoil[dm];//XL add for msoil
            break;
         }

    switch (dm) 
    {
    	
	    case 0: 
	    	{
	    		sthermal.hsnow19 = (soil.snowpack[CYCLE-1]+soil.snowpack[0])
			                         /2.0/1000.0; // satisfy soil thermal model
			    sthermal.hsnow29 = soil.snowpack[0]/1000.0;
			    sthermal.hsnow39 = (soil.snowpack[1]+soil.snowpack[0])/2.0/1000.0;
			  };break;
	    case 11:
	    	{
	    		sthermal.hsnow19 = (soil.snowpack[11]+soil.snowpack[10])/2.0/1000.0;
	 // satisfy soil thermal model
          sthermal.hsnow29 = soil.snowpack[11]/1000.0;
          sthermal.hsnow39 = (soil.snowpack[11]+soil.snowpack[0])/2.0/1000.0;
        }; break;
	      default:
	      	{
	      		sthermal.hsnow19 = (soil.snowpack[dm-1]+soil.snowpack[dm])
			                         /2.0/1000.0; // satisfy soil thermal model
            sthermal.hsnow29 = soil.snowpack[dm]/1000.0;
            sthermal.hsnow39 = (soil.snowpack[dm]+soil.snowpack[dm+1])/2.0/1000.0;
          };break;
        }
         
//   cout << airt19 << " " << airt29 << " " << airt39
//	<< " " << hsnow19 << " " << hsnow29 << " " << hsnow39 <<endl;

	  sthermal.is9 = 10;
    sthermal.ism19 = 9;
    sthermal.smass9 = 0.f;

   int i;

/*	for (i=1;i <=210; i++) {
    sthermal.t9[i] = 0;
  	sthermal.xfa9[i] = -1e10f;
	sthermal.xfb9[i] = 0.f;
	sthermal.water9[i] = 1.f;
	sthermal.dx9[i] = 0.f;
	sthermal.weight9[i] = 0.f; }*/ //Removed by J. Tang, 2007, April 18


   sthermal.calcu_cdsnow = 0.2;
   // constant for equilibrium run --- an assumpition
   sthermal.water2 = soil.pctp[dm]/100.0;
   // assume moss and organic have the same water content, 29/march/2001

/*  if (sthermal.soiltemp_(&sthermal.water2, &sthermal.calcu_cdsnow,
   &sthermal.airt19, &sthermal.airt29, &sthermal.airt39,
	&sthermal.hsnow19, &sthermal.hsnow29,
	&sthermal.hsnow39,sthermal.weight9, &sthermal.smass9,
	 &sthermal.is9, &sthermal.ism19, sthermal.x9, sthermal.dx9,
	 sthermal.xfb9, sthermal.xfa9,	sthermal.water9, sthermal.t9,
	  &sthermal.tsoil,&sthermal.frontd,&sthermal.thawbegin,
    &sthermal.thawend,sthermal.diffsoilt, veg.cmnt) !=0)
    { printf("bad tem"); 
//        getch();
};*/

   //      ++kswitch;
   //added for debugging, bu J. Tang, 2007, 4. 15
//    sthermal.Showsnowecd(veg.cmnt);
//    sthermal.Showsoillecd(veg.cmnt);
//    sthermal.Showsoiltecd(veg.cmnt);
//    sthermal.Inputwrite(veg.cmnt);
    
//    cout <<"soil thermal model stepyr"<<endl;
   
    if(sthermal.Soiltemp_(sthermal.water2, sthermal.calcu_cdsnow,
       sthermal.airt19, sthermal.airt29, sthermal.airt39,
	     sthermal.hsnow19, sthermal.hsnow29,
	     sthermal.hsnow39, sthermal.smass9,
	     sthermal.is9, sthermal.ism19, 
	     sthermal.tsoil, sthermal.frontd, sthermal.thawbegin,
       sthermal.thawend,sthermal.diffsoilt, veg.cmnt, sthermal.soilm) !=0)
    {
    	cout <<"bad tem"<<endl;
    }
//    cout <<"returned in stepyr"<<endl;


      atms.frontd[dm]=sthermal.frontd;
      atms.thawbe[dm]=sthermal.thawbegin;
      atms.thawend[dm]=sthermal.thawend;

      atms.tsoil[dm]=sthermal.tsoil;
      atms.dst5[dm]=sthermal.diffsoilt[0];
      atms.dst10[dm]=sthermal.diffsoilt[1];
      atms.dst20[dm]=sthermal.diffsoilt[2];
      atms.dst50[dm]=sthermal.diffsoilt[3];
      atms.dst100[dm]=sthermal.diffsoilt[4];
      atms.dst200[dm]=sthermal.diffsoilt[5];
      

//     cout <<sthermal.frontd<<" "<<sthermal.thawbegin<<" "<<sthermal.thawend<<" "<<sthermal.tsoil<<" "
//          <<sthermal.diffsoilt[0]<<" "<<sthermal.diffsoilt[1]<<" "<<sthermal.diffsoilt[2]<<" "
//          <<sthermal.diffsoilt[3]<<" "<<sthermal.diffsoilt[4]<<" "<<sthermal.diffsoilt[5]<<endl;
 /*** end of calling soil thermal model ***/
  } // stmflg ===1

 // * calculate the thawing-proportion in early month and later fall

  if (atms.dst5[dm] ==0.0) atms.dst5[dm] =0.001;
 if (atms.dst5[dm-1] ==0.0) atms.dst5[dm-1] =0.002;
 if (atms.dst5[dm+1] ==0.0) atms.dst5[dm+1] =0.003;

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] > 0.0)
  && (atms.dst5[dm+1] > 0.0))
	 veg.thawpercent[dm] = 1- (atms.dst5[dm] / (atms.dst5[dm]
	 - atms.dst5[dm-1]));

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] < 0.0)
 && (atms.dst5[dm+1] < 0.0))
     veg.thawpercent[dm] = 0.0;

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] < 0.0) && (atms.dst5[dm+1] >0.0))
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
	veg.thawpercent[dm] = ((1- (atms.dst5[dm] / (atms.dst5[dm]
	 - atms.dst5[dm-1])))
        + (atms.dst5[dm] / (atms.dst5[dm] - atms.dst5[dm+1])))/2.0;

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm] < 0.0)
  && (atms.dst5[dm+1] > 0.0))
   veg.thawpercent[dm] = ((atms.dst5[dm-1] / (atms.dst5[dm-1] - atms.dst5[dm]))
       + ( 1- (atms.dst5[dm+1] / (atms.dst5[dm+1] - atms.dst5[dm]))))/2.0;

 if (veg.thawpercent[dm] > 1.0) veg.thawpercent[dm] = 1.0;
 if (veg.thawpercent[dm] < 0.0) veg.thawpercent[dm] = 0.0;


/*
 if (atms.dst5[dm] ==0.0) atms.dst5[dm] =0.001;
 if (atms.dst5[dm-1] ==0.0) atms.dst5[dm-1] =0.002;
 if (atms.dst5[dm+1] ==0.0) atms.dst5[dm+1] =0.003;

 if (atms.dst5[dm] > 0.0)
     {
      if (atms.dst5[dm-1] > 0.0)
         {
          if (atms.dst5[dm+1] > 0.0)
            {
             veg.thawpercent[dm] = 1.0; // summer months
			 // To overcome wired JUNE phenomena,
			 Q. Z. Sept. 03 /2001, delete the following dm ==5 sentences
     //        if (dm == 5) {
	 //          if (atms.frontd[5] < 1.5)
	 veg.thawpercent[5] = abs(atms.frontd[5] / 1.5);

//              veg.thawpercent[4] = abs(atms.dst5[3]
 / (atms.dst5[4] - atms.dst5[3])); // may
//              veg.thawpercent[dm] = 1.1 * veg.thawpercent[4];
// Assume Jun thawing has not reached the rooting depth yet
      //        }
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
  // end of adding for thawing-frozen
  */
//    cout << "getenviron()"<<endl;


 // Min added for acclimation
 /*if (totyr >= 2001 && totyr <= 2010)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (12.36-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (12.36-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (12.36-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (12.36-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (12.36-10.9) + veg.topt;
 }
 else if (totyr >= 2011 && totyr <= 2020)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (12.9-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (12.9-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (12.9-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (12.9-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (12.9-10.9) + veg.topt;
 }
 else if (totyr >= 2021 && totyr <= 2030)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (13.61-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (13.61-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (13.61-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (13.61-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (13.61-10.9) + veg.topt;
 }
 else if (totyr >= 2031 && totyr <= 2040)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (14.60-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (14.60-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (14.60-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (14.60-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (14.60-10.9) + veg.topt;
 }
 else if (totyr >= 2041 && totyr <= 2050)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (16-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (16-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (16-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (16-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (16-10.9) + veg.topt;
 }
 else if (totyr >= 2051 && totyr <= 2060)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (17.51-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (17.51-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (17.51-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (17.51-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (17.51-10.9) + veg.topt;
 }
 else if (totyr >= 2061 && totyr <= 2070)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (19.37-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (19.37-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (19.37-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (19.37-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (19.37-10.9) + veg.topt;
 }
 else if (totyr >= 2071 && totyr <= 2080)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (20.91-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (20.91-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (20.91-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (20.91-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (20.91-10.9) + veg.topt;
 }
 else if (totyr >= 2081 && totyr <= 2090)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (22.23-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (22.23-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (22.23-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (22.23-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (22.23-10.9) + veg.topt;
 }
 else if (totyr >= 2091 && totyr <= 2100)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (23.35-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (23.35-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (23.35-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (23.35-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (23.35-10.9) + veg.topt;
 }*/

	getenviron(dm);
//	cout << "run getenviron()"<<endl;
//    cout <<"adpt()"<<endl;
	mintflag = adapt(NUMEQ,y,tol,dm);
//	cout <<"run adapt()"<<endl;
    if (mintflag == 1) { intflag = 1; }

    if (blackhol != 0) { qualcon[dyr][itype] = 10; }

    massbal(y,prevy);
    setPrevState(prevy,y);

    // Update carbon, nitrogen and water pools and fluxes from
    //  integrator results
    setMonth(dm, y);

	resetODEflux(y);
//	cout <<"run resetODEflux()"<<endl;
  }


  //  Update maximum EET, maximum PET, GPP optimum temperature (veg.topt),
  //  and maximum leaf cover (veg.prvleafmx) for the current year
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

  // Update optimum temperature parameters for GPP

  if (veg.topt > veg.toptmax[veg.cmnt]) { veg.topt = veg.toptmax[veg.cmnt]; }
  if (veg.topt < veg.toptmin[veg.cmnt]) { veg.topt = veg.toptmin[veg.cmnt]; }


// Determine vegetation C/N parameter as a function of vegetation type,
// annual PET, annual EET, CO2 concentration

  veg.updateC2N(veg.cmnt, atms.yreet, atms.yrpet, atms.co2[11], atms.initco2);

  soil.yravlh2o /= 12.0;
  soil.yrrgrndh2o /= 12.0;
  soil.yrsnowpack /= 12.0;
  soil.yrsgrndh2o /= 12.0;
  soil.yrsmoist /= 12.0;
  soil.yrpctp /= 12.0;
  soil.yrvsm /= 12.0;

  atms.yrtsoil /= 12.0; // for soil temperature
  atms.yrfrontd /= 12.0; // for soil temperature
  atms.yrthawbegin /= 12.0; // for soil temperature
  atms.yrthawend /= 12.0; // for soil temperature


  veg.yrcarbon  /= 12.0 ;
  veg.yrnitrogen  /= 12.0;
  veg.yrstructn /= 12.0;

  if (veg.yrstructn != 0.0)
  {
    veg.yrc2n  = veg.yrcarbon / veg.yrstructn;
  }
//  else { veg.yrc2n = MISSING; }

  soil.yrorgc /= 12.0;
  soil.yrorgn /= 12.0;

  if (soil.yrorgn != 0.0)
  {
    soil.yrc2n = soil.yrorgc / soil.yrorgn;
  }
//  else { soil.yrc2n = MISSING; }

  soil.yravln  /= 12.0;
  veg.yrstoren /= 12.0;
  veg.yrunleaf /= 12.0;
  veg.yrleaf /= 12.0;
  veg.yrfpc /= 12.0;


  // y[2] changed to y[I_SOLC] by DWK on 20000130 and
  // y[3] changed to y[I_SOLN] by DWK on 20000130
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

  if (endeq > 0) { ++endeq; }
//  cout <<"returned from stepyr"<<endl;

  return endeq;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */
// modified for Soil Thermal model

//int TTEM::transient(const int& dyr, const int& itype, double& tol)
int TTEM::transient(const int& dyr, const int& itype,
			 double& tol, const int& RTIME)
{

  endeq = 0;

  if (atms.tco2flag == 1) { totyr = atms.co2year[dyr]; }
  else if (atms.ttairflag == 1) { totyr = atms.tairyear[dyr]; }
  else if (atms.tprecflag == 1) { totyr = atms.precyear[dyr]; }

  if (atms.ttairflag != 0) { atms.mxtair = atms.mxttair[dyr]; }
  if (atms.tprecflag != 0) { atms.yrprec = atms.yrtprec[dyr]; }

  if (ag.tlulcflag == 1)
  {
    ag.state = ag.tstate[dyr];
    ag.RAP = ag.tRAP[dyr];
  }

 // stepyr(dyr,itype, intflag, tol);
   transtepyr(dyr,itype, intflag, tol, RTIME);


  // Update annual agricultural product pools and fluxes

  ag.updateyr(dyr);

  if (totyr == startyr) { wrtyr = 0;}
  if (totyr > startyr) {++wrtyr; }

  return wrtyr;

};

// addition in order to
  int dm;
  int mintflag; //use RTIME for soil thermal model

int TTEM::transtepyr(const int& dyr, const int& itype,
	 int& intflag, double& tol, const int& RTIME)
{


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


  // Reset annual fluxes to zero

  resetYrFlux();

  // Convert natural vegetation to agriculture

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

       //calculate CDM
     double CDM;
     CDM = 0.0;
     for (dm = 0; dm < CYCLE; dm++)
     CDM = CDM + (10.0 - atms.ttair[dyr][dm]);

  for (dm = 0; dm < CYCLE; dm++)
  {
    if (initFlag == 1)
    {
      if (atms.tco2flag != 0) { atms.co2[dm] = atms.tco2[dyr][dm]; }
      if (atms.ttairflag != 0) { atms.tair[dm] = atms.ttair[dyr][dm]; }
      if (atms.tprecflag != 0) { atms.prec[dm] = atms.tprec[dyr][dm]; }
      if (atms.tmsoilflag != 0) { atms.msoil[dm] = atms.tmsoil[dyr][dm]; }//XL ad for msoil

      if (ag.tlulcflag != 0)
      {
        if (ag.RAP0flag == 0)
        {
	  ag.potnpp[dm] = ag.tpotnpp[itype][dyr][dm];
        }
        else { ag.potnpp[dm] = 0.0; }
      }
    }


    // Get environmental conditions for month "dm"
  // calculate snowpack for the specific year

//  for (dm = 0; dm < CYCLE; dm++) {
 //    atms.pet[dm] = atms.petjh(atms.nirr[dm], atms.tair[dm], dm);
 //   if (atms.pet[dm] <= 0.0) { atms.pet[dm] = 0.001; }

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
	soil.snowinf[dm] = soil.snowmelt(elev, atms.tair[dm],
	 atms.prevtair, y[I_SNWPCK]);

    soil.snowpack[dm]= atms.snowfall[dm] - soil.snowinf[dm];
     if (soil.snowpack[dm] < 0.0) { soil.snowpack[dm] = 0.0; }


  if (tmflgs.stmflg==1) {
// call soil thermal subroutine
    switch (dm) {
     case 0: if (dyr==0) {
             sthermal.airt19= atms.ttair[dyr][0];; // satisfy soil thermal model
             sthermal.airt29= (atms.ttair[dyr][1]+atms.ttair[dyr][0])/2.0;
             sthermal.airt39= (atms.ttair[dyr][1]+atms.ttair[dyr][2])/2.0; }
             else {
			 sthermal.airt19= (atms.ttair[dyr-1][11]+atms.ttair[dyr][0])/2.0;
			  // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][0];
             sthermal.airt39= (atms.ttair[dyr][0]+atms.ttair[dyr][1])/2.0;
              }
             break;
     case 11:if ( dyr<RTIME-3) {
			 sthermal.airt19= (atms.ttair[dyr][CYCLE-3]
			 +atms.ttair[dyr][10])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][11];
             sthermal.airt39= (atms.ttair[dyr][11]+atms.ttair[dyr+1][0])/2.0; }
             else {
			 sthermal.airt19= (atms.ttair[dyr][CYCLE-2]
			 +atms.ttair[dyr][11])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][11];
             sthermal.airt39= (atms.ttair[dyr][11]+atms.ttair[dyr][0])/2.0; }
             break;
	 default: sthermal.airt19= (atms.ttair[dyr][dm-1]
	 	+atms.ttair[dyr][dm])/2.0; // satisfy soil thermal model
             sthermal.airt29= atms.ttair[dyr][dm];
             sthermal.airt39= (atms.ttair[dyr][dm]+atms.ttair[dyr][dm+1])/2.0;
             sthermal.soilm= atms.tmsoil[dyr][dm];//XL add for msoil
             break;
         }

	 // using soil.snowpack, which is calculated from the WBM,
	 // to drive soil thermal model
	 // convert mm to m for snow depth by dividing 1000,
	 //also assume the snow density is changed according to
	 //the snow classification, March,28, 2001, by Q. Z.
	 // According to the vegetation type determine the wind speed,
	 // calculate cooling degree month (CDM) and daily averaged
	 //precipitation to determine the snow classication
	 // snow classification determine the snow density, and
	 //therefore the snow depth, and snow thermal conductivity

     // determine the wind speed level
     int wind_speed; // high = 1; low =0.0

     switch (veg.cmnt) {
         case 1: wind_speed =1; break;
         case 2: wind_speed =1; break;
         case 3: wind_speed =0; break;
         case 4: wind_speed =0; break;
         case 5: wind_speed =0; break;
         case 6: wind_speed =1; break;
         case 7: wind_speed =1; break;
         case 8: wind_speed =1; break;
         case 9: wind_speed =0; break;
         case 10: wind_speed =1; break;
         case 11: wind_speed =0; break;
         case 12: wind_speed =1; break;
        }


     // determine the snow classes
     int snow_class;

     if (CDM < 50.0) snow_class = 6; // class Ephemeral
     else
       {
        if (CDM > 125.0) {
                         if (atms.prec[dm] / atms.daze[dm] >2.0)
                           {
                           if (wind_speed ==1) snow_class = 1;
                           else  snow_class = 2;
                           }
                         else
                            {
                           if (wind_speed ==0) snow_class = 2;
                           else snow_class = 1;
                            }
                         }

        else {
              if (atms.prec[dm] / atms.daze[dm] <2.0)
              {
               if (wind_speed ==1) snow_class = 4;
               else snow_class = 3;
              }
              else
               {
                if (wind_speed ==0) snow_class = 5;
                else snow_class = 5;
               }

            }
        }


// determine the snow denisty and snow thermal conductivity
  double snow_dens; // g mm-3

   switch (snow_class) {
   case 1: snow_dens = 0.28;  break;
   case 2: snow_dens = 0.225;  break;
   case 3: snow_dens = 0.25;  break;
   case 4: snow_dens = 0.25;  break;
   case 5: snow_dens = 0.30;  break;
   case 6: snow_dens = 0.35;  break;
   }

  sthermal.calcu_cdsnow = pow(10., (2.650 * snow_dens - 1.652));

  sthermal.water2 = soil.pctp[dm]/100.0;
  // assume moss and organic have the same water content, 29/march/2001

  // end of the program of snow classification

    switch (dm) {
	 case 0: sthermal.hsnow19=(soil.snowpack[CYCLE-1]
	 +soil.snowpack[0])/2.0/1000.0/snow_dens; // satisfy soil thermal model
             sthermal.hsnow29= soil.snowpack[0]/1000.0/snow_dens;
			 sthermal.hsnow39= (soil.snowpack[1]+soil.snowpack[0])
			 /2.0/1000.0/snow_dens;
             break;
	 case 11: sthermal.hsnow19=(soil.snowpack[11]+soil.snowpack[10])
	 /2.0/1000.0/snow_dens; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[11]/1000.0/snow_dens;
			  sthermal.hsnow39= (soil.snowpack[11]+soil.snowpack[0])
			  /2.0/1000.0/snow_dens;
            break;
	 default: sthermal.hsnow19=(soil.snowpack[dm-1]+soil.snowpack[dm])
	 /2.0/1000.0/snow_dens; // satisfy soil thermal model
              sthermal.hsnow29= soil.snowpack[dm]/1000.0/snow_dens;
			  sthermal.hsnow39= (soil.snowpack[dm]+soil.snowpack[dm+1])
			  /2.0/1000.0/snow_dens;
            break;
               }

    int i;

    sthermal.is9 = 10;
    sthermal.ism19 = 9;
    sthermal.smass9 = 0.f;

/*	for (i=1;i <=210; i++) {
    sthermal.t9[i] = 0;
  	sthermal.xfa9[i] = -1e10f;
	sthermal.xfb9[i] = 0.f;
	sthermal.water9[i] = 1.f;
	sthermal.dx9[i] = 0.f;
	sthermal.weight9[i] = 0.f; }*///Removed by J. Tang, 2007, April 18

// Need to change the soil moisture (water content for organic layer),
// snow thermal conductivity based on the snow classification ??


     //added for debugging, bu J. Tang, 2007, 4. 15
//    sthermal.Showsnowecd(veg.cmnt);
//    sthermal.Showsoillecd(veg.cmnt);
//    sthermal.Showsoiltecd(veg.cmnt);
//    sthermal.Inputwrite(veg.cmnt);
//        cout <<"Soil thermal model in transtepyr"<<endl;
     
	 if (sthermal.Soiltemp_(sthermal.water2, sthermal.calcu_cdsnow,
	     sthermal.airt19,  sthermal.airt29,  sthermal.airt39,
	     sthermal.hsnow19, sthermal.hsnow29, sthermal.hsnow39,
	     sthermal.smass9,  sthermal.is9, sthermal.ism19,sthermal.tsoil,         
		 sthermal.frontd, sthermal.thawbegin, sthermal.thawend,
		 sthermal.diffsoilt, veg.cmnt, sthermal.soilm) !=0)
      { printf("bad tem");
  //      getch();
        };

//        cout <<"returned in transtepyr"<<endl;

    //  ++ integer (kswitch);
      atms.frontd[dm]=sthermal.frontd;
      atms.thawbe[dm]=sthermal.thawbegin;
      atms.thawend[dm]=sthermal.thawend;

      atms.tsoil[dm]=sthermal.tsoil;
      atms.dst5[dm]=sthermal.diffsoilt[0];
      atms.dst10[dm]=sthermal.diffsoilt[1];
      atms.dst20[dm]=sthermal.diffsoilt[2];
      atms.dst50[dm]=sthermal.diffsoilt[3];
      atms.dst100[dm]=sthermal.diffsoilt[4];
      atms.dst200[dm]=sthermal.diffsoilt[5];

//     cout <<sthermal.frontd<<" "<<sthermal.thawbegin<<" "<<sthermal.thawend<<" "<<sthermal.tsoil<<" "
//          <<sthermal.diffsoilt[0]<<" "<<sthermal.diffsoilt[1]<<" "<<sthermal.diffsoilt[2]<<" "
//          <<sthermal.diffsoilt[3]<<" "<<sthermal.diffsoilt[4]<<" "<<sthermal.diffsoilt[5]<<endl;

     } //stmflg ==1
   /*** end of calling soil thermal model ***/

   // * calculate the thawing-proportion in early month and
   // later fall, 28/Nov/2001

 if (atms.dst5[dm] ==0.0) atms.dst5[dm] =0.001;
 if (atms.dst5[dm-1] ==0.0) atms.dst5[dm-1] =0.002;
 if (atms.dst5[dm+1] ==0.0) atms.dst5[dm+1] =0.003;

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] > 0.0)
 && (atms.dst5[dm+1] > 0.0))
	 veg.thawpercent[dm] = 1- (atms.dst5[dm] / (atms.dst5[dm]
	  - atms.dst5[dm-1]));

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] < 0.0)
  && (atms.dst5[dm+1] < 0.0))
     veg.thawpercent[dm] = 0.0;

 if ((atms.dst5[dm-1] < 0.0) && (atms.dst5[dm] < 0.0)
  && (atms.dst5[dm+1] > 0.0))
	 veg.thawpercent[dm] = 1- (atms.dst5[dm+1] / (atms.dst5[dm+1]
	  - atms.dst5[dm]));

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

 if ((atms.dst5[dm-1] > 0.0) && (atms.dst5[dm]
 < 0.0) && (atms.dst5[dm+1] > 0.0))
   veg.thawpercent[dm] = ((atms.dst5[dm-1]
    / (atms.dst5[dm-1] - atms.dst5[dm]))
       + ( 1- (atms.dst5[dm+1] / (atms.dst5[dm+1] - atms.dst5[dm]))))/2.0;

 if (veg.thawpercent[dm] > 1.0) veg.thawpercent[dm] = 1.0;
 if (veg.thawpercent[dm] < 0.0) veg.thawpercent[dm] = 0.0;

// mark following; added above new algorithm

/*
 if (atms.dst5[dm] ==0.0) atms.dst5[dm] =0.001;
 if (atms.dst5[dm-1] ==0.0) atms.dst5[dm-1] =0.002;
 if (atms.dst5[dm+1] ==0.0) atms.dst5[dm+1] =0.003;

 if (atms.dst5[dm] > 0.0)
     {
      if (atms.dst5[dm-1] > 0.0)
         {
          if (atms.dst5[dm+1] > 0.0)
            {
             veg.thawpercent[dm] = 1.0; // summer months
			 // To overcome wired JUNE phenomena, Q. Z. Sept. 03 /2001,
			 delete the following dm ==5 sentences
     //        if (dm == 5) {
	 //          if (atms.frontd[5] < 1.5) veg.thawpercent[5] =
	  abs(atms.frontd[5] / 1.5);

//              veg.thawpercent[4] = abs(atms.dst5[3]
 / (atms.dst5[4] - atms.dst5[3])); // may
//              veg.thawpercent[dm] = 1.1 * veg.thawpercent[4];
 // Assume Jun thawing has not reached the rooting depth yet
      //        }
            }
          else //atms.dst5[dm+1] < 0.0)
            {
		 //   veg.thawpercent[dm] = atms.dst5[dm]/(atms.dst5[dm]
		 - atms.dst5[dm+1]); // sept.
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

 // Min added for acclimation
 /*if (totyr >= 2001 && totyr <= 2010)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (12.36-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (12.36-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (12.36-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (12.36-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (12.36-10.9) + veg.topt;
 }
 else if (totyr >= 2011 && totyr <= 2020)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (12.9-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (12.9-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (12.9-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (12.9-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (12.9-10.9) + veg.topt;
 }
 else if (totyr >= 2021 && totyr <= 2030)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (13.61-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (13.61-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (13.61-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (13.61-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (13.61-10.9) + veg.topt;
 }
 else if (totyr >= 2031 && totyr <= 2040)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (14.60-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (14.60-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (14.60-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (14.60-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (14.60-10.9) + veg.topt;
 }
 else if (totyr >= 2041 && totyr <= 2050)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (16-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (16-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (16-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (16-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (16-10.9) + veg.topt;
 }
 else if (totyr >= 2051 && totyr <= 2060)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (17.51-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (17.51-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (17.51-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (17.51-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (17.51-10.9) + veg.topt;
 }
 else if (totyr >= 2061 && totyr <= 2070)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (19.37-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (19.37-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (19.37-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (19.37-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (19.37-10.9) + veg.topt;
 }
 else if (totyr >= 2071 && totyr <= 2080)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (20.91-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (20.91-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (20.91-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (20.91-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (20.91-10.9) + veg.topt;
 }
 else if (totyr >= 2081 && totyr <= 2090)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (22.23-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (22.23-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (22.23-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (22.23-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (22.23-10.9) + veg.topt;
 }
 else if (totyr >= 2091 && totyr <= 2100)
 {
	 veg.ac_tmax[veg.cmnt] = 0.38 * (23.35-10.9) + veg.tmax[veg.cmnt];
	 veg.ac_tmin[veg.cmnt] = 0.38 * (23.35-10.9) + veg.tmin[veg.cmnt];
	 veg.ac_toptmax[veg.cmnt] = 0.38 * (23.35-10.9) + veg.toptmax[veg.cmnt];
	 veg.ac_toptmin[veg.cmnt] = 0.38 * (23.35-10.9) + veg.toptmin[veg.cmnt];
	 veg.ac_topt = 0.38 * (23.35-10.9) + veg.topt;
 } */
    getenviron(dm);

    mintflag = adapt(NUMEQ,y,tol,dm);


    if (mintflag == 1) { intflag = 1; }

    if (blackhol != 0) { qualcon[dyr][itype] = 10; }

    massbal(y,prevy);

    setPrevState(prevy,y);

    // Update carbon, nitrogen and water pools and fluxes from
    //  integrator results

    setMonth(dm, y);
    resetODEflux(y);

  }

  //  Update maximum EET, maximum PET, GPP optimum temperature (veg.topt),
  //  and maximum leaf cover (veg.prvleafmx) for the current year
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

  // Update optimum temperature parameters for GPP

  if (veg.topt > veg.toptmax[veg.cmnt]) { veg.topt = veg.toptmax[veg.cmnt]; }
  if (veg.topt < veg.toptmin[veg.cmnt]) { veg.topt = veg.toptmin[veg.cmnt]; }


// Determine vegetation C/N parameter as a function of vegetation type,
// annual PET, annual EET, CO2 concentration

  veg.updateC2N(veg.cmnt, atms.yreet, atms.yrpet, atms.co2[11], atms.initco2);

  soil.yravlh2o /= 12.0;
  soil.yrrgrndh2o /= 12.0;
  soil.yrsnowpack /= 12.0;
  soil.yrsgrndh2o /= 12.0;
  soil.yrsmoist /= 12.0;
  soil.yrpctp /= 12.0;
  soil.yrvsm /= 12.0;

  atms.yrtsoil /= 12.0; // for soil temperature
  atms.yrfrontd /= 12.0; // for soil temperature
  atms.yrthawbegin /= 12.0; // for soil temperature
  atms.yrthawend /= 12.0; // for soil temperature

  veg.yrcarbon  /= 12.0;
  veg.yrnitrogen  /= 12.0;
  veg.yrstructn /= 12.0;

  if (veg.yrstructn != 0.0)
  {
    veg.yrc2n  = veg.yrcarbon / veg.yrstructn;
  }
//  else { veg.yrc2n = MISSING; }

  soil.yrorgc /= 12.0;
  soil.yrorgn /= 12.0;

  if (soil.yrorgn != 0.0)
  {
    soil.yrc2n = soil.yrorgc / soil.yrorgn;
  }
//  else { soil.yrc2n = MISSING; }

  soil.yravln  /= 12.0;
  veg.yrstoren /= 12.0;
  veg.yrunleaf /= 12.0;
  veg.yrleaf /= 12.0;
  veg.yrfpc /= 12.0;


  // y[2] changed to y[I_SOLC] by DWK on 20000130 and
  // y[3] changed to y[I_SOLN] by DWK on 20000130
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

  if (endeq > 0) { ++endeq; }
//  cout <<"returned from transtepyr"<<endl;
  return endeq;

};


