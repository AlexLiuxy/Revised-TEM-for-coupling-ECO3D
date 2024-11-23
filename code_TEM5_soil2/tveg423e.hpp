/* **************************************************************
*****************************************************************
TVEG423.HPP -  Vegetation characteristics used in the Terrestrial
               Ecosystem Model (TEM)
            -  Bug fix added by DWK on 19991028
            -  Additional modifications by DWK on 20000102
20000614 - DWK changes innpp to innpp[CYCLE] and ingpp to ingpp[CYCLE]
20010418 - Q. Z. changed for soil thermal model



*****************************************************************
************************************************************** */

#if !defined(BIOMS423_H)
  #include "bioms423.hpp"   // Tveg42 uses Biomass Class
#endif

#if !defined(TBIOME423_H)
  #include "tbiome423.cpp"   // Tveg42 inherits Biome Class
#endif

// Tveg42 also uses the global constants CYCLE and NUMVEG

#if !defined(TVEG423_H)
#define TVEG423_H

class Tveg42 : public Biome
{

  public:

     Tveg42();

/* **************************************************************
		 Public Functions
************************************************************** */

     double deltaleaf(const int& dcmnt, double& eet, double& prveetmx, double& prvleaf);
     void   getecd(ofstream& rflog1);
     void   getecd(char ecd[80]);
     void   getleafecd(ofstream& rflog1);
     void   getleafecd(char ecd[80]);
//     void   getvtype(ofstream& rflog1);
//     void   getvtype(char ecd[80]);
   // modified for thawing STM, hydrology model, adding thawpercent
     double gppxclm(int& dcmnt, double& co2, double& par, double& temp, double& gv, double& leaf, double& foliage,double& thawpercent);

//     double gppxclm(int& dcmnt, double& co2, double& par, double& temp, double& gv, double& leaf, double& foliage);
     void   leafinit(ofstream& rflog1);
     double nupxclm(int& dcmnt, double& soilh2o, double& availn, double& respq10, double& ksoil, double& foliage);
     double rmxclm(int& dcmnt, double& vegc, double& respq10);
     void   showecd(int& dcmnt);
     void   showleaf(int& dcmnt);
     void   updateC2N(const int& dcmnt, double& yreet, double& yrpet, double& currentco2, double& initco2);



/* **************************************************************
		 Public Variables
************************************************************** */

     Biomass plant[CYCLE];     // whole plant biomass (structural + labile)
     Biomass strctrl[CYCLE];   // structural plant biomass
     Biomass labile[CYCLE];    // labile plant biomass
     Biomass ltrfal[CYCLE];    // litterfall

     double nmobil[CYCLE];     // N mobilization by plants
     double nresorb[CYCLE];    // N resorption by plants
     double inuptake;          // initial N uptake by plants
     double nuptake[CYCLE];    // N uptake by plants
     double suptake[CYCLE];    // N uptake by plants for structural N
     double luptake[CYCLE];    // N uptake by plants for labile N
     double inprodcn;          // initial C/N of biomass production
     
     double ingpp[CYCLE];      // initial gross primary productivity
     double gpp[CYCLE];        // gross primary productivity (GPP)

     double innpp[CYCLE];      // initial net primary productivity
     double npp[CYCLE];        // net primary productivity (NPP)

     double rm[CYCLE];         // maintenance respiration
     double rg[CYCLE];         // growth respiration
     double gpr[CYCLE];        // gross plant respiration (rm + rg)

     double yrcarbon;          // annual sum of plant.carbon
     double yrnitrogen;        // annual sum of plant.nitrogen
     double yrstructn;         // annual sum of strctrl.nitrogen
     double yrstoren;          // annual sum of labile.nitrogen
     double yrc2n;             // ratio of yrcarbon to yrnitrogen
     double yrnpp;             // annual sum of npp
     double yrinnpp;           // annual sum of innpp
     double yrnup;             // annual sum of nuptake
     double yrnmobil;          // annual sum of nmobil
     double yrsup;             // annual sum of suptake
     double yrlup;             // annual sum of luptake
     double yrinnup;           // annual sum of innup
     double yrnrsorb;          // annual sum of nresorb
     double yrgpp;             // annual sum of gpp
     double yringpp;           // annual sum of ingpp
     double yrltrc;            // annual sum of ltrfal.carbon
     double yrltrn;            // annual sum of ltrfal.nitrogen
     double yrinpr;
     double yrprod;
     double yrunleaf;          // mean annual unnormalized leaf phenology
     double yrleaf;             // mean annual normalized leaf phenology

     double unnormleaf[CYCLE]; // monthly unnormalized leaf phenology
     double leaf[CYCLE];       // monthly normalized leaf phenology

     double lai[CYCLE];        // monthly leaf area index
     double fpc[CYCLE];        // monthly foliar projective cover
     double yrfpc;
     double alleaf;
     double foliage;

     double thawpercent[CYCLE]; // added for thaw-frown

    // Number of annual iterations for determining monthly phenology

     int leafyrs;

/* **************************************************************
		 Public Parameters
************************************************************** */

     // Biome-specific vegetation C/N parameters

     double cneven[MAXCMNT];
     double cnmin[MAXCMNT];
     double c2n;
     double c2na[MAXCMNT];
     double c2nb[MAXCMNT];
     double c2nmin[MAXCMNT];
     double initcneven[MAXCMNT];
     double dc2n;
     double adjc2n;

     // Biome-specific phenology parameters

     double minleaf[MAXCMNT];
     double aleaf[MAXCMNT];
     double bleaf[MAXCMNT];
     double cleaf[MAXCMNT];
     double initleafmx[MAXCMNT];  // Added by DWK on 19991028
     double prvleafmx[MAXCMNT];
     double unleaf12[MAXCMNT];

     // Biome-specific foliage projection cover parameters

     double fpcmax[MAXCMNT];
     double sla[MAXCMNT];
     double cov[MAXCMNT];

     // Biome-specific allocation parameters

     double kleafc[MAXCMNT];
     double leafmxc[MAXCMNT];

     // Biome-specific carbon uptake parameters for function gppxclm

     double cmax;
     double cmaxcut[MAXCMNT];
     double cmax1a[MAXCMNT];
     double cmax1b[MAXCMNT];
     double cmax2a[MAXCMNT];
     double cmax2b[MAXCMNT];

     // Biome-specific half saturation parameter for function gppxclm
     //   describing the effects of solar atmospheric carbon dioxide
     //   concentration on GPP

     double kc[MAXCMNT];

     // Biome-specific half saturation parameter for function gppxclm
     //   describing the effects of photosybtheically active radiation
     //   on GPP

     double ki[MAXCMNT];

    // Element-specific optimum temperature for GPP

     double topt;

     double tmin[MAXCMNT];
     double toptmin[MAXCMNT];
     double toptmax[MAXCMNT];
     double tmax[MAXCMNT];

	 // Min added for acclimation
	 double ac_tmin[MAXCMNT];
	 double ac_toptmin[MAXCMNT];
	 double ac_toptmax[MAXCMNT];
	 double ac_tmax[MAXCMNT];
	 double ac_topt;

     // Biome-specific parameter to describe the sensitivity of GPP
     //   to evapotranspiration

     double gva[MAXCMNT];

     // Biome-specific respiration parameters for function rmxclm

     double kr;
     double kra[MAXCMNT];
     double krb[MAXCMNT];

     // Biome-specific parameters for function rq10 to describe the
     //   effect of temperature on plant respiration

     double raq10a0[MAXCMNT];
     double raq10a1[MAXCMNT];
     double raq10a2[MAXCMNT];
     double raq10a3[MAXCMNT];

     // Biome-specific parameters to describe the effect of volumetric
     // soil mositure on GPP

     double vsmmin;
     double vsmmina[MAXCMNT];
     double vsmminb[MAXCMNT];

     // Biome-specific nitrogen uptake parameters for function nupxclm

     double nmax;
     double nmaxcut[MAXCMNT];
     double nmax1a[MAXCMNT];
     double nmax1b[MAXCMNT];
     double nmax2a[MAXCMNT];
     double nmax2b[MAXCMNT];

     double kn1[MAXCMNT];

     // Biome-specific proportion of vegetation lost as litterfall

     double cfall[MAXCMNT];  // proportion of vegetation carbon
     double nfall[MAXCMNT];  // proportion of vegetation nitrogen


     double labncon[MAXCMNT];

  private:

/* **************************************************************
		 Private Functions
************************************************************** */

     double rq10(int& dcmnt, double& tair);

};

#endif

