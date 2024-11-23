/* **************************************************************
*****************************************************************
TMCRB423.HPP - object describing characteristics of soil microbes
	       used in the Terrestrial Ecosystem Model (TEM)
*****************************************************************
************************************************************** */

// Tmicrobe4 uses the global constants CYCLE, NUMMSAC and NUMVEG

#if !defined(TMCRB423_H)
#define TMCRB423_H

class Tmicrobe4 {

   public:

/* **************************************************************
		 Public Functions
************************************************************** */

      void   getvegecd(ofstream& rflog1);
      void   getvegecd(char ecd[80]);
      double nminxclm(const int& dcmnt, const int& dm, double& soilh2o,
                      double& soilorgc, double& soilorgn, double& availn,
                      double& decay, double& rh, double& ksoil);
      double rhxclm(double& soilorgc, double& dq10, double& moist);
      void   showecd(const int& dcmnt);
      // ez changed to dcmnt in yrkd() by DWK on 20000130
      double yrkd(int nfeed, double& yrltrc, double& yrltrn, int dcmnt);


/* **************************************************************
		 Public Variables
************************************************************** */

      // monthly heterotrophic respiration (g C / (sq. meter * month))

      double rh[CYCLE];         

      // monthly nitrogen uptake by microbes (g N / (sq. meter * month))

      double nuptake[CYCLE];
      
      // monthly net nitrogen mineralization (g N / (sq. meter * month))

      double netnmin[CYCLE];    

      // annual heterotrophic respiration (g C / (sq. meter * year))

      double yrrh;

      // annual net nitrogen mineralization (g N / (sq. meter * year))

      double yrnmin;            


/* **************************************************************
		 Public Parameters
************************************************************** */

      // Parameter representing the quality of soil organic matter

      double decay;

      // Biome-specific decomposition parameters for function rhxclm

      double kd;
      double kda[MAXCMNT];
      double kdb[MAXCMNT];
      double kdc;

      double kdin[NUMMSAC];     // kd values read in from file
      double kdsave[NUMMSAC];   // kd values saved to a file

      double lcclnc[MAXCMNT];
      double propftos[MAXCMNT];

      // Biome-specific microbial nitrogen uptake parameters for function
      //   nminxclm

      double nup;
      double nupa[MAXCMNT];
      double nupb[MAXCMNT];

      // Biome-specific N fixation parameter

      double nfixpar[MAXCMNT];

      double cnsoil[MAXCMNT];


      double rhq10[MAXCMNT];

      //Biome-specific parameters describing the influence of soil moisture
      // on decomposition (i.e., moist)

      double moistmin[MAXCMNT];
      double moistopt[MAXCMNT];
      double moistmax[MAXCMNT];

      // Biome-specific half saturation parameter for function nminxclm
      //   describing the effect of available nitrogen on microbial
      //   nitrogen uptake

      double kn2[MAXCMNT];

};

#endif

