//--------------------------------
// qsoiltemp1.hpp
// headfile for Qsoiltemp1.cpp
//--------------------------------
#if !defined(QSOILTEMP1_H)
 #define QSOILTEMP1_H
 
   typedef long int int_ln;

   const int STM_DEBUG = 0;  //set stm_debug to 1 when debuging, otherwise set to zero

      
   #define ishafday  1
   
   #ifdef ishafday
      const int_ln NMAX0  = 3000;
      const int_ln NNMAX0 = 3000;
   #else
      const int_ln NMAX0  = 15000;
      const int_ln NNMAX0 = 150000;
   #endif
   

 class Qsoiltemp1{
   public:
	Qsoiltemp1();
	~Qsoiltemp1();

/*	int Soiltemp_(double& water2, double& calcu_cdsnow, double& airt199, double& airt299,
				  double& airt399, double& hsnow199, double& hsnow299, double& hsnow399,
				  double weight99[], double& smass9, int_ln& is9, int_ln& ism19, double x99[],
				  double dx99[], double xfb99[], double xfa99[], double water99[], double t99[], double& tsoil9,
				  double& frontd9, double& thawbegin9, double& thawend9,double diffsoilt9[10],
				  const int& cmnt); //, // this function will be used to communicate data with other TEM modules*/
	int Soiltemp_(double& water2, double& calcu_cdsnow9, double& airt199, double& airt299,
				  double& airt399, double& hsnow199, double& hsnow299, double& hsnow399,
				  double& smass99, int_ln& is99, int_ln& ism199,double& tsoil9,				  
				  double& frontd9, double& thawbegin9, double& thawend9,double diffsoilt9[10],
				  const int& cmnt, double& soilm);		// XL add for soilmoisture		  

    int Getsnowecd(ofstream& rflog1);
    int Getsnowecd(const char* ecd);
    int Getsoillecd(ofstream& rflog1);
    int Getsoillecd(const char* ecd);
    int Getsoiltecd(ofstream& rflog1);
    int Getsoiltecd(const char* ecd);
    
    void Showsoillecd(const int& cmnt);
    void Showsoiltecd(const int& cmnt);
    void Showsnowecd(const int& cmnt);
    void Inputwrite(const int& cmnt);
   
   //public variables used to transfer data with other classes/procedures
    int_ln ism19, is9;
    double airt19, airt29, airt39, tsoil, soilm, frontd, thawbegin, thawend, t9[211], x9[211]; // XL add for soilmoisture
    double water9[211], smass9, zz, hsnow19, hsnow29, hsnow39, dx9[211];
    double weight9[211], xfa9[211], xfb9[211], diffsoilt[10];
    
    double calcu_cdsnow; 
    double water1, water2, water3;
    double soilwater1, soilwater2, soilwater3; //blkwater_ struct in Zhuang's STM_TEM
    double hsnow1, hsnow2, hsnow3;
    double airt1, airt2, airt3;
    
                     
   private:
    
	int Lociso(double tiso[], int_ln liso, int_ln kiso);
        int Inthet(int kint[2], double source[]);
	int Temdep(double deptem[], int_ln lmax, int_ln ktemp, double diffsoilt[]);
	int Flxdep(double depflx[], int_ln& ldepf, int_ln& kflux);
	int Pram(int_ln& mat, double& ta, double& tb, double& dense, double& wat, double& sph, double& cnd);
	
	int Grid(int_ln& maxm, const int& cmnt, double& soilm); // XL add for soilmoisture	
	int Bctop(int_ln& maxm, double taer[], int& ktop, double& sflx1, double& sflx2, double& sflx3);
	int Bcbot(int_ln& maxm, double tber[]);
	int Data(double fitted[], const int index = 3, const int_ln np = -99, const double vtime = 1.0);
	
	int Tinitl(const int& cmnt);	
    int Snofal(double snow[], double sden[], double& comp, double& eta0, double& denfac, double& fcmelt, const int& cmnt);
	int Neige(double& acum, double& snoden, double weight[], double& comp, double& eta0, double& denfac, double& fcmelt, double& tdrive);
	int Schnee(double& acum, double& snoden, int_ln& ksnow);
	int Asmblo(int_ln& i1);
	int Asmabv(int_ln& i1);
	int Tabv(int_ln& i1);
	int Asmone(int_ln& i);
	int Asmtwo(int_ln& itzero, int_ln& ibzero);
	int Phase(int_ln& i, double& rlw, double cu, double rsu, double cd,
						double rsd, double& tt, double& td, double& dxx,
						double& gold, double& g, double& t1, double& t2);
	
	int Sethet(const int_ln& maxm, int kint[2], double heat[], double source[], const int& cmnt);
	int Intrpl(const double X[], const double Y[], const int_ln& NMAX,
			   const int& np, double xx[], double yy[],
			   const int_ln& NNMAX);
	int Cubsmh(const double X[], const double Y[],
			   double W[], const int_ln& M, const int& IWT,
			   const double& sval, const double& prec, const int& lim, double work[],
			   double coeff[], const int& idimsn, double& scalc);
	int Cubfit(const double X[], const double Y[], const int_ln& M, const int_ln& Konstl,
			   const double& constl, const int_ln& Konstr, const double& constr, double coeff[],
			   const int& idimsn);
	 
	int Cubval(const int& nprime, const double X[], const int_ln& M, double coeff[], const int& idimsn,
			   double xvalue[], double yvalue[], const int_ln& numpnt);


	int Surf(const int& ktop, const double& sflx1, const double& sflx2, const double& sflx3, double& tdrive);
	int Tblo(const int_ln& i1);

	int sign_t(const double& num);
	double sign_d(const double& a, const double& b);
	void msg_debug(const char* pstr);
    
                                  
     
    
     
    
    
     
	//private variables used to transfer data among the variables
	double dt, dtday, per, first, time, hlat, tf, sigma, emiss, htop, hbot;
	double tair, sflux, tbot, gflux, smass, dtfaz, final, theta, theta1;
	double spare[6], x[210], dx[210], xfa[210], xfb[210], xhet[210];
	double water[210], ddry[210], condt[30], condf[30], spht[30], sphf[30];
	double conx[210], capx[210], e[210], s[210], t[210], ht[210];
	double htold[210], caltim[NNMAX0];
	double eltim[NNMAX0];
    int_ln ism1, is, igm1, ig, imax1, imax, mmax, mater[210], nan;
    int_ln nst, nanmax, kallyr, knodes, kiso, ktemp, ksnow, kenvl, kflux; //read from input files
    int_ln liso, lmax; 
    double vdep;
    
    double blk_cdsnow;
    int ktop; 
    
//	char line_tmp[400];
//	char fmt0[200];
//	char fmt1[200];
//	char fmt[200];

	double my_first, my_final, my_per, my_dtday, my_theta, my_hlat, my_tf, my_gflux;
	int_ln my_ig, my_nst;
	int my_kiso, my_kallyr, my_knodes, my_ktemp, my_ksnow, my_kenvl, my_kflux;
	double my_cdsnow;
	double my_tiso[10];
	int_ln my_liso, my_lmax, my_nanmax, my_ldepf;
	double my_vdep, my_vdep1;
	double my_deptem[20], my_depflx[20];
//  for grid
	double my_thick[20], my_dxa[20], my_dxb[20], my_dense[20], my_wat[20];
	double my_vcond[20], my_vsph[20];
	int my_mat[20];
	double my_condt[30][20], my_spht[30][20];
	double my_condf[30][20], my_sphf[30][20];
	double my_vdens;
	double my_vspace, my_top;
	double my_epsmin;

//  for sethet
	int_ln my_kinti;
	double my_vdep2, my_dephet;
	
// for initl
  int_ln my_index, my_np;
  double my_vdepth;
  double my_tstart[211], my_xstart[211];	

// for snofal
	int_ln my_index2;
	double my_epssno, my_convrt, my_eta0, my_denfac, my_fcmelt, my_denmax;

//  read in parameters
	int_ln MAX[NUMVEG], NST[NUMVEG], KALLYR[NUMVEG], KNODES[NUMVEG], KISO[NUMVEG];
	int_ln KTEMP[NUMVEG], KSNOW[NUMVEG], KENVL[NUMVEG], KFLUX[NUMVEG], LISO[NUMVEG];
	double TISO[NUMVEG], VDEPTH[NUMVEG], DEPTEM1[NUMVEG], DEPTEM2[NUMVEG], DEPTEM3[NUMVEG];
	double DEPTEM4[NUMVEG], DEPTEM5[NUMVEG],VDEP[NUMVEG];
	int_ln LMAX[NUMVEG], NDEPF[NUMVEG], IG[NUMVEG];
	double VDEPTH1[NUMVEG], DEPFLX1[NUMVEG], DEPFLX2[NUMVEG], DEPFLX3[NUMVEG], DEPFLX4[NUMVEG], DEPFLX5[NUMVEG];
	double HLAT[NUMVEG], TF[NUMVEG], GFLUX[NUMVEG],CDSNOW[NUMVEG], FIRST[NUMVEG], FINAL[NUMVEG];
	double PER[NUMVEG], DTDAY[NUMVEG], THETA[NUMVEG], TOP[NUMVEG], EPSMIN[NUMVEG];
	double VSPACE[NUMVEG], VDEN[NUMVEG], VDEP1[NUMVEG], DEPHET[NUMVEG], EPSSNO[NUMVEG];
	double CONVRT[NUMVEG], ETAO[NUMVEG], DENFAC[NUMVEG], FCMELT[NUMVEG], DENMAX[NUMVEG];
	int    KINT[NUMVEG], SNOFAL[NUMVEG]; // for index

	int_ln MAT1[NUMVEG], MAT2[NUMVEG], MAT3[NUMVEG], MAT4[NUMVEG], MAT5[NUMVEG], MAT6[NUMVEG];
	double THICK1[NUMVEG], DXA1[NUMVEG], DXB1[NUMVEG], DENSE1[NUMVEG],
		WATER1[NUMVEG], VCOND1[NUMVEG], VSPH1[NUMVEG], COND1[NUMVEG], SPHT1[NUMVEG],
		CONDF1[NUMVEG], SPHF1[NUMVEG];
	double THICK2[NUMVEG], DXA2[NUMVEG], DXB2[NUMVEG], DENSE2[NUMVEG],
		WATER2[NUMVEG], VCOND2[NUMVEG], VSPH2[NUMVEG], COND2[NUMVEG], SPHT2[NUMVEG],
		CONDF2[NUMVEG], SPHF2[NUMVEG];
	double THICK3[NUMVEG], DXA3[NUMVEG], DXB3[NUMVEG], DENSE3[NUMVEG],
		WATER3[NUMVEG], VCOND3[NUMVEG], VSPH3[NUMVEG], COND3[NUMVEG], SPHT3[NUMVEG],
		CONDF3[NUMVEG], SPHF3[NUMVEG];
	double THICK4[NUMVEG], DXA4[NUMVEG], DXB4[NUMVEG], DENSE4[NUMVEG],
		WATER4[NUMVEG], VCOND4[NUMVEG], VSPH4[NUMVEG], COND4[NUMVEG], SPHT4[NUMVEG],
		CONDF4[NUMVEG], SPHF4[NUMVEG];
	double THICK5[NUMVEG], DXA5[NUMVEG], DXB5[NUMVEG], DENSE5[NUMVEG],
		WATER5[NUMVEG], VCOND5[NUMVEG], VSPH5[NUMVEG], COND5[NUMVEG], SPHT5[NUMVEG],
		CONDF5[NUMVEG], SPHF5[NUMVEG];
	double THICK6[NUMVEG], DXA6[NUMVEG], DXB6[NUMVEG], DENSE6[NUMVEG], WATER6[NUMVEG];

	double INDEX[NUMVEG], VDEPP[NUMVEG];
	int_ln NP[NUMVEG];

	double DEPTH[NUMVEG][25], TEMP[NUMVEG][25];
	
	double tsoil20dep[5];//depth to define the top 20cm averaged temperature
	double tsoil20[5];

    
//    ifstream snowdepth, initfile;
    
    
    
 };
 
#endif 
