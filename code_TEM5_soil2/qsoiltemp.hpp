#if !defined(QSOILTEMP_H)
  #define QSOILTEMP_H  // soil temperature class


typedef long int integer;
typedef char *address;
typedef short int shortint;
//typedef float double;
typedef double doubledouble;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;

#define TRUE_ (1)
#define FALSE_ (0)

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doubledouble)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doubledouble)min(a,b)
#define dmax(a,b) (doubledouble)max(a,b)

static integer c__4 = 4;
static integer c__1 = 1;
static double c_b127 = 1e-4f;
static integer c__25 = 25;
static doubledouble c_b358 = 10.;
static doubledouble c_b365 = .3333333;
static integer c__30 = 30;
static double c_b644 = 0.f;

class Soilthermal {

  public:

//  Soilthermal();

/* **************************************************************
		 Public Functions
************************************************************** */
 void   getsnowecd(ofstream& rflog1);
 void   getsnowecd(char ecd[80]);
 void   getsoillecd(ofstream& rflog1);
 void   getsoillecd(char ecd[80]);
 void   showsnowecd(int& cmnt);
 void   showsoillecd(int& cmnt);
 void   getsoiltecd(ofstream& rflog1);
 void   getsoiltecd(char ecd[80]);
 void   showsoiltecd(int& cmnt);

 int soiltemp_(double *, double *, double *, double *, double *, double *,
	    double *, double *, double *, double *, integer *, integer *,
	    double *, double *, double *, double *, double *, double *, double *, double *,double *, double *, double *, int& cmnt);
 int grid_(integer *, int& cmnt);
 int sethet_(integer *maxm, integer *kint, double *heat, double *
	source, int& cmnt);
 int tinitl_(int& cmnt);
 int snofal_(double *snow, double *sden, double *comp, double *eta0,
	double *denfac, double *fcmelt, int& cmnt);

 /* **************************************************************
		 Public Variables
************************************************************** */
 integer ism19;
 integer i__;
 double airt19, airt29, airt39, tsoil,frontd,thawbegin,thawend, t9[211], x9[211], water9[211],
    smass9, zz, hsnow19, hsnow29, hsnow39, dx9[211];
 integer is9;
 double weight9[10];
 double xfa9[211], xfb9[211];
 double diffsoilt[10];

 double calcu_cdsnow; // to calculate snow thermal conductivity

  // modified with linkage between STM and HYM

  double water1; // moss layer water
  double water2; // orgainc layer water
  double water3;


 //for basic data structure
 double my_first,my_final,my_per,my_dtday, my_theta,my_hlat,my_tf,my_gflux;
 integer my_ig,my_nst;
 double my1_cdsnow;
 integer mykiso,mykallyr, myknodes,myktemp, myksnow,mykenvl, mykflux;
 double mytiso[10];
 integer myliso, mylmax,mynanmax,myldepf;
 double myvdep,myvdep1;
 double mydeptem[20],mydepflx[20];

// for grid

 double mythick[20],mydxa[20], mydxb[20],mydense[20],mywat[20];
 double  myvcond[20],myvsph[20];
 integer mymat[20];
 double mycondt[30][20],myspht[30][20];
 double mycondf[30][20],mysphf[30][20];
 double myvdens;
 double myvspace,mytop;
 double  myepsmin;

//for sethet
 integer mykinti;
 double myvdep2,mydephet;

//for initl
 integer myindex, mynp;
 double myvdepth;
 double mytstart[211], myxstart[211];

//for snofal
 integer myindex2;
 double myepssno,myconvrt,myeta0,mydenfac,myfcmelt,mydenmax;

 // read in parameters
 integer MAX[NUMVEG], NST[NUMVEG],KALLYR[NUMVEG],KNODES[NUMVEG],KISO[NUMVEG];
 integer KTEMP[NUMVEG],KSNOW[NUMVEG],KENVL[NUMVEG],KFLUX[NUMVEG],LISO[NUMVEG];
 double TISO[NUMVEG],VDEPTH[NUMVEG],DEPTEM1[NUMVEG],DEPTEM2[NUMVEG],DEPTEM3[NUMVEG],DEPTEM4[NUMVEG],DEPTEM5[NUMVEG],VDEP[NUMVEG];
 integer LMAX[NUMVEG], NDEPF[NUMVEG],IG[NUMVEG];
 double VDEPTH1[NUMVEG],DEPFLX1[NUMVEG],DEPFLX2[NUMVEG],DEPFLX3[NUMVEG],DEPFLX4[NUMVEG],DEPFLX5[NUMVEG];
 double HLAT[NUMVEG],TF[NUMVEG],gflux[NUMVEG],cdsnow[NUMVEG],FIRST[NUMVEG],FINAL[NUMVEG];
 double PER[NUMVEG],DTDAY[NUMVEG],THETA[NUMVEG],TOP[NUMVEG],EPSMIN[NUMVEG];
 double VSPACE[NUMVEG],VDEN[NUMVEG],VDEP1[NUMVEG],DEPHET[NUMVEG],EPSSNO[NUMVEG];
 double CONVRT[NUMVEG],ETAO[NUMVEG],DENFAC[NUMVEG],FCMELT[NUMVEG],DENMAX[NUMVEG];
 integer kint[NUMVEG],SNOFAL[NUMVEG] ;// for index;

 integer MAT1[NUMVEG],MAT2[NUMVEG],MAT3[NUMVEG],MAT4[NUMVEG],MAT5[NUMVEG],MAT6[NUMVEG];
 double THICK1[NUMVEG],DXA1[NUMVEG],DXB1[NUMVEG],DENSE1[NUMVEG],
      WATER1[NUMVEG],vcond1[NUMVEG],vsph1[NUMVEG], cond1[NUMVEG],spht1[NUMVEG],
      condf1[NUMVEG],sphf1[NUMVEG];
 double THICK2[NUMVEG],DXA2[NUMVEG],DXB2[NUMVEG],DENSE2[NUMVEG],
      WATER2[NUMVEG],vcond2[NUMVEG],vsph2[NUMVEG], cond2[NUMVEG],spht2[NUMVEG],
      condf2[NUMVEG],sphf2[NUMVEG];
 double THICK3[NUMVEG],DXA3[NUMVEG],DXB3[NUMVEG],DENSE3[NUMVEG],
      WATER3[NUMVEG],vcond3[NUMVEG],vsph3[NUMVEG], cond3[NUMVEG],spht3[NUMVEG],
      condf3[NUMVEG],sphf3[NUMVEG];
 double THICK4[NUMVEG],DXA4[NUMVEG],DXB4[NUMVEG],DENSE4[NUMVEG],
      WATER4[NUMVEG],vcond4[NUMVEG],vsph4[NUMVEG], cond4[NUMVEG],spht4[NUMVEG],
      condf4[NUMVEG],sphf4[NUMVEG];
 double THICK5[NUMVEG],DXA5[NUMVEG],DXB5[NUMVEG],DENSE5[NUMVEG],
      WATER5[NUMVEG],vcond5[NUMVEG],vsph5[NUMVEG], cond5[NUMVEG],spht5[NUMVEG],
      condf5[NUMVEG],sphf5[NUMVEG];

 double THICK6[NUMVEG],DXA6[NUMVEG],DXB6[NUMVEG],DENSE6[NUMVEG],WATER6[NUMVEG];

 double INDEX[NUMVEG],VDEPP[NUMVEG];
 integer NP[NUMVEG];
 double DEPTH[NUMVEG][25],TEMP[NUMVEG][25];
};
#endif

