/* **************************************************************
*****************************************************************
ODEINT423.HPP - object contains functions to integrate ordinary
	       differential equations (ODE)
*****************************************************************

Modifications:

20000102 - DWK added compiler directives

References:

Runge - Kutta - Fehlberg (RKF):

Cheney, W., and D. Kincaid.  1985.  Numerical mathematics and
  computing.  Second edition.  Brooks/ Col Publishing Co.  Monterey,
  CA.  pp. 325-328.

*****************************************************************
************************************************************** */

// Odeint4 uses the global constants ACCEPT and NUMEQ

#if !defined(ODEINT423_H)
#define ODEINT423_H

class Odeint4 {
   
   public:

      Odeint4();

/* **************************************************************
		 Public Functions
************************************************************** */

      int adapt(const int& numeq, double pstate[], double& ptol, const int& pdm);
      void ask(ofstream& rflog1, Temstac& temstac);      
      void ask(Temstac& temstac);

      void rkf(const int& numeq, double pstate[], double& pdt, const int& pdm);

      // virtual functions for Odeint::adapt() and Odeint::rkf()  (Set to 0 by DWK on 19991102)

      virtual int boundcon(double ptstate[], double err[], double& ptol) = 0;
      virtual void delta(const int& pdm, double pstate[], double pdstate[]) = 0;


/* **************************************************************
		 Public Variables
************************************************************** */

      int blackhol;
      static double inittol;
      static int maxit;
      static long maxitmon;
      int retry;

      int syint;
      int test;

      double dum4[NUMEQ];
      double error[NUMEQ];

  private:

/* **************************************************************
		 Private Functions
************************************************************** */

      void step(const int& numeq, double pstate[], double pdstate[], double ptstate[], double& pdt);

/* **************************************************************
		 Private Variables
************************************************************** */


      double dum5[NUMEQ];
      double dumy[NUMEQ];

      double ydum[NUMEQ];
      double dy[NUMEQ];
      double yprime[NUMEQ];
      double rk45[NUMEQ];

      double f11[NUMEQ];
      double f2[NUMEQ];
      double f13[NUMEQ];
      double f3[NUMEQ];
      double f14[NUMEQ];
      double f4[NUMEQ];
      double f15[NUMEQ];
      double f5[NUMEQ];
      double f16[NUMEQ];
      double f6[NUMEQ];


      static double  a1;

      static double  a3;
      static double a31;
      static double a32;

      static double  a4;
      static double a41;
      static double a42;
      static double a43;

      static double  a5;
      static double a51;
      static double a52;
      static double a53;
      static double a54;

      static double  b1;
      static double  b3;
      static double  b4;
      static double  b5;

      static double  b6;
      static double b61;
      static double b62;
      static double b63;
      static double b64;
      static double b65;

      double dummy;

};

// Initialize static members

double Odeint4::inittol = 0.01;
int Odeint4::maxit = 20;
long Odeint4::maxitmon = 2000;

double Odeint4::a1 =   0.115740741;

double   Odeint4::a3 =   0.548927875;
double  Odeint4::a31 =   0.09375;
double  Odeint4::a32 =   0.28125;

double   Odeint4::a4 =   0.535331384;
double  Odeint4::a41 =   0.879380974;
double  Odeint4::a42 =  -3.277196177;
double  Odeint4::a43 =   3.320892126;

double   Odeint4::a5 =  -0.20;
double  Odeint4::a51 =   2.032407407;
double  Odeint4::a52 =  -8.0;
double  Odeint4::a53 =   7.173489279;
double  Odeint4::a54 =  -0.2058966866;

double   Odeint4::b1 =   0.118518519;
double   Odeint4::b3 =   0.518986355;
double   Odeint4::b4 =   0.50613149;
double   Odeint4::b5 =  -0.18;
double   Odeint4::b6 =   0.036363636;
double  Odeint4::b61 =  -0.296296296;
double  Odeint4::b62 =   2.0;
double  Odeint4::b63 =  -1.381676413;
double  Odeint4::b64 =   0.45297271;
double  Odeint4::b65 =  -0.275;



#endif
