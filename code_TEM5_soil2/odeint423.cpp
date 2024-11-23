/* **************************************************************
ODEINT423.CPP - object contains functions to integrate ordinary
	       differential equations (ODE)
*****************************************************************

Modifications:

19991102 - DWK deleted code for:
           virtual int boundcon(double ptstate[], double err[], double& ptol)
           virtual void delta(int& pdm, double pstate[], double pdstate[])

20000102 - DWK added compiler directives
************************************************************* */

#if !defined(ODEINT423_H)
  #include "odeint423.hpp"
#endif

/* *************************************************************
************************************************************* */

Odeint4::Odeint4()
{

  syint = 1;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Odeint4::adapt(const int& numeq, double pstate[], double& ptol,
			 const int& pdm)
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
	  }
	}    /* end rkf integrator */
  }      /* end time while */
//  cout <<" mflag "<< endl;
  return mflag;

};

/* *************************************************************
************************************************************* */
void Odeint4::ask(Temstac& temstac)
{
	inittol=temstac.inittol;
	maxit=temstac.maxit;
	maxitmon=temstac.maxitmon;
	
}


/* *************************************************************
************************************************************* */

void Odeint4::ask(ofstream& rflog1, Temstac& temstac)
{


/* **************************************************************
	      Parameters for Adaptive Integrator
************************************************************** */

  cout << endl << "Enter the proportional tolerance for the integrator: ";
//cin >> inittol;
  fpara >> inittol;
	temstac.inittol=inittol;
  rflog1 << endl << "Enter the proportional tolerance for the integrator: "
		 << inittol << endl;

  cout << "Enter the maximum number of iterations in the integrator: ";
//cin >> maxit;
  fpara >> maxit;
	temstac.maxit=maxit;
  rflog1 << "Enter the maximum number of iterations in the integrator: "
		<< maxit << endl;

  cout << "Enter the maximum number of times in a month that the" << endl;
  cout << "integrator can reach the maximum number of iterations: ";
//cin >> maxitmon;
  fpara >> maxitmon;
	temstac.maxitmon=maxitmon;
  rflog1 << "Enter the maximum number of times in a month that the" << endl;
  rflog1 << "integrator can reach the maximum number of iterations: "
  		<< maxitmon << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Odeint4::rkf(const int& numeq, double pstate[], double& pdt,
		 const int& pdm)
{

  int i;
  double ptdt = 0;

  for (i = 0; i < numeq;i++)
  {
    dum4[i] = dum5[i] = pstate[i];
    yprime[i] = rk45[i] = error[i] = 0.0;
  }

  ptdt = pdt * 0.25;

  delta(pdm,dum4,f11);
  step(numeq,yprime,f11,yprime,a1);
  step(numeq,rk45,f11,rk45,b1);
  step(numeq,dum4,f11,ydum,ptdt);
  delta(pdm,ydum,f2);
  for (i = 0; i < numeq; i++)
  {
    f13[i] = a31*f11[i] + a32*f2[i];
  }
  step(numeq,dum4,f13,ydum,pdt);
  delta(pdm,ydum,f3);
  step(numeq,yprime,f3,yprime,a3);
  step(numeq,rk45,f3,rk45,b3);
  for (i = 0; i < numeq; i++)
  {
    f14[i] = a41*f11[i] + a42*f2[i] + a43*f3[i];
  }
  step(numeq,dum4,f14,ydum,pdt);
  delta(pdm,ydum,f4);
  step(numeq,yprime,f4,yprime,a4);
  step(numeq,rk45,f4,rk45,b4);
  for (i = 0; i < numeq; i++)
  {
    f15[i] = a51*f11[i] + a52*f2[i] + a53*f3[i] + a54*f4[i];
  }
  step(numeq,dum4,f15,ydum,pdt);
  delta(pdm,ydum,f5);
  step(numeq,yprime,f5,yprime,a5);
  step(numeq,rk45,f5,rk45,b5);
  for (i = 0; i < numeq; i++)
  {
    f16[i] = b61*f11[i] + b62*f2[i] + b63*f3[i] + b64*f4[i] + b65*f5[i];
  }
  step(numeq,dum4,f16,ydum,pdt);
  delta(pdm,ydum,f6);
  step(numeq,rk45,f6,rk45,b6);
  step(numeq,dum4,yprime,dum4,pdt);
  step(numeq,dum5,rk45,dum5,pdt);
  for (i = 0; i < numeq; i++)
  {
    error[i] = fabs(dum4[i] - dum5[i]);
  }

};

/***************************************************************
 ***************************************************************/

void Odeint4::step(const int& numeq, double pstate[], double pdstate[],
		 double ptstate[],
			 double& pdt)
{

  for (int i = 0; i < numeq; i++)
  {
    ptstate[i] = pstate[i] + (pdt * pdstate[i]);
  }

};

