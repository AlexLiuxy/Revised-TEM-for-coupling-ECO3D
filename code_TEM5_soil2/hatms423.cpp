/* **************************************************************
*****************************************************************
HATMS423.CPP - uses the physical characteristics of the atmosphere
               as described by TEM, but uses transient CO2, specified
               every 6 months, from an ASCII file
*****************************************************************
************************************************************** */

#if !defined(HATMS423_H)
  #include "hatms423.hpp"
#endif


MPIatms::MPIatms() : Atmosphere() { };

/* *************************************************************
************************************************************** */

void MPIatms::loadfosfuel(ofstream& rflog1, const int& RTIME)
{

  int iyear;
  char ifilename[40];
  int testyear;
  FILE* ffossil;

  cout << endl << endl << "Enter the name of the file containing"
	   << " the fossil fuel data: ";
// cin >> ifilename;
  fpara >> ifilename;
  rflog1 << endl << endl << "Enter the name of the file"
  		 << " containing the fossil fuel data: ";
  rflog1 << ifilename << endl;
  ffossil = fopen(ifilename, "r");

  if (!ffossil)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

  cout << "What year should the data be normalized against?" << endl;
  cin >> testyear;


// Initialize historical fossil fuel emissions

  for (iyear = 0; iyear < RTIME; iyear++)
  {
    fscanf(ffossil, "%d %lf", ffuelyear+iyear, ffuel+iyear);
    if (ffuelyear[iyear] == testyear) { maxffuel = ffuel[iyear]; }
  }

  fclose(ffossil);

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MPIatms::loadmpitCO2(ofstream& rflog1, const int& totsptime,
	 const int& RTIME)
{

  int i;
  int j;
  int dm;
  char ifilename[40];
  CO2data co2dat;
  ifstream fco2;
  double mco2[MAXRTIME];         // CO2 concentration in July of year

  cout << endl << endl << "Enter the name of the file containing"
  	   << " the CO2 data: ";
//cin >> ifilename;
  fpara >> ifilename;
  
  rflog1 << endl << endl << "Enter the name of the file containing"
  	     << " the CO2 data: ";
  rflog1 << ifilename << endl;
  fco2.open(ifilename, ios::in);

  if (!fco2)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

// Initialize transient CO2 concentrations

  for (i = (totsptime+1); i < (RTIME+1); i++)
  {
    j = i - (totsptime + 1);
    co2dat.get(fco2);
    co2year[i] = (long) co2dat.year;
    mco2[j] = co2dat.mco2;
  }
  for (i = (totsptime+1); i < RTIME; i++)
  {
    j = i - (totsptime + 1);
    for (dm = 0; dm < CYCLE; dm++)
    {
      if (j == 0) { tco2[i][dm] = mco2[j]; }
      else
      {
        if (dm < 6)
        {
		  tco2[i][dm] = mco2[j-1] + ((dm + 6) * (mco2[j] - mco2[j-1])
			/ (double) CYCLE);
        }
        else
        {
		  tco2[i][dm] = mco2[j] + ((dm - 6) * (mco2[j+1] - mco2[j])
			/ (double) CYCLE);
        }
      }
    }
  }

  fco2.close();

  co2year[0] = co2year[totsptime+1] - totsptime - 1;
  for (dm = 0; dm < CYCLE; dm++) { tco2[0][dm] = co2level; }
  for (i = 1; i < (totsptime+1); i++)
  {
    co2year[i] = co2year[0] + i;
    for (dm = 0; dm < CYCLE; dm++) { tco2[i][dm] = mco2[0]; }
  }

};

/* *************************************************************
************************************************************** */


/* *************************************************************
************************************************************** */

void MPIatms::loadmpitCO2(char ifilename[80], const int& totsptime,
		 const int& RTIME)
{

  int i;
  int j;
  int dm;
  CO2data co2dat;
  ifstream fco2;
  double mco2[MAXRTIME];         // CO2 concentration in July of year

  fco2.open(ifilename, ios::in);

  if (!fco2)
  {
    cerr << "\nCannot open " << ifilename << " for data input" << endl;
    exit(-1);
  }

// Initialize transient CO2 concentrations

  for (i = (totsptime+1); i < (RTIME+1); i++)
  {
    j = i - (totsptime + 1);
    co2dat.get(fco2);
    co2year[i] = (long) co2dat.year;
    mco2[j] = co2dat.mco2;
  }
  for (i = (totsptime+1); i < RTIME; i++)
  {
    j = i - (totsptime + 1);
    for (dm = 0; dm < CYCLE; dm++)
    {
      if (j == 0) { tco2[i][dm] = mco2[j]; }
      else
      {
        if (dm < 6)
        {
		  tco2[i][dm] = mco2[j-1] + ((dm + 6) * (mco2[j] - mco2[j-1])
		  / (double) CYCLE);
		}
        else
        {
		  tco2[i][dm] = mco2[j] + ((dm - 6) * (mco2[j+1] - mco2[j])
		  / (double) CYCLE);
        }
      }
    }
  }

  fco2.close();

  co2year[0] = co2year[totsptime+1] - totsptime - 1;
  for (dm = 0; dm < CYCLE; dm++) { tco2[0][dm] = co2level; }
  for (i = 1; i < (totsptime+1); i++)
  {
    co2year[i] = co2year[0] + i;
    for (dm = 0; dm < CYCLE; dm++) { tco2[i][dm] = mco2[0]; }
  }

};



