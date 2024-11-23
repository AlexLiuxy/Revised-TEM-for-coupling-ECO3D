/* **************************************************************
TCO2DAT423.CPP - object to read and write the structure of transient
               atmospheric CO2 data from/to files used by the
               Terrestrial Ecosystem Model (TEM)

              20000102 - DWK added compiler directives
************************************************************** */

#if !defined(TCO2DAT423_H)
  #include "tco2dat423.hpp"
#endif

/* *************************************************************
************************************************************* */

CO2data::CO2data(void)
{

  co2end = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

int CO2data::get(ifstream& infile)
{

  lagpos = infile.tellg();

  infile >> year >> mco2;

  infile.seekg(0, ios::cur);
  curpos = infile.tellg();

  if (curpos < (lagpos + 5)) { co2end = -1; }

  return co2end;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int CO2data::getdel(FILE* infile)
{

  co2end = fscanf(infile, "%f,%lf", &year, &mco2);

  return co2end;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void CO2data::out(ofstream& ofile, float& year, double& mco2)
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << year << ' ';
  ofile << setprecision(4) << mco2;
  ofile << endl;


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void CO2data::outdel(ofstream& ofile, float& year, double& mco2)
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << year << ",";
  ofile << setprecision(4) << mco2;
  ofile << endl;

};

