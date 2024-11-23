/* **************************************************************
TELVDAT423.CPP - object to read and write the structure of elevation
                data from/to files used by the Water Balance Model
************************************************************** */

#if !defined(TELVDAT423_H)
  #include "telvdat423.hpp"
#endif

/* *********************************************************** */

Elevdata::Elevdata(void)
{

  elvend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Elevdata::get(ifstream& infile)
{

  lagpos = infile.tellg();

  infile >> col >> row;
  infile >> varname;
  infile >> carea;
  infile >> elev;
  infile >> contnent;

  infile.seekg(0, ios::cur);
  curpos = infile.tellg();

  if (curpos < (lagpos + 10)) { elvend = -1;}

  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Elevdata::getdel(FILE* infile)
{

  elvend = fscanf(infile,"%f,%f, %s ,%d,%lf, %s",
                  &col,&row,varname,&carea,&elev,contnent);

  return elvend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Elevdata::out(ofstream& ofile, float col, float row, char varname[9],
                   int carea, double elev, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << carea << ' ';
  ofile << setprecision(1) << elev << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Elevdata::outdel(ofstream& ofile, float col, float row, char varname[9],
                      int carea, double elev, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << carea << ",";
  ofile << setprecision(1) << elev << ", ";
  ofile << contnent;
  ofile << endl;

};

