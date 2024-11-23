/* **************************************************************
TVEGDAT423.CPP - object to read and write the structure of vegetation
               data from/to files used by the Water Balance Model
               (WBM) and the Terrestrial Ecosystem Model (TEM)
************************************************************** */

#if !defined(TVEGDAT423_H)
  #include "tvegdat423.hpp"  
#endif

Vegdata::Vegdata(void)
{

  vegend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Vegdata::get(ifstream& infile)
{

  lagpos = infile.tellg();

  infile >> col >> row;
  infile >> varname;
  infile >> year;
  infile >> temveg;
  infile >> source;
  infile >> contnent;

  infile.seekg(0, ios::cur);
  curpos = infile.tellg();

  if (curpos < (lagpos + 10)) { vegend = -1; }

  return vegend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Vegdata::getdel(FILE* infile)
{

  vegend = fscanf(infile, "%f,%f, %s ,%ld,%d, %s , %s",
                  &col, &row, varname, &year, &temveg, source, contnent);

  return vegend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Vegdata::out(ofstream& ofile, float col, float row, char varname[9],
                  long year, int temveg,char source[9],char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << year << ' ';
  ofile << temveg << ' ';
  ofile << source << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Vegdata::outdel(ofstream& ofile, float col, float row, char varname[9],
                     long year, int temveg,char source[9],char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << year << ",";
  ofile << temveg << ", ";
  ofile << source << " , ";
  ofile << contnent;
  ofile << endl;

};

