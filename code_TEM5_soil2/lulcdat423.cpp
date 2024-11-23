/* **************************************************************
LULCDAT423.CPP - object to read and write the structure of land use
               data from/to files
************************************************************** */

#if !defined(LULCDAT423_H)
  #include "lulcdat423.hpp"
#endif

/* *********************************************************** */
Lulcdata::Lulcdata(void)
{

  lulcend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Lulcdata::get(ifstream& infile)
{

  lagpos = infile.tellg();

  infile >> col >> row;
  infile >> varname;
  infile >> carea;
  infile >> year;
  infile >> agstate;
  infile >> pctag;
  infile >> RAP;
  infile >> contnent;

  infile.seekg(0, ios::cur);
  curpos = infile.tellg();

  if (curpos < (lagpos + 10)) { lulcend = -1;}

  return lulcend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Lulcdata::getdel(FILE* infile)
{

  lulcend = fscanf(infile,"%f,%f, %s ,%d,%d,%d,%lf,%lf, %s",
                   &col,&row,varname,&carea,&year,&agstate,&pctag,&RAP,
                   contnent);

  return lulcend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Lulcdata::out(ofstream& ofile, float col, float row, char varname[9],
                   int carea, int year, int agstate, double pctag, double RAP,
                   char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << carea << ' ';
  ofile << year << ' ';
  ofile << agstate << ' ';
  ofile << setprecision(3) << pctag << ' ';
  ofile << setprecision(2) << RAP << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Lulcdata::outdel(ofstream& ofile, float col, float row, char varname[9],
                      int carea, int year, int agstate, double pctag,
                      double RAP, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << carea << ",";
  ofile << year << ",";
  ofile << agstate << ",";
  ofile << setprecision(3) << pctag << ",";
  ofile << setprecision(2) << RAP << ", ";
  ofile << contnent;
  ofile << endl;

};

