/* **************************************************************
TKDDAT423.CPP - object to read and write the structure of spatially
               explicit values of Kd from/to files used by the 
               Terrestrial Ecosystem Model (TEM)
************************************************************** */

#if !defined(TKDDAT423_H)
  #include "tkddat423.hpp"
#endif

/* *********************************************************** */

KDdata::KDdata(void)
{

  kdend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int KDdata::get(ifstream& infile)
{

  lagpos = infile.tellg();

  infile >> col >> row;
  infile >> varname;
  infile >> tveg;
  infile >> subtveg;
  infile >> kd;
  infile >> contnent;

  infile.seekg(0, ios::cur);
  curpos = infile.tellg();

  if (curpos < (lagpos + 10)) { kdend = -1;}

  return kdend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int KDdata::getdel(FILE* infile)
{

  kdend = fscanf(infile,"%f,%f, %s ,%d,%d,%lf, %s",
                 &col,&row,varname,&tveg,&subtveg,&kd,contnent);

  return kdend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void KDdata::out(ofstream& ofile, float col, float row,
                 int tveg, int subtveg,
                 double kd, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << "KD" << ' ';
  ofile << setprecision(0) << tveg << ' ';
  ofile << subtveg << ' ';
  ofile << setprecision(6) << kd << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void KDdata::outdel(ofstream& ofile, float col, float row,
                    int tveg, int subtveg,
                    double kd, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row;
  ofile << ", KD ,";
  ofile << setprecision(0) << tveg << ",";
  ofile << subtveg << ",";
  ofile << setprecision(6) << kd << ", ";
  ofile << contnent;
  ofile << endl;

};

