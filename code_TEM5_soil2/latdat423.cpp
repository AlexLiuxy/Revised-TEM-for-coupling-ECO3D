/* **************************************************************
LATDAT423.CPP - object to read and write the structure of land use
               data from/to files
************************************************************** */

#if !defined(LATDAT423_H)
  #include "latdat423.hpp"
#endif

/* *********************************************************** */

Latdata::Latdata(void)
{

  latend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Latdata::get(ifstream& infile)
{

  lagpos = infile.tellg();

  infile >> col >> row;
  infile >> varname;
  infile >> lat;
  infile >> lon;
  infile >> contnent;

  infile.seekg(0, ios::cur);
  curpos = infile.tellg();

  if (curpos < (lagpos + 10)) { latend = -1;}

  return latend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Latdata::getdel(FILE* infile)
{

  latend = fscanf(infile,"%f,%f, %s ,%lf,%lf, %s",
                  &col,&row,varname,&lat,&lon,contnent);

  return latend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Latdata::out(ofstream& ofile, float col, float row, char varname[9],
                  double lat, double lon, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(7) << lat << ' ';
  ofile << lon << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Latdata::outdel(ofstream& ofile, float col, float row, char varname[9],
                     double lat, double lon, char contnent[9])
{

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(7) << lat << ",";
  ofile << lon << ", ";
  ofile << contnent;
  ofile << endl;

};

