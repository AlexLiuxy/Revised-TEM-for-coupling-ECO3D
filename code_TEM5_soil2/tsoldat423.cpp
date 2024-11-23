/* **************************************************************
TSOLDAT423.CPP - object to read and write the structure of soil
                texture data from/to files used by the Water Balance 
                Model (WBM) and the Terrestrial Ecosystem Model (TEM)
************************************************************** */

#if !defined(TSOLDAT423_H)
  #include "tsoldat423.hpp"  
#endif

Soildata::Soildata(void)
{

  soilend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */

int Soildata::get(ifstream& infile)
{

  lagpos = infile.tellg();

  infile >> col >> row;
  infile >> varname;
  infile >> carea;
  infile >> pctsand;
  infile >> pctsilt;
  infile >> pctclay;
  infile >> wsoil;
  infile  >> source;
  infile >> contnent;

  infile.seekg(0, ios::cur);
  curpos = infile.tellg();

  if (curpos < (lagpos + 10)) { soilend = -1; }

  return soilend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Soildata::getdel(FILE* infile)
{

  soilend = fscanf(infile,"%f,%f, %s ,%d,%lf,%lf,%lf,%d, %s , %s ",
                   &col,&row,varname,&carea,&pctsand,&pctsilt,&pctclay,
                   &wsoil,source,contnent);
  
  return soilend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soildata::out(ofstream& ofile, float col, float row, char varname[9],
                   int carea, double pctsand, double pctsilt, double pctclay,
                   int wsoil, char source[9], char contnent[9])
{

   ofile.setf(ios::fixed,ios::floatfield);
   ofile.setf(ios::showpoint);
   ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << carea << ' ';
  ofile << setprecision(2) << pctsand << ' ';
  ofile << pctsilt << ' ';
  ofile << pctclay << ' ';
  ofile << setprecision(0) << wsoil << ' ';
  ofile << source << ' ';
  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Soildata::outdel(ofstream& ofile, float col, float row, char varname[9],
                      int carea, double pctsand, double pctsilt, double pctclay,
                      int wsoil, char source[9], char contnent[9])
{

   ofile.setf(ios::fixed,ios::floatfield);
   ofile.setf(ios::showpoint);
   ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << carea << ",";
  ofile << setprecision(2) << pctsand << ",";
  ofile << pctsilt << ",";
  ofile << pctclay << ",";
  ofile << setprecision(0) << wsoil << ", ";
  ofile << source << " , ";
  ofile << contnent;
  ofile << endl;

};

