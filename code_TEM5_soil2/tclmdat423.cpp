/* *************************************************************
TCLMDAT423.CPP - object to read and write the structure of the
                climate data from files used by the Terrestrial
                Ecosystem Model (TEM)
************************************************************* */
#if !defined(TCLMDAT423_H)
  #include "tclmdat423.hpp"
#endif

Clmdata::Clmdata(void)
{

  clmend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

int Clmdata::get(ifstream& ifile)
{

  lagpos = ifile.tellg();

  ifile >> col >> row;
  ifile >> varname;
  ifile >> carea;
  ifile >> year;

  ifile >> total;
  ifile >> max;
  ifile >> ave;
  ifile >> min;

  for (int i = 0; i < CYCLE; i++) { ifile >> mon[i]; }

  ifile >> contnent;

  ifile.seekg(0, ios::cur);
  curpos = ifile.tellg();

  if (curpos < (lagpos + 10)) { clmend = -1; }

  return clmend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Clmdata::getdel(FILE* ifile)
{

  clmend = fscanf(ifile,"%f,%f, %s ,%d,%ld,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,\
				  %lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf, %s",
                  &col, &row, varname, &carea, &year, &total, &max, &ave, &min,
                  mon+0, mon+1, mon+2, mon+3, mon+4, mon+5,
                  mon+6, mon+7, mon+8, mon+9, mon+10, mon+11, contnent);

  return clmend;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Clmdata::out(ofstream& ofile, float col, float row, char varname[9],
                  int carea, long year, double mon[CYCLE], char contnent[9])
{

  int j;
  double predtotl;
  double predmax;
  double predave;
  double predmin;

  predtotl = 0.0;
  predmax   = -900000.0;
  predmin   =  900000.0;

  for (j = 0; j < CYCLE; j++) {
    if (mon[j] > predmax) { predmax = mon[j]; }
    if (mon[j] < predmin) { predmin = mon[j]; }
    predtotl += mon[j];
  }

  predave = predtotl / (double) CYCLE;
  if (predtotl <= (MISSING*CYCLE)) { predtotl = MISSING; }

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << carea << ' ';
  ofile << year << ' ';

  ofile << setprecision(1) << predtotl << ' ';
  ofile << predmax << ' ';
  ofile << setprecision(2) << predave << ' ';
  ofile << setprecision(1) << predmin << ' ';

  for (int i = 0; i < CYCLE; i++)  { ofile << mon[i] << ' '; }

  ofile << contnent << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Clmdata::outdel(ofstream& ofile, float col, float row, char varname[9],
                     int carea, long year, double mon[CYCLE], char contnent[9])
{

  int j;
  double predtotl;
  double predmax;
  double predave;
  double predmin;

  predtotl = 0.0;
  predmax   = -900000.0;
  predmin   =  900000.0;

  for (j = 0; j < CYCLE; j++) {
    if (mon[j] > predmax) { predmax = mon[j]; }
    if (mon[j] < predmin) { predmin = mon[j]; }
    predtotl += mon[j];
  }

  predave = predtotl / (double) CYCLE;
  if (predtotl <= (MISSING*CYCLE)) { predtotl = MISSING; }

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << carea << ",";
  ofile << year << ",";
  ofile << setprecision(1) << predtotl << ",";
  ofile << predmax << ",";
  ofile << setprecision(2) << predave << ",";
  ofile << setprecision(1) << predmin << ",";

  for (int i = 0; i < CYCLE; i++)  { ofile << mon[i] << ","; }

  ofile << " " << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Clmdata::pctout(ofstream& ofile, float col, float row, char varname[9],
                     int carea, long year, double mon[CYCLE], char contnent[9])
{

  int j;
  double predtotl;
  double predmax;
  double predave;
  double predmin;

  predtotl = 0.0;
  predmax   = -900000.0;
  predmin   =  900000.0;

  for (j = 0; j < CYCLE; j++) {
    if (mon[j] > predmax) predmax = mon[j];
    if (mon[j] < predmin) predmin = mon[j];
    predtotl += mon[j];
  }

  predave = predtotl / (double) CYCLE;
  if (predtotl <= (MISSING*CYCLE)) predtotl = MISSING;

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << ' ' << row << ' ';
  ofile << varname << ' ';
  ofile << setprecision(0) << carea << ' ';
  ofile << year << ' ';

  ofile << setprecision(2) << predtotl << ' ';
  ofile << predmax << ' ';
  ofile << predave << ' ';
  ofile << predmin << ' ';

  for (int i = 0; i < CYCLE; i++)  { ofile << mon[i] << ' '; }

  ofile << contnent;
  ofile << endl;

};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Clmdata::poutdel(ofstream& ofile, float col, float row, char varname[9],
                      int carea, long year, double mon[CYCLE], char contnent[9])
{

  int j;
  double predtotl;
  double predmax;
  double predave;
  double predmin;

  predtotl = 0.0;
  predmax   = -900000.0;
  predmin   =  900000.0;

  for (j = 0; j < CYCLE; j++) {
    if (mon[j] > predmax) { predmax = mon[j]; }
    if (mon[j] < predmin) { predmin = mon[j]; }
    predtotl += mon[j];
  }

  predave = predtotl / (double) CYCLE;
  if (predtotl <= (MISSING*CYCLE)) predtotl = MISSING;

  ofile.setf(ios::fixed,ios::floatfield);
  ofile.setf(ios::showpoint);
  ofile.precision(1);

  ofile << col << "," << row << ", ";
  ofile << varname << " ,";
  ofile << setprecision(0) << carea << ",";
  ofile << year << ",";

  ofile << setprecision(2) << predtotl << ",";
  ofile << predmax << ",";
  ofile << predave << ",";
  ofile << predmin << ",";

  for (int i = 0; i < CYCLE; i++) { ofile << mon[i] << ","; }

  ofile << " " << contnent;
  ofile << endl;

};
