/* *************************************************************
ELMNT423.CPP - Contains functions to manage elements used in GIS
************************************************************** */

#if !defined(ELMNT423_H)
  #include "elmnt423.hpp"
#endif

/* ********************************************************** */

Elmnt::Elmnt(void) {

  count = 0;
  numskip  = 0;
  numgrids = 64000;
  grdcnt  = 64000;
  totyr    = -99;
  endflag  = 1;
  stopflag = 0;

  col = row = -999.9;
  carea = -99;
  elev  = -999.9;

  return;
};

/* *************************************************************
************************************************************* */

void Elmnt::ask(ofstream& rflog1) {

  cout << endl << endl << "For how many grid cells would"
	   << " you like a timestamp?" << endl;
  cout << "(After 'x' number of grid cells)" << endl;
  cout << "(Enter 0 for no timestamp):  ";
//  cin >> grdcnt;
  fpara >> grdcnt;
  if (grdcnt == 0) {
    grdcnt = -99;
    cout << grdcnt << endl;
  }

  rflog1 << endl << endl << "For how many grid cells would you"
		 << " like a timestamp?" << endl;
  rflog1 << "After 'x' number of grid cells" << endl;
  rflog1 << "Enter 0 for no timestamp): " << grdcnt << endl;

  cout << endl << "Do you want to start at the beginning"
  	   << " of the GIS files?" << endl;
  cout << "  Enter 1 for YES" << endl;
  cout << "  Enter 0 for NO:  ";
//cin >> strtflag;
  fpara >> strtflag;

  rflog1 << endl << endl << "Do you want to start at the beginning"
  		 << " of the GIS files?" << endl;
  rflog1 << "  Enter 1 for YES" << endl;
  rflog1 << "  Enter 0 for NO: " << strtflag << endl << endl << endl;

  if (strtflag == 0) {
    cout << endl << "How many records would you like to skip? ";
//  cin >> numskip;
	fpara >> numskip;

	rflog1 << endl << "How many records would you like to skip? "
	       << numskip << endl << endl << endl;
  }

  cout << endl << "Do you want to finish at the end of the GIS files?" << endl;
  cout << "  Enter 1 for YES" << endl;
  cout << "  Enter 0 for NO:  ";
//cin >> endflag;
  fpara >> endflag;

  rflog1 << endl << endl << "Do you want to finish at the end of"
  	     << " the GIS files?" << endl;
  rflog1 << "  Enter 1 for YES" << endl;
  rflog1 << "  Enter 0 for NO: " << endflag << endl << endl << endl;

  if (endflag == 0) {
    cout << endl << "For how many records would you like to make estimates? ";
//  cin >> numgrids;
	fpara >> numgrids;
	rflog1 << endl << endl << "For how many records would you like"
		   << " to make estimates" << numgrids << endl << endl << endl;
  }

};

/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

void Elmnt::show(ofstream& rflog1, const float& col, const float& row) {

  time_t timer;

  timer = time(NULL);

   cout.setf(ios::fixed,ios::floatfield);
   cout.setf(ios::showpoint);
   rflog1.setf(ios::fixed,ios::floatfield);
   rflog1.setf(ios::showpoint);


  if (count == 0 || (grdcnt != -99 && count <= grdcnt)) {
	cout << "Finished cell " << (count+numskip) << " ("
		<< setprecision(1) << col << " , " << row << ") " << ctime(&timer);
  }

	if (count == grdcnt || (count < grdcnt && endflag == 0
		&& count == numgrids)) cout << "Finished printing to the screen."
			<< "  GOOD LUCK!!!!!\n  ";

	rflog1 << "Finished cell " << setprecision(0) << (count+numskip)
		<< " (" << setprecision(1) << col << " , " << row << ") "
		<< ctime(&timer);

    if (endflag == 0 && count == numgrids) stopflag = 1;

    ++count;


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

void Elmnt::show(ofstream& rflog1, const float& col,
			 const float& row, const long& totyr, const double& tol)
{

  time_t timer;

  timer = time(NULL);

  cout.setf(ios::fixed,ios::floatfield);
  cout.setf(ios::showpoint);
  rflog1.setf(ios::fixed,ios::floatfield);
  rflog1.setf(ios::showpoint);

   if (count == 0 || (grdcnt != -99 && count <= grdcnt)) {
	 cout << "Finished cell " << (count+numskip) << " ("
	 << setprecision(1) << col << " , " << row << ")  TOTYR = "
	 << setprecision(0) << totyr << " TOL = " << setprecision(6)
	 << tol << " " << ctime(&timer);
   }

   if (count == grdcnt || (count < grdcnt && endflag == 0
   && count == numgrids)) cout << "Finished printing to the screen."
   		<< "  GOOD LUCK!!!!!\n  ";

   rflog1 << "Finished cell " << setprecision(0) << (count+numskip)
		<< " (" << setprecision(1) << col << " , " << row
		<< ")  TOTYR = " << setprecision(0) << totyr
		<< " TOL = " << setprecision(6) << tol << " " << ctime(&timer);

   if (endflag == 0 && count == numgrids) stopflag = 1;

   ++count;


};

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Elmnt::coregerr(ofstream& rflog1, char varname1[9],
	 const float& col1, const float& row1, char varname2[9],
	 const float& col2, const float& row2)
{

  int fatalerr = 0;

  if (col1 != col2 || row1 != row2) {
	fatalerr = 1;

	cout << "ERROR:  " << varname1 << " data and "
	<< varname2 << "data are not coregistered.\n";
	cout << "COL = " << col1 << " and ROW = " << row1
	<< " in " << varname1 << " data" << endl;
	cout << "COL = " << col2 << " and ROW = " << row2
	<< " in " << varname2 << " data" << endl;


	rflog1 << "ERROR:  " << varname1 << " data and "
		<< varname2 << "data are not coregistered.\n";
	rflog1 << "COL = " << col1 << " and ROW = " << row1 << " in "
		<< varname1 << " data" << endl;
	rflog1 << "COL = " << col2 << " and ROW = " << row2 << " in "
		<< varname2 << " data" << endl;

  }

  return fatalerr;
};

