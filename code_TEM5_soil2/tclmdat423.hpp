/* *************************************************************
TCLMDAT423.HPP - object to read and write the structure of the
                climate data from files used by the Terrestrial
                Ecosystem Model (TEM)
************************************************************** */

#if !defined(TCLMDAT423_H)
#define TCLMDAT423_H

class Clmdata {
  
  public:
    
     Clmdata(void);

/* **************************************************************
		      Public Functions
************************************************************** */

// read data structure.
     int get(ifstream& ifile);
     int getdel(FILE* ifile);
//write data structure.
     void out(ofstream& ofile, float col, float row, char varname[9], int carea, long year, double mon[], char contnent[9]);
     void outdel(ofstream& ofile, float col, float row, char varname[9], int carea, long year, double mon[CYCLE], char contnent[9]);
     void pctout(ofstream& ofile, float col, float row, char varname[9], int carea, long year, double mon[], char contnent[9]);
     void poutdel(ofstream& ofile, float col, float row, char varname[9], int carea, long year, double mon[CYCLE], char contnent[9]);


/* **************************************************************
		     Public Variables
************************************************************** */

/* NOTE:  substitute "daily" for "monthly" if CYCLE = 365 rather than
			 CYCLE = 12 */
	  
     float col;          // column or longitude of grid cell (degrees)
     float row;          // row or latitude of grid cell (degrees)
     char varname[9];    // climate variable name
     int carea;          // area covered by grid cell (sq. km)
     long year;          // date (year) of data
     double total;       // annual sum of monthly data for grid cell
     double max;         // maximum monthly value for grid cell
     double ave;         // mean annual value for grid cell
     double min;         // minimum monthly value for grid cell
     double mon[CYCLE];  // monthly values for the grid cell
     char contnent[9];   // name of continent containing grid cell


  private:

/* **************************************************************
		      Private Variables
************************************************************** */

     int clmend;
     long curpos;
     long lagpos;

};

#endif
