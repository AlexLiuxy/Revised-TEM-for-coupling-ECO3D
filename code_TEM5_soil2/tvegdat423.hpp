/* **************************************************************
TVEGDAT423.HPP - object to read and write the structure of vegetation
	       data from/to files used by the Water Balance Model
	       (WBM) and the Terrestrial Ecosystem Model (TEM)
************************************************************** */

#if !defined(TVEGDAT423_H)
#define TVEGDAT423_H

class Vegdata {
   
  public:
      
     Vegdata(void);
     
/* **************************************************************
		    Public Functions
************************************************************** */
     
     int get(ifstream& ifile);
     int getdel(FILE* infile);
     void out(ofstream& ofile, float col, float row, char varname[9], long year, int temveg, char source[9], char contnent[9]);
     void outdel(ofstream& ofile, float col, float row, char varname[9], long year, int temveg, char source[9], char contnent[9]);

/* **************************************************************
		     Public Variables
************************************************************** */
	  
     float col;          // column or longitude of grid cell (degrees)
     float row;          // row or latitude of grid cell (degrees)
     char varname[9];    // name of vegetation classification
     long year;          // date (year) of data (transient version of TEM)
			 // (year = -111 for equilibrium version of TEM)
     int temveg;         // vegetation type of grid cell (categorical data)
     char source[9];     // source of vegetation data 
     char contnent[9];   // name of continent containing grid cell


  private:

/* **************************************************************
		      Private Variables
************************************************************** */

     int vegend;
     long curpos;
     long lagpos;


};

#endif

