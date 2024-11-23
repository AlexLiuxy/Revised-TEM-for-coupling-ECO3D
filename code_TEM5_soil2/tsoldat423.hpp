/* **************************************************************
TSOLDAT423.HPP - object to read and write the structure of soil
                texture data from/to files used by the Water Balance 
                Model (WBM) and the Terrestrial Ecosystem Model (TEM)
************************************************************** */

#if !defined(TSOLDAT423_H)
#define TSOLDAT423_H

class Soildata {
  
  public:
	 
     Soildata(void);

/* **************************************************************
		    Public Functions
************************************************************** */
     
     int get(ifstream& infile);
     int getdel(FILE* infile);
     void out(ofstream& ofile, float col, float row, char varname[9], int carea, double pctsand, double pctsilt, double pctclay, int wsoil, char source[9], char contnent[9]);
     void outdel(ofstream& ofile, float col, float row, char varname[9], int carea, double pctsand, double pctsilt, double pctclay, int wsoil, char source[9], char contnent[9]);

/* **************************************************************
		     Public Variables
************************************************************** */
	 
     float col;           // column or longitude of grid cell (degrees)
     float row;           // row or latitude of grid cell (degrees)
     char varname[9];     // "TEXTURE"
     int carea;           // area covered by grid cell (sq. km)
     double pctsand;      // percent sand of grid cell's soil texture
     double pctsilt;      // percent silt of grid cell's soil texture
     double pctclay;      // percent clay of grid cell's soil texture
     int wsoil;           // wetland soil type designation (categorical data)
     char source[9];      // reference to data source
     char contnent[9];    // name of continent containing grid cell


  private:

/* **************************************************************
		      Private Variables
************************************************************** */

     int soilend;
     long curpos;
     long lagpos;

};

#endif

