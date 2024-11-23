/* **************************************************************
TELVDAT423.HPP - object to read and write the structure of elevation
               data from/to files used by the Water Balance Model
************************************************************** */

#if !defined(TELVDAT423_H)
#define TELVDAT423_H

class Elevdata {
   
  public:
          
     Elevdata(void);

/* **************************************************************
                      Public Functions
************************************************************** */

// read data structure.
     int get(ifstream& infile);
     int getdel(FILE* infile);
//write data structure.
     void out(ofstream& ofile, float col, float row, char varname[9], int carea, double elev, char contnent[9]);
     void outdel(ofstream& ofile, float col, float row, char varname[9], int carea, double elev, char contnent[9]);

          
/* **************************************************************
                     Public Variables
************************************************************** */

     float col;          // column or longitude of grid cell (degrees)
     float row;          // row or latitude of grid cell (degrees)
     char varname[9];    // "ELEV"
     int carea;          // area covered by grid cell (sq. km)
     double elev;        // elevation of grid cell (m)
     char contnent[9];   // name of continent containing grid cell


  private:

/* **************************************************************
                      Private Variables
************************************************************** */

     int elvend;
     long curpos;
     long lagpos;

};

#endif

