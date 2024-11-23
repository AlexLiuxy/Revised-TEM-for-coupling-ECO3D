/* **************************************************************
LATDAT423.HPP - object to read and write the structure of latitude
               and longitude data from/to files
************************************************************** */

#if !defined(LATDAT423_H)
#define LATDAT423_H

class Latdata {

  public:

     Latdata(void);

/* **************************************************************
                      Public Functions
************************************************************** */

// read data structure.
     int get(ifstream& infile);
     int getdel(FILE* infile);
//write data structure.
     void out(ofstream& ofile, float col, float row, char varname[], double lat, double lon, char contnent[9]);
     void outdel(ofstream& ofile, float col, float row, char varname[], double lat, double lon, char contnent[9]);


/* **************************************************************
                     Public Variables
************************************************************** */

     float col;          // column of grid cell
     float row;          // row or of grid cell
     char varname[9];    // "LATITUDE?"
     double lat;          // latitude of grid cell (degrees)
     double lon;          // longitude of grid cell (degrees)
     char contnent[9];   // name of continent containing grid cell


  private:

/* **************************************************************
                      Private Variables
************************************************************** */

     int latend;
     long curpos;
     long lagpos;

};

#endif

