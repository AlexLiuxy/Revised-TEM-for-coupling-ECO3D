/* **************************************************************
TKDDAT423.HPP - object to read and write the structure of spatially
              values of Kd from/to files used by the Terrestrial
              Ecosystem Model (TEM)
************************************************************** */

#if !defined(TKDDAT423_H)
#define TKDDAT423_H

class KDdata {
   
  public:
          
     KDdata(void);

/* **************************************************************
                      Public Functions
************************************************************** */

// read data structure.
     int get(ifstream& infile);
     int getdel(FILE* infile);
//write data structure.
     void out(ofstream& ofile, float col, float row, int tveg, int subtveg, double kd, char contnent[9]);
     void outdel(ofstream& ofile, float col, float row, int tveg, int subtveg, double kd, char contnent[9]);

          
/* **************************************************************
                     Public Variables
************************************************************** */

     float col;          // column or longitude of grid cell (degrees)
     float row;          // row or latitude of grid cell (degrees)
     char varname[9];    // "KD"
     int tveg;           // Vegetation type
     int subtveg;        // Plant functional type
     double kd;          // the parameter Kd
     char contnent[9];   // name of continent containing grid cell


  private:

/* **************************************************************
                      Private Variables
************************************************************** */

     int kdend;
     long curpos;
     long lagpos;

};

#endif

