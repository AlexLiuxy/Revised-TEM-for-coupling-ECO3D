/* **************************************************************
LULCDAT423.HPP - object to read and write the structure of land use
               data from/to files
************************************************************** */

#if !defined(LULCDAT423_H)
#define LULCDAT423_H

class Lulcdata {

  public:

     Lulcdata(void);

/* **************************************************************
                      Public Functions
************************************************************** */

// read data structure.
     int get(ifstream& infile);
     int getdel(FILE* infile);
//write data structure.
     void out(ofstream& ofile, float col, float row, char varname[], int carea, int year, int agstate, double pctag, double RAP, char contnent[9]);
     void outdel(ofstream& ofile, float col, float row, char varname[], int carea, int year, int agstate, double pctag, double RAP, char contnent[9]);


/* **************************************************************
                     Public Variables
************************************************************** */

     float col;          // column or longitude of grid cell (degrees)
     float row;          // row or latitude of grid cell (degrees)
     char varname[9];    // "AGRICUL"
     int carea;          // area of a grid cell
     int year;
     int agstate;        // flag whether or not grid cell is cultivated
     double pctag;       // percent of grid cell covered with agriculture
     double RAP;        // relative agricultural production
     char contnent[9];   // name of continent containing grid cell


  private:

/* **************************************************************
                      Private Variables
************************************************************** */

     int lulcend;
     long curpos;
     long lagpos;

};

#endif

