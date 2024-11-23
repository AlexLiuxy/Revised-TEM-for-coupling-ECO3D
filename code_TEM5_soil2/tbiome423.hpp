/* **************************************************************
*****************************************************************
Tbiome423.hpp - object describing general characteristics of
                vegetation mosaic used in the Terrestrial Ecosystem
                Model (TEM)

              - modified from Veg4.hpp by DWK 20000102
*****************************************************************
************************************************************** */

// Global constants

const int NUMVEG = 35;
const int MAXCMNT = 13;
const int NUMMSAC = 5;

#if !defined(TBIOME423_H)
#define TBIOME423_H
class Biome {

  public:

     Biome(void);

 /* **************************************************************
		 Public Functions
************************************************************** */

     void   getvtype(ofstream& rflog1);
     void   getvtype(char ecd[80]);

/* **************************************************************
		 Public Variables
************************************************************** */

     // biome type or ecozone (categorical data)
     int temveg;

     // vegetation community type (categorical data)
     int cmnt;

     // number of community types in a vegetation mosaic
     int numtype[NUMVEG];

     // community types of a vegetation mosaic
     int subtype[NUMVEG][5];

     // percent coverage of a community type in a vegetation
     //   mosaic
     double pcttype[NUMVEG][5];

     //Description of a vegetation community type

     // string changed to char[80] by DWK on 20000210
     char cmnt_name[80];
};

#endif

