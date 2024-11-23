/* **************************************************************
*****************************************************************
TBIOME423.CPP - object describing general characteristics of
                vegetation mosaic used in the Terrestrial
	        Ecosystem Model (TEM)

              - modified from Veg4 by DWK 20000102
*****************************************************************
************************************************************** */

#if !defined(TBIOME423_H)
  #include "tbiome423.hpp"
#endif

/* **************************************************************
************************************************************** */

Biome::Biome(void) { temveg = -99; }

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Biome::getvtype(ofstream& rflog1)
{

  char ecd[80];

  cout << endl << "Enter name of the file prescribing vegetation"
	   << " mosaics (.ECD):" << endl;
//cin >> ecd;
  fpara >> ecd;

  rflog1 << endl << "Enter name of the file prescribing vegetation"
         << " mosaics (.ECD):" << endl << ecd << endl << endl;

  getvtype(ecd);

};

/* **************************************************************
************************************************************** */


/* **************************************************************
************************************************************** */

void Biome::getvtype(char ecd[80])
{

  const int NUMVAR = 14;
  char dummy[NUMVAR][8];
  ifstream infile;
  int i;
  int j;

  int subvegid[NUMVEG];
  char subvgname[NUMVEG][31];
  int update[NUMVEG];

  infile.open(ecd, ios::in);

  if (!infile) {
    cerr << "\nCannot open " << ecd << " for data input" << endl;
    exit(-1);
  }

  for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
  for (i = 0; i < NUMVEG; i++)
  {
    infile >> subvegid[i] >> subvgname[i];
    infile >> numtype[i];
    for (j = 0; j < 5; j++) { infile >> subtype[i][j] >> pcttype[i][j]; }
    infile >> update[i];
  }

  infile.close();

};

