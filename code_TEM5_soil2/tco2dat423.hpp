/* **************************************************************
TCO2DAT423.HPP - object to read and write the structure of transient
	       CO2 data from/to files used by the
	       Terrestrial Ecosystem Model (TEM)

               20000102 - DWK added compiler directives
************************************************************** */

#if !defined(TCO2DAT423_H)
#define TCO2DAT423_H

class CO2data {
   
  public:
      
     CO2data(void);
     
/* **************************************************************
		    Public Functions
************************************************************** */
     
     int get(ifstream& ifile);
     int getdel(FILE* infile);
     void out(ofstream& ofile, float& year, double& bco2);
     void outdel(ofstream& ofile, float& year, double& bco2);

/* **************************************************************
		     Public Variables
************************************************************** */
	  
     float year;        // year at beginning of year
     double mco2;        // atmospheric CO2 concentration in July (ppmv)


  private:

/* **************************************************************
		      Private Variables
************************************************************** */

     int co2end;
     long curpos;
     long lagpos;


};

#endif

