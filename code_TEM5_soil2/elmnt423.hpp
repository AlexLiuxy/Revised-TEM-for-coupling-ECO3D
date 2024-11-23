/* *************************************************************
ELMNT423.HPP - Contains functions to manage elements used in GIS
************************************************************** */

#if !defined(ELMNT423_H)
#define ELMNT423_H

class Elmnt {
   
   public:

     Elmnt(void);

/* **************************************************************
                 Public Functions
************************************************************** */

     void ask(ofstream& rflog1);
	 int coregerr(ofstream& rflog1, char varname1[9], const float& col1, \
		 const float& row1, char varname2[9], const float& col2, \
		 const float& row2);
	 void show(ofstream& rflog1, const float& col, const float& row);
	 void show(ofstream& rflog1, const float& col, const float& row, \
		  const long& totyr, const double& tol);

/* **************************************************************
                 Public Variables
************************************************************** */
         
     int end;
     long count;
     long grdcnt;
     long numskip;
     int strtflag;
     int stopflag;
     int  totyr;

     float col;
     float row;
     int carea;
     double elev;
     char contnent[9];

   private:

/* **************************************************************
                 Private Variables
************************************************************** */

     int endflag;
     long numgrids;
};

#endif

