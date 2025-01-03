#if !defined(QSOILTEMP1_H)
  #include "qsoiltemp1.hpp"
#endif
  
Qsoiltemp1 basicf;

 Qsoiltemp1::Qsoiltemp1()
 {
	 int i;
	 for(i = 0; i < 5; i ++)
		 tsoil20dep[i] = 0.0 + double(i) * 0.05;

    return;
 };

 Qsoiltemp1::~Qsoiltemp1()
 {
 };
 

/* int Qsoiltemp1::Soiltemp_(double& water2, double& calcu_cdsnow9, double& airt199, double& airt299,
				  double& airt399, double& hsnow199, double& hsnow299, double& hsnow399,
				  double weight99[], double& smass99, int_ln& is99, int_ln& ism199, double x99[],
				  double dx99[], double xfb99[], double xfa99[], double water99[], double t99[], double& tsoil9,
				  double& frontd9, double& thawbegin9, double& thawend9,double diffsoilt9[10],
				  const int& cmnt)*/
int Qsoiltemp1::Soiltemp_(double& water2, double& calcu_cdsnow9, double& airt199, double& airt299,
				  double& airt399, double& hsnow199, double& hsnow299, double& hsnow399,
				  double& smass99, int_ln& is99, int_ln& ism199,double& tsoil9,				  
				  double& frontd9, double& thawbegin9, double& thawend9,double diffsoilt9[10],
				  const int& cmnt, double& soilm)		// XL add for soilmoisture
 {
//------------------------------------------------
// History:
// Feb. 22, 2007
//------------------------------------------------ 	  
// Usage:
//
//------------------------------------------------
// subprocedures called:
//------------------------------------------------
// input variables
//
//------------------------------------------------
// return value:
//--------------------------------------------------

  msg_debug("Running soil thermal model in soiltemp_");
  
  double taer[NNMAX0], tber[NNMAX0], snow[NNMAX0], sden[NNMAX0];
  double heat[2 * NNMAX0], heatt[NNMAX0];
  double thalt[210], adif[210], tanee[210], thi[210], tlo[210];
  double weight[210];
  double deptem[20], tiso[10], depflx[20];
  double phase_count;
  int kint[2];
  int_ln ldepf;
  int_ln maxm = NNMAX0;
  int_ln jjj = 98;

  double sflx1, sflx2, sflx3;
  double comp, eta0, denfac, fcmelt;
 
  

//  cout <<" snow[0] = "<<snow[0]<<" sden[0] = "<<sden[0]<<endl;  

  
/*  --t99;
  --water99;
  --xfa99;
  --xfb99;
  --dx99;
  --x99;
  --weight99;*/
  
  blk_cdsnow = calcu_cdsnow9;
  airt1 = airt199;
  airt2 = airt299;
  airt3 = airt399;
  hsnow1 = hsnow199;
  hsnow2 = hsnow299;
  hsnow3 = hsnow399;
  
  soilwater1 = water2;
  soilwater2 = water2;
//  cout <<"water2 = "<<water<<endl;
  
  
  
/*  for(int_ln i = 1; i <= 210; i ++)
  {
    x[i - 1] = x99[i];
    dx[i - 1] = dx99[i];
    xfa[i - 1] = xfa99[i];
    xfb[i - 1] = xfb99[i];
    water[i - 1] = water99[i];
    t[i - 1] = t99[i];
    	
  } 
  for(int_ln i = 1; i <= 10; i ++)
   weight[i - 1] = weight99[i];*/
   
   nanmax = MAX[cmnt]; 
   nst = NST[cmnt];
   kallyr = KALLYR[cmnt];
   kiso = KISO[cmnt];
   knodes = KNODES[cmnt];
   ktemp = KTEMP[cmnt];
   ksnow = KSNOW[cmnt];
   kenvl = KENVL[cmnt];
   kflux = KFLUX[cmnt];
   
   my_nanmax = nanmax;
   my_kallyr = kallyr;
   my_knodes = knodes;
   my_kiso = kiso;
   my_ktemp = ktemp;
   my_ksnow = ksnow;
   my_kenvl = kenvl;
   my_kflux = kflux;
   my_nst = nst;
   
  int_ln kall = kallyr; 
  liso = LISO[cmnt];
  my_liso = liso;

  for(int_ln i = 1; i <= liso; i ++)
  {
    tiso[i - 1] =TISO[cmnt];
    my_tiso[i - 1] = tiso[i - 1];	
  }

  
  lmax = LMAX[cmnt];
  vdep = VDEPTH[cmnt];
  my_lmax = lmax;
  my_vdep = vdep;
  
  deptem[0] = 0.05;//5cm
  deptem[1] = 0.1;//10cm
  deptem[2] = 0.2;//20cm
  deptem[3] = 0.5;//50cm
  deptem[4] = 1.0;//100cm
  deptem[5] = 2.0;//200cm
  
  for(int_ln i = 1; i <= lmax; i ++)
  {
  	my_deptem[i - 1] = deptem[i - 1];
  	deptem[i - 1] *= vdep;
  	  	
  }

  
  ldepf = NDEPF[cmnt];
  vdep = VDEP[cmnt];
  my_ldepf = ldepf; 
  my_vdep1 = vdep;
  
  
  depflx[0] = DEPFLX1[cmnt];
  depflx[1] = DEPFLX2[cmnt];
  depflx[2] = DEPFLX3[cmnt];
  depflx[3] = DEPFLX4[cmnt];
  depflx[4] = DEPFLX5[cmnt];
  
  for(int_ln i = 1; i <= ldepf; i ++)
  {
  	my_depflx[i - 1] = depflx[i - 1];
  	depflx[i - 1] *= vdep;
//  	cout <<" "<<depflx[i - 1];
  }

  time = 0.0;
  sigma = 0.0;
  double rad = 0.0;
  double turb = 0.0;
  msg_debug("read hlat, tf");  

  hlat = HLAT[cmnt];
  tf = TF[cmnt];
  gflux = GFLUX[cmnt];
  
  my_hlat = hlat;
  my_tf = tf;
  my_gflux = gflux;
  

  msg_debug("call function Grid");
  
  msg_debug("Entering function Grid from Soiltemp_");
  Grid(maxm, cmnt, soilm);		// XL add for soilmoisture
//  for(int i0 = 0; i0 < 210; i0 ++)
//	  cout <<"water["<<i0<<"] = "<<water[i0]<<endl;
//  cin.get();
  msg_debug("Exit function Grid to Soiltemp_");
//  cout <<"after Grid mmax = "<<mmax<<endl;
//
  int_ln imaxp = imax + 1;
  double dep[18];
  dep[0] = x[9];
  dep[1] = x[11];
  dep[9] = x[24];
  dep[10] = x[59];
  dep[11] = x[64];
  dep[12] = x[69];
  dep[13] = x[74];
  dep[14] = x[79];
  dep[15] = x[84];
  dep[16] = x[89];
  dep[17] = x[94];
  
  
  for(int i = 2; i < 9; i ++)  
      dep[i] = x[12 + i];

// for(int i = 0; i < 18; i ++)  
//  cout <<"dep["<<i<<"]="<<dep[i]<<endl;  
//  cin.get();
      
// Bctop sets surface B.C. Type
// Reads data as appropriate   
  msg_debug("Entering function Bctop from Soiltemp_");  
  Bctop(maxm, taer, ktop, sflx1, sflx2, sflx3);
//  for(int i = 0; i < mmax; i ++)
//	  cout <<"taer["<<i<<"] = "<<taer[i]<<endl;
//  cin.get();
  msg_debug("Exit function Bctop to Soiltemp_");

// call Sethet

  msg_debug("Entering function Sethet from Soiltemp_");
  Sethet(maxm, kint, heat, heatt, cmnt);
  msg_debug("Exit function Sethet to Soiltemp_");
   

// call Bcbot
  msg_debug("Entering function Bcbot from Soiltemp_");
  Bcbot(maxm, tber);
  msg_debug("Exit function Bcbot to Soiltemp_");

// call Snofal
  msg_debug("Entering function Snofal from Soiltemp_");
  Snofal(snow, sden, comp, eta0, denfac, fcmelt, cmnt);
  msg_debug("Exit function Snofal to  Soiltemp_");
// call tinitl

  msg_debug("Entering function Tinitl from Soiltemp_");
  Tinitl(cmnt);
  msg_debug("Exit function Tinitl to Soiltemp_");

// read cdsnow

// Initialization
   is = is99;  
   ism1 = ism199;
   smass = smass99;
   double dxx, tu, td, g; 
   if(kswitch == 0)     	     
   {		
   	is = ig;
   	ism1 = igm1;
//   	cout <<"is = "<<is <<" ism1 = "<<ism1<<endl;
//   	cin.get();
   	
   	smass = 0.0;
   	
   for(int_ln i = 0; i < igm1; i ++)
   {
   	 xfa[i] = -1.0e10;
   	 xfb[i] = 0.0;
   	 water[i] = 1.0e0;
   	 dx[i] = 0.0e0;
   	 thalt[i] = -999.9; //attention
   	 weight[i] = 0.0;
   	 ht[i] = 0.0;
   	 htold[i] = 0.0;   	 
   }
   weight[ig - 1] = 0.0;  

// set xfa value for initial phase boundary
   dx[imax - 1] = 0.0; 
     
   for(int_ln i = ig; i <= imax; i ++)
   {
	   thalt[i - 1] = -999.9;
	   htold[i - 1] = 0.0;
	   ht[i - 1] = 0.0;
	   dxx = dx[i - 1];
	   tu = t[i - 1] - tf;
	   if( i < imax ) td = t[i] - tf;
	   xfa[i - 1] = -1.0e10;
	   xfb[i - 1] = 0.0;
	   if(tu * td >= 0.0) continue;
	   double dbl_ = tu - td;
	   if(dbl_ == 0.0) dbl_ = 0.01; //to protect
       g = dxx * tu / dbl_;
       if(hlat != 0.0) xfa[i - 1] = g;		
   } 

  }

// calculate air freezing & thawing index
// calcualte maximum & minimum modified air temperatures
   double frza = 0.0;
   double thwa = 0.0;
   double tairlo = 1.0e3;
   double tairhi = -1.0e3;
   for(int_ln i = 0; i < mmax; i++)
   {
   	 tair = taer[i];
   	 if(tair >= tairhi) tairhi = tair;
   	 if(tair <= tairlo) tairlo = tair;
   	 if(tair < 0)
   	 {
   	 	 frza -= tair;
   	 	 continue;
   	 }

   	 thwa += tair;		   	 		   	 		
   }
   thwa *= dtday;
   frza *= dtday;


   msg_debug("==================begin annual cycle # nan===================");
   nan = 0;   
   nanmax -= 1;   //nanmax is read from outside files

   while(1)   
   {
     nan ++;
     if(nan > nanmax) kallyr = 1;
     double tslo = 1.0e3;	
     double tshi = -1.0e3;
     double tgrlo = 1.0e3;
     double tgrhi = -1.0e3;
     double thws = 0.0;
     double frzs = 0.0;
     double thwg = 0.0;
     double frzg = 0.0;
     double tairan = 0.0;
     double tgran = 0.0;
     int_ln nsum = 0;
     for(int_ln i = ig - 1; i < imax; i ++)
     {
   	   thi[i] = - 1.0e3;
   	   tlo[i] = 1.0e3;
   	   tanee[i] = 0.0;
     }
//     cout <<"================ begin time steping ====================="<<endl;
     msg_debug("================ begin time steping =====================");     
     int_ln kkk = int_ln(first);
//     cout <<"kkk = "<<kkk<<endl;
//     cin.get();

     time = 0.0;
     int_ln m = 0;

     double calend, acum, topold, topnew, tdrive, snoden;
     int_ln isp1;
//     cout <<" mmax = "<<mmax<<endl;
 
     while(m < mmax)
     {
//   	    cout <<"dtday = "<<dtday<<" mmax = "<<mmax<<endl;
//		cin.get();
		 time += dtday;   //elapsed time in days at end of step

        m ++;
        calend = first + time;
        if(calend > per) 
        {

          calend -= per;
        }


    
   //evaluate new air temperature (tair)
 
   	    tair = taer[m - 1];
   	    tdrive = tair;
   	
//  evaluate new surface conditions
//		for(int ii = 0; ii < imax; ii ++)
//			cout <<"xfa["<<ii<<"] = "<<xfa[ii]<<endl; 
        msg_debug("Entering function Surf from Soiltemp_");
        Surf(ktop, sflx1, sflx2, sflx3, tdrive);   	
        msg_debug("Exit function Surf to Soiltemp_");
           
        dtfaz = dt;
        
//  snow cover grid
        acum = snow[m - 1];
        snoden = sden[m - 1];  
//		for(int ii = 0; ii < imax; ii ++)
//			cout <<"xfa["<<ii<<"] = "<<xfa[ii]<<endl;
   
        msg_debug("Entering function Neige from Soiltemp_");
		Neige(acum, snoden, weight, comp, eta0, denfac, fcmelt, tdrive); 
        msg_debug("Exit function Neige to Soiltemp_");
          
//		cin.get();
        topold = t[is - 1];
           
        ism1 = is - 1;
        isp1 = is + 1;
        t[ism1 - 1] = tdrive;
        e[ism1 - 1] = tdrive;
        s[ism1 - 1] = 0.0;
    
//  update internal heat sources
        for(int iheat = 1; iheat <= 2; iheat ++)
          heatt[iheat - 1] = heat[(iheat - 1) * NNMAX0 + m - 1];
        
        msg_debug("Entering functio Inthet from Soiltemp_");    
        Inthet(kint, heatt);    
        msg_debug("Exit function Inthet to Soiltemp_");
     
//  evaluate new bottom conditions     
        tbot = tber[m - 1];
        t[imaxp - 1] = tbot;
    
//  evaluate thermal properties
        double tt = t[is -1];
		double sph, cnd;

        for(int_ln i = is - 1; i < imax1; i ++)
        {
    	     dxx = dx[i];
    	     td = t[i + 1];
    	     
    	     msg_debug("Entering function Pram from Soiltemp_");
    	     Pram(mater[i], tt, td, ddry[i], water[i], sph, cnd);
//			 cout <<"i = "<<i <<" cnd = "<<cnd <<" sph = "<<sph<<endl;
    	     msg_debug("Exit function Pram to Soiltemp_");
    	     
		     tt = td;
             if(dxx == 0.0) dxx = 0.01; //to protect

    	     conx[i] = cnd / dxx;
    	     capx[i] = ddry[i] * sph * dxx;
        }
    
        conx[ism1 - 1] = htop;
        capx[ism1 - 1] = 0.0;
     
        int_ln itzero, jj, ibzero;
        if(hlat > 0.0) 
        {

    	      for(int_ln i = is; i <= imax; i ++)
    	      {
    		      itzero = i;
//				  cout <<"xfa["<<i - 1<<"] = "<<xfa[i - 1]<<endl;

    		      if(xfa[i - 1] >= 0.0)
    		      {

    		      	break;
    		      }
    	      }
    	
    	      jj = ism1 + imax1;

	          for(int_ln j = ism1; j <= imax1; j ++)
	          {
    		      int_ln i = jj - j;
    		      ibzero = i;
    		      if(xfa[i - 1] >= 0.0) 
    		      {

    		      	break;
    		      }
    	      }	      
	          int int_= sign_t(ibzero - itzero);

//			  cout <<"ibzero = "<<ibzero<<" itzero = "<<itzero<<endl;
//			  cin.get();
            int l60 = 0;
	          if(int_ == -1)
	          {
//33	    	 	
                 l60 = 0;
                 msg_debug("Entering function Asmblo from Soiltemp_");
	  	         Asmblo(is);
	  	         msg_debug("Exit function Asmblo to Soiltemp_");
	  	 
	   	         topnew = tdrive * s[is - 1] + e[is - 1];
	   	         if((topnew - tf) * (topold - tf) < 0.0)
	   	         {
	   	  // start of single phase plane at surface

	   	            itzero = is;
	   	            ibzero = is;
	   	            xfa[itzero - 1] = 0.0;
	   	            if(per > 0.0)
	   	            {
						double dbl_ = topold - topnew;
	   		            if(dbl_ == 0.0) dbl_ = 0.01; //to protect
						 dtfaz = dt * (topold - tf) / dbl_;
	   		             int_ln itz1 = itzero + 1;
	   		             t[itzero - 1] = tf;
	   		  
	   		             msg_debug("Entering function Asmblo from Soiltemp_");
	   		             Asmblo(itz1);
	   		             msg_debug("Exit function Asmblo to Soiltemp_");
	   		  
	   		             msg_debug("Entering function Tblo from Soiltemp_");
	   		             Tblo(itz1);
	   		             msg_debug("Exit function Tblo to Soiltemp_");
	   		  
	   		             topold = tf + sign_d(1.0e-4, topnew - tf);
	   		             t[itzero - 1] = topold;
	   		             dtfaz = dt - dtfaz;	    					  
	   	            }
	    	        else
	    	        {
	    		  //Transient problem only
	    		         topold = topnew;
	    		         t[itzero - 1] = topold;
	    			  	
	    	        }		
	    	     }
	    	     else
	    	     {

	    	     	     msg_debug("Entering function Tblo from Soiltemp_");
	    		         Tblo(is);
	    		         msg_debug("Exit function Tblo to Soiltemp_");
	    		
	    		         l60 = 1;
	    	     }	
	    	    }
// 34	    
 		    
	    	    if(!l60)
	    	    {
	    	       if(int_ != 1)
	    	       {

	    	          if(xfb[itzero - 1] <= 0.0)	    		    		    		    		    		    	
	    	          {

	    		           msg_debug("Entering function Asmone from Soiltemp_");
				           Asmone(itzero);
	    		           msg_debug("Exit function Asmone to Soiltemp_");
	    		 
						   topnew = t[is - 1];


	    		           if((topnew - tf) * (topold - tf) >= 0.0) 
	    		           {

	    		           	 l60 = 1;
	    		           	
	    		           }
	    		           else
	    		           {

					           double tset = tf - sign_d(1.0e-4, topnew - tf);
	    		               for(int i = is - 1; i < itzero; i ++)
	    		               {
					             e[i] = 0.0;
					             s[i] = 0.0;
	    		  	             t[i] = tset;
	    		               }	
							   double dbl_ = topnew - topold;
	    		               if(dbl_ == 0.0) dbl_ = 0.01; //to protect
							   dtfaz = dt * (topnew - tf) / dbl_;
					           ibzero = itzero;

//   To ensure 2nd call to asmtwo igores lower phase plane
//   In case previously and multiple phase planes
 				               if(ibzero > itzero) ibzero = 1;

						       itzero = is;

				               topold = tf + sign_d(1.0e-4, topnew - tf);

					           t[is - 1] = topold;

						       if(xfa[is - 1] > 0.0) xfb[is - 1] = xfa[is - 1];
						       xfa[is - 1] = 0.0;


	    		           }
	
	    		        }		
	             }	   
// More than one phase plane, 
// 35	    	
//		for(int ii = 0; ii < imax; ii ++)
//			cout <<"xfa["<<ii<<"] = "<<xfa[ii]<<endl; 
//		cin.get(); 

    
	             while(!l60)
	             {	 

		              msg_debug("Entering function Asmtwo from Soiltemp_");	  	
//		              cout <<"itzero ="<<itzero<<endl;    	      	    		    
	    	          Asmtwo(itzero, ibzero);
	    	          msg_debug("Exit function Asmtwo to Soiltemp_");
	    		     
	    	          topnew = t[is - 1];
//	    	          cout <<"topnew ="<<topnew<<" is = "<<is<<endl;
//	    	          cout <<"topold ="<<topold<<" tf = "<<tf<<endl;
	    	          if((topnew - tf) * (topold - tf) >= 0.0)
	    	          {

	    	            l60 = 1;
//	    	            cout <<"l60 ="<<l60<<endl;
	    	          }
	    	          else
	    	          {
	    	         	
				            double tset = tf - sign_d(1.0e-4, topnew - tf);
//				            cout <<"tset = "<<tset<<endl;
//				            cout <<"itzero = "<<itzero<<endl;
				            for(int i = is - 1; i < itzero; i ++)
				            {
//				            	 cout <<"i = "<<i<<endl;
				               e[i] = 0.0;
				               s[i] = 0.0;
	    	 	               t[i] = tset;
	    	                }
//	    	                cout <<"loop ok"<<endl;	
							double dbl_ = topnew - topold;
	    	                if(dbl_ == 0.0) dbl_ = 0.01; //to protect
							dtfaz = dt * (topnew - tf) / dbl_;
	    	                ibzero = itzero;
       msg_debug("To ensure 2nd call to asmtwo igores lower phase plane");
//   In case previously and multiple phase planes
                            if(ibzero > itzero) ibzero = 1;
                            itzero = is;
                            topold = tf + sign_d(1.0e-4, topnew - tf);	             	    
                            t[is - 1] = topold;
                            if(xfa[is - 1] > 0.0) xfb[is - 1] = xfa[is - 1];
                            xfa[is - 1] = 0.0;
	    	          }

	    	       }

//		for(int ii = 0; ii < imax; ii ++)
//			cout <<"xfa["<<ii<<"] = "<<xfa[ii]<<endl;
	    	    }	
   			  	    	  	
					             	    	    	
    	
        }	
        else
        {

       	    msg_debug("Entering function Asmblo from Soiltemp_");
    	    Asmblo(is);
    	    msg_debug("Exit function Asmblo to Soiltemp_");
    	  
    	    msg_debug("Entering function Tblo from Soiltemp_");    	  
    	    Tblo(is);    	
    	    msg_debug("Exit function Tblo to Soiltemp_");
 

         }

    
// Screen output monitoring    
//60 

   
       msg_debug(" calculate surface freezing and thawing indices");
        double tsurf = t[is - 1]; 

        if(tsurf < 0)
        {
    	    frzs -= tsurf;
        }
        else
        {
    	    thws += tsurf;
        }
      msg_debug("Calculate ground freezing and thawing indices");
        double tgrnd = t[ig - 1];

        if(tgrnd < 0)
        {
    	    frzg -= tgrnd;
        }	
        else
        {
    	    thwg += tgrnd;
        }
       	
//  Calculate maximum and minimum temperatures
      msg_debug("Calculate mean annual ground temperatures");

        if(tsurf >= tshi) tshi = tsurf;
        if(tsurf <= tslo) tslo = tsurf;
        if(tgrnd >= tgrhi) tgrhi = tgrnd;
        if(tgrnd <= tgrlo) tgrlo = tgrnd;
        if(mmax == 0) mmax =5; //to protect        
        tairan += taer[m - 1] / double(mmax);
        tgran += tgrnd / double(mmax);

        for(int_ln i = ig - 1; i < imax; i ++)
        {
    	    if(t[i] >= thi[i]) thi[i] = t[i];
    	    if(t[i] <= tlo[i]) tlo[i] = t[i];
    	    tanee[i] += t[i] / double(mmax);	
//			cout <<"tanee["<<i<<"] = "<<tanee[i]<<endl;

        }  
//  output        
//        cout <<"kallyr = "<<kallyr <<endl;
//        cout <<"nst = "<<nst <<endl;
//        cin.get();

        if(kallyr != 0)
        {
    	    nsum ++;
//			cout <<"nsum = "<<nsum<<endl;

    	    if(nsum >= nst)
    	    {

    	    	nsum = 0;
    	    	kkk ++;
//    	    	cout <<kkk<<" ";
//                for(int i = 12; i < 20; i ++)
//                  cout <<setw(7)<<setprecision(2)<<t[i];
                
//    	    	cin.get();



//         Calculae depth of isotherms
//                msg_debug("Entering function Lociso from Soiltemp_");
//                Lociso(tiso, liso, kiso);
//                msg_debug("Exit function Lociso to Soiltemp_");
//         Print air, surface & interpolated ground temperatures                        	 
//                msg_debug("Entering function Temdep from Soiltemp_");
//                Temdep(deptem, lmax, ktemp); 	    	
//                msg_debug("Exit function Temdep to Soiltemp_");
//         Compute heat flux at selected levels xflux[lflux]
//                msg_debug("Entering function Flxdep from Soiltemp_");
//                Flxdep(depflx, ldepf, kflux);            
//                msg_debug("Exit function Flxdep to Soiltemp_");
//         pirnt snowcover information
//                msg_debug("Entering function Schnee from Soiltemp_");
//                Schnee(acum, snoden, ksnow);             
//                msg_debug("Exit function Schnee to Soiltemp_");
            
    	    }	
   	    
        }     	                   
  				     		  	
     }


 
	 //calculate temperature at specific depths
	 Temdep(deptem, 6, 1, diffsoilt);

	 Temdep(tsoil20dep, 5, 1, tsoil20);
	 tsoil9 = 0.0;
	 for(int i = 0; i < 5; i ++)
		 tsoil9 += tsoil20[i];
	 tsoil9 /= 5.0;
//--------------------------------------------------------   	
//  annual cycle nan complete   
     thws *= dtday;
     frzs *= dtday;
     thwg *= dtday;
     frzg *= dtday;

     double tdif = tgran - tairan;	
    
/*     if(kenvl != 0) 
     {
//    Print out is for all years    	
	
        if(kallyr != 0 || nan >= nanmax)
        {

        }	 
     }*/
 
     if(nan > nanmax)
     {
//    Print out and stop
       int first_phase = 0;
	   double dratio;
       for(int_ln i = 10; i < 37; i ++)
       {

       	 if(tanee[i] < tf)
       	 {
			 dratio = (tanee[i - 1] - tf) / (tanee[i - 1] - tanee[i]);
			 if(dratio <= 1.0)
       			frontd = (- 1.0) * (x[i - 1] + (x[i] - x[i - 1]) * dratio);	
			 else
				 frontd = 0.0;
       	 
       	   first_phase = i;
       	   if(frontd > 0 || fabs(frontd) < 0.05) frontd = 0.0;
       	   break;
       	 }
       	 
       }
       int second_phase = 0;
       thawbegin = 0.0;
       if(first_phase != 0)
       {
       	  for(int_ln i = first_phase + 1; i < 37; i ++)
       	  {
       	  	if(tanee[i] > tf)
       	  	{
       	  	  thawbegin = x[i];
       	  	  second_phase = i;
       	  	  break;
       	  	}
       	  }
       
       
          thawend = 0.0;
          if(second_phase != 0)
          {
       	    for(int_ln i = second_phase + 1; i < 37; i ++)
       	    {
       	 	  if(tanee[i] < tf)
       	 	  {
       	 		thawend = x[i];
       	 		break;
       	 	  }
       	    }
          }
       }
/*
       switch(cmnt)
       {
       	case 2:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11])
       	          / (x[9] + x[10] + x[11]);  
       	  diffsoilt9[0] = tanee[10];        	
       	}
       	break;
       	case 3:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11])
       	          / (x[9] + x[10] + x[11]);
       	  diffsoilt9[0] = tanee[10];        	
       	}
       	break;
       	case 4:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11]
       	           + tanee[12] * x[12] + tanee[13] * x[13]) / (x[9] + x[10]
       	           + x[11] + x[12] + x[13]);
       	  diffsoilt9[0] = tanee[11];         	
       	}
       	break;
       	case 5:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11])
       	  	   / (x[9] + x[10] + x[11]);
       	  diffsoilt9[0] = tanee[10];	   
       	}
       	break;
       	case 6:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11])
       	           / (x[9] + x[10] + x[11]);
       	  diffsoilt9[0] = tanee[10];         	
       	}
       	break;
       	case 7:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10])
       	           / (x[9] + x[10]);	
       	  diffsoilt9[0] = tanee[10];         
       	}
       	break;
       	case 8:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11]
       	            + tanee[12] * x[12]) / (x[9] + x[10] + x[11] + x[12]);	
       	  diffsoilt9[0] = tanee[11];          
       	}
       	break;
       	case 9:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11]
       	           + tanee[12] * x[12]) / (x[9] + x[10] + x[11] + x[12]);
       	  diffsoilt9[0] = tanee[10];         	
       	}
       	break;
       	case 10:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11] 
       	          + tanee[12] * x[12]) / (x[9] + x[10] + x[11] + x[12]);
       	  diffsoilt9[0] = tanee[10];        	
       	}
       	break;
       	case 11:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11]) 
       	          / (x[9] + x[10] + x[11]);	
       	  diffsoilt9[0] = tanee[10];        
       	}
       	break;
       	case 12:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11]
       	         + tanee[12] * x[12]) / (x[9] + x[10] + x[11] + x[12]);	
       	  diffsoilt9[0] = tanee[10];       
       	}
       	break;
       	default:
       	{
       	  tsoil9 = (tanee[9] * x[9] + tanee[10] * x[10] + tanee[11] * x[11])
       	          / (x[9] + x[10] + x[11]);	
       	  diffsoilt9[0] = tanee[10];        
        }
       	break;
       }
       
       diffsoilt9[1] = tanee[10]; // depth .12m
       diffsoilt9[2] = tanee[11]; // depth .22m
       diffsoilt9[3] = tanee[12]; // depth .32m
       diffsoilt9[4] = tanee[13]; // depth .42m
       diffsoilt9[5] = tanee[14]; // depth .52m
//     diffsoilt9[6] = tanee[15]; // depth .62m
//     diffsoilt9[7] = tanee[16]; // depth .72m  
*/
       is99 = is;
       ism199 = ism1;
       smass99 = smass;
/*       for(int_ln i = 1; i <= 210; i ++)
       {
       	 x99[i] = x[i - 1];
       	 dx99[i] = dx[i - 1];
       	 xfa99[i] = xfa[i - 1];
       	 xfb99[i] = xfb[i - 1];
       	 water99[i] = water[i - 1];
       	 t99[i] = t[i - 1];
       }
       for(int_ln i = 1; i <= 10; i ++)
       {
       	 weight99[i] = weight[i - 1];
       }*///Removed on April, 18, 2007

       return 0;    	
     }

     double thalto, dif;
     for(int_ln i = ig - 1; i < imax; i ++)
     {
    	thalto = thalt[i];
    	thalt[i] = t[i];
	    dif = fabs(thalt[i]) - fabs(thalto);
	    adif[i] = fabs(dif);
     }
     double adifmx = 0.0;
     for(int_ln i = ig - 1; i < imax; i ++)
     {
    	if(adif[i] > adifmx) adifmx = adif[i];
     }
     if(adifmx <= 5.0e-3)
     {
        nanmax = 0; 	
     }
     else
     {
    	if(nan > nanmax) return 0;
    		
     }	
   
  }
			
 
};				   
				  
				  
 

//============================================
// End of main function
//============================================


int Qsoiltemp1::Lociso(double tiso[], int_ln liso, int_ln kiso)
{
//----------------------------
// History:
// March 22, 2007
//---------------------------
// Usuage:
//  finds position of isotherms tiso(l)
// 0-isotherm got directly from stored xfa[i], xfb[i] values
//  all other isotherms by linear interpolation
//
//-----------------------
// subprocedure called:
//-------------------------
// input variables
//----------------------------
// return value
// integer
//-------------------------
   double p[10];
   
   
   if(kiso != 1) return 0;
   double tmp, tu, td, elap;
   int kont;
   int_ln isx, i;
   for(int_ln l = 0; l < liso; l ++)
   {
     tmp = tiso[l];
     if(tmp == tf)
     {
		   if(hlat != 0.0)
     	 {
//  finds position of moving interface
           for(int_ln k = 0; k < 10; k ++)
           {
             p[k - 1] = 1.0E6;	
           }    	  	
           kont = 1;
           isx = ig + imax1;
           for(int_ln j = ig; j <= imax1; j ++)
           {
             i = isx - j;
             if(xfa[i - 1] == -1.0E10)continue;
             if(xfb[i - 1] > 0.0)
             {
             	p[kont - 1] = x[i - 1] + xfa[i - 1] + xfb[i - 1];
             	kont ++;
             } 	
             p[kont - 1] = x[i - 1] + xfa[i - 1];
             kont ++;
             
           }
     	 }
     	
     }
//   finds position of other than 
//   freezing isotherm by linear interpolation     
     kont = 1;
     for(int_ln k = 0; k < 10; k ++)
     {
     	 p[k] = 1.0E6;
     }
     isx = is + imax1;
     for(int_ln j = is; j <= imax1; j ++)
     {
       i = isx - j;
       tu = t[i - 1] - tmp;
       td = t[i] - tmp;
       if(tu * td < 0.0) 
       {
       	 p[kont - 1] = x[i - 1] + dx[i - 1] * tu / (tu - td);
       	 kont ++;
       }
       
     	  	
     }
     elap = 1.0 * (nan - 1) + time / 365.0;
     
     
     
     	
   }
   return 0;
   
     

   	
};




int Qsoiltemp1::Inthet(int kint[2], double source[])
{
//-----------------------------
// History:
// March 22, 2007
//----------------------------
// Usuage:
// updates internal heat source vector
// returns if kint # 1 ( all sources switched off)
//   elements having value for xhet () other than 999.9 have heat source
//  numbered 0... nheat - 1, in in creasing order of depth
// formulation allows for point sources at depths xhet()
//---------------------------------
// Subprocedure called
//-----------------------------
// input variables
//--------------------------
// return value
// integer
//---------------------------------

// Initialzie ht() & update htold()
   for(int_ln i = is - 1; i < imax; i ++)
   {
   	 htold[i] = ht[i];
   	 ht[i] = 0.0;
   }
// assemble new nodal heating vector
   int_ln nheat = 0, ip1;
   double hti, htip1;
   for(int_ln i = ig; i <= imax1; i ++)   
   {
     ip1 = i + 1;	
     if(xhet[i - 1] == 999.9) continue;
//   element contains a point heat source
     nheat ++;
     if(nheat > 2) return 0;     
     if(kint[nheat - 1] != 1) continue;
     if(t[i - 1] < -10.0 && t[i] < -10.0) source[nheat - 1] = 0.0;
     if(dx[i - 1] == 0.0) dx[i - 1] = 0.01;	
     hti = source[nheat - 1] * (x[ip1 - 1] - xhet[i - 1]) / dx[i - 1];
     htip1 = source[nheat - 1] * (xhet[i - 1] - x[i - 1]) / dx[i - 1];
     ht[i - 1] += hti;
     ht[ip1 - 1] += htip1;
      
     
   }
   return 0; 
  	
};


int Qsoiltemp1::Temdep(double deptem[], int_ln lmax, int_ln ktemp, double diffsoilt[])
{
//-----------------------------------
// History:
// March 22, 2007
// modifed Dec, 22, 2007
//----------------------------------------
// Usuage:
// Evaluates temperature at any depths deptem[l], l = 0, lmax - 1
// Linear interpolation but including phase boundaries
// prints out air, surface and ground temperatures on unit 13
//-----------------------------------------
// subprocedure called:
//---------------------------------------
// input variables
//------------------------------------
// return value
// integer
//-----------------------------------
   double temp[15];
   if(ktemp == 0) return 0;
   int_ln i;
   double xd, xu, tu, td, gu, gd;
   double depxu, dep;
   for(int_ln l = 0; l < lmax; l ++)
   {
     dep = deptem[l];
//  locate nodes i, i + 1 on either side of deptem[l]    	
     i = ism1;
     for(int_ln j = is - 1; j < imax1; j ++)
     {
       i ++;
       xu = x[i - 1];
       xd = x[i];
       if(dep >= xu && dep <= xd)break;
       	
     }
//   evaluate temperature at deptem[l]
     tu = t[i - 1];
     td = t[i];
	   gu = xfa[i - 1];
	   double dxx;
     if(gu != -1.0E10)
     {
// element contains phase planes[s]     	
       dxx = dx[i - 1];
       gd = gu + xfb[i - 1];
       depxu = dep - xu;
       temp[l]= 0.0;      	
       if(depxu < gu) temp[l] = tu * (gu - depxu) / gu;
       if(depxu > gd) temp[l] = td * (depxu - gd) / (dxx -gd);
       
     }
     else
     {
       temp[l] = tu + (td - tu) * (dep - xu) / dx[i - 1];
     }
	 diffsoilt[l] = temp[l];
          
   } 
//  double elap = 1.0 * double(nan - 1) + time / 365.0;

   
   return 0;
   

   
   
   
	
};


int Qsoiltemp1::Flxdep(double depflx[], int_ln& ldepf, int_ln& kflux)
{
//---------------------------
// History:
// March 22, 2007
//------------------------
// Usuage:
// Evaluates fluxes at any depths depflx[k],k = 0, ldepf - 1
//---------------------------
// subprocedure called:
//-------------------------
// Input variables
//-------------------------
// return value
// integer
//---------------------------
   double flux[21];
   if(kflux == 0) return 0;	
   double dep, xu, xd, tu, td;
   int_ln i;
   for(int_ln k = 0; k < ldepf; k ++)
   {
     dep = depflx[k];
//locate nodes i, i + 1, on either side of depflx[k]
     i = ism1;
     for(int_ln j = is - 1; j < imax1; j ++)
     {
     	 i ++;
     	 xu = x[i - 1];
     	 xd = x[i];
     	 if(dep >= xu && dep <= xd)break;
     }
//  
//   evaulate flux for element containing depflx[k]
     tu = t[i - 1];
     td = t[i];
     flux[k] = conx[i - 1] * (tu - td);
               
     	
   }
   double elap = 1.0 * double(nan - 1) + time / 365.0;
//   cout <<nan <<"  "<<time<<"  "<<elap<<"  ";
//   for(int_ln k = 0; k < ldepf; k ++)
//   {
//      cout << flux[k] <<"  ";
//      if(!((k + 1) % 5)) cout <<endl; 	
//   }
   return 0;
   
};


int Qsoiltemp1::Pram(int_ln& mat, double& ta, double& tb, double& dense, double& wat, double& sph, double& cnd)
{
//--------------------------------------------------------------
// History:
// March 22, 2007
//----------------------------------------------------------
// Note:
// all units are S.I.U.
//------------------------------------------
// Usuage:
//------------------------------------------------
// subprocedure called
//-----------------------------------------
// input variables
//-----------------------------------------
// return value
// integer
//----------------------------------------------

  if(mat <= 10) 
  {
    if(ta <= tf) 
    {
    	cnd = condf[mat - 1];
    	sph = sphf[mat - 1];
    	return 0;
    }
    cnd = condt[mat - 1];
    sph = spht[mat - 1];
    return 0;
  }
  
  int_ln mat10;
  double temp = 0.5 * (ta + tb);
  if(mat >= 31)return 0;    	       

  
  
  mat10 = mat - 10;
  if(mat <=20 && time > dtday)
  {
      if(ta <= tf) 
      {
    	  cnd = condf[mat - 1];
    	  sph = sphf[mat - 1];
    	  return 0;
      }
      cnd = condt[mat - 1];
      sph = spht[mat - 1];
      return 0;     	
  }
     
  
  switch(mat10)
  {

	  case 1:
      {
// type # 11
// polystyrene insulation

          condt[mat - 1] = 0.03;
	      spht[mat - 1] = 4186.0 * 0.3;
	      condf[mat - 1] = condt[mat - 1];
          sphf[mat - 1] = spht[mat - 1];       	
      };break;
	  case 2:
      {
//  type # 12    	
//   constant

          condt[mat - 1] = 0.0;
	      if(condt[mat - 1] == 0.0) return 0;
           
      };break;	
	  case 3:
      { 
///  type # 13
//   gravel or sand according to kersten

//   frozen

	      sphf[mat - 1] = 745.1 + 2093.0 * wat;
	      condf[mat - 1] = 1.09E-2 * pow(10.0, 8.12E-4 * dense)
                     + 4.6E-1 * wat * pow(10.0, 9.12E-4 * dense);
//   thawed
	      spht[mat - 1] = 745.1 + 4186.0 * wat;
              if(wat <= 0.0) wat = 0.01;//to protect
          condt[mat - 1] = (0.44E-1 * log(wat) + 0.26) * pow(10.0, 6.25E-4 * dense);                           
                       	
      };break;
	  case 4:
     {
// type # 14
//    clay or silt soil according to kersten
//    frozen
	     sphf[mat - 1] = 745.1 + 2093.0 * wat;
	     condf[mat - 1] = 1.43E-3 * pow(10.0, 1.375E-3 * dense)
                     + 1.22 * wat * pow(10.0, 5.00E-4 * dense);
//    thawed
       spht[mat - 1] = 745.1 + 4186.0 * wat;
	     condt[mat - 1] = (0.565E-1 * log(wat) + 0.2304) * pow(10.0,6.25E-4 * dense);
                      	
    	
     };break;
     case 5:
     {

       
         flog1 <<" **************** Material index "<< mat <<" undefined ***********"<<endl;	
         return 0;    	
     };break;
	 case 6:
	 {		   	    

	     flog1 <<" **************** Material index "<< mat <<" undefined ***********"<<endl;
	     return 0;
	 };break;
	 case 7:
	 {

	    
	     flog1 <<" **************** Material index "<< mat <<" undefined ***********"<<endl;
	     return 0;
	 };break;
	 case 8:
	 {	    
	     flog1 <<" **************** Material index "<< mat <<" undefined ***********"<<endl;
	     return 0;
	 };break;
	 case 9:
	 {	     
	     flog1 <<" **************** Material index "<< mat <<" undefined ***********"<<endl;
	     return 0;
	 };break;
	 case 10:
	 {	     
	     flog1 <<" **************** Material index "<< mat <<" undefined ***********"<<endl;
	     return 0;
	 };break;

//..... Note ..... materials 21-30 reserved for cases where material
//  properties must be re-evaluated at each time-step
	 case 11:
	 {
//  type # 21
//    snow conductivity modified woodside formula

	     sph = 2093.0;
	     cnd = blk_cdsnow;
//  due to pitman & zuckerman
//  accounts for bonds between spherical ice particles
//  vapor diffusion as per formulae of leg (thesis)
//  assumes vapor diffusion coeff  d = 0.22 C.G.S. (in pores)
	     if(ta < 0.0) return 0;
//   thawed
	     sph = 4086.0;
	     cnd = 1.0E1;
	     return 0;
	  };break;
	  case 12:
	  {
//    type # 22
//    snow devaux formula for conductivity
	     cnd = 1.0E1;
	     sph = 4086.0;
	     if(ta >= 0.0) return 0;
//   frozen
	     sph = 2093.0;
	     cnd = 2.9E-2 + 2.9E-6 * dense * dense;
	     return 0;
	  };break;
	  case 13:
	  {
//   type # 23
//  theoretical gravel or sand with variable unfrozen water content
//  thermal conductivity - euken theoretical formula (de vries)
//  assumes air is continous medium
//  unfrozen water content data guessed (based on anderson for clay)
	    double unf = wat;
	    double du = 0.0;
		double dbl_;
	    if(temp < tf)
	    {

	      double um1 = 0.15 * wat;
	      double um4 = 0.01 * wat;
		  dbl_ = wat - um1;
		  if(dbl_ == 0.0) dbl_ = 0.01;  //to protect
	      double exbm1 = (um1 - um4) / dbl_;
		  dbl_ = 1.0 - exbm1;
		  if(dbl_ == 0.0) dbl_ = 0.01;  // to protect
	      double a = (wat - um1) / dbl_;
		  if(exbm1 <= 0.0 ) exbm1 = 0.01; // to protect
	      double b = -log(exbm1);
	      double c = wat - a;
		  if(temp >= 0.0) temp = -0.01;// to protect
	      double sqt = sqrt(-temp);
	      unf = a * exp(-b * sqt) + c;
	      du = 0.5 * b * (unf - c) / sqt;
	      if(unf - c <= 1.0E-5)
	      {

		      unf = c;
		      du = 0.0E0;
	      }
	    }
// apprant heat capacity per unit mass
	    sph = 711.62 + 2093.0 * (wat + unf) + 334.0e3 * du;
	    double rs = 2700.0;
	    double rl = 1000.0;
	    double ri = 917.0;
	    double xs = dense / rs;
	    double xl = dense * unf / rl;
	    double xi = dense * (wat - unf) / ri;
	    double xw = xl + xi;
	    double xwsat = 1.0 - xs;
	    double xa = xwsat - xw;
	    if(xw > xwsat)
	    {
		 
		    xa = 0.0;
		    xs = 1- xw;
	    }
	    double cs = 7.0;
	    double cl = 0.56;
	    double ci = 2.10;
		dbl_ = temp + 273.2;
//		if(dbl_ == 0.0) dbl_ = -0.01; // to protect, usually doesn't happen
	    double ca = 0.0242 + 5.96 * exp(0.082 * temp) / dbl_;
	    double ga = 0.10;
	    double gc = 0.10;
		if(ca == 0.0) ca = 0.01; //to protect
	    double csca = cs / ca;
	    dbl_ = 1.0 + ga * (csca - 1.0);
	    if(dbl_ == 0.0) dbl_ = 0.01;
	    double fs = 0.666667 / dbl_;
	    dbl_ = 1.0 + gc * (csca - 1.0);
	    if(dbl_ == 0.0) dbl_ = 0.01;
	    fs += 0.333333 / dbl_;
	    double cica = ci / ca;
	    dbl_ = 1.0 + ga * (cica - 1.0);
	    if(dbl_ == 0.0) dbl_ = 0.01; 
	    double fi = 0.666667 / dbl_;
	    dbl_ = 1.0 + gc * (cica - 1.0);
	    if(dbl_ == 0.0) dbl_ = 0.01;
	    fi += 0.333333 / dbl_;
	    double clca = cl / ca;
	    dbl_ = 1.0 + ga * (clca - 1.0);
	    if(dbl_ == 0.0) dbl_ = 0.01;
	    double fl = 0.666667 / dbl_;
	    dbl_ = 1.0 + gc * (clca - 1.0);
	    if(dbl_ == 0.0) dbl_ = 0.01;
	    fl += 0.333333 / dbl_;
	    double hs = xs * fs;
	    double hi = xi * fi;
	    double hl = xl * fl;
		dbl_ = xa + hs + hi + hl;
		if(dbl_ == 0.0) dbl_ = 0.01; //to protect
	    cnd = (xa * ca + hs * cs + hi * ci + hl * cl) / dbl_;
	    return 0;



	  };break;
	  case 14:
	  {
//   type # 24
//   theoretical saturated clay soil ... with variable unfrozen water
//   thermal conductivity - kuken theoretical formula (de vries)
//   unfrozen water content curve adapted from D. Anderson
	    double unf = wat;
	    double du = 0.0;
	    double um1, um4, exbm1, a, b, c, sqt;
		double dbl_;
	    if(temp < tf)
	    {
		   um1 = 0.35 * wat;
		   um4 = 0.17 * wat;
		   dbl_ = wat - um1;
		   if(dbl_ == 0.0) dbl_ = 0.01;
		   exbm1 = (um1 - um4) / dbl_;
		   dbl_ = 1.0 - exbm1;
		   if(dbl_ == 0.0) dbl_ = 0.01;
		   a = (wat - um1) / dbl_;
		   if(exbm1 <= 0.0) exbm1 = 0.01; // to protect
		   b = -log(exbm1);
		   c = wat - a;
		   if(temp >= 0.0) temp = -0.01; // to protect
		   sqt = sqrt(-temp);
		   unf = a * exp( -b * sqt) + c;
		   du = 0.5 * b * (unf - c) / sqt;
		   if(unf - c <= 1.0E-5)
		   {
		     unf = c;
		     du = 0.0E0;
		   }

	   }
// apprant heat capacity per unit mass
	   sph = 711.62 + 2093.0 * (wat + unf) + 334.0E3 * du;
	   double rs = 2700.0;
	   double rl = 1000.0;
	   double ri = 917.0;
	   double xs = dense / rs;
	   double xl = dense * unf / rl;
	   double xi = dense * (wat - unf) / ri;
	   double xw = xl + xi;
	   double xwsat = 1.0 - xs;
	   double xa = xwsat - xw;
	   if(xw > xwsat)
	   {
		    xa = 0.0;
		    xs = 1 - xw;

	   }

	   double cs = 2.90;
	   double cl = 0.56;
	   double ci = 2.10;
	   double ca = 0.0;
	   double ga = 0.125;
	   double gc = 0.750;
	   double cscl = cs / cl;
	   dbl_ = 1.0 + ga * (cscl - 1.0);
	   if(dbl_ == 0.0) dbl_ = 0.01;// to protect
	   double fs = 0.666667 / dbl_;
	   dbl_ = 1.0 + gc * (cscl - 1.0); 
	   if(dbl_ == 0.0) dbl_ = 0.01; // to protect
	   fs += 0.333333 / dbl_;
	   double cicl = ci / cl;
	   dbl_ = 1.0 + ga * (cicl - 1.0);
	   if(dbl_ == 0.0) dbl_ = 0.01; // to protect
	   double fi = 0.666667 / dbl_;
	   dbl_ = 1.0 + gc * (cicl - 1.0); 
	   if(dbl_ == 0.0) dbl_ = 0.01; // to protect
	   fi += 0.333333 / dbl_;
	   double fa = 0.0;
	   double hs = xs * fs;
	   double hi = xi * fi;
	   double ha = xa * fa;
	   dbl_ = xl + hs + hi + ha;
	   if(dbl_ == 0.0) dbl_ = 0.01;
	   cnd = (xl * cl + hs * cs + hi * ci + ha * ca) / dbl_;
	   return 0;

	 };break;
	 case 15:
	 {
//   type # 25
//   theoretical peat soil - with variable unfrozen water content
//   unfrozen water content from emprical formula
//   thermal conductivity - Kuken theortical formula (de Vries)
	   double unf = wat;
	   double du = 0.0;
	   double tab = temp + 273.15;
	   double dbl_;
//   vapor disffusion by evaporation
	   double v = 5.40E3;
//	   if(tab == 0.0) tab = 0.01;//to protect, usually not happen
	   double ca = 0.0242 + 0.110295E-2 * (v / tab - 1.0) * exp(v * (1.0 / 273.15 - 1.0 / tab));
	   double um1, um2, exbm1, a, b, c;
	   if(temp < 0.0)
	   {
//    vapour diffusion by sublimation
	      v = 6.12E3;
//	      if(tab == 0.0) tab = 0.01;//to protect, usually not happen
	      ca = 0.0242 + 0.110295E-2 * (v / tab - 1.0) * exp(v * (1.0 / 273.15 - 1.0 / tab));
	      um1 = 0.579980;
	      um2 = 0.557970;
		  dbl_ = wat - um1;
		  if(dbl_ == 0.0) dbl_ = 0.01; // to protect
	      exbm1 = (um1 - um2) / dbl_;
		  dbl_ = 1.0 - exbm1;
		  if(dbl_ == 0.0) dbl_ = 0.01;
	      a = (wat - um1) / dbl_;
		  if(exbm1 <= 0.0) exbm1 = 0.01; // to protect
	      b = -log(exbm1);
	      c = wat - a;
	      unf = a * exp(b * temp) + c;
	      du = b * (unf - c);
//
	 }
//   heat capacity per unit mass
	      sph = 1925.6 + 2093.0 * (wat + unf) + hlat * du;
	      double rs = 1550.0;
	      double rl = 1000.0;
	      double ri = 917.0;
	      double ra = 1.25;
	      double xs = dense /rs;
	      double xl = dense * unf / rl;
	      double xi = dense * (wat - unf) / ri;
	      double xw = xl + xi;
	      double xp = xs + xw;
	      double xa = 1.0 - xp;
	      double xwsat = 1.0 - xs;
	      if(xw > xwsat)
	      {
			  flog1 <<setw(10)<<""<<setfill('?')<<setw(10)<<""<<setfill(' ')
				          <<" MAT #25--IMPOSSIBLE COMBINATION DENSE,WAT"<<endl;
//		           cout << " mat #25 impossible combination dense, wat in Pram"<<endl;
		       
		      return 0;
	      }
	      double cs = 0.25;
	      double cl = 0.56;
	      double ci = 2.10;
	      double fs = 1.25;
	      double fi = 0.53;
	      double hs = xs * fs;
	      double hi = xi * fi;
	      double fl, fa, ha;
	      if(xa < 0.5)
	      {
//   use wat as continous phase with peat, ice, & air inclusions
	         fl = 1.0;
			 dbl_ = cl + ca;
			 if(dbl_ == 0.0) dbl_ = 0.01;
	         fa = (5.0 * cl + ca) / (3.0 * dbl_);
	         ha = xa * fa;
			 dbl_ = xl + hs + hi + ha;
			 if(dbl_ == 0.0) dbl_ = 0.01;
	         cnd = (xl * cl + hs * cs + hi * ci + ha * ca) / dbl_;
	         return 0;
	      }
//   Use air as continous phase with inclusion of " moist peat"
//   "moist peat" -- use continous liquid phase with ice & peat inclusion
		  dbl_ = xl + hs + hi;		  
		  if(dbl_ == 0.0) dbl_ = 0.01;
	      double cp = (xl * cl + hs * cs + hi * ci) / dbl_;
		  dbl_ = ca + cp;
		  if(dbl_ == 0.0) dbl_ = 0.01;
	      double fp = (5.0 * ca + cp) / (3.0 * dbl_);
	      double hp = xp * fp;
		  dbl_ = xa + hp;
		  if(dbl_ == 0.0) dbl_ = 0.01;
	      cnd = (xa * ca + hp * cp) / dbl_;
	      return 0;


	 };break;
	 case 16:
	 {
//  type # 26
//   second sand soil model .. with variable unfrozen water content
//   thermal conductivity - Kuken theoretical formula (de Vries)
//   Assumes continuous water phase
//   Unfrozen water content guessed
	   double unf = wat;
	   double du = 0.0;
	   double um1, um4, exbm1, a, b, c, sqt;
	   double dbl_;
	   if(temp < tf)
	   {
		   um1 = 0.15 * wat;
		   um4 = 0.01 * wat;
		   dbl_ = wat - um1;
		   if(dbl_ == 0.0) dbl_ = 0.01;
		   exbm1 = (um1 - um4) / dbl_;
		   dbl_ = 1.0 - exbm1;
		   if(dbl_ == 0.0) dbl_ = 0.01;
		   a = (wat - um1) / dbl_;
		   if(exbm1 <= 0.0) exbm1 = 0.01;
		   b = -log(exbm1);
		   c = wat - a;
		   if(temp >= 0.0) temp = -0.01;
		   sqt = sqrt(-temp);
		   unf = a * exp(-b * sqt) + c;
		   du = 0.5 * b * (unf - c) / sqt;
		   if(unf - c <= 1.0E-5)
		   {
		     unf = c;
		     du = 0.0E0;
		   }
	   }
//   apparant  heat capacity per unit mass
	   sph = 711.62 + 2093.0 * (wat + unf) + 334.0E3 * du;
	   double rs = 2700.0;
	   double rl = 1000.0;
	   double ri = 917.0;
	   double xs = dense / rs;
	   double xl = dense * unf / rl;
	   double xi = dense * (wat - unf) / ri;
	   double xw = xl + xi;
	   double xwsat = 1.0 - xs;
	   double xa = xwsat - xw;

	   double cs = 7.00;
	   double cl = 0.56;
	   double ci = 2.10;
	   double ca = 0.0242 + 5.96 * exp(0.082 * temp) / (273.2 + temp);

	   double ga = 0.100;
	   double gc = 0.100;
	   double cscl = cs / cl;
	   dbl_ = 1.0 + ga * (cscl - 1.0);
	   if(dbl_ == 0.0) dbl_ = 0.01; // to protect
	   double fs = 0.666667 / dbl_;
	   dbl_ = 1.0 + gc * (cscl - 1.0);
	   if(dbl_ == 0.0) dbl_ = 0.01;
	   fs += 0.333333 / dbl_;
	   double cicl = ci / cl;
	   dbl_ = 1.0 + ga * (cicl - 1.0);
	   if(dbl_ == 0.0) dbl_ = 0.01;
	   double fi = 0.666667 / dbl_;
	   dbl_ = 1.0 + gc * (cicl - 1.0);
	   if(dbl_ == 0.0) dbl_ = 0.01;
	   fi += 0.333333 / dbl_;
	   double cacl = ca / cl;
	   dbl_ = 1.0 + ga * (cacl - 1.0);
	   if(dbl_ == 0.0) dbl_ = 0.01;
	   double fa = 0.666667 / dbl_;
	   dbl_ = 1.0 + gc * (cacl - 1.0);
	   if(dbl_ == 0.0) dbl_ = 0.01;
	   fa += 0.333333 / dbl_;
	   double hs = xs * fs;
	   double hi = xi * fi;
	   double ha = xa * fa;
	   dbl_ = xl + hs + hi + xa;
	   if(dbl_ == 0.0) dbl_ = 0.01;
	   cnd = (xl * cl + hs * cs + hi * ci + ha * ca) / dbl_;
	   return 0;

	 };break;
	 case 17:
	 {
//  type # 27
//  third sand soil type ... with variable unfrozen water content
//  thermal conductivity --Kersten
//  Assumes continuous water phase
//  unfrozen water content guessed
	   double unf = wat;
	   double du = 0.0;
	   double um1, um4, exbm1, a, b, c, sqt;
	   double dbl_;
	   if(temp < tf)
	   {
		    um1 = 0.15 * wat;
		    um4 = 0.01 * wat;
			dbl_ = wat - um1;
			if(dbl_ == 0.0) dbl_ = 0.01;
		    exbm1 = (um1 - um4) / dbl_;
			dbl_ = 1.0 - exbm1;
			if(dbl_ == 0.0) dbl_ = 0.01;
		    a = (wat - um1) / dbl_;
			if(exbm1 <= 0.0) exbm1 = 0.01;
		    b = -log(exbm1);
		    c = wat - a;
			if(temp >= 0.0) temp = -0.01;
		    sqt = sqrt(-temp);
		    unf = a * exp(-b * sqt) + c;
		    du = 0.5 * b * (unf - c) / sqt;
		    if(unf - c <= 1.0E-5)
		    {
		       unf = c;
		       du = 0.0E0;
		    }
	   }
//   apparant heat capacity per unit mass
	   sph = 711.62 + 2093.0 * (wat + unf) + 334.0E3 * du;
	   double rs = 2700.0;
	   double rl = 1000.0;
	   double ri = 917.0;
	   double xs = dense / rs;
	   double xl = dense * unf / rl;
	   double xi = dense * (wat - unf) / ri;
	   double xw = xl + xi;
	   double xwsat = 1.0 - xs;
	   double xa = xwsat - xw;
           if(wat <= 0.0) wat = 0.01;//to protect


	   if(temp < tf)
		   cnd = 1.09E-2 * pow(10.00, 8.12E-4 * dense)
			    + 4.6E-1 * wat * pow(10.00, 9.12E-4 * dense);
	   else
		   cnd = (0.44E-1 * log(wat) + 0.26) * pow(10.00, 6.25E-4 * dense);
	   return 0;





	 };break;
	 case 18:
	 {
	   
	   	 flog1 <<" **************** Material index "<< mat <<" undefined ***********"<<endl;
	     return 0;

	 };break;
	 case 19:
	 {
	   
	   	 flog1 <<" **************** Material index "<< mat <<" undefined ***********"<<endl;
	     return 0;
	 };break;

	 case 20:
     {
     
     	 flog1 <<" **************** Material index "<< mat <<" undefined ***********"<<endl;	
         return 0;     	
     };break;
     default:
     {
     
     	 flog1 <<"Bad assignment to mat in Pram "<<endl;
         return 0;
     };break;    
    
  }
  
  if(ta <= tf) 
  {
    cnd = condf[mat - 1];
    sph = sphf[mat - 1];
    return 0;
  }
  cnd = condt[mat - 1];
  sph = spht[mat - 1];
  return 0;            	  
  
     	
};


int Qsoiltemp1::Grid(int_ln& maxm, const int& cmnt, double& soilm) // XL add for soilmoisture
{
//--------------------------------------------
// History:
// April 3, 2007
//--------------------------------------------
// Usuage:
// Set grid spacing, node depths, material type, density, water	content
// and explicit thermal properties is mat <= 10
// either constant or geometric progression spacing for each layer
// --------------------------------------
// subprocedure called:
// Pram()
//-----------------------------------------
// input variables
// none
// -----------------------------------------
// return value
// integer
//------------------------------------

   double tmpb[20], capb[20], conb[20];
   double top, epsmin, vspace, vdens, total, tgiv, xx, dxx, dtbig, thick, dxa, dxb, dense, wat;   
   double temp, sph, cnd, eps, dxagiv, vcond, vsph, dtmax, thmdxa, err, df; 
   double usedxa, ratio, rwas, rnew, fwas, power, r, f, dxmin;
   int_ln i, mat, layer, nodes, nn;
   char comment[25];

   int_ln lmimax = 210;
   lmimax --;
// data assignment
   if(kswitch == 0)
   {
	   first = FIRST[cmnt];
	   final = FINAL[cmnt];
	   per = PER[cmnt];
	   dtday = DTDAY[cmnt];
	   theta = THETA[cmnt];

	   top = TOP[cmnt];
	   ig = IG[cmnt];
	   epsmin = EPSMIN[cmnt];
	   vspace = VSPACE[cmnt];
	   vdens = VDEN[cmnt];

	   my_first = first;
	   my_final = final;
	   my_per = per;
	   my_dtday = dtday;
	   my_theta = theta;


	   my_ig = ig;
	   my_top = top;
	   my_epsmin = epsmin;
	   my_vspace = vspace;
	   my_vdens = vspace;
   }
   else
   {
	   top = my_top;
	   epsmin = my_epsmin;
	   vspace = my_vspace;
	   vdens = my_vspace;

	   first = my_first;
	   final = my_final;
	   per = my_per;
	   dtday = my_dtday;
	   theta = my_theta;
	   ig = my_ig;
   }


// data assignment

   if(theta <= 0.5) theta = 0.5;
   theta *= 2.0;
   theta1 = 2.0 - theta;
   igm1 = ig - 1;
//   cout <<endl;
//   cout <<"final = "<<final<<endl;
//   cout <<endl;
   if(per <= 0.0 && final <= first)
   {
     return 0;
   }
   if(per > 0.0 && final <= 0.0) final = first + per;
   total = final - first;

   if(dtday == 0.0) dtday = 0.5; //to protect
   mmax = int_ln(ceil(total / dtday));
//   cout <<endl<<" in grid"<<endl;
//   cout <<"total = "<<total<<endl;
//   cout <<"mmax = "<<mmax <<endl;
//   cout <<endl<<endl;

   dtday = total / double(mmax);
   double dtgiv = dtday;
   xx = top * vspace;
   x[ig - 1] = xx;
   i = ig - 1;
   
   layer = 0;
   dtbig = 1.0E10;

   while(1)
   {
	   layer ++;
//	   cout <<"layer = "<<layer<<" vspace = "<<vspace<<endl;
//	   cout <<"mmax = "<<mmax<<endl;
//	   cout <<"kswitch = "<<kswitch<<endl;
	   if(kswitch == 0)
	   {
		   switch(layer)
		   {
		   case 1:
			   {
				   thick = THICK1[cmnt];
				   dxa = DXA1[cmnt];
				   dxb = DXB1[cmnt];
				   mat = MAT1[cmnt];
				   dense = DENSE1[cmnt];
				   wat = soilm; // XL add for soilmoisture
			   }
			   break;
		   case 2:
			   {
				   thick = THICK2[cmnt];
				   dxa = DXA2[cmnt];
				   dxb = DXB2[cmnt];
				   mat = MAT2[cmnt];
				   dense = DENSE2[cmnt];
				   wat = WATER2[cmnt];

			   }
			   break;
		   case 3:
			   {
				   thick = THICK3[cmnt];
				   dxa = DXA3[cmnt];
				   dxb = DXB3[cmnt];
				   mat = MAT3[cmnt];
				   dense = DENSE3[cmnt];
				   wat = WATER3[cmnt];
			   }
			   break;
		   case 4:
			   {
				   thick = THICK4[cmnt];
				   dxa = DXA4[cmnt];
				   dxb = DXB4[cmnt];
				   mat = MAT4[cmnt];
				   dense = DENSE4[cmnt];
				   wat = WATER4[cmnt];

			   }
			   break;
		   case 5:
			   {
				   thick = THICK5[cmnt];
				   dxa = DXA5[cmnt];
				   dxb = DXB5[cmnt];
				   mat = MAT5[cmnt];
				   dense = DENSE5[cmnt];
				   wat = WATER5[cmnt];
			   }
			   break;
		   case 6:
			   {
				   thick = THICK6[cmnt];
				   dxa = DXA6[cmnt];
				   dxb = DXB6[cmnt];
				   mat = MAT6[cmnt];
				   dense = DENSE6[cmnt];
				   wat = WATER6[cmnt];
			   }
		   
			   break;
		   default:		   	 		   				   		
		   				   	
			   break;
		   }
		   my_thick[layer - 1] = thick;
		   my_dxa[layer - 1] = dxa;
		   my_dxb[layer - 1] = dxb;
		   my_mat[layer - 1] = mat;
		   my_dense[layer - 1] = dense;
		   my_wat[layer - 1] = wat;
	   }
	   else
	   {
		   thick = my_thick[layer - 1];
		   dxa = my_dxa[layer - 1];
		   dxb = my_dxb[layer - 1];
		   mat = my_mat[layer - 1];
		   dense = my_dense[layer - 1];
/*		   if(layer == 3)
			   wat = soilwater2;
		   else
			   wat = my_wat[layer - 1];*/
		     wat = my_wat[layer - 1];		   		   

	   }
//	   cout <<"before return mmax = "<<mmax<<endl;
	   if(thick <= 0.0)
	   {
	   	
//	   	cout <<"return here"<<endl;

     	   imax1 = i;
     	   imax = i + 1;
		     dt = 86400.0 * dtday;
/*		   for(i = ig - 1; i < imax1; i ++)
		   {
			   if( i <= lmimax) continue;

		   }*/
		   x[imax - 1] = xx;
//		   cout <<"mmax = "<<mmax<<endl;
		   return 0;



	   }
       thick *= vspace;
       dxa *= vspace;
       dxb *= vspace;
       dense *= vdens;
       
	   if(thick < dxa)  return 0;
	   if(mat > 10)
	   {
		   if(mat > 20)
		   {
			   if(mat > 30)return 0;
			   temp = tf + 1.0;
			   for(int_ln n = 0; n < 10; n ++)
			   {
     	           msg_debug("1: Entering function Pram from function Grid");	
     	           Pram(mat, temp, temp, dense, wat, sph, cnd);
     	           msg_debug("1: Exit function Pram to function Grid");
				   temp -= 1.0;
			   }


		   }
		   else
		   {

     	       msg_debug("2: Entering function Pram from function Grid");	
     	       Pram(mat, tf, tf, dense, wat, sph, cnd);
     	       msg_debug("2: Exit function Pram to function Grid");
		   }
	   }
	   else
	   {
		   if(kswitch == 0)
		   {
			   switch(layer)
			   {
			   case 1:
				   {
					   vcond = VCOND1[cmnt];
					   vsph = VSPH1[cmnt];
					   condt[mat - 1] = COND1[cmnt];
					   spht[mat - 1] = SPHT1[cmnt];
					   condf[mat - 1] = CONDF1[cmnt];
					   sphf[mat - 1] = SPHF1[cmnt];
				   }
				   break;
			   case 2:
				   {
					   vcond = VCOND2[cmnt];
					   vsph = VSPH2[cmnt];
					   condt[mat - 1] = COND2[cmnt];
					   spht[mat - 1] = SPHT2[cmnt];
					   condf[mat - 1] = CONDF2[cmnt];
					   sphf[mat - 1] = SPHF2[cmnt];

				   }
				   break;
			   case 3:
				   {
					   vcond = VCOND3[cmnt];
					   vsph = VSPH3[cmnt];
					   condt[mat - 1] = COND3[cmnt];
					   spht[mat - 1] = SPHT3[cmnt];
					   condf[mat - 1] = CONDF3[cmnt];
					   sphf[mat - 1] = SPHF3[cmnt];

				   }
				   break;
			   case 4:
				   {
					   vcond = VCOND4[cmnt];
					   vsph = VSPH4[cmnt];
					   condt[mat - 1] = COND4[cmnt];
					   spht[mat - 1] = SPHT4[cmnt];
					   condf[mat - 1] = CONDF4[cmnt];
					   sphf[mat - 1] = SPHF4[cmnt];
				   }
				   break;
			   case 5:
				   {
					   vcond = VCOND5[cmnt];
					   vsph = VSPH5[cmnt];
					   condt[mat - 1] = COND5[cmnt];
					   spht[mat - 1] = SPHT5[cmnt];
					   condf[mat - 1] = CONDF5[cmnt];
					   sphf[mat - 1] = SPHF5[cmnt];
				   }
				   break;
			   default:
				   break;
			   }
			   my_vcond[layer - 1] = vcond;
			   my_vsph[layer - 1] = vsph;
			   my_condt[mat - 1][layer - 1] = condt[mat - 1];
			   my_spht[mat - 1][layer - 1] = spht[mat - 1];
			   my_condf[mat - 1][layer - 1] = condf[mat - 1];
			   my_sphf[mat - 1][layer - 1] = sphf[mat - 1];

		   }
		   else
		   {
			   vcond = my_vcond[layer - 1];
			   vsph = my_vsph[layer - 1];
			   condt[mat - 1] = my_condt[mat - 1][layer - 1];
			   spht[mat - 1] = my_spht[mat - 1][layer - 1];
			   condf[mat - 1] = my_condf[mat - 1][layer - 1];
			   sphf[mat - 1] = my_sphf[mat - 1][layer - 1];

		   }

           if(vcond == 0.0) vcond = 1.0;
           if(vsph == 0.0) vsph = 1.0;
           condt[mat - 1] *= vcond;
           spht[mat - 1] *= vsph;
           condf[mat - 1] *= vcond;
	       sphf[mat - 1] *= vsph;
           cnd = condf[mat - 1];
           sph = sphf[mat  - 1];

	   }
	   double dbl_ = sph * dense;
	   if(dbl_ == 0.0) dbl_ = 0.01; // to protect
	   eps = 86400.0 * epsmin * cnd / dbl_;
	   dxagiv = dxa;
//	   cout <<"MMAX = "<<mmax<<endl;
//	   cout <<"LAYER = "<<layer<<endl;
	   while(1)
	   {
		   if(eps * dtday < 0.0) eps = 0.01; //to protect
		   dxmin = sqrt(eps * dtday);
		   
//		       cout <<"dxmin = "<<dxmin<<" dxa = "<<dxa<<endl;
//		       cin.get();  
           if(dxa >= dxmin) break;
           dxa = dxmin;
//           cout <<"dxa = "<<dxa<<"  thick = "<<thick<<endl;
           if(dxa < thick) break;
    
           dxa = dxagiv;
           dtday *= 0.5;
		       mmax *= 2;

	   }
//	   cout <<"mmax == "<<mmax<<" vspace ="<<vspace<<endl;
	   if(eps == 0.0) eps = 0.01; //to protect
           dtmax = dxa * dxa / eps;
	   if(dtmax < dtbig) dtbig = dtmax; 
	   if(dxa < 0.5 * thick && dxb > dxa)
	   {
			if(vspace == 0.0) vspace = 0.01; // to protect
		    usedxa = dxa / vspace;
			dbl_ = thick - dxb;
			if(dbl_ == 0.0) dbl_ = 0.01; // to protect						
			r = (thick - dxa ) / dbl_;
                        if(dxa == 0.0) dxa = 0.1; // to protect
			ratio = dxb / dxa;
			nn = int_ln(10.0 * log(dxb / dxa) / log(r) + 5.0);
			nn /= 10;
                        if(nn < 0) flog1 <<"nn < 0 in function Grid()"<<endl;
			nodes = nn + 1;

            thmdxa = thick - dxa;
            power = 1.0 / double(nn);
            rwas = pow(ratio, power);         
            fwas = dxa * pow(rwas, double(nodes)) - rwas * thick + thmdxa;
		    err = 1.0E-6;
		    for(int_ln k = 0; k < 20; k ++)
			{
                f = dxa * pow(r, double(nodes)) - rwas * thick + thmdxa;
		        if(fabs(f) <= err) break;
		        df = f - fwas;
				if(df == 0.0) break;
                rnew = r - (r - rwas) * f / df;
                rwas = r;
                fwas = f;
                r = rnew;
			}
			r = rnew;
			dxx = dxa;
			for(int_ln j = 0; j < nodes; j ++)
			{
				i ++;
				if( i > lmimax)
				{

      	            imax1 = i;
      	            imax = i + 1;
					dt = 86400.0 * dtday;
/*				    for(i = ig - 1; i < imax1; i ++)
				    {
					    if( i <= lmimax) continue;

		            }*/
		            x[imax - 1] = xx;
		            return 0;
				}
	            x[i - 1] = xx;
                dx[i - 1] = dxx;
                xx += dxx;
                dxx *= r;
                mater[i - 1] = mat;
                ddry[i - 1] = dense;
                water[i - 1] = wat;
			}

	   }
	   else
	   {
		   if(dxa == 0.0) dxa = 0.01; //to protect
		   if(nodes == 0) nodes = 1;
		   if(vspace == 0.0) vspace = 0.01;

		   nodes = int_ln(thick / dxa + 0.01);
		   dxa = thick / double(nodes);
           usedxa = dxa / vspace;
           for(int_ln j = 0; j < nodes; j ++)
		   {
			   i ++;
      		   if(i > lmimax) 
      	       {

      	           imax1 = i;
      	           imax = i + 1;
			       dt = 86400.0 * dtday;
				   for(i = ig - 1; i < imax1; i ++)
/*				   {
					   if( i <= lmimax) continue;

		           }*/
		           x[imax - 1] = xx;
		           return 0;				   
			   }
			   else
			   {

			       x[i - 1] = xx;
			       dx[i - 1] = dxa;
			       xx += dxa;
			       mater[i - 1] = mat;
			       ddry[i - 1] = dense;
			       water[i- 1] = wat;
			   }

		   }

	   }
//	   cout <<"next layer"<<"  vspace = "<<vspace<<endl;






   }



};





int Qsoiltemp1::Bctop(int_ln& maxm, double taer[], int& ktop, double& sflx1, double& sflx2, double& sflx3)
{
  //-----------------------------------------------------
  // History:
  // March, 1, 2007
  //-----------------------------------------------------
  // Usage:	
	//determines type of B.C. at top surface .. level x[is]
	//reads data according to B.C. type
	//ktop = 1 ... prescribed temperature
	//ktop = 2 ... prescribed heat flux (positive downward)
	//ktop = 3 ... linearized heat flux b.c. (not complete)
	//-----------------------------------------------------
	//subprocedure called:
	// Data()
	//-----------------------------------------------------
	// Input variable
	//
	//-----------------------------------------------------
	// return value
	// integer
	//-----------------------------------------------------

	
	

//	cout <<"ktop = "<<ktop<<endl;
//	system("pause");
	ktop = 1;
	switch(ktop)
	{
	   case 1:
	   {

//----- Prescribed temperature ... taer(m)				
	       htop = 1.0e20;
	       sflux = 0.0;
		  

		   msg_debug("Entering function Data from Bctop");
		   Data(taer, 3, -99, 1.0);
		   msg_debug("Exit function Data to Bctop");
//		  system("pause");
// n-factor multiplication
//   thnfac... thawing n-factor
//   fznfac... freezing n-factor
//   both are set to 1.0 if no value is given
	       double thnfac = 1.0;
//   read thnfac
		   double fznfac = 1.0;

		   double tint, fac;

		   for(int_ln m = 0; m < mmax; m ++)
		   {
			   tint = taer[m];
			   fac = thnfac;
			   if(tint <= 0.0) fac = fznfac;
	   	       taer[m] = fac * tint;	
	     }	   					
	   }
	   break;
	   case 2:
	   {
//---- prescribed heat flux
		   htop = 0.0;
		   for(int_ln m = 0; m < mmax; m ++)
		   {
		      taer[m] = 0.0;
		   }


					
						
		 }
		 break;
		 default:
		 {

			sigma = 56.7e-9;

//      read emiss here
            msg_debug("Entering function Data from Bctop");
			Data(taer);
			msg_debug("Exit function Data to Bctop");
		 }break;		
	}
	return 0;
};

int Qsoiltemp1::Bcbot(int_ln& maxm, double tber[])
{
//-----------------------------------
// History:
// March 19, 2007
//----------------------------------
// Usage:
// Determines type of B.C. at level x[imax]	
// Reads data according to B.C. type
// kbot = 1 ........ prescribed temperature
// kbot = 2 ........ constant heat flux (positive downward)
// kbot = 3 ........ Radiation B.C. (not used)
//-----------------------------------
// input variables
//
//-------------------------------------
// return variable
// integer
//-------------------------------------
   int kbot;

//   cout <<"kbot = "<<kbot<<endl;
//   system("pause");
   kbot = 2;
   switch(kbot)
   {
      case 1:
      {
   	      hbot = 1.0e20;
	      gflux = 0.0;

	      msg_debug("Entering function Data from Bctop");
   	      Data(tber);
   	      msg_debug("Exit function Data to Bctop");
   	      return 0;
      }
      break;
      case 2:
      {
   	      hbot = 0.0;
   	
//  read glfux
          for(int_ln i = 0; i < mmax; i ++)
   	      {
   		      tber[i] = 0.0;
   	      }
   	   	   	  
   	
   	      return 0;
      }
      break;	
      case 3:
      {
//-- radiation B.C.
//........ Not operational   	
   	      gflux = 0.0;
//......... hbot to be set
//......... tber[m] to be set
          return 0;   	   	
      }
      break;
      default:
      {
		  flog1 <<"Bad assignment to ktop in function Bcbot"<<endl;
   	      return 0;
      }
      break;
 }
 
 return 0;
   
};


int Qsoiltemp1::Data(double fitted[], const int index, const int_ln np, const double vtime)
{
// ----------------------------------------
// History:
// March 19, 2007
// ----------------------------------------
// Usage:
// Creates temperature values for each time step m
// allows for 4 different functions temperature vs time
// -----------------------------------------
// subprocedure called:
// ----------------------------------------
// input variables
// ----------------------------------------
// output variables
// integer
// -----------------------------------------
   double elraw[NMAX0], raw[NMAX0];
   char comment[30];	
// definition of input parameters read
//  index = 1 ... constant value
//  index = 2 ... sinusoidal variation
//  index = 3 ... interpolation from tabulated data(nonperiodic)
//  index = 4 ... interpolation from tabulated data (periodic) 
//  vtime ... convert time to days
//  np .... spline smoothing parameter
//  int index = 3;
//	int_ln np = -99;

//  double vtime = 1.0;
    int_ln nmax;
//    cout <<"Entering here"<<endl;

//	system("pause");
//  set up time arrays
    double cal;
    for(int_ln m = 0; m < mmax; m ++)
    {
    	eltim[m] = double(m + 1) * dtday;
    	cal = eltim[m] + first;
    	if(cal > per) cal -= per;
    	caltim[m] = cal;
    	
    }  
	double value;
//	cout <<"index = "<<index<<endl;
//	system("pause");
    switch(index)
    {
    	case 1:
    	{
//----- type # 1
//----- constant value



        
            for(int_ln m = 0; m < mmax; m ++)
                fitted[m] = value;
            return 0;               		
    	}
    	break;
    	case 2:
    	{
//----- type # 2
//------ sinusoidal variation about annual mean
            double tmean, tampsn, tampcs, phas;
            // need input here
            phas = 0.0;
            phas *= vtime;

            const double pie = 3.141593;
            const double omega = 2.0 * pie / per;              		
            double omt;
            for(int_ln m = 0; m < mmax; m ++)
            {
                omt = omega * (eltim[m] + phas);
                fitted[m] = tmean + tampsn * sin(omt) + tampcs * cos(omt);
            	
            }
//         cout <<" interpolated values "<<endl;
//         cout <<" index " << " time "<<" calendar " << " value "<<endl;
//            for(int_ln i = 0; i < mmax; i += nst)
//         cout << i <<" " <<eltim[i] <<" "<<caltim[i] <<" "<< fitted[i]<<endl;
    	    return 0;
    	}
    	break;
    	case 3:
    	{
//------ type # 3
//------ interpolation of non periodic data

            
            elraw[0] = 0.0;
			raw[0] = airt1;
			elraw[1] = 15.0;
		    raw[1] = airt2;
		    elraw[2] = 30.0;
	        raw[2] = airt3;
//                  cout <<setw(7)<<setprecision(2)<<airt1<<setw(7)<<setprecision(2)<<airt2<<setw(7)<<setprecision(2)<<airt3<<endl;
		    nmax = 3;
//	        ofsnore <<setw(8)<<nmax<<setw(8)<<mmax<<endl;	  
            msg_debug("Entering function Intrpl from Data");
            Intrpl(elraw, raw, nmax, np, eltim, fitted, mmax);
            msg_debug("Exit function Intrpl to Data");
    	    return 0;          
        
    	}
    	break;
    	case 4:
    	{
// ----- type # 4
// interpolation of periodic data
// creates three successive data cycles before interpolation
//   need input from external file
            int_ln n = 0;
            double day = -1.0e1; 
            double delta = 0.0;
            double dayo;
            while(1)
            {

                dayo = day;    

                if(day < 0.0) break;       
                if(day < dayo) delta = per;
                n ++;
                elraw[n - 1] = vtime * day - first + delta;
                raw[n - 1] = value;
            

                elraw[n - 1] -= per;
            
                         
           }
           int_ln nn = n;
           int_ln nn2 = 2 * nn;
           nn ++;
           int_ln no = 0;
           
           for(n = nn - 1; n < nn2; n ++)
           {

           	   no ++;
           	   elraw[n] = elraw[no - 1] + per;
           	   raw[n] = raw[no - 1];
           	   nmax = nn2 + no;
           	   elraw[nmax - 1] = elraw[no - 1] + 2.0 * per;
           	   raw[nmax - 1] = raw[no - 1];
           }
//    interpolate
           msg_debug("Entering function Intrpl from function Data()"); 
           Intrpl(elraw, raw, nmax, np, eltim, fitted, mmax);
           msg_debug("Exit function Intprl to function Data()");
          
    	   return 0;                 

    	}
    	break;
    	default:
    	{
    	
    	  	flog1 <<"Bad assignment to index in Data"<<endl;
    	    return 0;	
    	}
    	
    	
	}
	return 0;
}

int Qsoiltemp1::Tinitl(const int& cmnt)
{
//------------------------------------
// History:
// April 3, 2007
//-----------------------------------------
// Usage:
// Creates initial temperatures for nodes on calculation grid
// used for nodes between ig & imax only
// vdepth = conversion factor depth units to meters
//  np = decimal precision of temperature data
//
// -----------------------------------
// subprocedure called:
// -------------------------------------
// input variables
// ------------------------------------
// output variables
// ------------------------------------

    double xstart[210], tstart[210], xx[210], tt[210];
    int_ln index, np, num, nn1, i, j, jmax, n, nmax, nn;
    double vdepth, temp, dxstr, dtemp, x0, xbot;
	if(kswitch == 0)
	{
		index = int_ln(INDEX[cmnt]);
		vdepth = VDEPP[cmnt];
		np = NP[cmnt];
		my_index = index;
		my_vdepth = vdepth;
		my_np = np;

	}
	else
	{
		index = my_index;
		vdepth = my_vdepth;
		np = my_np;
	}
//	cout <<"index = "<<index<<" in Tinitl"<<endl;
	if(index==1)return 0;
	else
	{
		n = 0;		
		while(1)
		{
			n ++;
			if(kswitch == 0)
			{
				xstart[n - 1] = DEPTH[cmnt][n - 1];
				tstart[n - 1] = TEMP[cmnt][n - 1];

				my_xstart[n - 1] = xstart[n - 1];
				my_tstart[n - 1] = tstart[n - 1];
			}
			else
			{
				xstart[n - 1] = my_xstart[n - 1];
				tstart[n - 1] = my_tstart[n - 1];

			}
			xstart[n - 1] *= vdepth;
			if(xstart[n - 1] <= -1.00 * vdepth)break;
		}
		nmax = n - 1;
//		cout <<"nmax = "<<nmax <<" in Tinitl"<<endl;
		if(nmax > 2)
		{
			nn = n - 1;
			xbot = x[imax - 1];
			if(xstart[nn - 1] < xbot)
			{
				dxstr = xstart[nn - 1] - xstart[nn - 2];
				if(dxstr == 0.0) dxstr = 0.01;//to protect
				num = int_ln((xbot - xstart[nn - 1]) / dxstr);
				if(num <= 1 ) num = 1;
	            dxstr = (xbot - xstart[nn - 1]) / double(num);
                dtemp = (tstart[nn] - tstart[nn - 1]) / double(num);
                nn1 = nn + 1;
                nmax = num + nn1;
                for(int_ln n = nn1 - 1; n < nmax; n ++)
                {

    	            xstart[n] = xstart[n - 1] + dxstr;
    	            tstart[n] = tstart[n - 1] + dtemp;
                }

			}
			jmax = imax - igm1;
	        for(int_ln j = 1; j <= jmax; j ++)
	        {
		        i = igm1 + j;
     	        xx[j - 1] = x[i - 1];   
//				cout <<"xx["<<j-1<<"] = "<<xx[j - 1]<<endl;
            }
            msg_debug("Entering function Intrpl from function Tinitl");
            Intrpl(xstart, tstart, nmax, np, xx, tt, jmax);
            msg_debug("Exit function Intrpl to function Tinitl");
     
	        for(int_ln j = 1; j <= jmax; j ++)
	        {
		        i = igm1 + j;
		        t[i - 1] = tt[j - 1];
//				cout <<"t["<<i-1<<"] = "<<t[i - 1]<<endl;
	        }
			return 0;
		}
		else
		{
			double dbl_ = xstart[1] - xstart[0];
			if(dbl_ == 0.0) dbl_ = 0.01;//to protect
			dtemp = (tstart[1]  - tstart[0]) / dbl_;
			temp = tstart[0];   	
	        x0 = xstart[0];
			for(int_ln i = ig - 1; i < imax; i ++)
			{
				t[i] = temp + dtemp * (x[i] - x0);

			}
			return 0;

		}




	}


};

int Qsoiltemp1::Snofal(double snow[], double sden[], double& comp, double& eta0, double& denfac, double& fcmelt, const int& cmnt)
{
//-----------------------------------
// History:
// April 3, 2007
//-----------------------------------
// Usuage:
// Snow accumulation rate for each full time step
//---------------------------------------
// subprocedure called:
//	
//-----------------------------------------
// input variables
//-------------------------------------
// return value
//----------------------------------------
	double elraw[365],arate[365],dnew[365],height[365];
	double cal, epssno, convrt, denmax, abel, cnd, sphsno, daton, datmax, snoden, depmax; 
	double cummax, rate, day, elap, dend, rat;
	int index;
	int_ln mdat1, mdat2, m, m1, m2, nmax, nm1;

	if(kswitch == 0)
	{
		index = SNOFAL[cmnt];
		my_index2 = index;
	}
	else
	{
		index = my_index2;
	}
	if(index <= 0) return 0;

	if(kswitch == 0)
	{
		epssno = EPSSNO[cmnt];
		convrt = CONVRT[cmnt];
		eta0 = ETAO[cmnt];
		denfac = DENFAC[cmnt];
		fcmelt = FCMELT[cmnt];
		denmax = DENMAX[cmnt];

		my_epssno = epssno;
		my_convrt = convrt;
		my_eta0 = eta0;
		my_denfac = denfac;
		my_fcmelt = fcmelt;
		my_denmax = denmax;

	}
	else
	{
		epssno = my_epssno;
		convrt = my_convrt;
		eta0 = my_eta0;
		denfac = my_denfac;
		fcmelt = my_fcmelt;
		denmax = my_denmax;
	}

    abel = 2.9E-6;
    cnd = abel * denmax * denmax;
    sphsno = 2000.0;
    comp = sqrt(cnd * denmax * dt * epssno / sphsno);
    if(epssno < 0.1) comp = 1.0E10;
	switch(index)
	{
	case 1:
		{
//-------  Type # 1
// snow thickness starts at daton & increases linearly to depmax
// snow thickness constant at depmax from time datmax until melting

//  daton = day of onset of snow
//  datmax = day when snow reaches maximum depth (after linear increase
//            in thickness from time daton)
//  snoden = snow density (homogenous snow-cover)
//  depmax = maximum snow thickness
               //need input from external files

			if(datmax < daton) datmax += per;          
            cummax = snoden * depmax;
            double dbl_ = datmax - daton + 1.0;
            if(dbl_ == 0.0) dbl_ = 0.01;
            rate = cummax / dbl_;
            for(int_ln m = mdat1 - 1; m < mdat2; m ++)
            {			
				
          	    sden[m] = snoden;
          	    snow[m] = rate;
            }
            return 0;
		}
		break;
	case 2:
		{			
// type #2
// accumulation rate & new snow density given function of time

// end of data when day is negative
		}
		break;
	case 3:
		{

//----- Type # 3
//----- total snow height and mean density given
//----- spreads accumulation linearly between given dates
			nmax = mmax / 2 - 1;
			if(mmax == 0) mmax = 1; //to protect
//                  cout <<setw(7)<<setprecision(2)<<hsnow1<<setw(7)<<setprecision(2)<<hsnow2<<setw(7)<<setprecision(2)<<hsnow3<<endl;
			for(m = 1; m <= nmax; m ++)
			{
				rat = (hsnow2 - hsnow1) * 4.0 / double(mmax);
				sden[m - 1] = 250.0;
				snow[m - 1] = rat * sden[m - 1];
//				cout <<setw(7)<<setprecision(2)<<snow[m - 1]<<setw(7)<<setprecision(2)<<sden[m - 1]<<endl;
			}
			nmax ++;
			for(m = nmax; m <= mmax; m ++)
			{
				rat = (hsnow3 - hsnow2) * 4.0 / double(mmax);				
				sden[m - 1] = 250.0;
				snow[m - 1] = rat * sden[m - 1];
//				cout <<setw(7)<<setprecision(2)<<snow[m - 1]<<setw(7)<<setprecision(2)<<sden[m - 1]<<endl;				
			}
//      cout <<"mmax = "<<mmax<<endl;
//      cout <<"snow[28] = "<<snow[28]<<" sden[28] = "<<sden[28]<<endl;
//      cout <<"snow[29] = "<<snow[29]<<" sden[29] = "<<sden[29]<<endl;      
			return 0;
		}
		break;
	default:
		{
		}
		break;
	}
   
    return 0;   










};
 

 
 
int Qsoiltemp1::Neige(double& acum, double& snoden, double weight[], double& comp, double& eta0, double& denfac, double& fcmelt, double& tdrive)
{
//----------------------------------------
// History:
// March 20, 2007
//----------------------------------------
// Usage:
// This subroutine calculates the density & thickness at time m of layers
// is, is + 1, ... igm1, where is = top of top layer, igm1 = top of bottom layer
// also calculates compaction of surface layer is
// number of layers changes only when sufficient snow has fallen
// (i.e. a layer may correspond to several time steps)
// comp = initial mass of a full layer
// includes accumulation, mechanical ablation & and melting
// densification by successive iteration centered difference
// all units are S.I.U. except acum & fcmelt (days)
// -----------------------------------------	
// Subprocedure called:
//
//------------------------------------------
// Input variables
//------------------------------------------
// return value
//-----------------------------------------
//--- check for viable values of IS

    if(is <= 2 || is > ig) 
    {
		flog1 <<"Warning: Snow node counter out of range"<<endl;		    	
    	return 0;
	}
	int_ln jj;

	double rit, dxold, ratio, xx;
    double dsmass = acum * dtday;
    if(smass <= 0.0 && dsmass <= 0.0) return 0;
    smass += dsmass;
    if(dsmass <= 0.0) 
    {
//50     mechanical ablation


        double ablat = - dsmass;
        while(1)
        {
              if(ablat < weight[is - 1]) break;
              if(is == ig)
              {

       	          smass = 0.0;
       	          weight[igm1 - 1] = 0.0;
       	          xfa[igm1 - 1] = - 1.0E10;
       	          return 0;
              }
              ablat -= weight[is - 1];
              weight[is - 1] = 0.0;
              dx[is - 1] = 0.0;
              xfa[is - 1] = -1.0E10;
              is ++;
        }
        weight[is - 1] -= ablat;    	
	}
	else
	{
// accumulation
	    double tsno = t[is - 1];
        if(is == ig) is = igm1;
        double wtis = weight[is - 1] + dsmass;
        double wtisd, remain;
        if(wtis > comp)
	    {
//40
		    wtisd = comp - weight[is - 1];
			if(snoden == 0.0) snoden = 0.01; //to protect
			double dbl_;
			dbl_ = dx[is - 1] + wtisd / snoden;
			if(dbl_ == 0.0) dbl_ = 0.01;//to protect
    	    ddry[is - 1] = comp / dbl_;
    	    weight[is - 1] = comp;
    	    remain = dsmass - wtisd;
    	    t[is - 1] = tsno;
    	    while(1)
    	    {


    	        is --;
    	        weight[is - 1] = comp;
    	        ddry[is - 1] = snoden;
    	        t[is - 1] = tsno;
    	        xfa[is - 1] = -1.0E10;
    	        if(remain <= comp) break;
    	        remain -= comp;
            } 
            weight[is - 1] = remain;
        }
        else
	    {
// add only a fraction of a layer

			if(snoden == 0.0) snoden = 0.01; //to protect
			double dbl_ = dx[is - 1] + dsmass / snoden;
			if(dbl_ == 0.0) dbl_ = 0.01; // to protect
            ddry[is - 1] = wtis / dbl_;
            if(weight[is - 1] <= 0.0) xfa[is - 1] = -1.0E10;
            weight[is - 1] = wtis;
            t[is - 1]= tsno;    	
	    }
      }
	  if(tdrive > 0.0)
	  {
//   melt ablation

	      double fonte = fcmelt * dtday * tair;
	      smass -= fonte;
	      if(smass > 0.0) tdrive = 1.0E-2;
	      while(1)
	      {


      	      if(fonte < weight[is  - 1])break; 
     	      if(is == ig)
     	      {

// complete disapperance of snow cover     		

     	          smass = 0.0;
     	          weight[igm1 - 1] = 0.0;
     	          xfa[igm1 - 1] = -1.0E10;	
     	          return 0;
     	      }
     	      fonte -= weight[is - 1];
     	      weight[is - 1] = 0.0;
     	      dx[is - 1] = 0.0;
     	      xfa[is  - 1] = -1.0E10;
     	      is ++;
     	
          }
          weight[is - 1] -= fonte;
// Compaction of all layers including surface layer is
// assumes mass of each layer constant
// sload includes 1/2 weight of node i
         	
     }
    
     for(int_ln i = is - 1; i < igm1; i ++)
     {

//---- calculate layer thickness ... assume mass constant     	
      	 mater[i] = 21;
     	 rit = ddry[i];
     	 dxold = dx[i];
		 if(rit == 0.0) rit = 0.01; //to protect
     	 dx[i] = weight[i] / rit;
     	 if(xfa[i] < 0.0) continue;
//----- update phase boundaries for changed dx[i]
		 if(dxold == 0.0) dxold = 0.01; // to protect
         ratio = dx[i] / dxold;
         if(xfb[i] != 0.0)                  	                
         {
             xfb[i] += (ratio - 1.0) * (dx[i] - xfa[i]);
             if(xfb[i] < 0.0) xfb[i] = 0.01 * dx[i];	
         }
         xfa[i] *= ratio; 
     }
// update layer heights
     xx = x[ig];     
	 jj = is + igm1;
     int_ln i;
     for(int_ln j = is; j <= igm1; j ++)
     {

      	 i = jj - j;
      	 xx -= dx[i - 1];
      	 x[i - 1] = xx;
     }
     for(i = is - 1; i < ig; i ++)
     {
      	if(t[i] > 0.0) t[i] = 1.0E-2;
      
     }
	 return 0;
     
    
}; 


int Qsoiltemp1::Schnee(double& acum, double& snoden, int_ln& ksnow)
{
//--------------------------------------
// History:
// March 20, 2007
//--------------------------------------
// Usuage:
//-- Print out snowcover data (layer height, density, temperature, etc.)
//--------------------------------------
// Subprocedure called:
//--------------------------------------
// Input variables
//--------------------------------------
// return value
// integer
//--------------------------------------
   double xmid[100];
   if(ksnow == 0) return 0;
   for(int_ln i = is - 1; i < igm1; i ++)
   {
     xmid[i] = - x[i] - 0.5 * dx[i];	
   }
//   cout <<" Time = " << time <<" smass = " << smass <<" " <<is<<endl;
   if(smass <= 0.0) return 0;
//   cout <<" accum rate = " << acum <<" surf density = " << snoden
//        <<" air temperature = "<<tair <<" surf temperature = "<<t[is - 1]
//        <<" base temperature = "<<t[ig - 1]<<endl;
//   cout <<" Layer top    " << " Thickness " <<" density "<<" T(i) " <<" x-mid "<<endl;
   for(int_ln i = is - 1; i < igm1; i ++)
   {
//   	cout <<i + 1<<" "<<x[i]<<" "<<dx[i]<<" "<<ddry[i]<<" "<<t[i] <<" "<<xmid[i]<<endl;
   }     
   return 0;

   
};


int Qsoiltemp1::Asmblo(int_ln& i1)
{
//---------------------------
// History:
// March 20, 2007
//---------------------------
// Usuage:
// Assembles e[i], s[i] for ordinary node
// for upsweep e, s (imax to i1)
//----------------------------	
// subprocedure called:
//----------------------------
// input variables
//----------------------------
// return value
//----------------------------
// bottom node ordinary implicit
   double con = conx[imax1 - 1];
   if(dt == 0.0) dt = 0.01; // to protect
   double rc = 0.5 * capx[imax1 - 1] / dt;
   double denm = rc + con + hbot;
   if(denm == 0.0) denm = 30.0; // to protect
   s[imax - 1] = con /denm;
   e[imax - 1] = (rc * t[imax - 1] - gflux + hbot * tbot + ht[imax - 1]) / denm;
   if(i1 >= imax) return 0;
//
// Theat - method
   int_ln i2 = i1;   
   if(i1 == is) i2 = is + 1;
   int_ln jj = i2 + imax1;
   int_ln im1, ip1;
   double rhs, conut1, condt1;
   double r, conuth, condth;
   int_ln i;

   for(int_ln j = i2; j <= imax1; j ++)
   {
	    i = jj - j;
	    im1 = i - 1;
	    ip1 = i + 1;
		if(dt == 0.0) dt = 0.01; // to protect
	    rc = (capx[i - 1] + capx[im1 - 1]) / dt;
	    conut1 = conx[im1 - 1] * theta1;
	    condt1 = conx[i - 1] * theta1;
	    rhs = (rc - conut1 - condt1) * t[i - 1] + conut1 * t[im1 - 1] + condt1 * t[ip1 - 1];
	    rhs += theta * ht[i - 1] + theta1 * htold[i - 1];
	    conuth = conx[im1 - 1] * theta;
   	    condth = conx[i - 1] * theta;
   	    r = conuth + rc + condth;
   	    denm = r - condth * s[ip1 - 1];
		if(denm == 0.0) denm = 30.0; // to protect
   	    s[i - 1] = conuth / denm;
   	    e[i - 1] = (rhs + condth * e[ip1 - 1]) / denm;
   }
   if(i1 > is) return 0;
// 
// Surface node ordinary implicit - prescribed flux permitted
   if(dt == 0.0) dt = 0.01; // to protect
   rc = 0.5 * capx[is - 1] / dt;
   double conu = conx[ism1 - 1];
   double cond = conx[is - 1];
   rhs = rc * t[is - 1] + sflux + ht[i - 1];
   r = conu + rc + cond;
   denm = r - cond * s[i2 - 1];
   s[is - 1] = conu / denm;
   if(denm == 0.0) denm = 30.0; // to protect
   e[is - 1] = (rhs + cond * e[i2 - 1]) /denm;
   return 0;
   
   

};



int Qsoiltemp1::Asmabv(int_ln& i1)
{
//---------------------------------
// History:
// March 20, 2007
//--------------------------------
// Usuage:
// Assembles e[i], s[i] for ordinary node
// for downsweep e, s (is to i1)
//----------------------------------
// Subprocedure called:
//
//---------------------------------
// Input variables:
//----------------------------------
// Return values
// Integer
// ------------------------------------
//
// Surface node ordinary implicit -prescribed flux permitted
   if(i1 < is) return 0;
   if(dt == 0.0) dt = 0.01; // to protect
   double rc = 0.5 * capx[is - 1] / dt;
   double conu = conx[ism1 - 1];
   double cond = conx[is - 1];
   double rhs = rc * t[is - 1] + sflux + ht[is - 1];
   double r = conu + rc + cond;
   double denm = r - conu * s[ism1 - 1];
   if(denm == 0.0) denm = 30.0; // to protect
   s[is - 1] = cond / denm;
   e[is - 1] = (rhs + conu * e[ism1 - 1]) / denm;
   if(i1 == is) return 0;

// Theta method
   int_ln isp1 = is + 1;   
   int_ln im1, ip1;
   double conut1, condt1;
   double conuth, condth;
   for(int_ln i = isp1; i <= i1; i ++)
   {
   	   im1 = i - 1;
   	   ip1 = i + 1;
   	   rc = (capx[i - 1] + capx[im1 - 1]) / dt;
   	   conu = conx[im1 - 1];
   	   cond = conx[i - 1];
   	   conut1 = theta1 * conu;
   	   condt1 = theta1 * cond;
   	   rhs = (rc - conut1 - condt1) * t[i - 1] 
   	      + conut1 * t[im1 - 1] + condt1 * t[ip1 - 1];
   	   rhs += theta * ht[i - 1] + theta1 * htold[i - 1];
   	   conuth = theta * conu;
   	   condth = theta * cond;
   	   r = conuth + rc + condth;
   	   denm = r - conuth * s[im1 - 1];
	   if(denm == 0.0) denm = 30.0; //to protect
   	   s[i - 1] = condth / denm;
   	   e[i - 1] = (rhs + conuth * e[im1 - 1]) /denm;   	        	      
   	
   }
   return 0;

	
};

int Qsoiltemp1::Tabv(int_ln& i1)
{
//-------------------------------
// History:
// March 20, 2007
//--------------------------------
// Usuage:
// Evaluates new node temperatures between is & i1
// backward substitution (for asmabv)
//-----------------------------------
// Subprocedure called:
//---------------------------------
// Input variables
//-----------------------
// return value
// Integer
//----------------------------
//
   if(i1 < is) return 0;
   int_ln jj = is + i1;
   int_ln  i;
   for(int_ln j = is; j <= i1; j ++)
   {
       i = jj - j;
       t[i - 1] = s[i - 1] * t[i] + e[i - 1];	
   }
   return 0;
	
};


int Qsoiltemp1::Asmone(int_ln& i)
{
//-----------------------------------
// History:
// March 20, 2007
//-------------------------------------
// Usuage:
// Calcualtes new position of single phase plane
// Calculates new temperatures at all nodes
//--------------------------------------
// Subprocedure called:
// Asmabv(), Asmblo()
// Pram(), Phase()
// Tblo(), Tabv()
//--------------------------------------
// Input variables:
//---------------------------------------
// return variables
//---------------------------------------
   double dtwas = dtfaz;
   double gold = xfa[i - 1];
   int_ln counter =0;
   int_ln im1 = i - 1;
   int_ln ip1 = i + 1;
   int_ln ip2 = i + 2;
   double dxx, wat, dens, tt, cnd, sph;
   double cu, rsu, td, cd, rsd, rlw, t1, t2, g;
   int llc_n=0;
   //llc
//   printf("1,\n");
   while(1)   
   {
	//  cout <<"counter = "<<counter ++ <<endl;
     llc_n++;
      //llc
     Asmabv(im1);
     Asmblo(ip2);
     dxx = dx[i - 1];
     wat = water[i - 1];
//	 cout <<"water["<<i-1<<"] = "<<water[i - 1]<<endl;
     dens = ddry[i - 1];
     tt = t[i - 1];
       
     msg_debug("1: Entering function Pram from function Asmone");     
     Pram(mater[i - 1], tt, tf, ddry[i - 1], water[i - 1], sph, cnd);
     msg_debug("1: Exit function Pram to function Asmone");
     
     cu = cnd;
     rsu = dens * sph;
     td = t[ip1 - 1];
     
     msg_debug("2: Entering function Pram from function Asmone");
     Pram(mater[i - 1], tf, td, ddry[i - 1], water[i - 1], sph, cnd);
     msg_debug("2: Exit function Pram to function Asmone");
     
     cd = cnd;
     rsd = dens * sph;
     rlw = dens * hlat * wat;
  
//     cout <<"dens = "<<dens<<" hlat = "<<hlat<<" wat ="<<wat<<endl;
     msg_debug("1: Entering function Phase from function Asmone");
     Phase(i, rlw, cu, rsu, cd, rsd, tt, td, dxx, gold, g, t1, t2);
     msg_debug("1: Exit function Phase from function Asmone");

// Evaluate whether g goes out of element
     if( g >= 0.0)
     {
       if(g > dxx)
       {
// out below       	

       	   capx[i - 1] = rsu * dxx;
		   if(dxx == 0.0) dxx = 0.01; // to protect
       	   conx[i - 1] = cu / dxx;
       	   xfa[i - 1] = -1.0E10;
       	   if(i == imax1) return 0;
       	   i ++;
       	   im1 = i - 1;
       	   ip1 = i + 1;
       	   ip2 = i + 2;
       	   t[i - 1] = tf - sign_d(1.0E-4, td - tf);
       	
       	   msg_debug("Entering function Asmabv from Asmone");
       	   Asmabv(im1);
       	   msg_debug("Exit function Asmabv to Asmone");
       	
       	   msg_debug("Entering function Tabv from Asmone");
       	   Tabv(im1);
       	   msg_debug("Exit function Tabv to Asmone");
       	
       	   msg_debug("Entering function Asmblo to Asmone");
       	   Asmblo(ip1);
       	   msg_debug("Exit function Asmblo to Asmone");	
       	
       	   msg_debug("Entering function Tblo from Asmone");
       	   Tblo(ip1);
       	   msg_debug("Exit function Tblo to Asmone");
       	
       	   dtfaz = dtwas - dtfaz;
       	   dtwas = dtfaz;
       	   gold = 0.0;

       	   continue;   	
       }
//     g remains within element i     
       t[i - 1] = t1;
       msg_debug("Entering function Tabv from Asmone");
       Tabv(im1);
       msg_debug("Exit function Tabv to Asmone");
   
       t[ip1 - 1] = t2;
   
       msg_debug("Entering function Tblo from Asmone");
       Tblo(ip2);
       msg_debug("Exit function Tblo to Asmone");
   
       xfa[i - 1] = g;  

       return 0;  	
     }
// out above
     capx[i - 1] = rsd * dxx;
	 if(dxx == 0.0) dxx = 0.01; // to protect
     conx[i - 1] = cd / dxx;
     xfa[i - 1] = -1.0E10;
     if(i == is) return 0;
     i --;
     im1 = i - 1;
     ip1 = i + 1;
     ip2 = i + 2;
     t[ip1 - 1] = tf - sign_d(1.0E-4,tt - tf);   
     
     msg_debug("1: Entering function Asmabv from function Asmone");
     Asmabv(i);
     msg_debug("1: Exit function Asmabv from function Asmone");
     
     msg_debug("1: Entering function Tabv from function Asmone");
     Tabv(i);
     msg_debug("1: Exit function Tabv to function Asmone");
     
     msg_debug("1: Entering function Asmblo from function Asmone");
     Asmblo(ip2);
     msg_debug("1: Exit function Asmblo to function Asmone");
     
     msg_debug("1: Entering function Tblo from function Asmone");
     Tblo(ip2);
     msg_debug("1: Exit function Tblo to function Tblo");
     
     dtfaz = dtwas - dtfaz;
     dtwas = dtfaz;
     gold = dx[i - 1];

  
   if(llc_n>1000) {
 //   printf("llcn=%i,g=%f,i=%i,is=%i,\n",llc_n,g,i,is);
   return 0;
    }
   }   


   return 0;
   	    
   	
};

int Qsoiltemp1::Asmtwo(int_ln& itzero, int_ln& ibzero)
{
//----------------------------------
// History:
// March 20, 2007
//----------------------------------
// Usuage:
// Calculates new positions of more than one phase planes
// Calcualtes new temperatures at all nodes
// Phase planes are uncoupled
// Temperature = 0.0 between phase planes
// only exterior phase planes are mobile
//-------------------------------------------
// Subprocedure called:
//--------------------------------------------
// Input variables
//-------------------------------------------
// Output variables
//-------------------------------------------
   double gaold = xfa[itzero - 1];
   double gbold = xfa[ibzero - 1] + xfb[ibzero - 1];
//   cout <<"gaold = "<<gaold<<" gbold = "<<gbold<<endl;
//   cin.get();
// section for upper phase plane
   int_ln i = itzero;
   double dtwas = dtfaz;
   double v = xfb[i - 1];
   int_ln im1 = i - 1;
   int_ln ip1 = i + 1;
   int_ln ip2 = i + 2;
   double dxx, wat;
   double dens, tt, td, cu, cnd, rsu, rlw, sph;
   double t1, t2, ga, gv;
   double gdn, vdn;
   while(1)
   {
     msg_debug("1: Entering function Asmabv from function Asmtwo");	
     Asmabv(im1);
     msg_debug("1: Exit function Asmabv to function Asmtwo");
     
     dxx  = dx[i - 1];
     wat = water[i - 1];
     dens = ddry[i - 1];
     tt = t[i - 1];
     td = t[ip1 - 1]; 
     
     msg_debug("1: Entering function Pram from function Asmtwo");
     Pram(mater[i - 1], tt, tf, ddry[i - 1], water[i - 1], sph, cnd);
     msg_debug("1: Exit function Pram to function Asmtwo");
//     cout <<"cnd = "<<cnd <<" sph = "<<sph<<" wat = "<<wat<<endl;
     cu = cnd;
     rsu = dens * sph;
     rlw = dens * hlat * wat;
//     cout <<i<<" "<<rlw<<" "<<cu<<" "<<rsu<<" "<<tt<<" "<<tf<<" "<<dxx<<" "<<gaold<<" "<<ga<<endl;
     msg_debug("Entering function Phase from function Asmtwo");
     Phase(i, rlw, cu, rsu, 0.0, 0.0, tt, tf, dxx, gaold, ga, t1, t2); 
     msg_debug("Exit function Phase to function Asmtwo");
//     cout <<"t1 = "<<t1<<" t2 = "<<t2<<endl;
//	 cin.get();
//   Evaluate whether ga goes out of element
     if(v != 0.0)
     {
// Second phase plane within elemnt exists     	
     	t[i - 1] = t1;
     	
     	msg_debug("1: Entering function Tabv from function Asmtwo");	
     	Tabv(im1);
     	msg_debug("1: Exit function Tabv to function Asmtwo");
     	
     	gv = v + gaold;
     	xfb[i - 1] = gv - ga;
     	xfa[i - 1] = ga;
     	if(gv <= ga)
     	{
     	  xfa[i - 1] = -1.0E10;
     	  xfb[i - 1] = 0.0;
     	  if(i >= ibzero)
     	  {
     	    xfa[i - 1] = ga;
     	    xfb[i - 1] = gv - ga;	
     	  }	
     	}
     	break;

     }
           
     else
     {
     	if(ga > dxx)
     	{
//        out below     	   
     	   gdn = xfa[ip1 - 1];
     	   capx[i - 1] = rsu * dxx;
		   if(dxx == 0.0) dxx = 0.01; // to protect
     	   conx[i - 1] = cu / dxx;
     	   xfa[i - 1] = -1.0E10;
     	   if(i == imax1) return 1;
     	   i ++;
     	   im1 = i - 1;
     	   ip1 = i + 1;
     	   ip2 = i + 2;
     	   dtfaz = dtwas - dtfaz;
     	   dtwas = dtfaz;
     	   t[i - 1] = tf - sign_d(1.0E-4, td - tf);
     	   
     	   msg_debug("2: Entering function Asmabv from function Asmtwo");
     	   Asmabv(im1);
     	   msg_debug("2: Exit function Asmabv to function Asmtwo");
     	   
     	   msg_debug("2: Entering function Tabv from function Asmtwo"); 
     	   Tabv(im1);
     	   msg_debug("2: Exit function Tabv to function Asmtwo");
     	   
     	   gaold = 0.0;
     	   if(gdn == -1.0E10)continue;	
     	   vdn = xfb[i - 1];	
     	   if(vdn != 0.0)
     	   {
     	     gaold = vdn + gdn;
     	     v = 0.0;
     	     xfb[i - 1] = 0.0;	
     	   }
     	   else
     	   {
     	     v = gdn;	
     	   }
     	}
     	else
     	{
// ga remains within element i
           t[i - 1] = t1;
           
           msg_debug("3: Entering function Tabv from function Asmtwo"); 
           Tabv(im1);
           msg_debug("3: Exit function Tabv to function Asmtwo");
           
           xfa[i - 1] = ga;
   	       break;	
     	}
     }

   }
   
   
   if(ibzero < itzero) return 0;
     	
//  section for lower phase plane     	
    i = ibzero;
    dtwas = dt;
    dtfaz = dt;
    double gold = xfa[i - 1];
    v  = xfb[i - 1];
    im1 = i - 1;
    ip1 = i + 1;
	ip2 = i + 2;
	double gb, rsd, cd, gup, vup;
    while(1)
    {       
       
       msg_debug("1: Entering function Asmblo from function Asmtwo"); 
       Asmblo(ip2);
       msg_debug("1: Exit function Asmblo to function Asmtwo");
       
       dxx = dx[i - 1];
       wat = water[i - 1];
       dens = ddry[i - 1];
       tt = t[i - 1];
       td = t[ip1 - 1];
       
       msg_debug("2: Entering function Pram from function Asmtwo");
       Pram(mater[i - 1], td, tf, ddry[i - 1], water[i - 1], sph, cnd);
       msg_debug("2: Exit function Pram to function Asmtwo");
       
       cd = cnd;
       rsd = dens * sph;
       rlw = dens * hlat * wat;
       
	     msg_debug("2: Entering function Phase from function Asmtwo");
       Phase(i, rlw, 0.0, 0.0, cd, rsd, tf, td, dxx, gbold, gb, t1, t2);        
       msg_debug("2: Exit function Phase to function Asmtwo");
       
//    Evaulate whether g goes out of element
       if(v != 0.0)
       {
//    second phase plane within element exists
          t[ip1 - 1] = t2;
          
          msg_debug("1: Entering function Tblo from function Asmtwo");
          Tblo(ip2);
          msg_debug("1: Exit function Tblo to function Asmtwo");
          
          xfa[i - 1] = gold;
          xfb[i - 1] = gb - gold;
          if(gb > gold) return 0;
          xfa[i - 1] = -1.0E10;
          xfb[i - 1] = 0.0;
          return 0;	
        }
        if(gb >= 0.0)
        {
//   gb remains within element i
            t[ip1 - 1] = t2;
          
            msg_debug("2: Entering function Tblo from function Asmtwo");
            Tblo(ip2);
            msg_debug("2: Exit function Tblo to function Asmtwo");
          
            xfa[i - 1] = gb;
            return 0;            	
        }
//     out above
        gup = xfa[im1 - 1];
        capx[i - 1] = rsd * dxx;
		if(dxx == 0.0) dxx = 0.01;// to protect
        conx[i - 1] = cd / dxx;
        xfa[i - 1] = -1.0E10;
        if(i == is) return 0;
        i --;
        im1 = i - 1;
        ip1 = i + 1;
        ip2 = i + 2;
        dtfaz = dtwas - dtfaz;
        dtwas = dtfaz;
        t[ip1 - 1] = tf - sign_d(1.0E-4, tt- tf);
        
        msg_debug("2: Entering function Asmblo from function Asmtwo");
        Asmblo(ip2);
        msg_debug("2: Exit function Asmblo to function Asmtwo");
        
        msg_debug("2: Entering function Tblo from function Asmtwo");
        Tblo(ip2);
        msg_debug("2: Exit function Tblo to function Asmtwo");
        
        gbold = dx[i - 1];
        v = 0.0;
        if(gup == -1.0E10) continue;
        gold = gup;
        vup = xfb[i - 1];
        if(vup != 0.0)
        {
          gbold = gup;
          v = 0.0;
          xfb[i - 1] = 0.0;
         
        	
         }
         else
         {
           v = gbold - gup;          	
         	
         }                
    }   
          
   	
};



int Qsoiltemp1::Phase(int_ln& i, double& rlw, double cu, double rsu, double cd,
                                double rsd, double& tt, double& td,  double& dxx, 
                                double& gold,double& g,  double& t1,  double& t2)
{
//-----------------------
// History:
// March 21, 2007
//-----------------------
// Usuage:
// Secant solution for phase plane or dry node equation
// uses time average heat fluxes for phase element
// formulation is ordinary implicit if flux = 0 on one side of plane
// also if g is closer than star to element boundary
//-------------------------	
// Subprocedure called:
//---------------------
// Input variables
//-------------------------
// retrun value
// integer
//---------------------------
	double dbl_;
   double star = 5.0E-2 * dxx;
   double hold = dxx - gold;
   int_ln im1 = i - 1;
   int_ln ip1 = i + 1;
   int_ln ip2 = i + 2;
//   cout <<"dxx = "<<dxx<<"star = "<<star<<" hold = "<<hold<<" rlw = "<<rlw<<endl;
//   cin.get();
   double conu = conx[im1 - 1];
   if(dtfaz == 0.0) dtfaz = dtday * 86400.0;//to protect
   double capu = 0.5 * capx[im1 - 1] / dtfaz;
   double cond = conx[ip1 - 1];
   double capd = 0.5 * capx[ip1 - 1] / dtfaz;
   double rsu4 = 0.25 * rsu;
   double d1 = rsu4 / dtfaz;
   double c1 = conu * (1.0 - s[im1 - 1]) + capu + d1 * gold;
   double b1 = d1 * tt;
   double a1 = (capu + d1 * gold) * tt + 1.00001 * conu * e[im1 - 1] + ht[i - 1];
//   cout <<a1<<" "<<b1<<" "<<c1<<" "<<d1<<" "<<sflux<<" i = "<<i<<endl;
//   cin.get();
   if(i == is) a1 += sflux;
   double cutf = cu * tf;   
   double a1p = a1 - tf * c1;
   double b1p = b1 - tf * d1;
   
   double rsd4 = 0.25 * rsd;
   double d2 = rsd4 / dtfaz;
   double c2 = cond * (1.0 - s[ip2 - 1]) + capd + d2 * hold;
   double b2 = d2 * td;
   double a2 = (capd + d2 * hold) * td + cond * e[ip2 - 1] + ht[ip1 - 1];
   double cdtf = cd * tf;
   double a2p = a2 - tf * c2;
   double b2p = b2 - tf * d2;
//   cout <<a1p<<" "<<b1p<<" "<<a2p<<" "<<b2p<<endl;
//   cin.get();
   if(a1p * cu < 0.0)
   {
   	rlw = -rlw;
   }
   else
   {

   	   if(a1p * cu == 0.0)
   	   {
   		   if(a2p * cd == 0.0) return 0;
   	       else
   	       {

   	          if(a2p * cd > 0.0) rlw = -rlw;
   	       }
       }
   }
//   cout <<"rlw = "<<rlw<<" cu = "<<cu<<endl;
//   cin.get();
   
   double rlwt = rlw / dtfaz; 
   double alf1 = cu * (a1p + b1p * gold);
   double gam1 = c1 + d1 * gold;
   double bet1 = cu + gam1 * gold;
   double alf2 = cd * (a2p + b2p * hold);
   double gam2 = c2 + d2 * hold;
   double bet2 = cd + gam2 * hold;
   dbl_ = rlwt * bet1 * bet2 + alf1 * gam2 - alf2 * gam1;
   if(dbl_ == 0.0) dbl_ = 0.01;
   double dg = (alf2 * bet1 + alf1 * bet2) / dbl_;
//   cout <<"dg = "<<setprecision(7)<<dg<<" dxx = "<<setprecision(7)<<dxx<<endl;
   if(dg >= dxx) dg = 0.99 * dxx;
   if(dg <= 0.0) dg = 0.1 * dxx;
   g = gold + 0.5 * dg;
//   cout <<"g = "<<g<<" gold = "<<gold<<" dg = "<<dg<<endl;
   double h = dxx - g;
   double ft = rlwt * (g - gold);
   dbl_ = g * (c1 + d1 * g) + cu;
   if(dbl_ == 0.0) dbl_ = 0.01;
   double flu = cu * (a1p + b1p * g) / dbl_;
   dbl_ = h * (c2 + d2 * h) + cd;
   if(dbl_ == 0.0) dbl_ = 0.01;
   double fld = cd * (a2p + b2p * h ) / dbl_;
   double fwas = ft - flu - fld;
   double gwas = g;
//   cout <<"flu = "<<flu<<" fld = "<<fld<<" ft = "<<ft<<endl;
   g = gold + dg;
//   cout <<"g = "<<g<<"fwas = "<<fwas<<" gwas = "<<gwas<<endl;

   double f;
   double err;
   double gnew;
   if(rlw == 0.0)
   {
//   no latent heat
//   uses ordinary implicit formulation   	
       err = 1.0E-5;
	   int_ln n;
	   double df;
       for(n = 0; n < 50; n ++)
       {

     	   h = dxx - g;
     	   flu = 0.0;
		   dbl_ = cu + (c1 + d1 * g) * g;
		   if(dbl_ == 0.0) dbl_ = 0.01;
     	   if(cu > 0.0) flu = cu * (a1p + b1p * g) / dbl_;
     	   fld = 0.0;
		   dbl_ = cd + (c2 + d2 * h) * h;
		   if(dbl_ == 0.0) dbl_ = 0.001;
     	   if(cd > 0.0) fld = cd * (a2p + b2p * h) / dbl_;
     	   f = -flu - fld;
     	   if(fabs(f) <= err) break;
     	   df = f - fwas;
     	   if(df == 0.0) break;
     	   gnew = g - (g - gwas) * f / df;
     	   gwas = g;
     	   fwas = f;
     	   g = gnew;     	
      }
	  if(n == 50) return 0;
     
      	
   }
   else
   {
     err = 1.0E-5;
     int_ln n;
     for(n = 0; n < 50; n ++)     
     {

      	 h = dxx - g;
     	 ft = rlwt * (g - gold);
		 dbl_ = g * (c1 + d1 * g) + cu;
		 if(dbl_ == 0.0) dbl_ = 3.01;

     	   flu = cu * (a1p + b1p * g) / dbl_;
		   dbl_ = h * (c2 + d2 * h) + cd;
		   if(dbl_ == 0.0) dbl_ = 3.001;
     	   fld = cd * (a2p + b2p * h) / dbl_;
     	   if(gold > star && hold > star) 
     	   {
			   if(gold == 0.0) gold = 0.01;
			   if(hold == 0.0) hold = 0.01;
     	     flu = 0.5 * (theta * flu + theta1 * cu * (tt - tf) / gold);
     	     fld = 0.5 * (theta * fld + theta1 * cd * (td - tf) / hold);
         }
     	   f = ft - flu - fld;
         if(fabs(f) <= err)break;
         double df = f - fwas;
         if(df == 0.0) break;
         gnew = g - (g - gwas) * f / df;
         gwas = g;
         fwas = f;
         g = gnew;        
     	
     }
	 if(n == 50)return 0;

     
 
     	
   }
   dbl_ = g * (c1 + d1 * g) + cu;
   if(dbl_ == 0.0) dbl_ = 0.01; // to protect
   t1 = (g * (a1 + b1 * g) + cutf) / dbl_;
   dbl_ = h * (c2 + d2 * h) + cd;
   if(dbl_ == 0.0) dbl_ = 0.001; // to protect
   t2 = (h * (a2 + b2 * h) + cdtf) / dbl_;
   if(g >= 0.0) 
   {
     if(g <= dxx) return 0;
//   out at lower side
     dbl_ = g - gold;
     if(dbl_ == 0.0) 
	     dtfaz = dtday * 86400.0;
     else
         dtfaz *= hold / dbl_;
     t2 = td;
     return 0;   	
   }
   dbl_ = gold - g;
   if(dbl_ == 0.0) 
	   dtfaz = dtday * 86400.0;
   else

// out ot upper side
       dtfaz *= gold / dbl_;
   t1 = tt;
   return 0;   
          
   
      	
	
};

int Qsoiltemp1::Sethet(const int_ln& maxm, int kint[2], double heat[], double source[], const int& cmnt)
{
//----------------------------------------------
// History:
// April 3, 2007
//----------------------------------------------
//Usage:
//Sets internal point heat source functions	 
// Sets kint = 0 if source is inactive
// 
// individual sources can be turned off independently: kint(iheat) = 0
// for each successive source from top to bottom
// reads kint, vdep, source depth and time variation type
// Then reads appropriate strength and time data
//--------------------------------------------------------------------
// subprocedure called:
//---------------------------------------------------
// input variables
//-------------------------------------------------
// return value
// integer
//-------------------------------------------------
    double dephet;
    int_ln kinti;
	int iheat;
	


	kint[0] = 1;
	kint[1] = 1;
	iheat = 0;
	iheat ++;

	if(kswitch == 0)
	{
		kinti = KINT[cmnt];
		vdep = VDEP1[cmnt];
		dephet = DEPHET[cmnt];

		my_kinti = kinti;
		my_vdep2 = vdep;
		my_dephet = dephet;

	}
	else
	{
		kinti = my_kinti;
		vdep = my_vdep2;
		dephet = my_dephet;
	}
//        cout <<"kinti = "<<kinti<<" in Sethet"<<endl;
	if( kinti == 99999 ) return 0;

    dephet = vdep * dephet + x[ig - 1];   
    kint[iheat - 1] = kinti;
    for(int_ln i = ig - 1; i < imax1; i ++)
	{
		xhet[i] = 999.9;
		if(x[i] <= dephet && dephet <= x[i + 1]) xhet[i] = dephet;


	}

    msg_debug("Entering function Data from Sethet");
    Data(source); 
    msg_debug("Exit function Data to Sethet");

    for(int_ln m = 0; m < mmax; m ++)
    {
        heat[maxm * (iheat - 1) + m] = source[m];	
    }
	return 0;




};




	


//remeber to test the const float pointer, Feb/20/2007,
int Qsoiltemp1::Intrpl(const double X[], const double Y[], const int_ln& NMAX,
            const int& np, double xx[], double yy[], 
            const int_ln& NNMAX)
//---------------------------------------------------------------------------
// History:
// J. Tang adpated from f77 code by Goodrich
//---------------------------------------------------------------------------
// Usage
// procedure for interpolation:
// NP <= - 99     : linear interpolation
// NP >= 99       : ordinary spline interpolation
// NP = otherwise : smoooth spline interpolation, np ctrl the degree of smooth
//----------------------------------------------------------------------------
// subprocedure called: 
// int Cubsmh()
// int Cubfit()
// int Cubval()
//---------------------------------------------------------------------------- 
// input variables:
// (source points) X, (values)Y at X, (target points)xx
// (decimal precision) np,
// (log file object) flog;
//----------------------------------------------------------------------------
// output variables:
// yy 
//----------------------------------------------------------------------------
{
//  Temporary variables
    double XD[100], YD[100], weight[1000], work[3000];
    double  C[4 * 100];
    double XXD[NNMAX0],	YYD[NNMAX0];
    double sv, prec, sc, conl, conr;
//  Exit when bad input encountered
    if(NMAX > NMAX0 || NNMAX > NNMAX0)
    {

   		flog1 <<"Bad input of NMAX or NNMAX in Intrpl"<<endl;
//    	cout <<"Bad input of NMAX or NNMAX in Intrpl"<<endl;
    	exit(-1);
    }
    xx[NNMAX] = 1.0e10;	
//    if(np <= -99)
//    {
    	// linear interpolation
    	int m = 1;
    	double x1, y1, x2, slope; 
	    for(int_ln n = 2; n < NMAX + 1; n ++)
	    {
	      y1 = Y[n - 2];
	      x1 = X[n - 2];
	      x2 = X[n - 1];
	      slope = (Y[n - 1] - y1) / (x2 - x1);
	      int lg = 1;
	      while(lg)
	      {
	         yy[m - 1] = y1 + slope * (xx[m - 1] - x1);
	         if(xx[m - 1] > x2) break;
           m ++ ;
           if(m > NNMAX)lg = 0; 
	      }
	    }

	    return 0;
//	}
    //input arrays to double precision for spline routines

	for(int_ln n = 0; n < NMAX; n ++)
	{
	  XD[n] = X[n];

	  YD[n] = Y[n];
	}

    for(int_ln nn = 0; nn < NNMAX; nn ++)
    XXD[nn] = xx[nn];
    
	  int IERR, IFAIL;
    if(np < 99)
    {
    	//smooth spline interpolation
	
    	double stndev = 0.5e0 / (sqrt(3.0e0) * pow(10.,double(np))); 
    	double var = stndev * stndev;
    	sv = var * (NMAX + sqrt(2.0e0 * NMAX));
      prec = 0.5e-10;
		  int lim = 10;
		
		  IERR = 1;
      while(IERR > 0)
      {
           	
      	   IERR = Cubsmh(XD, YD, weight, NMAX, 1, sv, prec, lim, work, C, 4, sc);
      
	         if(IERR < 0)	
      	   {      	      
      			  flog1 <<" "<<"SUBROUTINE INTRPL ... SPLINEFIT FAILED,IERR<0"<<endl;
      			  cout << " Failed in Intrpl ... spline interpolation, IERR < 0"<<endl;  
//      			flog <<" Failed in Intrpl ... spline interpolation, IERR < 0"<<endl;
      			  exit(-1);
      	   }	
      		
      	   if(IERR == 0)
      	   {
      		   IFAIL = Cubval(0,XD,NMAX,C,4,XXD,YYD,NNMAX);
      		   if( IFAIL >= 0)
      		   {
      		   	 for(int_ln nn = 1 ; nn < NNMAX + 1; nn ++)
 		   	       yy[nn - 1] = YYD[nn - 1];
 			         return 0;
 		         } 			    
 		         else
 		         {
 			          cout <<"Intrpl ... Cubval failed"<<endl;
// 			flog <<"Intrpl ... Cubval failed"<<endl;
 			          exit(-1);
 		         }      			      				
      	    }
          prec *= 10.0e0;
      		
       }	
    }
//  Ordinary spline interpolation 
	   
//    flog <<'Ordinary spline interpolation   np = '<<np<<endl;  				
	   conl = (YD[1] - YD[0]) / (XD[1] - XD[0]);
	   conr = (YD[NMAX - 1] - YD[NMAX - 2]) / (XD[NMAX - 1] - XD[NMAX - 2]);
	   msg_debug("Entering function Cubfit from Intrpl");
     IERR = Cubfit(XD, YD, NMAX, 1, conl, 1, conr, C, 4);    	
     msg_debug("Exit function Cubfit to Intrpl");
     
     msg_debug("Entering function Cubval from Intrpl");
     IFAIL = Cubval(0, XD, NMAX, C, 4, XXD, YYD, NNMAX);
     msg_debug("Exit function Cubval to Intrpl");
    if( IFAIL >= 0)
    {
      	for(int_ln nn = 1 ; nn < NNMAX + 1; nn ++)
 		    yy[nn - 1] = YYD[nn - 1];
 		    return 0;
    } 			    
    else
    {
 	    flog1 <<" "<<"SUBROUTINE INTRPL ... CUBVAL FAILED"<<endl;
 	    cout <<"Intrpl ... Cubval failed"<<endl;
// 	flog <<"Intrpl ... Cubval failed"<<endl;
 	    exit(-1);
	}
	return 0;           	
};


int Qsoiltemp1::Cubsmh(const double X[], const double Y[], double W[], const int_ln& M, const int& IWT,
           const double& sval, const double& prec, const int& lim, double work[],
           double coeff[], const int& idimsn, double& scalc)
{
//------------------------------------------------------
//    History: 
//    J. Tang, adapted from F77 code by Goodrich
//    Feb. 21, 2007
//------------------------------------------------------
//    Usage:
//------------------------------------------------------	
//    Subprocedure called
//------------------------------------------------------
//    About the input variables
//    X, Y, W, are of size M
//    Caution: Problem exist with vecotr W, check before use
//    coeff is of size idimsn * M
//---------------------------------------------------------
//    output put variables
//    return value: IERR
//----------------------------------------------------------
      double XVALUE[2001],YVALUE[2001];
      int IERR = 0;                   // This is the return value
      double realmx = 8.0 * 16.0 * 16.0 - 1.0;
      if(idimsn < 4)
      {
      	return -4;
      }
      if( M <= 2)
      {
      	return -1;
      }	
      
      if(!(IWT == 0 || IWT == 1))
      {
      	return -2;
      }	
      if(sval <= 0.0)
      {
      	return -3;
      }	
      double accrcy = 0.5e-13;
      if( prec > 0.0 ) accrcy = prec;
      int kount = 15;
      if( lim > 0 ) kount = lim;
      int_ln M1 = 7 * M - 6;		
      int_ln M2 = M - 1;
	    work[0] = 0.0;
	    work[5] = 0.0;
	    work[7] = 0.0;
	    work[12] = 0.0;
	    work[M1 + 7] = 0.0;
	    work[M1 + 8] = 0.0;
	    work[M1 + 11] = 0.0;
	    work[M1 + 15] = 0.0;
      work[M1 + 18] = 0.0;
	    double p = 0.0;
	    double h = X[1] - X[0];
	    double xsign = h;
	    double g, gg;
	    double hh = 1.0 / h;
	    double f = (Y[1] - Y[0]) * hh;
	    double e;
      int_ln j = 8;
      for (int_ln i = 2; i < M2 + 1; i ++)
      {
      	j += 7;
      	g = h;
      	h = X[i] - X[i -1];
      	gg = hh;
      	hh = 1.0e0 / h;
      	e = f;
		    f = (Y[i] - Y[i - 1]) * hh;
		    coeff[i - 1] = f - e;
		    work[j - 1] = hh;
		    work[j] = - (gg + hh);
		    work[j + 1] = gg;
		    work[j + 2] = 2.0e0 * (g + h) / 3.0e0;
		    work[j + 3] = h / 3.0e0;
      }
      double wtm1 = 1.0e0;
      double wt = 1.0e0;
      double wtp1 = 1.0e0;
      if( IWT != 1)
      {
		     wt = 1.0e0 / W[0];
      	 wtp1 = 1.0e0 / W[1];
      }	
      j = 1;
      double r, r1, r2;
      for (int_ln i = 15; i < M1 + 1; i += 7)
      {
      	j ++;
      	if( IWT != 1) 
      	{
      		wtm1 = wt;
      	  wt = wtp1;
      	  wtp1 = 1.0e0 / W[j];
      	}
		    r = work[i - 1];
		    r1 = work[i];
      	coeff[M + j - 1] = wtp1 * r * r 
							  + wt * r1 * r1 + wtm1 * work[i + 1]
							  * work[i + 1];
		    coeff[2 * M + j - 1] = wtp1 * r * work[i + 7] + wt * r1 * work[i + 8];
      	coeff[3 * M + j - 1] = wtp1 * r * work[i + 15];                       
      }
      double f2;
      f2 = realmx;
      int_ln iter;
      for(iter = 1; iter < kount + 1; iter ++)
      {
      	j = 1;
      	for(int_ln i = 15; i < M1 + 1; i += 7)
      	{
      		j ++;
			    work[i - 7] = f *  work[i - 8];
			    work[i - 13] = g * work[i - 15];
			    r1 = work[i - 7];
			    r2 = work[i - 13];
			    work[i - 1] = 1.0e0 / (coeff[M + j -1] + p * work[i + 2] - f * r1 - g * r2);
			    work[i + 4] = coeff[j - 1] - r1 * work[i - 3] - r2 * work[i - 10];
			    f = coeff[2 * M + j - 1] + p * work[i + 3] - h * r1;
			    g = h;
			    h = coeff[3 * M + j - 1];
		    }
      	int_ln M3 = M1 + 15;
      	for(int_ln i = 15; i < M1 + 1; i += 7)
      	{
      		j = M3 - i;
			    work[j + 4] = work[j - 1] * work[j + 4] - work[j] * work[j + 11]
      		                 - work[j + 1] * work[j + 18];
      		                 
      	}
      	e = 0.0e0;
      	h = 0.0e0;
      	j = 1;
      	for(int_ln i = 8; i < M1 + 1; i += 7)
      	{
      		g = h;
			    h = (work[i + 11] - work[i + 4]) / (X[j] - X[j - 1]);
			    if(IWT == 0) wt = W[j - 1];
			    work[i + 5] = (h - g) / wt;
			    e += work[i + 5] * (h -g);
			    j ++;
		    }
		    if(IWT == 0) wt = W[M - 1];
		    g = - h / wt;
		    work[M1 + 12] = g;
		    e -= g * h;
		    g = f2;
		    f2 = e;
		    double darg = sval / f2;
		    darg = sqrt(darg);
		    if(darg >= 1.0e0 || f2 >= g)
		    {
			    if(darg < 1.0e0 || iter != 1)
			    {
				    darg = (1.0e0 - darg) / darg;
				    darg = fabs(darg);
				    if(darg > accrcy)
			      {
				      IERR ++;
			      }
      		}
          
      		scalc = f2;
			    f = p * work[12];
			    g = Y[0] - work[13];
			    j = 1;
			    for(int_ln i = 1; i < M2 + 1; i ++)
			    {
			      j += 7;
			      coeff[i - 1] = g;
			      gg = g;
			      coeff[2 * M + i - 1] = f;
			      e = f;
      		  h = X[i] - X[i - 1];
			      g = Y[i] - work[j + 12];
			      f = p * work[j + 11];
			      coeff[3 * M + i - 1]= (f - e) / (3.0e0 * h);
			      coeff[M + i - 1] = (g - gg) / h - (f + 2.0e0 * e) * h / 3.0e0;
			    }
			    coeff[M - 1] = g;
      		coeff[3 * M - 1] = f;
      		return IERR;
      		  
      	}	
      	
      	darg = (1.0e0 - darg) / darg;
      	if(darg <= accrcy)
        {
           scalc = f2;
		       f = p * work[12];
		       g = Y[0] - work[13];
      	   j = 1;		
      	   for(int_ln i = 1; i < M2 + 1; i ++)
      	   {
      		   j += 7;
			       coeff[i - 1] = g;
			       gg = g;
      		   coeff[2 * M + i - 1] = f;
      		   e = f;
			       h = X[i] - X[i - 1];
			       g = Y[i] - work[j + 12];
			       f = p * work[j + 11];
			       coeff[3 * M + i - 1]= (f - e) / (3.0e0 * h);
			       coeff[M + i - 1] = (g - gg) / h - (f + 2.0e0 * e) * h / 3.0e0;
		       }
		       coeff[M - 1] = g;
      	   coeff[3 * M - 1] = f;
      	   return IERR;
         }
         h = 0.0e0;
         for(int_ln i = 15; i < M1 + 1; i += 7)
         {
			      f = work[i - 4] * work[i - 3] + work[i + 2]
			        * work[i + 4] + work[i + 3] * work[i + 11];
			      g = f - work[i - 7] * work[i - 8] - work[i - 13] * work[i - 15];
			      f *= work[i + 4];
			      h = h - f + p * g * g * work[i - 1];
            work[i - 1] = g;    
         }
         if((h * xsign) >= 0.0e0)
         {
           IERR ++;
           scalc = f2;
	         f = p * work[12];
	         g = Y[0] - work[13];
	         j = 1;
	         for(int_ln i = 1; i < M2 + 1; i ++)
	         {
		          j += 7;
		          coeff[i - 1] = g;
		          gg = g;
      	      coeff[2 * M + i - 1] = f;
      	      e = f;
		          h = X[i] - X[i - 1];
		          g = Y[i] - work[j + 12];
		          f = p * work[j + 11];
		          coeff[3 * M + i - 1]= (f - e) / (3.0e0 * h);
		          coeff[M + i - 1] = (g - gg) / h - (f + 2.0e0 * e) * h / 3.0e0;
	         }
	         coeff[M - 1] = g;
           coeff[3 * M - 1] = f;
           return IERR;               	 
         }
         p -= f2 * darg / h;        
             	      			      	
         
      }
      IERR = 1;
      IERR ++;
      scalc = f2;
	    f = p * work[12];
	    g = Y[0] - work[13];
	    j = 1;
	    for(int_ln i = 1; i < M2 + 1; i ++)
	   {
		    j += 7;
		    coeff[i - 1] = g;
		    gg = g;
      	coeff[2 * M + i - 1] = f;
      	e = f;
		    h = X[i] - X[i - 1];
		    g = Y[i] - work[j + 12];
		    f = p * work[j + 11];
		    coeff[3 * M + i - 1]= (f - e) / (3.0e0 * h);
		    coeff[M + i - 1] = (g - gg) / h - (f + 2.0e0 * e) * h / 3.0e0;
	   }
	    coeff[M - 1] = g;
      coeff[3 * M - 1] = f;
      return IERR;            
};

int Qsoiltemp1::Cubfit(const double X[], const double Y[], const int_ln& M, const int_ln& Konstl, const double& constl,
		   const int_ln& Konstr, const double& constr, double coeff[], const int& idimsn)
{
//------------------------------------------------------
//      History:
//      J. Tang, adapted from F77 code by Goodrich
//      Feb. 21, 2007
//------------------------------------------------------
//      Usage:
//------------------------------------------------------
//     Subprocedure called:
//------------------------------------------------------
//     input variables:
//------------------------------------------------------
//     return value: ierr	
//------------------------------------------------------
	int ierr = 0;
	if(M <= 2)return -1;
	if(idimsn < 4)return -2;
	int_ln M2 = M - 1;
	double h = X[1] - X[0];
	double f = (Y[1] - Y[0]) / h;
	switch(Konstl)
	{
		case 1:
		{
			coeff[M] = 0.5e0;
			coeff[M * 2] = 1.5e0 * (f - constl) / h;
		}
		break;
		case 2:
		{
			coeff[M] = 0.0e0;
			coeff[M * 2] = constl / 2.0e0;
		}
		break;
		case 3:
		{
			if(constl > 0.0e0 && constl <= 1.0e0)
			{
				coeff[M] = - constl;
				coeff[M * 2] = 0.0e0;
			}	
			else
			{
				return -4;
			}	 
		}
		break;
		default:
		return -3;			
	}
	double g, hh, e, gg;
	for(int_ln i = 2; i < M2 + 1; i ++)
	{
		g = h;
		h = X[i] - X[i - 1];
		hh = 2.0e0 * (h + g);
		e = f;
		f = (Y[i] - Y[i - 1]) / h;
		gg = hh - g * coeff[M + i - 2];
		coeff[M + i - 1] = h / gg;
		coeff[2 * M + i - 1] = (3.0e0 * (f - e)
								  - g * coeff[M * 2 + i - 2]) / gg;
					
	}
	switch(Konstr)
	{
		case 1:
		{
			coeff[M * 3 - 1] = (3.0e0 * (constr - f) / h
							 - coeff[2 * M + M2 - 1])
							 / (2.0e0 - coeff[M + M2 - 1]);
		}
		break;
		case 2:
		{
			coeff[M * 3 - 1] = constr / 2.0e0;
		}
		break;
		case 3:
		{
			if(constr > 0.0 && constr <= 1.0e0)
			{
				coeff[3 * M - 1] = constr * coeff[2 * M + M2 - 1]
				                     / (1.0e0 + constr * coeff[M + M2 - 1]);
			}
			else
			{
				return -6;
			}	
		}
		break;
		default:
		return -5;			
	}
	g = coeff[3 * M - 1];
	int_ln j;
	for(int_ln i = 1; i < M2 + 1; i ++)
	{
		j = M - i;
		h = X[j] - X[j - 1];
		coeff[j - 1] = Y[j - 1];
		f = coeff[2 * M + j - 1] - g * coeff[M + j - 1];
		coeff[3 * M + j - 1] = (g - f) / (3.0e0 * h);
		coeff[M + j - 1] = (Y[j] - Y[j - 1]) / h
							  - (h * coeff[M * 3 + j - 1] + f) * h;
		coeff[2 * M + j - 1] = f;
		g = f;
	}
	coeff[M - 1] = Y[M - 1];
	return ierr;	 	
}


int  Qsoiltemp1::Cubval(const int& nprime, const double X[], const int_ln& M, double coeff[], const int& idimsn, 
           double xvalue[], double yvalue[], const int_ln& numpnt)
{
//-------------------------------------------------------------------
//History:
// J.Tang, adpated from F77 code by Goodrich
// Feb. 21, 2007
//-------------------------------------------------------------------
// Usage:
//------------------------------------------------------------------
// subprocedure called:
//-----------------------------------------------------------------
// input variables:
//------------------------------------------------------------------
// return value: ifail
//------------------------------------------------------------------	
  int ifail = 0;          // this is the return value
  if(numpnt < 1) ifail = - 1;
  if(nprime < 0) ifail = - 2;
  if(nprime > 2) ifail = - 3;
  if(idimsn < 4) ifail = - 4;
  if(M < 2) ifail = - 5;
  if(ifail < 0) return ifail;
  int_ln low = 1;
  int_ln iupper = M;
  double xval, h;
  int_ln middle;
  for(int i = 1; i < numpnt + 1; i ++)
  {
  	xval = xvalue[i - 1];
	  if((xval - X[0]) * (xval - X[M - 1]) > 0.0)
	  {
		  ifail ++;
		  return ifail;
	  }
	  if((xval - X[0]) * (xval - X[low - 1]) <= 0.0) low = 1;
	  if((xval - X[M - 1]) * (xval - X[iupper - 1]) <= 0.0) iupper = M;
	  while(iupper != (low + 1))
	  {
	  	middle = (low + iupper) / 2;
		  if((xval - X[middle - 1]) * (xval - X[iupper - 1]) > 0.0)
		  {
			   iupper = middle;
		  }
		  else
		  {
			  low = middle;
		  }
	  }
	  h = xval - X[low - 1];
	  if(nprime > 0)
	  {
		 if(nprime > 1)
		 {
			  yvalue[i - 1] = 2.0e0 * (3.0e0 * coeff[3 * M + low - 1] * h + coeff[2 * M + low - 1]);
		 }
		 else
		 {
			  yvalue[i - 1] = (3.0e0 * coeff[3 * M + low - 1] * h + 2.0e0 * coeff[2 * M + low - 1])
			    *h + coeff[M + low - 1];
		 }
	 }
	 else
	 {
		  yvalue[i - 1] = ((coeff[3 * M + low - 1] * h + coeff[2 * M + low - 1]) * h + coeff[M + low - 1])
			  *h + coeff[low - 1];
   }
  	
  	
  }
  return ifail;	
}

int Qsoiltemp1::Surf(const int& ktop, const double& sflx1, const double& sflx2, const double& sflx3, double& tdrive)
{
//------------------------------------
// History:
// March 21, 2007
//--------------------------------
// Usuage:
// Evaluates new surface conditions according to B.C. type ktop
//---------------------------------
// subprocedure called:
//-------------------------------
// input variables
//--------------------------------
// return value
// integer
//-----------------------------------
  double onoff2, onoff3;
  
//  ofsnore <<"ktop="<<ktop<<endl;
  switch(ktop)
  {
    case 1:
    {
    	return 0;
    }
    break;	
    case 2:
    {
//    prescribed flux
       onoff2 = 0.0; 
       onoff3 = 0.0;
      
      if(time - first < per / 2.0)
       onoff2 = 1.0;
      else
       onoff3 = 1.0; 	
       
       sflux = sflx1 + onoff2 * sflx2 - onoff3 * sflx3;
       tdrive = tf + onoff2 - onoff3;
       return 0; 
    }
    break;
    case 3:
    {
//    linearized heat balance treatment
//    evaluate net solar radiation rad
//   dummy version indicates usauage only
//   actually not used
       double time0 = 0.0;
       double rad = 150.0 + 100.0 * cos(6.28 * (time - time0) / per);
       emiss = 0.70;	
//   Evaluate new surface heat flux
       double tairab = 273.15 + tair;
       double sigta3 = sigma * tairab * tairab * tairab;
	   double rlwfix = sigta3 * (emiss - 1.0) * tairab;
       sflux = rad + rlwfix;
//   evaluate heat trasnfer coefficient trub
       double turb = 6.0;
	   htop = turb + 4.0 * sigta3;
       return 0;          
    }
    break;
    default:
    {

    }
	break;
  }
  return 0;
	
};


int Qsoiltemp1::Tblo(const int_ln& i1)
{
//------------------------------------------
// History:
// March 22, 2007
//-------------------------------------------
// Usuage:	
// evaluates new node temperature between i1 & imax
// forward substitution (for asmblo)
// ---------------------------------------
// subroutine called:
//-------------------------------------
// input variables
//-----------------------------------
// return value
// integer
//---------------------------------
   if(i1 > imax) return 0;
   for(int_ln i = i1 - 1; i < imax; i ++)
   {
     t[i] = s[i] * t[i - 1] + e[i];	
   }
   return 0;




	
};




//determine sign of a given integer
int Qsoiltemp1::sign_t(const double& num)
{
	if(num > 0) return 1;
	if(num == 0) return 0;
	if(num < 0 ) return -1;
	return 0;
}

double Qsoiltemp1::sign_d(const double& a, const double& b)
{
   if(b >=0) return fabs(a);
   return -fabs(a);
   	
}

void Qsoiltemp1::msg_debug(const char* pstr)
{
  if(STM_DEBUG)cout <<endl<<endl<<pstr<<endl<<endl;	
  return;
};

int Qsoiltemp1::Getsnowecd(ofstream& rflog1)
{
	
    char ecd[80];
    cout << "Enter name of the snow (.ECD) data file with the"
	     << " parameter values:" << endl;
//cin >> ecd;
    fpara >> ecd;

    rflog1 << "Enter name of the snow (.ECD) data file with"
		 << " the parameter values:" << ecd << endl << endl;

	return Getsnowecd(ecd);
 
	
};
int Qsoiltemp1::Getsnowecd(const char* ecd)
{
    const int NUMVAR = 52;
    char dummy[NUMVAR][10];
    ifstream infile;
    int i;

    int vegid[NUMVEG], update[NUMVEG];
    char vegname[NUMVEG][31];

    infile.open(ecd, ios::in);
    if (!infile) 
	{

        cerr << endl << "Cannot open " << ecd << " for data input" << endl;
        exit(-1);
    }

    for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
    for (i = 1; i < NUMVEG+1; i++) 
	{
//    infile >> vegid[i] >> vegname[i] >> kc[i] >> ki[i] >> gva[i]
//	>> tmin[i] >> toptmin[i] >> toptmax[i] >> tmax[i] >> raq10a0[i]
//	>> raq10a1[i] >> raq10a2[i] >> raq10a3[i] >> kn1[i] >> labncon[i]
//	>> leafmxc[i] >> kleafc[i] >> sla[i] >> cov[i] >> fpcmax[i] >> update[i];
	    infile >>vegid[i] >> vegname[i]>> MAX[i] >> NST[i]
	           >> KALLYR[i]>> KNODES[i] >> KISO[i]>> KTEMP[i]
	           >> KSNOW[i] >> KENVL[i]>> KFLUX[i]
	           >> LISO[i]>> TISO[i]>> LMAX[i]>> VDEPTH[i]
	           >> DEPTEM1[i]>> DEPTEM2[i]>> DEPTEM3[i]
	           >> DEPTEM4[i]>> DEPTEM5[i]>> NDEPF[i]
	           >> VDEP[i]>> DEPFLX1[i]>> DEPFLX2[i]
	           >> DEPFLX3[i]>> DEPFLX4[i]
	           >> DEPFLX5[i]>> HLAT[i]>> TF[i]>> GFLUX[i]
	           >> CDSNOW[i]>> FIRST[i]
	           >> FINAL[i]>> PER[i]>> DTDAY[i]>> THETA[i]
	           >> TOP[i]>> IG[i]>> EPSMIN[i]
	           >> VSPACE[i]>> VDEN[i]>> KINT[i]>> VDEP1[i]
	           >> DEPHET[i]>> SNOFAL[i]
	           >> EPSSNO[i]>> CONVRT[i]>> ETAO[i]>> DENFAC[i]
	           >> FCMELT[i]>> DENMAX[i]>>update[i];
     }
     infile.close();
	 return 1;
};


int Qsoiltemp1::Getsoillecd(ofstream& rflog1)
{

    char ecd[80];
    cout << "Enter name of soil layer (.ECD) data file with"
  	     << " parameter values:" << endl;
// cin >> ecd;
    fpara >> ecd;

    rflog1 << "Enter name of soil layer (.ECD) data file with parameter"
  	  	   << " values:" << ecd << endl << endl;
	return Getsoillecd(ecd);
};
int Qsoiltemp1::Getsoillecd(const char* ecd)
{
    const int NUMVAR = 69;
    char dummy[NUMVAR][8];
    ifstream infile;
    int i;

    int lfvegid[NUMVEG], update[NUMVEG];
    char lfvegname[NUMVEG][31];

    infile.open(ecd, ios::in);

    if (!infile) 
	{
        cerr << "\nCannot open " << ecd << " for data input" << endl;
        exit(-1);
    }

    for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
    for (i = 1; i < NUMVEG+1; i++)
    {
//    infile >> lfvegid[i] >> lfvegname[i] >> minleaf[i]
//	>> aleaf[i] >> bleaf[i] >> cleaf[i] >> update[i];
        infile >> lfvegid[i] >> lfvegname[i]>> THICK1[i]>> DXA1[i]
               >> DXB1[i]>> MAT1[i]>> DENSE1[i]>> WATER1[i]
               >> VCOND1[i]>> VSPH1[i]>> COND1[i]>> SPHT1[i]
               >> CONDF1[i]>> SPHF1[i]>> THICK2[i]>> DXA2[i]
               >> DXB2[i]>> MAT2[i]>> DENSE2[i]>> WATER2[i]
			   >> VCOND2[i]>> VSPH2[i]>> COND2[i]>> SPHT2[i]
               >> CONDF2[i]>> SPHF2[i]>> THICK3[i]>> DXA3[i]
               >> DXB3[i]>> MAT3[i]>> DENSE3[i]>> WATER3[i]
               >> VCOND3[i]>> VSPH3[i]>> COND3[i]>> SPHT3[i]
               >> CONDF3[i]>> SPHF3[i]>> THICK4[i]>> DXA4[i]
               >> DXB4[i]>> MAT4[i]>> DENSE4[i]>> WATER4[i]
			   >> VCOND4[i]>> VSPH4[i]>> COND4[i]>> SPHT4[i]
			   >> CONDF4[i]>> SPHF4[i]>> THICK5[i]>> DXA5[i]
			   >> DXB5[i]>> MAT5[i]>> DENSE5[i]>> WATER5[i]
               >> VCOND5[i]>> VSPH5[i]>> COND5[i]>> SPHT5[i]
			   >> CONDF5[i]>> SPHF5[i]>> THICK6[i]>> DXA6[i]
			   >> DXB6[i]>> MAT6[i]>> DENSE6[i]
			   >> WATER6[i]>>update[i];
    }
    infile.close();
	return 1;
};
int Qsoiltemp1::Getsoiltecd(ofstream& rflog1)
{

    char ecd[80];
    cout << "Enter name of the soil temperature initial (.ECD) data"
	     << " file with the parameter values:" << endl;
//cin >> ecd;
    fpara >> ecd;

    rflog1 << "Enter name of the soil temperature initial (.ECD) data"
		   << " file with the parameter values:" << ecd << endl << endl;	
	return Getsoiltecd(ecd);
};
int Qsoiltemp1::Getsoiltecd(const char* ecd)
{
    const int NUMVAR = 56;
    char dummy[NUMVAR][10];
    ifstream infile;
    int i,j;

    int vegid[NUMVEG], update[NUMVEG];
    char vegname[NUMVEG][31];

    infile.open(ecd, ios::in);

    if (!infile) 
	{
        cerr << endl << "Cannot open " << ecd << " for data input" << endl;
        exit(-1);
    }

    for (i = 0; i < NUMVAR; i++) { infile >> dummy[i]; }
    for (i = 1; i < NUMVEG+1; i++) 
	{
	    infile >>vegid[i] >> vegname[i]>> INDEX[i]
			>> VDEPP[i]>> NP[i];
        for (j=0;j<25; j++)
            infile >> DEPTH[i][j] >> TEMP[i][j];
        infile >>update[i];
    }
    infile.close();
	return 1;
};


void Qsoiltemp1::Showsnowecd(const int& cmnt) {

 // clrscr();
    cout << endl << "             SNOW PARAMETERS FOR VEGETATION TYPE STANDS "
		     << endl << endl;
    cout <<endl;
    cout <<endl;

  printf("MAX     = %d     NST    = %d       ",
		 MAX[cmnt], NST[cmnt]);
  printf("KALLYR = %d     KNODES = %d	   KISO   = %d\n",
		 KALLYR[cmnt],  KNODES[cmnt], KISO[cmnt]);
  printf("KTEMP   = %5.1lf KSNOW  = %5.1lf KENVL  = %5.1lf KFLUX  = %5.1lf\n",
		 KTEMP[cmnt],   KSNOW[cmnt],
		 KENVL[cmnt],   KFLUX[cmnt]);
  printf("LISO    = %5.1lf TISO   = %5.1lf LMAX   = %d     VDEPTH = %5.1lf\n",
		 LISO[cmnt],    TISO[cmnt],  LMAX[cmnt],
		 VDEPTH[cmnt]);
  printf("DEPTEM1 = %5.1lf DEPTEM2= %5.1lf", DEPTEM1[cmnt],
		 DEPTEM2[cmnt]);
  printf(" DEPTEM3= %5.1lf DEPTEM4= %5.1lf DEPTEM5= %5.1lf\n",
	   DEPTEM3[cmnt], DEPTEM4[cmnt], DEPTEM5[cmnt]);
  printf("NDEPF   = %5.1lf VDEP   = %5.1lf",  NDEPF[cmnt],
		 VDEP[cmnt]);
  printf(" DEPFLX1= %5.1lf DEPFLX2= %5.1lf DEPFLX3= %5.1lf\n",
		 DEPFLX1[cmnt], DEPFLX2[cmnt], DEPFLX3[cmnt]);
  printf("DEPFLX4 = %5.1lf DEPFLX5= %5.1lf HLAT   = %5.1lf ",
		 DEPFLX4[cmnt], DEPFLX5[cmnt], HLAT[cmnt]);
  printf("TF     = %5.1lf GFLUX  = %5.1lf\n", TF[cmnt],
		 GFLUX[cmnt]);
  printf("CDSNOW  = %5.1lf FIRST  = %5.1lf FINAL  = %5.1lf",
		 CDSNOW[cmnt], FIRST[cmnt], FINAL[cmnt]);
  printf(" PER    = %5.1lf DTDAY  = %5.1lf\n", PER[cmnt],
		 DTDAY[cmnt]);
  printf("THETA   = %5.1lf TOP    = %5.1lf IG     = %5.1lf ",
		 THETA[cmnt],  TOP[cmnt],  IG[cmnt]);
  printf("EPSMIN = %5.1lf VSPACE = %5.1lf\n",
		 EPSMIN[cmnt],  VSPACE[cmnt]);
  printf("VDEN    = %5.1lf KINT   = %5.1lf VDEP1  = %5.1lf ",
		 VDEN[cmnt],  KINT[cmnt], VDEP1[cmnt]);
  printf("DEPHET = %5.1lf SNOFAL = %5d\n",
		  DEPHET[cmnt],  SNOFAL[cmnt]);
  printf("EPSSNO  = %5.1lf CONVRT = %5.1lf ", EPSSNO[cmnt],
		 CONVRT[cmnt]);
  printf("ETAO   = %5.1lf DENFAC = %5.1lf FCMELT = %5.1lf\n",
		 ETAO[cmnt],  DENFAC[cmnt],  FCMELT[cmnt]);
  printf("DENMAX  = %5.1lf \n",  DENMAX[cmnt]);
  cout << endl << "Press any key to continue . . ." << endl;
  cin.get();
 // while (getch() == '\0');

};

/* **************************************************************
************************************************************** */

void Qsoiltemp1::Showsoillecd(const int& cmnt) {
  // clrscr();
  cout << endl << "         PARAMETERS FOR SOIL LAYERS " << endl << endl;
  printf("THICK1 =%5.1f  DXA1 =%5.1f DXB1 =%5.1f MAT1 =%d", THICK1[cmnt],
		 DXA1[cmnt],  DXB1[cmnt],  MAT1[cmnt]);
  printf("    DENSE1=%5.1f WATER1= %5.3f\n",
		 DENSE1[cmnt], WATER1[cmnt]);
  printf("VCOND1 =%5.1f  VSPH1=%5.1f COND1=%5.1f SPHT1=%5.1f",
		 VCOND1[cmnt], VSPH1[cmnt],
		 COND1[cmnt],  SPHT1[cmnt]);
  printf(" CONDF1=%5.1f SPHF1 = %5.1f\n",  CONDF1[cmnt],
		 SPHF1[cmnt]);
  printf("THICK2 =%5.1f  DXA2 =%5.1f DXB2 =%5.1f MAT2 =%d",
		 THICK2[cmnt], DXA2[cmnt],
		 DXB2[cmnt],   MAT2[cmnt]);
  printf("    DENSE2=%5.1f WATER2= %5.3f\n", DENSE2[cmnt],
		 WATER2[cmnt]);
  printf("VCOND2 =%5.1f  VSPH2=%5.1f COND2=%5.1f SPHT2=%5.1f",
		 VCOND2[cmnt], VSPH2[cmnt],
		 COND2[cmnt],  SPHT2[cmnt]);
  printf(" CONDF2=%5.1f SPHF2 = %5.1f\n", CONDF2[cmnt],
  	 SPHF2[cmnt]);
  printf("THICK3 =%5.1f  DXA3 =%5.1f DXB3 =%5.1f MAT3 =%d",
		 THICK3[cmnt], DXA3[cmnt],
		 DXB3[cmnt],   MAT3[cmnt]);
  printf("    DENSE3=%5.1f WATER3= %5.3f\n",
		 DENSE3[cmnt],  WATER3[cmnt]);
  printf("VCOND3 =%5.1f  VSPH3=%5.1f COND3=%5.1f SPHT3=%5.1f",
		 VCOND3[cmnt], VSPH3[cmnt],
		 COND3[cmnt],  SPHT3[cmnt]);
  printf(" CONDF3=%5.1f SPHF3 = %5.1f\n",
		 CONDF3[cmnt], SPHF3[cmnt]);
  printf("THICK4 =%5.1f  DXA4 =%5.1f DXB4 =%5.1f MAT4 =%d",
		 THICK4[cmnt], DXA4[cmnt],
		 DXB4[cmnt],  MAT4[cmnt]);
  printf("    DENSE4=%5.1f WATER4= %5.3f\n", DENSE4[cmnt],
		 WATER4[cmnt]);
  printf("VCOND4 =%5.1f  VSPH4=%5.1f COND4=%5.1f SPHT4=%5.1f",
		 VCOND4[cmnt], VSPH4[cmnt],
		 COND4[cmnt],  SPHT4[cmnt]);
  printf(" CONDF4=%5.1f SPHF4 = %5.1f\n", CONDF4[cmnt],
		 SPHF4[cmnt]);
  printf("THICK5 =%5.1f  DXA5 =%5.1f DXB5 =%5.1f MAT5 =%d",
		 THICK2[cmnt], DXA2[cmnt],
		 DXB2[cmnt],   MAT2[cmnt]);
  printf("    DENSE5=%5.1f WATER5= %5.3f\n",
		 DENSE2[cmnt], WATER2[cmnt]);
  printf("VCOND5 =%5.1f  VSPH5=%5.1f COND5=%5.1f SPHT5=%5.1f",
		 VCOND2[cmnt], VSPH2[cmnt],
		 COND2[cmnt],  SPHT2[cmnt]);
  printf(" CONDF5=%5.1f SPHF5 = %5.1f\n",
		 CONDF2[cmnt], SPHF2[cmnt]);
  printf("THICK6 =%5.1f  DXA6 =%5.1f DXB6 =%5.1f MAT6 =%d",
		 THICK6[cmnt], DXA6[cmnt],
		 DXB6[cmnt],   MAT6[cmnt]);
  printf("    DENSE6=%5.1f WATER6= %5.3f\n",
		 DENSE6[cmnt], WATER6[cmnt]);
  cout << endl << "Press any key to continue . . ." << endl;
  cin.get();
//  while (getch() == '\0');

};	

void Qsoiltemp1::Showsoiltecd(const int& cmnt) {
  int j;
//  clrscr();
  cout << "INITIALIZATION FOR SOIL TEMPERATURES:" << endl;
  printf(" INDEX = %5.1f  VDEP = %5.1f  NP = %d \n",  INDEX[cmnt],
  		VDEPP[cmnt],  NP[cmnt]);
  printf("DEPTH     TEMP \n");
    for (j=0;j<15; j++)
     { cout << DEPTH[cmnt][j] << "    " << TEMP[cmnt][j]<< endl;
     }
  cout << endl << "Press any key to continue . . ." << endl;
//  while (getch() == '\0');
  cin.get();
  printf("DEPTH     TEMP \n");
    for (j=15;j<25; j++)
     { cout << DEPTH[cmnt][j] << "    " << TEMP[cmnt][j]<< endl;
    }
  cout << endl << "Press any key to continue . . ." << endl;
  cin.get();
//  while (getch() == '\0');


};


void Qsoiltemp1::Inputwrite(const int& cmnt)
{
	FILE* pSnoda;
	pSnoda = fopen("snoda_t.dat","w");
	if(pSnoda == NULL)
	{	
		cout <<"Failed to create file snoda.dat for output"<<endl;
		exit(-1);
	}
	
	fprintf(pSnoda, "     MAX     NST  KALLYR  KNODES    KISO   KTEMP   KSNOW   KENVL   KFLUX\n");	
	printf("     MAX     NST  KALLYR  KNODES    KISO   KTEMP   KSNOW   KENVL   KFLUX\n");	
	
	fprintf(pSnoda," %7d %7d %7d %7d %7d %7d %7d %7d %7d\n",MAX[cmnt], NST[cmnt],
	        KALLYR[cmnt], KNODES[cmnt], KISO[cmnt], KTEMP[cmnt], KSNOW[cmnt],
	        KENVL[cmnt], KFLUX[cmnt]);
	printf(" %7d %7d %7d %7d %7d %7d %7d %7d %7d\n",MAX[cmnt], NST[cmnt],
	        KALLYR[cmnt], KNODES[cmnt], KISO[cmnt], KTEMP[cmnt], KSNOW[cmnt],
	        KENVL[cmnt], KFLUX[cmnt]);	        
	        
	fprintf(pSnoda," LISO      TISO      TISO      TISO      TISO      TISO\n");
	printf(" LISO      TISO      TISO      TISO      TISO      TISO\n");	
	
	fprintf(pSnoda,"%5d", LISO[cmnt]);
	printf("%5d", LISO[cmnt]);	
	
	fprintf(pSnoda,"%10.3f                                         \n",TISO[cmnt]);
	printf("%10.3f                                         \n",TISO[cmnt]);
		
	fprintf(pSnoda," LMAX    VDEPTH    DEPTEM    DEPTEM    DEPTEM    DEPTEM    DEPTEM\n");
	printf(" LMAX    VDEPTH    DEPTEM    DEPTEM    DEPTEM    DEPTEM    DEPTEM\n");	
	
	fprintf(pSnoda,"%5d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               LMAX[cmnt], VDEPTH[cmnt], DEPTEM1[cmnt], DEPTEM2[cmnt], DEPTEM3[cmnt],
	               DEPTEM4[cmnt], DEPTEM5[cmnt]);        

	printf("%5d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               LMAX[cmnt], VDEPTH[cmnt], DEPTEM1[cmnt], DEPTEM2[cmnt], DEPTEM3[cmnt],
	               DEPTEM4[cmnt], DEPTEM5[cmnt]);
	               	               
	fprintf(pSnoda,"NDEPF      VDEP    DEPFLX    DEPFLX    DEPFLX    DEPFLX    DEPFLX\n");
	printf("NDEPF      VDEP    DEPFLX    DEPFLX    DEPFLX    DEPFLX    DEPFLX\n");
		
	fprintf(pSnoda,"%5d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               NDEPF[cmnt], VDEP[cmnt], DEPFLX1[cmnt], DEPFLX2[cmnt],
	               DEPFLX3[cmnt], DEPFLX4[cmnt], DEPFLX5[cmnt]);
	printf("%5d %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               NDEPF[cmnt], VDEP[cmnt], DEPFLX1[cmnt], DEPFLX2[cmnt],
	               DEPFLX3[cmnt], DEPFLX4[cmnt], DEPFLX5[cmnt]);
	               	               
	fprintf(pSnoda,"      HLAT        TF\n");   
	printf("      HLAT        TF\n");
		
	fprintf(pSnoda," %9.3f %9.3f\n",HLAT[cmnt], TF[cmnt]);           
	printf(" %9.3f %9.3f\n",HLAT[cmnt], TF[cmnt]);
		                
	fprintf(pSnoda,".... $GRID--\n");
	printf(".... $GRID--\n");
		
	fprintf(pSnoda,"     FIRST     FINAL       PER     DTDAY     THETA\n");
	printf("     FIRST     FINAL       PER     DTDAY     THETA\n");
		
	fprintf(pSnoda," %9.3f %9.3f %9.3f %9.3f %9.3f\n", FIRST[cmnt],
	              FINAL[cmnt], PER[cmnt], DTDAY[cmnt], THETA[cmnt]);
	printf(" %9.3f %9.3f %9.3f %9.3f %9.3f\n", FIRST[cmnt],
	              FINAL[cmnt], PER[cmnt], DTDAY[cmnt], THETA[cmnt]);
	              	              
	fprintf(pSnoda,"       TOP   IG    EPSMIN    VSPACE      VDEN\n");       
	printf("       TOP   IG    EPSMIN    VSPACE      VDEN\n");
		       
	fprintf(pSnoda," %9.3f %4d %9.2f %9.3f %9.3f\n",
	               TOP[cmnt], IG[cmnt], EPSMIN[cmnt],VSPACE[cmnt], VDEN[cmnt]);
	printf(" %9.3f %4d %9.2f %9.3f %9.3f\n",
	               TOP[cmnt], IG[cmnt], EPSMIN[cmnt],VSPACE[cmnt], VDEN[cmnt]);
	               	               
	fprintf(pSnoda,"     THICK       DXA       DXB  MAT     DENSE     WATER\n");
	printf("     THICK       DXA       DXB  MAT     DENSE     WATER\n");
	
	
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %4ld %9.2f %10.3f %25s\n",
	               THICK1[cmnt], DXA1[cmnt], DXB1[cmnt], MAT1[cmnt], DENSE1[cmnt], WATER1[cmnt]," ");            
	printf("%9.3f %9.3f %9.3f %4ld %9.2f %10.3f %25s\n",
	       THICK1[cmnt], DXA1[cmnt], DXB1[cmnt], MAT1[cmnt], DENSE1[cmnt], WATER1[cmnt]," "); 
	               	                  
	fprintf(pSnoda,"   vcond      vsph       cond   spht         cond   sphf\n");          
	printf("   vcond      vsph       cond   spht         cond   sphf\n");
		     
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               VCOND1[cmnt],VSPH1[cmnt], COND1[cmnt], SPHT1[cmnt], CONDF1[cmnt], SPHF1[cmnt]);
	printf("%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               VCOND1[cmnt],VSPH1[cmnt], COND1[cmnt], SPHT1[cmnt], CONDF1[cmnt], SPHF1[cmnt]);
	               	               
	fprintf(pSnoda,"     THICK       DXA       DXB  MAT     DENSE     WATER\n");
	printf("     THICK       DXA       DXB  MAT     DENSE     WATER\n");
		
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %4d %9.2f %10.3f %25s\n",
	               THICK2[cmnt], DXA2[cmnt], DXB2[cmnt], MAT2[cmnt], DENSE2[cmnt], WATER2[cmnt]," ");          
	printf("%9.3f %9.3f %9.3f %4d %9.2f %10.3f %25s\n",
	               THICK2[cmnt], DXA2[cmnt], DXB2[cmnt], MAT2[cmnt], DENSE2[cmnt], WATER2[cmnt]," ");
	               	                    
	fprintf(pSnoda,"   vcond      vsph       cond   spht         cond   sphf\n");   
	printf("   vcond      vsph       cond   spht         cond   sphf\n");
		            
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               VCOND2[cmnt],VSPH2[cmnt], COND2[cmnt], SPHT2[cmnt], CONDF2[cmnt], SPHF2[cmnt]);
	printf("%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               VCOND2[cmnt],VSPH2[cmnt], COND2[cmnt], SPHT2[cmnt], CONDF2[cmnt], SPHF2[cmnt]);
	               	               
	fprintf(pSnoda,"     THICK       DXA       DXB  MAT     DENSE     WATER\n");
	printf("     THICK       DXA       DXB  MAT     DENSE     WATER\n");
		
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %4d %9.2f %10.3f %25s\n",
	               THICK3[cmnt], DXA3[cmnt], DXB3[cmnt], MAT3[cmnt], DENSE3[cmnt], WATER3[cmnt], " ");   
	printf("%9.3f %9.3f %9.3f %4d %9.2f %10.3f %25s\n",
	               THICK3[cmnt], DXA3[cmnt], DXB3[cmnt], MAT3[cmnt], DENSE3[cmnt], WATER3[cmnt], " ");
	               	                           
	fprintf(pSnoda,"   vcond      vsph       cond   spht         cond   sphf\n");     
	printf("   vcond      vsph       cond   spht         cond   sphf\n");
		          
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               VCOND3[cmnt],VSPH3[cmnt], COND3[cmnt], SPHT3[cmnt], CONDF3[cmnt], SPHF3[cmnt]);
	printf("%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               VCOND3[cmnt],VSPH3[cmnt], COND3[cmnt], SPHT3[cmnt], CONDF3[cmnt], SPHF3[cmnt]);
	               	               
	fprintf(pSnoda,"     THICK       DXA       DXB  MAT     DENSE     WATER\n");
	printf("     THICK       DXA       DXB  MAT     DENSE     WATER\n");
	
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %4d %9.2f %10.3f %25s\n",
	               THICK4[cmnt], DXA4[cmnt], DXB4[cmnt], MAT4[cmnt], DENSE4[cmnt], WATER4[cmnt]," ");  
	printf("%9.3f %9.3f %9.3f %4d %9.2f %10.3f %25s\n",
	               THICK4[cmnt], DXA4[cmnt], DXB4[cmnt], MAT4[cmnt], DENSE4[cmnt], WATER4[cmnt]," ");
	               	                            
	fprintf(pSnoda,"   vcond      vsph       cond   spht         cond   sphf\n");        
	printf("   vcond      vsph       cond   spht         cond   sphf\n");
		       
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               VCOND4[cmnt],VSPH4[cmnt], COND4[cmnt], SPHT4[cmnt], CONDF4[cmnt], SPHF4[cmnt]);
	printf("%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               VCOND4[cmnt],VSPH4[cmnt], COND4[cmnt], SPHT4[cmnt], CONDF4[cmnt], SPHF4[cmnt]);
	               	               
	fprintf(pSnoda,"     THICK       DXA       DXB  MAT     DENSE     WATER\n");
	printf("     THICK       DXA       DXB  MAT     DENSE     WATER\n");
		
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %4d %9.2f %10.3f %25s\n",
	               THICK5[cmnt], DXA5[cmnt], DXB5[cmnt], MAT5[cmnt], DENSE5[cmnt], WATER5[cmnt], " ");       
	printf("%9.3f %9.3f %9.3f %4d %9.2f %10.3f %25s\n",
	               THICK5[cmnt], DXA5[cmnt], DXB5[cmnt], MAT5[cmnt], DENSE5[cmnt], WATER5[cmnt], " ");
	               	                       
	fprintf(pSnoda,"   vcond      vsph       cond   spht         cond   sphf\n");     
	printf("   vcond      vsph       cond   spht         cond   sphf\n");          
	
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               VCOND5[cmnt],VSPH5[cmnt], COND5[cmnt], SPHT5[cmnt], CONDF5[cmnt], SPHF5[cmnt]);
	printf("%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	               VCOND5[cmnt],VSPH5[cmnt], COND5[cmnt], SPHT5[cmnt], CONDF5[cmnt], SPHF5[cmnt]);
	               	               
	fprintf(pSnoda,"     THICK       DXA       DXB  MAT     DENSE     WATER\n");
	printf("     THICK       DXA       DXB  MAT     DENSE     WATER\n");
		
	fprintf(pSnoda,"%9.3f %9.3f %9.3f %4d %9.2f %10.3f %25s\n",
	               THICK6[cmnt], DXA6[cmnt], DXB6[cmnt], MAT6[cmnt], DENSE6[cmnt], WATER6[cmnt], " ");   
	printf("%9.3f %9.3f %9.3f %4d %9.2f %10.3f %25s\n",
	               THICK6[cmnt], DXA6[cmnt], DXB6[cmnt], MAT6[cmnt], DENSE6[cmnt], WATER6[cmnt], " ");
	               	                           

	fprintf(pSnoda,"....$BCTOP....\n");
	printf("....$BCTOP....\n");
		
	fprintf(pSnoda,"    1     ---- Prescribed Temperature\n");          
	printf("    1     ---- Prescribed Temperature\n");  
		     
	fprintf(pSnoda,"index     vtime   np\n");
	printf("index     vtime   np\n");
		
	fprintf(pSnoda,"    3     1.000  -99\n");
	printf("    3     1.000  -99\n");
		
	fprintf(pSnoda,"....$sethet ------ internal heat sources 0\n");
	printf("....$sethet ------ internal heat sources 0\n");
		
	fprintf(pSnoda," kint      vdep    dephet\n");
	printf(" kint      vdep    dephet\n");
		
	fprintf(pSnoda,"%5d %9.3f %9.3f\n",KINT[cmnt], VDEP[cmnt], DEPHET[cmnt]);
	printf("%5d %9.3f %9.3f\n",KINT[cmnt], VDEP[cmnt], DEPHET[cmnt]);
		
	fprintf(pSnoda,"....$BCBOT....\n");
	printf("....$BCBOT....\n");
		
	fprintf(pSnoda,"    2     ----- Constant heat flux\n");
	printf("    2     ----- Constant heat flux\n");
		
	fprintf(pSnoda,"Flux value\n");
	printf("Flux value\n");
		
	fprintf(pSnoda,"%10.3f\n", GFLUX[cmnt]);
	printf("%10.3f\n", GFLUX[cmnt]);
		
	fprintf(pSnoda,"....$SNOFAL....  -- SEASONAL SNOWCOVER\n");
	printf("....$SNOFAL....  -- SEASONAL SNOWCOVER\n");
		
	fprintf(pSnoda,"%5d               snow\n",SNOFAL[cmnt]);
	printf("%5d               snow\n",SNOFAL[cmnt]);
		
	fprintf(pSnoda,"    epssno    convrt      etao    denfac    fcmelt   denmax\n");
	printf("    epssno    convrt      etao    denfac    fcmelt   denmax\n");
		
	fprintf(pSnoda,"%10.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",EPSSNO[cmnt],
	              CONVRT[cmnt], ETAO[cmnt], DENFAC[cmnt], FCMELT[cmnt],
	              DENMAX[cmnt]);
	printf("%10.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",EPSSNO[cmnt],
	              CONVRT[cmnt], ETAO[cmnt], DENFAC[cmnt], FCMELT[cmnt],
	              DENMAX[cmnt]);
	              	              
	fprintf(pSnoda,"....$TINITL....\n");
	printf("....$TINITL....\n");
		
	fprintf(pSnoda,"index      vdep   np\n");        
	printf("index      vdep   np\n"); 
		      
	fprintf(pSnoda,"%5.1f %9.3f %4d    --- Linear interpolation\n",INDEX[cmnt], VDEPP[cmnt], NP[cmnt]);
	printf("%5.1f %9.3f %4d    --- Linear interpolation\n",INDEX[cmnt], VDEPP[cmnt], NP[cmnt]);
		
	fprintf(pSnoda,"     Depth      Temp\n");
	printf("     Depth      Temp\n");
		
	fprintf(pSnoda,"_______________________\n");
	printf("_______________________\n");
		

	for(int j = 0; j < 25; j ++)
	{  
		fprintf(pSnoda,"%9.2f %9.2f\n",DEPTH[cmnt][j],TEMP[cmnt][j]);
		printf("%9.2f %9.2f\n",DEPTH[cmnt][j],TEMP[cmnt][j]);		
	}
        fprintf(pSnoda,"%9.3f\n",CDSNOW[cmnt]);
        printf("%9.3f\n",CDSNOW[cmnt]);
	                 	               	               	               	               	                              
	fclose(pSnoda);	

	cout <<"boundary condition for temperature interpolation"<<endl;

	FILE* pFile;
	pFile = fopen("bound_t.dat","w");
	if(pFile == NULL)
	{
		cout <<"Failed to create file bound_t.dat "<<endl;
		exit(-1);
	}
	fprintf(pFile,"%9.2f %9.2f\n",0.0, airt19);
	printf("%9.2f %9.2f\n",0.0, airt19);
	fprintf(pFile,"%9.2f %9.2f\n",15.0, airt29);
	printf("%9.2f %9.2f\n",15.0, airt29);
	fprintf(pFile,"%9.2f %9.2f\n",30.0, airt39);
	printf("%9.2f %9.2f\n",30.0, airt39);
	fprintf(pFile,"%9.2f %9.2f\n",-100.0, 0.0);
	fclose(pFile);
	cout <<"snow depth data"<<endl;

	pFile = fopen("sndep_t.dat","w");
	if(pFile == NULL)
	{
		cout <<"Failed to create file sndep_t.dat"<<endl;
		exit(-1);
	}
	fprintf(pFile,"%9.2f %9.2f %9.2f\n",0.0, hsnow19, 250.0);
	printf("%9.2f %9.2f %9.2f\n",0.0, hsnow19, 250.0);
	fprintf(pFile,"%9.2f %9.2f %9.2f\n",15.0, hsnow29, 250.0);
	printf("%9.2f %9.2f %9.2f\n",15.0, hsnow29, 250.0);
	fprintf(pFile,"%9.2f %9.2f %9.2f\n",30.0, hsnow39, 250.0);
	printf("%9.2f %9.2f %9.2f\n",30.0,hsnow39, 250.0);
	fprintf(pFile,"%9.2f %9.2f %9.2f\n",-100.0, 0.0, 0.0);
	fclose(pFile);
	cin.get();
	return;




};
