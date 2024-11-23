// Peng should modify this file to couple the new ideas.

#include "lpj.h"

#define LAG_CONV (NDAYYEAR*0.5*M_1_PI)  /* conversion factor for oscillation 
lag from angular units to days (=365/(2*PI))*/

Real soiltemp(const Soil *soil, 
			  const Climbuf *climbuf, 
			  Real Tlastday, 
			  Real Tsurface, 
			  Real Tinitial[], 
			  Real theta_init[], 
			  Real theta_i_init[], 
			  Real hinitial[], 
			  int year,
			  Real Tsoil[])
{
	const int Nlayer = 6;

	const int Nz = 120;

	const int Nt = 24 * 2;	// the time step is half hour.

//	for tundra.
	Real dz[121] = {
	0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
	0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
	0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 
	0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 
	0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 
	0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
	0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
	0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 
	0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
	0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
	0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
	0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

	int layer[121] =  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
		2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
		5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

	Real multiple1[6] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};

	Real multiple2[6] = {12350.61, 7054.782, 15085.07, 37852.83, 17702.5, 22526.8};

	Real M = 0.018015;	// The molecular weight of water, kg/mol.
	Real g = 9.81;		// The gravitational acceleration, ms-2.
	Real R = 8.314;		// universal gas constant. J/mol/K.
	Real theta_r = 0.05;	// the residual water contents. L3L-3.
	Real theta_s = 0.535;	// the saturated water content. L3L-3.
	Real K_s = 0.0000032;	// the saturated hydraulic conductivity. LT-1, m/s.
	Real n = 1.48;
	Real m = 0.2;
	Real l = 0.5;			// empirical parameters.
	Real alpha = 1.11;	// empirical parameter. L-1, m-1.
	Real alpha_v = 0.66;
	Real rou_i = 931.0;			// density of ice. ML-3, 931 kg/m3.
	Real rou_w = 1000.0;			// density of liquid water. ML-3, 1000 kg/m3.
	Real Lf = 3.34 * 100000.0;	// the latent heat of freezing (J/kg, L2T-2).  
	Real C1 = 0.55;
	Real C2 = 0.8;
	Real C33 = 3.07;
	Real C44 = 0.13;
	Real C5 = 4.0;	// Wm-1K-1.
	Real F1 = 13.05;
	Real F2 = 1.06;
	Real G_wT = 10.0;	// Gradient gain factor.
	Real fc = 0.02;	// The mass fraction of clay in the soil (unitless).
	Real dt = 300.0;	// Half hour resolution, the unit is s.

	int i, j, k, ii, kk;

	double theta[Nz+1];		// volumetric liquid water content. L3L-3.
	double theta_i[Nz+1];	// volumetric ice content, L3L-3.

	double T_upperB[Nt+1];
	double hi[Nz+1];
	double theta_0[Nz+1][Nt+1];
	double theta_i0[Nz+1][Nt+1];

	double h[Nz+1][Nt+1];	// the pressure heat. L.
	double T[Nz+1][Nt+1];	// temperature. K.

	double hlast[Nz+1][Nt+1];
	double Tlast[Nz+1][Nt+1];

	double K_Lh[Nz+1];	// pressure head gradient. LT-1. for liquid flow.
	double K_LT[Nz+1];	// temperature gradient. L2T-1K-1. for liquid flow.
	double K_vh[Nz+1];	// pressure head gradient. LT-1. for vapor flow.
	double K_vT[Nz+1];	// temperature gradient. L2T-1K-1. for vapor flow.
	double lamda[Nz+1];
	double L0[Nz+1];	// the volumetric latent heat of vaporization of water (Jm-3, ML-1T-2).
	double Ca[Nz+1];	// the apparant volumetric heat capacoty, (Jm-3K-1, ML-1T-2K-1). 
	double C[Nz+1];		// C is the hydraulic (or soil moisture) capacity, T-1, in other words, is the slope of the retention curve.

	double Lw;		// The latent heat of vaporization of water. (J/kg).
	double Hr;		// Relative Humidity.
	double D;		// Vapor diffusitivy in soil. m2s-1.
	double Da;		// diffusivity of water vapor in air. m2s-1.
	double etha;	// enhancement factor.

	double a_ij1, a_ij2, b_ij1, b_ij2;

	double delta2 = 999.9;

	double theta_w;		// Volumetric liquid water contents.
	double rou_vs;		// the empirical formulation of water vapor density.
	double Se;			// sink term. T-1, accounting for root water uptake.
	double F;			// Account for the difference between the thermal conductivities of ice and water in soils.


	Real soilt25 = 0.0;

	Real a,b,temp_lag;

	if(year <= 1958)
	{
		if(soil->w[0]==0)
			return climbuf->temp[NDAYS-1];

		linreg(&a,&b,climbuf->temp,NDAYS);

		temp_lag=a+b*(NDAYS-1-soil->alag*LAG_CONV);

		return climbuf->atemp_mean+soil->amp*(temp_lag-climbuf->atemp_mean);
	}

	for(i = 0; i < Nz + 1; i++)
	{
		Tsoil[i] = 0.0;

		for(j = 0; j < Nt + 1; j++)
		{
			theta_0[i][j] = 0.0;
			theta_i0[i][j] = 0.0;

			h[i][j] = 0.0;
			T[i][j] = 0.0;

			hlast[i][j] = 0.0;
			Tlast[i][j] = 0.0;
		}
	}

	// Calculate the functions.

	// Initialization.
	for(i = 0; i < Nz + 1; i++)
	{
		hi[i] = 0.0;
	}

	for(i = 0; i <= Nz; i++)	// initial condition for the first time step.
	{
		theta[i] = theta_init[i];

		theta_i[i] = theta_i_init[i];

		T[i][0] = Tinitial[i];

		h[i][0] = hinitial[i];

		hlast[i][0] = hinitial[i];
	}

	for(k = 0; k <= Nt; k++)	// upper boundary condition.
	{
		T_upperB[k] = (273.15 + Tlastday)  + (Tsurface - Tlastday) * (double)k / 48.0;
	}

	for(k = 0; k < Nt; k++)	// Time loops.
	{
		for(i = 0; i < Nz; i++)
		{
			theta_0[i][k] = theta[i];

			theta_i0[i][k] = theta_i[i];
		}

		//	Upper boundary conditions.

		T[0][k] = T_upperB[k];

		Tlast[0][k] = T_upperB[k];

		h[0][k] = -fabs(pow((pow(( 0.35 - theta_r) / (theta_s - theta_r), -1.0/m) - 1.0), 1.0/n) / alpha);

		hlast[0][k] = -fabs(pow((pow(( 0.35 - theta_r) / (theta_s - theta_r), -1.0/m) - 1.0), 1.0/n) / alpha);

		for(i = 1; i < Nz; i++)
		{
			j = k;

			kk = 0;

			if(i > 1)
			{				
				T[i-1][j] = Tlast[i-1][j+1];	// use the value of new time level (i-1).

				h[i-1][j] = hlast[i-1][j+1];
			}

			while(1)
			{
				kk++;

				for(ii = i - 1; ii <= i + 1; ii++)
				{
					theta[ii] = theta_r + (theta_s - theta_r) * pow(1.0 + pow(fabs(alpha * h[ii][j]), n), -m);

					theta_w = theta[ii];

					C[ii] = 1.0;

					theta_i[ii] = 0.0;	// temporally set the ice content to 0, should change.

					Ca[ii] = fabs(C[ii] * Lf * Lf * rou_i / g / T[ii][j] + theta_w * 4200.0 + theta_i[ii] * 2100.0);

					F =	1.0 + F1 * pow(theta_i[ii], F2);	// checked.

					Se = pow(1.0 + pow(fabs(alpha * h[ii][j]), n), -m);	// checked.

					K_Lh[ii] = K_s * pow(Se, l) * pow(1.0 - pow(1.0 - pow(Se, 1.0/m), m), 2.0) / 2.5;

					K_LT[ii] = K_Lh[ii] * (fabs(h[ii][j]) * G_wT * (-0.00198 - 6.62 * 0.000001 * T[ii][j])) * multiple1[layer[i]];

					Hr = exp(-fabs(h[ii][j]) * M * g / R / T[ii][j]);

					rou_vs = exp(130.8868 - 6790.4985 / T[ii][j] - 6.02808 * 3.391 * log(T[ii][j]) - 0.0376 * (T[ii][j] - 257.0)) * 1000.0;	

					Da = 2.12 * 0.00001 * pow(T[ii][j] / 273.15, 2.0);

					D = Da * alpha_v * (theta_s - theta[ii]);

					K_vh[ii] = D / rou_w * rou_vs * M * g / R / T[ii][j] * Hr;	// checked. saturated flow.

					etha = 9.5 + 3.0 * theta[ii] / theta_s - 8.5 * exp(-pow(((1.0 + 2.6 / sqrt(fc)) * theta[ii] / theta_s), 4.0));

					K_vT[ii] = D / rou_w * rou_vs * Hr * etha * (6790.4985 / T[ii][j] / T[ii][j] - 6.02808 * 3.391 / T[ii][j] + 0.052) * multiple2[layer[i]];	// checked.

					lamda[ii] = C1 + C2 * (theta_w + F * theta_i[ii]) - (C1 - C44) * exp(-pow(C33 * (theta_w + F * theta_i[ii]), C5));	// The function has been checked, from the Hansson et al, 2004.

					Lw = 2.501 * 1000000.0 - 2369.2 * T[ii][j];	// checked, I think the equation in Saito et al, 2006 is not correct for its unit.

					L0[ii] = Lw * rou_w;	// checked, from Saito et al, 2006.
				}

				ii = i;

				a_ij1 = ((K_Lh[ii] + K_Lh[ii-1]) / 2.0 + (K_vh[ii] + K_vh[ii-1]) / 2.0) * (h[ii][j] - h[ii-1][j]) / dz[i] + (K_Lh[ii] + K_Lh[ii-1]) / 2.0 + ((K_LT[ii] + K_LT[ii-1]) / 2.0 + (K_vT[ii] + K_vT[ii-1]) / 2.0) * (T[ii][j] - T[ii-1][j]) / dz[i];

				a_ij2 = ((K_Lh[ii+1] + K_Lh[ii]) / 2.0 + (K_vh[ii+1] + K_vh[ii]) / 2.0) * (h[ii+1][j] - h[ii][j]) / dz[i] + (K_Lh[ii+1] + K_Lh[ii]) / 2.0 + ((K_LT[ii] + K_LT[ii+1]) / 2.0 + (K_vT[ii+1] + K_vT[ii]) / 2.0) * (T[ii+1][j] - T[ii][j]) / dz[i];	// use the previous t

				b_ij2 = (h[ii+1][j] - h[ii][j]) / dz[i] * ((4200.0 * (K_Lh[ii] + K_Lh[ii+1]) / 2.0 + 1.88 * (K_vh[ii+1] + K_vh[ii]) / 2.0) * (T[ii+1][j] + T[ii][j]) / 2.0 + L0[ii] * (K_vh[ii] + K_vh[ii+1]) / 2.0) + (T[ii+1][j] - T[ii][j]) / dz[i] * ((lamda[ii] + lamda[ii+1]) / 2.0 + L0[ii] * (K_vT[ii] + K_vT[ii+1]) / 2.0 + 1.88 * (K_vT[ii] + K_vT[ii+1]) / 2.0 * (T[ii+1][j] + T[ii][j]) / 2.0);

				b_ij1 = (h[ii][j] - h[ii-1][j]) / dz[i] * ((4200.0 * (K_Lh[ii] + K_Lh[ii-1]) / 2.0 + 1.88 * (K_vh[ii-1] + K_vh[ii]) / 2.0) * (T[ii-1][j] + T[ii][j]) / 2.0 + L0[ii] * (K_vh[ii] + K_vh[ii-1]) / 2.0) + (T[ii][j] - T[ii-1][j]) / dz[i] * ((lamda[ii] + lamda[ii-1]) / 2.0 + L0[ii] * (K_vT[ii] + K_vT[ii-1]) / 2.0 + 1.88 * (K_vT[ii] + K_vT[ii-1]) / 2.0 * (T[ii][j] + T[ii-1][j]) / 2.0);

				hlast[ii][j] = -fabs(h[ii][j] + ((a_ij2 - a_ij1) * (dt / dz[i]) - (theta[ii] - theta_0[ii][k])) / 4.0); // the water transport.

				Tlast[ii][j] = T[ii][j] + (dt * (b_ij2 - b_ij1) / dz[i] - (theta_i[ii] - theta_i0[ii][k]) * Lf * 931.0) / Ca[ii]; // the heat transport.

				delta2 = Tlast[i][j] - T[i][j];

				h[ii][j] = hlast[ii][j];
				T[ii][j] = Tlast[ii][j];

				h[Nz][j] = h[Nz-1][j];	// Bottom boundary conditions, suppose no flux exists.
				T[Nz][j] = T[Nz-1][j];

				if(kk == 1)
				{
					for(ii = i - 1; ii <= i + 1; ii++)
					{
						h[ii][j+1] = h[ii][j];

						T[ii][j+1] = T[ii][j];

						hlast[ii][j+1] = h[ii][j];

						Tlast[ii][j+1] = T[ii][j];

						if(fabs(h[ii][j]) < 1.0)
						{
							h[ii][j+1] = -1.01;
							hlast[ii][j+1] = -1.01;
						}
					}

					j = j + 1;
				}

				if(kk > 1)
				{
					if(fabs(delta2) < 0.05)
					{
						break;
					}
				}
			}
		}

		for(i = 0; i < Nz; i++)
		{	
			Tsoil[i] += (Tlast[i][k] - 273.15) / 48.0;
		}
	}

	for(i = 0; i < Nz; i++)
	{
		Tinitial[i] = T[i][k-1];
		hinitial[i] = h[i][k-1];

		theta_init[i] = theta[i];
		theta_i_init[i] = theta_i[i];
	}

	for(i = 1; i <= 25; i++)
	{
		soilt25 += Tsoil[i] / 25.0;
	}

	return soilt25;

}

