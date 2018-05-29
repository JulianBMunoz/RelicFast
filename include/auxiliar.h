
#ifndef __AUX_H_INCLUDED__
#define __AUX_H_INCLUDED__

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"





//array functions:

double *allocate_1D_array(int n1);

int *allocate_1D_array_int(int n1);


double **allocate_2D_array(int n1, int n2);

void free_2D_array(double **array, int n1);


double ***allocate_3D_array(int n1, int n2, int n3);

void free_3D_array(double ***array, int n1,int n2);

void free_cosmology(Cosmology *cosmo);

void do_check(int var_to_check);


long find_value(long length, double xtab[], double x);

long find_value_reverse(long length, double xtab[], double x);

void reverse(double *list,int N);


//mathemtical functions:

double interpol(double data[], double xtab[],int length, double x);

double interpol_2D(double **data, double xtab[], int lengthx, double ytab[], int lengthy, double x, double y);

double nintegrate(double data[],double xtab[], int length);

int triangle(double l1,double l2,double l3);

double interpol_cubic(double x0, double dx, double *ytab, int Nx, double x);

void solve_syst(double  **A, double  *X, double  *B, int N);

double det(double **M, int N);


//for progressbar:
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
void printProgress (double percentage);


//H/H0 and its log derivative for LCDM only
double E_LCDM(Cosmology *cosmo, double z);
double dlogE_dz_LCDM(Cosmology *cosmo, double z);


//time as a function of z with 1 additional rho.
void find_time(Cosmology *cosmo, double *zlist_Fri, long Nz_Fri, double *tlist_Fri,
				 double zmin_EoS, double dz_EoS, long Nz_EoS, double *rholist_EoS);


//Compute growth factor (7.73 in Dodelson's book)
void find_growth(Cosmology *cosmo, double *z_growth_array, double *Dp_growth_array, long N_growth);

double getsigma_M(Cosmology *cosmo, double *k_transfer, double *transfer_function);
  //calculates sigma(M)
double getsigma_8(Cosmology *cosmo, double *k_transfer, double *transfer_function);
	//for R=8 Mpc/h. It's sigma_8=0.83 for Planck2015 best-fit

double getsigma_M_z(Cosmology *cosmo, double z, double *z_transfer,  double *k_transfer,  double **transfer_function);
//calculates sigma(M) @z.


double getsigma_8_z(Cosmology *cosmo, double z, double *z_transfer, double *k_transfer, double **transfer_function);
//calculates sigma8 @z.


void findeta(Cosmology *cosmo, double *etalist, double zmin_EoS, double dz_EoS, long Nz_EoS,
            double* rholist_nu1_EoS, double* rholist_nu2_EoS, double* rholist_extra_EoS);
  //finds superconformal time etalist, at z_EoS list. Given LCDM + neutrinos and extra species.

#endif
