#ifndef __COLL_H_INCLUDED__
#define __COLL_H_INCLUDED__


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "auxiliar.h"
#include "boltzmann_Calls.h"
#include "pressure.h"

#ifdef _OPENMP
#include <omp.h>    // For nested loops.
#endif

//\delta_crit as a function of \delta_long
int collapse(Cosmology *cosmo, double *zlist_transfer);

//these are the functions that find z_collapse as a function of the initial
//conditions: for nothing massive (nu or extra)
double find_z_collapse_nothing(
    Cosmology *cosmo, 
    double Ri, 
    double Rpi, 
    double delta_long, 
    double *zlist_transfer,
    double Tfm_klong, 
    double *transfer_gamma_klong, 
    double *transfer_nu_klong
);

//for one massive species (nu or extra)
double find_z_collapse_1nu(
    Cosmology *cosmo, 
    double Ri, 
    double Rpi, 
    double delta_long, 
    double *zlist_transfer,
    double Tfm_klong, 
    double *transfer_gamma_klong, 
    double *transfer_nu_klong, 
    double *transfer_nu_massive_klong,
    double zmin_EoS, 
    double dz_EoS, 
    long Nz_EoS, 
    double *rholist_EoS, 
    double *plist_EoS,
    long Nz_solution, 
    double *R_solution, 
    double *Mnu_solution
);

//for 2 massive species (nu1+nu2 OR nu1+extra)
double find_z_collapse_2nu(
    Cosmology *cosmo, 
    double Ri, 
    double Rpi, 
    double delta_long, 
    double *zlist_transfer,
    double Tfm_klong, 
    double *transfer_gamma_klong, 
    double *transfer_nu_klong, 
    double *transfer_nu1_klong, 
    double *transfer_nu2_klong,
    double zmin_EoS, 
    double dz_EoS, 
    long Nz_EoS, 
    double *rholist1_EoS, 
    double *plist1_EoS, 
    double *rholist2_EoS, 
    double *plist2_EoS,
    long Nz_solution, 
    double *R_solution, 
    double *Mnu1_solution, 
    double *Mnu2_solution
);

//for 2 massive neutrinos + extra
double find_z_collapse_3nu(
    Cosmology *cosmo, 
    double Ri, 
    double Rpi, 
    double delta_long, 
    double *zlist_transfer,
    double Tfm_klong, 
    double *transfer_gamma_klong, 
    double *transfer_nu_klong, 
    double *transfer_nu1_klong, 
    double *transfer_nu2_klong, 
    double *transfer_extra_klong,
    double zmin_EoS, 
    double dz_EoS, 
    long Nz_EoS, 
    double *rholist1_EoS, 
    double *plist1_EoS, 
    double *rholist2_EoS, 
    double *plist2_EoS, 
    double *rholist_extra_EoS, 
    double *plist_extra_EoS,
    long Nz_solution, 
    double *R_solution, 
    double *Mnu1_solution, 
    double *Mnu2_solution, 
    double *Mnu3_solution
);

//For massless neutrinos plus axion
double find_z_collapse_masslessnu_axion(
    Cosmology *cosmo,
    double Ri,
    double Rpi,
    double delta_long,
    double *zlist_transfer,
    double Tfm_klong,
    double *transfer_gamma_klong,
    double *transfer_nu_massless_klong,
    double *transfer_axion_klong,
    double zmin_EoS,
    double dz_EoS,
    long Nz_EoS,
    double *rholist_axion_EoS,
    double *plist_axion_EoS,
    long Nz_solution,
    double *R_solution,
    double *axion_solution
);

//for 2 massive neutrinos + extra + axion
double find_z_collapse_3nu_axion(
    Cosmology *cosmo, 
    double Ri, 
    double Rpi, 
    double delta_long, 
    double *zlist_transfer,
    double Tfm_klong, 
    double *transfer_gamma_klong, 
    double *transfer_nu_klong, 
    double *transfer_nu1_klong, 
    double *transfer_nu2_klong, 
    double *transfer_extra_klong,
    double *transfer_axion_klong,
    double zmin_EoS, 
    double dz_EoS, 
    long Nz_EoS, 
    double *rholist1_EoS, 
    double *plist1_EoS, 
    double *rholist2_EoS, 
    double *plist2_EoS, 
    double *rholist_extra_EoS, 
    double *plist_extra_EoS,
    double *rholist_axion_EoS,  
    double *plist_axion_EoS,
    long Nz_solution, 
    double *R_solution, 
    double *Mnu1_solution, 
    double *Mnu2_solution, 
    double *Mnu3_solution,
    double *axion_solution
);

//this gets the solution for R(t) and finds the Mnu(t) within Rhalo. One for
//each species (could be made more efficient by solving together, but it is 
//irrelevant for neutrinos, so we leave it as is)
int findmass_1nu(
    Cosmology *cosmo, 
    double mass_0, 
    double temp_0, 
    double zmin_EoS, 
    double dz_EoS, 
    long Nz_EoS,
    double* rho1list_EoS, 
    double* rho2list_EoS, 
    double* rho3list_EoS, 
    double* etalist,
    double logz_solution_min, 
    double dlogz_solution, 
    long Nz_solution, 
    double* R_solution, 
    double* Mnu_solution
);


#endif
