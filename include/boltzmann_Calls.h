#ifndef __BOLTZ_H_INCLUDED__
#define __BOLTZ_H_INCLUDED__


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include "common.h"
#include "auxiliar.h"

//tests
int tests_pre_run(Cosmology *cosmo);

//Read the transfer functions from the chosen Boltzmann code
int gettransfer_matter(Cosmology *cosmo, char *filename, double *kgrid, double *TF);

int gettransfer_gamma(Cosmology *cosmo, char *filename, double *kgrid, double *TF);

int gettransfer_nu_massless(Cosmology *cosmo, char *filename, double *kgrid, double *TF);

int gettransfer_nu1(Cosmology *cosmo, char *filename, double *kgrid, double *TF);

int gettransfer_nu2(Cosmology *cosmo, char *filename, double *kgrid, double *TF);

int gettransfer_extra(Cosmology *cosmo, char *filename, double *kgrid, double *TF);

int gettransfer_axion(Cosmology *cosmo, char *filename, double *kgrid, double *TF);

//save the transfer functions to then read them:
int save_transfers_camb(
    Cosmology *cosmo, 
    double zlist_transfer[], 
    double *kgrid,
    double **TFm, 
    double **TFgamma, 
    double **TFnu_massless, 
    double **TFnu1, 
    double **TFnu2, 
    double **TFextra
);

int save_transfers_axioncamb(
    Cosmology *cosmo,
    double zlist_transfer[],
    double *kgrid,
    double **TFm,
    double **TFgamma,
    double **TFnu_massless,
    double **TFnu1,
    double **TFnu2,
    double **TFextra,
    double **TFaxion
);

int save_transfers_class(
    Cosmology *cosmo, 
    double zlist_transfer[], 
    double *kgrid,
    double **TFm, 
    double **TFgamma, 
    double **TFnu_massless, 
    double **TFnu1, 
    double **TFnu2, 
    double **TFextra
);



//// Run Boltzmann over the redshift list zlist_transfer
int run_camb(Cosmology *cosmo, double zlist_transfer[]);

int run_axioncamb(Cosmology *cosmo, double zlist_transfer[]);

int run_class(Cosmology *cosmo, double zlist_transfer[]);

int boltzmann(Cosmology *cosmo, double *zlist_transfer); 
//this calls the Boltzmann solver

//// Read the input file, and store the values:
int read_input_file (char *filename, double parameter_values[Ninput]);

double read_param (std::string buffer, std::string parameter_name); 
//returns parameter value from buffer

int prepare_cosmology (Cosmology *cosmo, double *parameter_values);

#endif
