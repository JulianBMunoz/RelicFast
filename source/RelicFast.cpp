////////////////////////////////////////////////////////////////////////////////////////
//
//	Code to find the bias of DM haloes of different masses, at different redshifts--
//   --including scale-dependent effects from relics, such as neutrinos.
//	By Julian B Munoz (05/2018 @Harvard)
//
////////////////////////////////////////////////////////////////////////////////////////

#include "RelicFast.h"

int main(int argc, char** filenameinput){

    //we start by reading the ini file:
    if(filenameinput[1]==NULL){
        printf("Error, no input file specified. Use ./relicfast FILENAME \n");
        return 0;
    }

    double *parameter_values;//input parameter array.
    parameter_values = (double *) calloc(Ninput, sizeof(double));

    //we read the input file and save the parameter values
    int read_ini_check;
    read_ini_check = read_input_file(filenameinput[1], parameter_values);

    //and we prepare the cosmological parameters (read + derived)
    struct Cosmology *cosmo = (Cosmology*) malloc(sizeof(Cosmology));

    //here we get all parameters from the inputs read
    int cosmo_check;
    cosmo_check = prepare_cosmology(cosmo, parameter_values);

    //we check that everything is well defined
    if(read_ini_check==0){
        printf("Error reading .ini file \n");
    }
    if(cosmo_check==0){
        printf("Error making cosmology \n");
    }

    //we save Mhalos and/or z_collapses if there are more than one.
    int lengthname=200;
    FILE *fp;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));

    //we create a folder to store our results
    lengthname=sprintf(filename,"mkdir output/result-%d", cosmo->file_tag);
    system(filename);

    //and copy the parameter file for future reference.
    lengthname=sprintf(
        filename,
        "cp %s output/result-%d/NM_%d-Nz_%d-Nk_%d-%s", 
        filenameinput[1], 
        cosmo->file_tag,
	cosmo->N_Mhalo, 
        cosmo->N_zcoll, 
        cosmo->N_klong, 
        filenameinput[1]
    );

    system(filename);

    ///We run the boltzmann code
    double *zlist_transfer;
    zlist_transfer=allocate_1D_array(Nz_transfer);

    int boltzmann_check = boltzmann(cosmo, zlist_transfer); //run CLASS/CAMB.

    int collapse_check=0;
    int bias_check=0;
    int tests = tests_pre_run(cosmo); //to make sure everything is well defined.
    if(tests<=0){
        printf("Something is wrong in tests_pre_run \n");
    }

    int iM, iz;

    //we now import the axion background equation of state
    int j;
    int axion_N=5000; 

    cosmo->axion_N = &axion_N; 
    cosmo->axion_a = allocate_1D_array(2*axion_N);
    cosmo->axion_z = allocate_1D_array(2*axion_N);
    cosmo->axion_w = allocate_1D_array(2*axion_N);
    cosmo->axion_rho = allocate_1D_array(2*axion_N);
    cosmo->axion_p = allocate_1D_array(2*axion_N);

    double axion_a[2*axion_N];
    double axion_w[2*axion_N];
    double axion_rho[2*axion_N]; 
    double temp;    
    double axion_osc; 

    FILE *fp3=fopen("/Users/nicholasdeporzio/Downloads/axion_aosc.dat", "r"); 

    fscanf(fp3, "%le", &axion_osc); 
    cosmo->axion_osc = &axion_osc;
    printf("Axion a_{osc} = %le \n", axion_osc); 
    printf("Axion a_{osc} = %le \n", cosmo->axion_osc);  

    FILE *fp2=fopen("/Users/nicholasdeporzio/Downloads/axion_background.dat", "r"); 

    for(
        j=0; 
        fscanf(
            fp2, 
            "%le %le %le %le",
            &axion_a[j], 
            &axion_w[j], 
            &temp, 
            &axion_rho[j]
        )==4; 
        ++j
    ){
        printf("%le \t %le \t %le \n", axion_a[j], axion_w[j], axion_rho[j]);
        cosmo->axion_a[j]=axion_a[j]; 
        cosmo->axion_z[j]=((1./axion_a[j])-1.); 
        cosmo->axion_w[j]=axion_w[j]; 
        cosmo->axion_rho[j]=axion_rho[j];
        cosmo->axion_p[j]=(cosmo->axion_w[j])*(cosmo->axion_rho[j]);  
    }; 

    //Fill in values after a_osc
    double rho_osc; 
    double a_osc; 
    double z_osc, logz_osc; 
    double dz_late, dlogz_early; 
    int linear_idx=2*axion_N-1000; 

    for(j=0; j<(2*axion_N); ++j){
        if (j==4999){
            rho_osc = cosmo->axion_rho[j]; 
            a_osc = cosmo->axion_a[j]; 
            z_osc = (1./a_osc)-1.; 
            dlogz_early = (log10(1.)-log10(z_osc))/4000.; 
            printf("OSCILLATION BEGINS HERE\n"); 
        }; 
        if (j>=4999){
            if (j<linear_idx){
                cosmo->axion_w[j]=0.; 
                cosmo->axion_p[j]=0.;  
                cosmo->axion_z[j]=pow(10., log10(z_osc)+((j-4998)*dlogz_early)); 
                cosmo->axion_a[j]=1./(cosmo->axion_z[j] + 1.);
                cosmo->axion_rho[j]=rho_osc * pow((1+cosmo->axion_z[j])/(1+z_osc), 3.); 
            } 
            else{
                if (j==linear_idx){
                    dz_late = cosmo->axion_z[j-1]/1000.;
                    printf("SWITCHING TO LINEAR Z SPACING... dz = %le \n", dz_late); 
                }
                cosmo->axion_w[j]=0.; 
                cosmo->axion_p[j]=0.;  
                cosmo->axion_z[j]=cosmo->axion_z[linear_idx-1]-((j+1-linear_idx)*dz_late); 
                cosmo->axion_a[j]=1./(cosmo->axion_z[j] + 1.);
                cosmo->axion_rho[j]=rho_osc * pow((1+cosmo->axion_z[j])/(1+z_osc), 3.); 
            };  
            printf("%le \t %le \t %.10e \n", cosmo->axion_a[j], cosmo->axion_w[j], cosmo->axion_rho[j]); 
        };
    }; 

    //we now solve for the collapse and calculate the biases at each z
    for(iz=0;iz<cosmo->N_zcoll;iz++){
        printf("TEST1\n");
        cosmo->z_collapse = cosmo->z_collapse_array[iz];
        if(debug_mode>=0) printf("z_coll = %le \n", cosmo->z_collapse);

        printf("TEST2\n"); 

        collapse_check = collapse(cosmo, zlist_transfer); 
        //solve for the collapse, and save delta_crit as a function of delta_L

        printf("TEST3\n"); 

        bias_check = get_bias(cosmo, zlist_transfer);    
        //find Lagrangian and Eulerian biases from the saved files.

        printf("TEST4\n"); 
    }

    printf("TEST5\n"); 

    //this tells the user which of the things have been done.
    printf("Boltzmann? %d Collapse? %d Bias? %d \n",boltzmann_check,collapse_check,bias_check);

    ////////////////////////////////////////////////////////
    ////  we save the masses and redshifts in a file  /////
    //////////////////////////////////////////////////////

    lengthname=sprintf(
        filename,
        "output/result-%d/info_-NM_%d-Nz_%d-Nk_%d.dat",
        cosmo->file_tag, 
        cosmo->N_Mhalo, 
        cosmo->N_zcoll, 
        cosmo->N_klong
    );

    fp=fopen(filename,"w");

    if(print_headers==1){
        fprintf(fp,"log10(Mhalo/Msun) \t");
    }

    for(iM=0;iM<cosmo->N_Mhalo;iM++){
        fprintf(fp,"%.2f \t", log10(cosmo->Mhalo_array[iM]));
    }

    fprintf(fp,"\n");

    if(print_headers==1){
        fprintf(fp,"z_collapse \t");
    }
    for(iz=0;iz<cosmo->N_zcoll;iz++){
        fprintf(fp,"%.2f \t", cosmo->z_collapse_array[iz]);
    }

    fclose(fp);

    free(filename);
    free(parameter_values);

    free_cosmology(cosmo);

    return 0;
}
