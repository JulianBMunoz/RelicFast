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
		lengthname=sprintf(filename,"cp %s output/result-%d/NM_%d-Nz_%d-Nk_%d-%s", filenameinput[1], cosmo->file_tag,
		 										cosmo->N_Mhalo, cosmo->N_zcoll, cosmo->N_klong, filenameinput[1]);
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


//we now solve for the collapse and calculate the biases at each z
	for(iz=0;iz<cosmo->N_zcoll;iz++){

		cosmo->z_collapse = cosmo->z_collapse_array[iz];

		if(debug_mode>=0) printf("z_coll = %le \n", cosmo->z_collapse);

		collapse_check = collapse(cosmo, zlist_transfer); //solve for the collapse, and save delta_crit as a function of delta_L
		bias_check = get_bias(cosmo, zlist_transfer);	//find Lagrangian and Eulerian biases from the saved files.

	}




//this tells the user which of the things have been done.
	 printf("Boltzmann? %d Collapse? %d Bias? %d \n",boltzmann_check,collapse_check,bias_check);




////////////////////////////////////////////////////////
////  we save the masses and redshifts in a file  /////
//////////////////////////////////////////////////////

	lengthname=sprintf(filename,"output/result-%d/info_-NM_%d-Nz_%d-Nk_%d.dat",cosmo->file_tag, cosmo->N_Mhalo, cosmo->N_zcoll, cosmo->N_klong);

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
