////////////////////////////////////////////////////////////////////////////////////////
//
//	Code to find the delta_crit necessary for halo of mass Mhalo
//	to collapse by some redshift z_collapse.
//	By Julian B Munoz (05/2018 @Harvard)
//
////////////////////////////////////////////////////////////////////////////////////////

#include "collapse.h"


#define option_initial_velocity 0 //	0 for sigma_dot/sigma, 1 for T_dot/T (@k_star=2/Mpc). They agree, but it's a good sanity check.


int collapse(Cosmology *cosmo, double *zlist_transfer){
//one of the main functions, computes \delta_crit at collapse for Mhalo at zcollapse.
//zlist_transfer is the redshift list we ran Boltzmann code for (since CLASS does not save zs, although CAMB does)


//list of k_longs for delta_long (to compute derivative)

	double *klong_list;//in 1/Mpc. Input and closest CAMB values.
	klong_list = allocate_1D_array(cosmo->N_klong);


	int i,j;



	int lengthname=200;
 	char *filename; //To open files.
 	filename=(char *)malloc((lengthname+1)*sizeof(char));
	FILE *fp;




//we stablish if for the light relic and Mhalo we have to do clustering or not. Threshold given by mnu_min_clustering @ Mhalo=Mhalo_min_clustering.
	int i_clustering_index;
	int do_clustering;



/////////////////////////////////////////////////////////////////////////////////////////
////    we READ the boltzmann output to get the transfer functions between z=zi and z=0.	////
///////////////////////////////////////////////////////////////////////////////////////

	double *k_transfer_array; // k for transfer functions (common for all) (NOT over h)

	double **transfer_matter,**transfer_nu_massless,**transfer_nu1,**transfer_gamma;// transfer functions (not over k^2), for matter, massless nu, and photons.


	k_transfer_array = allocate_1D_array(length_transfer);
	transfer_matter = allocate_2D_array(Nz_transfer,length_transfer);
	transfer_gamma = allocate_2D_array(Nz_transfer,length_transfer);
	transfer_nu_massless = allocate_2D_array(Nz_transfer,length_transfer);
	transfer_nu1 = allocate_2D_array(Nz_transfer,length_transfer);

	//have to create these arrays anyway, to read transfer files. They can be nonsense if mnu2=0 or no extra species, is ok.
	double **transfer_nu2;
	transfer_nu2 = allocate_2D_array(Nz_transfer,length_transfer);

	double **transfer_extra;
	transfer_extra = allocate_2D_array(Nz_transfer,length_transfer);




	int transfer_check;



	if(boltzmann_tag == _CLASS_){ //CLASS
		transfer_check=save_transfers_class(cosmo, zlist_transfer, k_transfer_array, transfer_matter, transfer_gamma,
						transfer_nu_massless, transfer_nu1, transfer_nu2, transfer_extra);
	}
	else if(boltzmann_tag == _CAMB_){ //CAMB
		transfer_check=save_transfers_camb(cosmo, zlist_transfer, k_transfer_array, transfer_matter, transfer_gamma,
						transfer_nu_massless, transfer_nu1, transfer_nu2, transfer_extra);
	}
	else{
		printf("ERROR: Select CLASS or CAMB \n");
		return 0;
	}
	do_check(transfer_check);

	const double kmin=k_transfer_array[0];
	const double kmax=k_transfer_array[length_transfer-1];	//	printf("kmin=%.1le, kmax=%.1le \n",kmin, kmax);




//find the closest value to k_input in k_transfer_array to avoid unnecesary interpolation.
	long jk;


	for(i=0;i<cosmo->N_klong;i++){
		jk=find_value(length_transfer, k_transfer_array, cosmo->klong_list_input[i]);
		klong_list[i]=k_transfer_array[jk];//	printf("j=%ld, out=%le, in=%le \n", jk, klong_list[i],cosmo->klong_list_input[i]);
	}


//We compute T_dot to have two different initial conditions for comparison.
	double *d_transfer_array_zi; //dT/dz. The matter(c+b) transfer function and k/h obtained from CAMB code
	d_transfer_array_zi = allocate_1D_array(length_transfer);

	const int index_zip1 = Nz_transfer-1;
	const int index_zip2 = Nz_transfer-3;
	const int index_zi = Nz_transfer-2;
	const double zip1=zlist_transfer[index_zip1]; //we save zs before and after zi
	const double zip2=zlist_transfer[index_zip2];


	for(j=0;j<length_transfer;j++){
		d_transfer_array_zi[j] = (transfer_matter[index_zip1][j]-transfer_matter[index_zip2][j])/(zip1-zip2); //dT/dz@zi
	}


	double z; //redshift






	double *transfer_array_z0; //T*k^2, k: The matter(c+b) transfer function and k/h obtained from CAMB code, and k (NOT over h)


	transfer_array_z0 = allocate_1D_array(length_transfer);
	if(boltzmann_tag == _CLASS_ ){//CLASS
		lengthname=sprintf(filename,"Boltzmann_%d/transfer_files_%d/z%d_tk.dat",boltzmann_tag,cosmo->file_tag, Nz_transfer); //zlist reversed, starts at 1
	}
	else {
		lengthname=sprintf(filename,"Boltzmann_%d/transfer_files_%d/_transfer_out_z%.3f",boltzmann_tag,cosmo->file_tag,0.); //We reuse the same filename variable name.
	}

	transfer_check=gettransfer_matter(cosmo, filename,  k_transfer_array, transfer_array_z0);

	do_check(transfer_check);

	double *transfer_array_zi, *transfer_array_zip1, *transfer_array_zip2; //for CDM+b, derivative
	double *transfer_array_gamma_zi, *transfer_array_nu_massless_zi, *transfer_array_nu1_zi; //for photons, massless neutrinos, and massive neutrino1

	transfer_array_zi = allocate_1D_array(length_transfer);
	transfer_array_gamma_zi = allocate_1D_array(length_transfer);
	transfer_array_nu_massless_zi = allocate_1D_array(length_transfer);
	transfer_array_nu1_zi = allocate_1D_array(length_transfer);

	double *transfer_array_nu2_zi, *transfer_array_extra_zi; //for mnu2 and extra species
	transfer_array_nu2_zi = allocate_1D_array(length_transfer);
	transfer_array_extra_zi = allocate_1D_array(length_transfer);

	transfer_array_zip1 = allocate_1D_array(length_transfer);
	transfer_array_zip2 = allocate_1D_array(length_transfer);

	if(boltzmann_tag == _CLASS_){//CLASS
		lengthname=sprintf(filename,"Boltzmann_%d/transfer_files_%d/z%d_tk.dat",boltzmann_tag,cosmo->file_tag, 1); //zlist reversed, starts at 1
	}
	else { //CAMB
		lengthname=sprintf(filename,"Boltzmann_%d/transfer_files_%d/_transfer_out_z%.3f",boltzmann_tag,cosmo->file_tag,zi); //We reuse the same filename variable name.
	}
	transfer_check=gettransfer_matter(cosmo, filename,  k_transfer_array, transfer_array_zi);
	do_check(transfer_check);

	//we have a different function to extract photon transfer function
	transfer_check=gettransfer_gamma(cosmo, filename,  k_transfer_array, transfer_array_gamma_zi);
	do_check(transfer_check);
	//same for massless neutrinos
	transfer_check=gettransfer_nu_massless(cosmo, filename,  k_transfer_array, transfer_array_nu_massless_zi);
	do_check(transfer_check);
//and for massive neutrino 1
	transfer_check=gettransfer_nu1(cosmo, filename,  k_transfer_array, transfer_array_nu1_zi);
	do_check(transfer_check);




	//and for massive neutrino 2
	if(cosmo->mnu2>0){
		transfer_check=gettransfer_nu2(cosmo, filename,  k_transfer_array, transfer_array_nu2_zi);
		do_check(transfer_check);
	}
	//and for extra species
	if(cosmo->Omega_extra > 0){
		transfer_check=gettransfer_extra(cosmo, filename,  k_transfer_array, transfer_array_extra_zi);
		do_check(transfer_check);
}



//for zip1 and zip2, to take a better derivative
	if(boltzmann_tag == _CLASS_){//CLASS
		lengthname=sprintf(filename,"Boltzmann_%d/transfer_files_%d/z%d_tk.dat",boltzmann_tag,cosmo->file_tag, index_zip1); //zlist reversed, starts at 1
	}
	else { //CAMB
		lengthname=sprintf(filename,"Boltzmann_%d/transfer_files_%d/_transfer_out_z%.3f",boltzmann_tag,cosmo->file_tag,zip1); //We reuse the same filename variable name.
	}
	transfer_check=gettransfer_matter(cosmo, filename,  k_transfer_array, transfer_array_zip1);
	do_check(transfer_check);

	if(boltzmann_tag == _CLASS_){//CLASS
		lengthname=sprintf(filename,"Boltzmann_%d/transfer_files_%d/z%d_tk.dat",boltzmann_tag,cosmo->file_tag, index_zip2); //zlist reversed, starts at 1
	}
	else { //CAMB
		lengthname=sprintf(filename,"Boltzmann_%d/transfer_files_%d/_transfer_out_z%.3f",boltzmann_tag,cosmo->file_tag,zip2); //We reuse the same filename variable name.
	}
	transfer_check=gettransfer_matter(cosmo, filename,  k_transfer_array, transfer_array_zip2);
	do_check(transfer_check);






//we also read the transfer function at zcollapse
	double *transfer_array_z_collapse;

	transfer_array_z_collapse = allocate_1D_array(length_transfer);

	const int number_collapse=find_value(cosmo->N_zcoll,cosmo->z_collapse_array,cosmo->z_collapse);
	const int j_collapse = cosmo->j_collapse_array[number_collapse];

	if (debug_mode>0){
		printf("number_collapse=%d, j_collapse=%d \n",number_collapse,j_collapse);
	}




	if(boltzmann_tag == 0){//CLASS
		lengthname=sprintf(filename,"Boltzmann_%d/transfer_files_%d/z%d_tk.dat",boltzmann_tag,cosmo->file_tag, j_collapse+1); //zlist reversed, starts at 1
	}
	else {
		lengthname=sprintf(filename,"Boltzmann_%d/transfer_files_%d/_transfer_out_z%.3f",boltzmann_tag,cosmo->file_tag,cosmo->z_collapse); //We reuse the same filename variable name.
	}
	transfer_check=gettransfer_matter(cosmo, filename, k_transfer_array, transfer_array_z_collapse);
	do_check(transfer_check);




/////////////////////////////////////////////////////////////////////////////
////    we create an interpolation table for the equations  of state	  ////
///////////////////////////////////////////////////////////////////////////

	const int Nz_EoS=1000; //how many redshifts we take for the EoS. 1000 is fine for z

//first for nu1
	double *zlist_EoS;
	double *plist_nu1_EoS, *rholist_nu1_EoS;
	double Temp_nu;

	zlist_EoS=allocate_1D_array(Nz_EoS);
	plist_nu1_EoS=allocate_1D_array(Nz_EoS);
	rholist_nu1_EoS=allocate_1D_array(Nz_EoS);

	const double zmin_EoS=0.;//we go down to z=0, needed to take ratios with Omega.
	const double dz_EoS=(zi-zmin_EoS)/(Nz_EoS-1);


	if(debug_mode > 0){
		lengthname=sprintf(filename,"tests/rho_nu1_%d.dat",cosmo->file_tag);
		fp=fopen(filename,"w");
		if(print_headers!=0){
			fprintf(fp, "z\t rho(z)/rho(0) w(z) \n");
		}
}

	for(i=0;i<Nz_EoS;i++){
		zlist_EoS[i]=zmin_EoS+i*dz_EoS; //z linearly from 0 to zi EQUALLY SPACED! Since we will use interpol_cubic.
		Temp_nu = cosmo->T0_nu * (1.+zlist_EoS[i]); //neutrino temperature in K
		rholist_nu1_EoS[i] = density_WDM(cosmo->mnu1,Temp_nu);
		plist_nu1_EoS[i] = pressure_WDM(cosmo->mnu1,Temp_nu);		//						printf("z=%.1le,  w=%.1le, Tnu/mnu=%.1le \n",zlist_EoS[i],plist_EoS[i]/rholist_EoS[i],Temp_nu*KtoeV/(mnu1));//
		if(debug_mode > 0){
			fprintf(fp, "%le %le %le \n",zlist_EoS[i], rholist_nu1_EoS[i]/rholist_nu1_EoS[0], plist_nu1_EoS[i]/rholist_nu1_EoS[i]);
		}
	}
	if(debug_mode > 0) fclose(fp);


	Temp_nu = cosmo->T0_nu * (1.+zi);
	double rho_nu1_zi=density_WDM(cosmo->mnu1,Temp_nu);
	Temp_nu = cosmo->T0_nu;
	double rho_nu1_z0=density_WDM(cosmo->mnu1,Temp_nu);

	double rho_nu1_ratio_zi=0.;
	if(cosmo->mnu1>0){
		rho_nu1_ratio_zi=rho_nu1_zi/rho_nu1_z0;
	}

	if(debug_mode > 0){
		printf("kmin=%.3le and kmax=%.3le \n",kmin,kmax);
	}


//now for nu2.
	double *plist_nu2_EoS, *rholist_nu2_EoS;
	double rho_nu2_zi,rho_nu2_z0,rho_nu2_ratio_zi=0;

	plist_nu2_EoS=allocate_1D_array(Nz_EoS);
	rholist_nu2_EoS=allocate_1D_array(Nz_EoS);


	if(debug_mode > 0){
		lengthname=sprintf(filename,"tests/rho_nu2_%d.dat",cosmo->file_tag);
		fp=fopen(filename,"w");
		if(print_headers!=0){
			fprintf(fp, "z\t rho(z)/rho(0) w(z) \n");
		}
	}

	if(cosmo->Omeganu2>0){
		for(i=0;i<Nz_EoS;i++){
			zlist_EoS[i]=zmin_EoS+i*dz_EoS; //z linearly from 0 to z_collapse EQUALLY SPACED! Since we will use interpol_cubic.
			Temp_nu = cosmo->T0_nu * (1.+zlist_EoS[i]); //neutrino temperature in K
			rholist_nu2_EoS[i] = density_WDM(cosmo->mnu2,Temp_nu);
			plist_nu2_EoS[i] = pressure_WDM(cosmo->mnu2,Temp_nu);		//						printf("z=%.1le,  w=%.1le, Tnu/mnu=%.1le \n",zlist_EoS[i],plist_EoS[i]/rholist_EoS[i],Temp_nu*KtoeV/(mnu1));

			if(debug_mode > 0){
				fprintf(fp, "%le %le %le \n",zlist_EoS[i], rholist_nu2_EoS[i]/rholist_nu2_EoS[0], plist_nu2_EoS[i]/rholist_nu2_EoS[i]);
			}
		}

		Temp_nu = cosmo->T0_nu * (1.+zi);
		rho_nu2_zi=density_WDM(cosmo->mnu2,Temp_nu);
		Temp_nu = cosmo->T0_nu;
		rho_nu2_z0=density_WDM(cosmo->mnu2,Temp_nu);

		rho_nu2_ratio_zi=rho_nu2_zi/rho_nu2_z0;

	}
	if(debug_mode > 0) fclose(fp);




//and for the extra species

	double *plist_extra_EoS, *rholist_extra_EoS;
	double rho_extra_zi,rho_extra_z0,rho_extra_ratio_zi=0;

	plist_extra_EoS=allocate_1D_array(Nz_EoS);
	rholist_extra_EoS=allocate_1D_array(Nz_EoS);

	double T0_extra, Temp_extra, m_extra;

	if(cosmo->tag_sterile_nu == _TRUE_){
		T0_extra=cosmo->T0_nu;
		m_extra=cosmo->m_SN;
	}
	else if(cosmo->tag_thermal_relic == _TRUE_){
		T0_extra=cosmo->T0_TR;
		m_extra=cosmo->m_TR;
	}
	else{
		m_extra=1.0;//just to make sure there are no divergencies.
		T0_extra=0.0;
	}

	if(debug_mode > 0){
		lengthname=sprintf(filename,"tests/rho_extra_%d.dat",cosmo->file_tag);
		fp=fopen(filename,"w");
		if(print_headers!=0){
			fprintf(fp, "z\t rho(z)/rho(0) w(z) \n");
		}
	}

	if(cosmo->Omega_extra>0){
		for(i=0;i<Nz_EoS;i++){
			zlist_EoS[i]=zmin_EoS+i*dz_EoS; //z linearly from 0 to z_collapse EQUALLY SPACED! Since we will use interpol_cubic.
			Temp_extra = cosmo->T0_extra * (1.+zlist_EoS[i]); //neutrino temperature in K
			rholist_extra_EoS[i] = density_WDM(cosmo->m_extra,Temp_extra);
			plist_extra_EoS[i] = pressure_WDM(cosmo->m_extra,Temp_extra);		//						printf("z=%.1le,  w=%.1le, Tnu/mnu=%.1le \n",zlist_EoS[i],plist_EoS[i]/rholist_EoS[i],Temp_nu*KtoeV/(mnu1));

			if(debug_mode > 0){
				fprintf(fp, "%le %le %le \n",zlist_EoS[i], rholist_extra_EoS[i]/rholist_extra_EoS[0], plist_extra_EoS[i]/rholist_extra_EoS[i]);
			}
		}

		Temp_extra = T0_extra * (1.+zi);
		rho_extra_zi=density_WDM(m_extra,Temp_extra);
		Temp_extra = T0_extra;
		rho_extra_z0=density_WDM(m_extra,Temp_extra);

		rho_extra_ratio_zi=rho_extra_zi/rho_extra_z0;

	}
	if(debug_mode > 0) fclose(fp);




	//if we ask for clustering do it only if at least one of the species is heavy enough to matter.
	if(do_clustering_tag>0){
		if((cosmo->mnu1 > mnu_min_clustering*sqrt(Mhalo_min_clustering/cosmo->Mhalo)) ){
			do_clustering=1;
		}
		else if((cosmo->m_extra > cosmo->T0_extra/cosmo->T0_nu * mnu_min_clustering*sqrt(Mhalo_min_clustering/cosmo->Mhalo))){
			do_clustering=1;
		}
		else{
			do_clustering=0;
		}
	}
	else{
		do_clustering=0;
	}


//we find the superconformal time eta for clustering. If required.
	double *etalist;
	etalist=allocate_1D_array(Nz_EoS); //we find eta in the same z range we find EoS, where it matters.

	if(do_clustering_tag>0){
		findeta(cosmo, etalist, zmin_EoS, dz_EoS, Nz_EoS,
		            rholist_nu1_EoS, rholist_nu2_EoS, rholist_extra_EoS);
	}//creates a cubic interpolator for eta at z_EoS list.




///////////////////////////////////////////////////////////////////////////////////
////    we set the R(z) array that we will use to find extra species clustering////
///////////////////////////////////////////////////////////////////////////////////

	const int Nz_solution=400;//number of z in solution of R(t). LOGSPACED.
	const double logz_solution_min=log(zi);
	const double logz_solution_max=log(cosmo->z_collapse/2);//as in collapse code.
	const double dlogz_solution=(logz_solution_max-logz_solution_min)/(Nz_solution-1);

	double *Rhalo_solution; //radius of halo at different zs.
	Rhalo_solution=allocate_1D_array(Nz_solution);
	double *Mnu1_solution; //neutrino 1 mass within R_halo, needed to find collapse.
	Mnu1_solution=allocate_1D_array(Nz_solution);
	double *Mnu2_solution; //neutrino 1 mass within R_halo, needed to find collapse.
	Mnu2_solution=allocate_1D_array(Nz_solution);
	double *Mextra_solution; //neutrino 1 mass within R_halo, needed to find collapse.
	Mextra_solution=allocate_1D_array(Nz_solution);

	//z array logspaced between zi and zcollapse/2, to match collapse code.




//////////////////////////////////////////////////////
////    we set initial conditions									////
/////////////////////////////////////////////////////


	double OmL_i, OmM_i, OmR_i, Omnu1_i, Omnu2_i, Omextra_i;
	OmL_i=cosmo->OmegaL;
	OmM_i=cosmo->OmegaM*pow(1.+zi,3.);
	OmR_i=cosmo->OmegaR*pow(1.+zi,4.);
	Omnu1_i=cosmo->Omeganu1*rho_nu1_ratio_zi;
	Omnu2_i=cosmo->Omeganu2*rho_nu2_ratio_zi;
	Omextra_i=cosmo->Omega_extra*rho_extra_ratio_zi;
	const double Hi=cosmo->H0_Mpc*sqrt(OmL_i + OmM_i + OmR_i +
				Omnu1_i + Omnu2_i + Omextra_i);//H(zi) in Mpc-1

	if(debug_mode>0){
		printf("Hi=%.2le \n",Hi);
		printf("OmM_i=%.2le, OmR_i=%.2le, Omnu1=%.2le \n",OmM_i,OmR_i,Omnu1_i);
		printf("Omnu2_i=%.2le, Omextra=%.2le \n",Omnu2_i,Omextra_i);
	}


	const double sigma8_0 = getsigma_8(cosmo, k_transfer_array, transfer_matter[0]);
	int coll_index = index_zi-j_collapse;
	if(debug_mode>0) printf("coll_index=%d, z_coll = %lf \n", coll_index, zlist_transfer[coll_index]);
	const double sigma8_collapse = getsigma_8(cosmo, k_transfer_array, transfer_matter[coll_index]);
	const double sigma8_i = getsigma_8(cosmo, k_transfer_array, transfer_matter[index_zi]);





//we set delta_short.
 	const double delta_short_min_constant=delta_short_constant_centroid/(2.0+boost_initial_conditions); //initial conditions, to do bipartition. these are the benchmark.
 	const double delta_short_max_constant=delta_short_constant_centroid*(2.0+boost_initial_conditions);



	const double tolerance=0.0002/fmax(fmin(precision_scale,10.),1.); //relative error before we call it converged. 3e-4 should guarantee 0.1% in bias.
	const double tolerance_ini= tolerance + (do_clustering>0)*tolerance*4.0; //first try does not have to be as precise (if we do clustering), since it is used to find Mnu collapse only.
	const double tolerance_z = cosmo->z_collapse * 0.1; //this is just to make sure that we are not artificially converging to a "bisection" \delta_crit if our initial conditions are bad.



//we set delta_long. We assume Lambda (DE) does not cluster.
//for matter first, rest will depend on z.
	double *Ti_klong; //transfer function of CDM+b @ zi
	double *dTi_klong; //derivative of transfer function of CDM+b @ zi, wrt to z
	Ti_klong = allocate_1D_array(cosmo->N_klong);
	dTi_klong = allocate_1D_array(cosmo->N_klong);





//for usual species:
	double **transfer_gamma_klong;//transfer evaluated at klong, for quick interpolation.
	double **transfer_nu_massless_klong;
	double **transfer_nu1_klong;
	transfer_gamma_klong = allocate_2D_array(cosmo->N_klong,Nz_transfer);
	transfer_nu_massless_klong = allocate_2D_array(cosmo->N_klong,Nz_transfer);
	transfer_nu1_klong = allocate_2D_array(cosmo->N_klong,Nz_transfer);


//we set the initial cdm+b transfer functions, and the interpolators for all fluids.
		for(int i_klong=0;i_klong<cosmo->N_klong;i_klong++){ //we find transfer function at initial redshift
			double k_long = klong_list[i_klong];
			Ti_klong[i_klong]=interpol_2D(transfer_matter, zlist_transfer, Nz_transfer, k_transfer_array, length_transfer, zi, k_long);

			dTi_klong[i_klong]=(interpol_2D(transfer_matter, zlist_transfer, Nz_transfer, k_transfer_array, length_transfer, zip2, k_long)
			-interpol_2D(transfer_matter, zlist_transfer, Nz_transfer, k_transfer_array, length_transfer, zip1, k_long))/(zip2-zip1);

			for(i=0;i<Nz_transfer;i++){
				z=zlist_transfer[i];
				transfer_gamma_klong[i_klong][i]=interpol_2D(transfer_gamma, zlist_transfer, Nz_transfer, k_transfer_array, length_transfer, z, k_long);
				transfer_nu_massless_klong[i_klong][i]=interpol_2D(transfer_nu_massless, zlist_transfer, Nz_transfer, k_transfer_array, length_transfer, z, k_long);	//			printf("i_klong=%d, i=%ld, z=%.3le, Tf_gamma=%.1le, Tf_nu=%.1le \n",i_klong, i, z,transfer_gamma_klong[i_klong][i],transfer_nu_massless_klong[i_klong][i]);
				transfer_nu1_klong[i_klong][i]=interpol_2D(transfer_nu1, zlist_transfer, Nz_transfer, k_transfer_array, length_transfer, z, k_long);	//			printf("i_klong=%d, i=%ld, z=%.3le, Tf_gamma=%.1le, Tf_nu=%.1le \n",i_klong, i, z,transfer_gamma_klong[i_klong][i],transfer_nu1_klong[i_klong][i]);//			printf("k_long=%.1le, z=%.1le , T_g=%.1le \n",k_long, z, transfer_gamma_klong[i_klong][i]);
			}
		}


	// and for extra species. No need to add clause: if(mnu2>0), since it will just be zeroes if so.
		double **transfer_nu2_klong;
		double **transfer_extra_klong;
		transfer_nu2_klong = allocate_2D_array(cosmo->N_klong,Nz_transfer);
		transfer_extra_klong = allocate_2D_array(cosmo->N_klong,Nz_transfer);

		for(int i_klong=0;i_klong<cosmo->N_klong;i_klong++){
			double k_long = klong_list[i_klong];
			for(i=0;i<Nz_transfer;i++){
				z=zlist_transfer[i];
				transfer_nu2_klong[i_klong][i]=interpol_2D(transfer_nu2, zlist_transfer, Nz_transfer, k_transfer_array, length_transfer, z, k_long);	//			printf("i_klong=%d, i=%ld, z=%.3le, Tf_gamma=%.1le, Tf_nu=%.1le \n",i_klong, i, z,transfer_gamma_klong[i_klong][i],transfer_nu1_klong[i_klong][i]);
				transfer_extra_klong[i_klong][i]=interpol_2D(transfer_extra, zlist_transfer, Nz_transfer, k_transfer_array, length_transfer, z, k_long);	//
			}
		}




	///////////////////////////////////////////////////////////////////////////////////
	////    we iterate over klong array defined above, as well as Mhalos					////
	/////////////////////////////////////////////////////////////////////////////////

		int iM;
		const int N_klong_calc [N_delta_long] = {1, cosmo->N_klong}; //what k_longs we calculate over. For delta_long=0 only one is needed, otherwise one per kmode.
		//N_delta_long is defined in common.h

		double *delta_long_list;
		delta_long_list=allocate_1D_array(N_delta_long);

		const double delta_long_min=0.; //the maximum value is specified in common.h for convenience to change it.
		const double delta_long_step=(delta_long_max-delta_long_min)/(N_delta_long-1.);

		for(i=0;i<N_delta_long;i++){
			delta_long_list[i]=delta_long_min+delta_long_step*i;
		}

		double **delta_short_collapse; //delta_short (initial) to collapse for each delta_long and k_long.

		delta_short_collapse=allocate_2D_array(N_delta_long,cosmo->N_klong);



		int nuclustering_check;



//these variables are to save the delta_L and delta crit to files
	double **delta_short_crit; //delta_short_critical, extrapolated to z_collapse.
	double **delta_long_collapse; //delta_long, extrapolated to z_collapse.

	delta_short_crit=allocate_2D_array(N_delta_long,cosmo->N_klong);
	delta_long_collapse=allocate_2D_array(N_delta_long,cosmo->N_klong);

	double T_z_collapse_klong; //transfer function T(k_long) at zcollapse, to extrapolate \delta_long.




for(iM=0;iM<cosmo->N_Mhalo;iM++){


	cosmo->Mhalo = cosmo->Mhalo_array[iM];
	cosmo->Mhalo_Mpc = cosmo->Mhalo * MsuntoMpc;//Mhalo in Mpc (*G)

	if(debug_mode>=0)  printf("M_halo = %.1le Msun \n", cosmo->Mhalo);




//average R_i
	double Ribar= pow(cosmo->OmegaM*pow(1.+zi,3.)/2.*cosmo->H0_Mpc*cosmo->H0_Mpc/cosmo->Mhalo_Mpc,-1./3); //R in Mpc.


	//sigma(M) at zi, z)collapse, and derivative (and z=0 just in case), for initial conditions of ODE. (option_initial_velocity=0)
	double sigmaM_0 = getsigma_M(cosmo,  k_transfer_array, transfer_matter[0]);
	double sigmaM_collapse = getsigma_M(cosmo,  k_transfer_array, transfer_matter[coll_index]);
	double sigmaM_i = getsigma_M(cosmo,  k_transfer_array, transfer_matter[index_zi]);
	double sigmaM_zip1 = getsigma_M(cosmo,  k_transfer_array, transfer_matter[index_zip1]);
	double sigmaM_zip2 = getsigma_M(cosmo,  k_transfer_array, transfer_matter[index_zip2]);
	double d_sigmaM_i = (sigmaM_zip2 - sigmaM_zip1)/(zip2-zip1);



	//We also use the transfer functions (option_initial_velocity=1)
	double kstar = 1.0/pow(cosmo->OmegaM/2.*cosmo->H0_Mpc*cosmo->H0_Mpc/cosmo->Mhalo_Mpc,-1./3); //Mpc-1; where we evaluate Tdot/T.

//we find the k that is closest to kstar, since that is easier than interpolating for no reason.
	int indexkstar = find_value(length_transfer, k_transfer_array, kstar);
	double Ti_kstar=transfer_array_zi[indexkstar];
	double dTi_kstar=d_transfer_array_zi[indexkstar];
	double T_z_collapse_kstar=transfer_array_z_collapse[indexkstar]; //transfer function T(k_pivot) at zcollapse, to extrapolate \delta_crit.






	if(debug_mode>0){
		printf("IC_1:%.3le  IC_2:%.3le \n", d_sigmaM_i/sigmaM_i, dTi_kstar/Ti_kstar);
		printf("Ribar=%.3le \n", Ribar);

		printf("We will run some checks: \n");
		printf("@z=%.1le \t sigma_M= %.3le \n",0., sigmaM_0);
		printf("@z=%.1le \t sigma_8= %.3le \n",0., sigma8_0);
		printf("@z=%.1le \t sigma_M= %.3le \n", cosmo->z_collapse, sigmaM_collapse);
		printf("@z=%.1le \t sigma_8= %.3le \n", cosmo->z_collapse, sigma8_collapse);
		printf("@z=%.1le \t sigma_M= %.3le \n",zi, sigmaM_i);
		printf("@z=%.1le  d sigma_M/dz= %.3le \n",zi, d_sigmaM_i);
		printf("@z=%.1le \t sigma_8= %.3le \n",zi, sigma8_i);
	}


	#ifdef _OPENMP
		omp_set_nested(1); //we need to nest the two parallel loops since the limits on the inner for() depend on the outer one, so the scheduler does not know how to do collapse(2).
	#endif

	#pragma omp parallel for
	for(int i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
		#pragma omp parallel for
		for(int i_klong=0;i_klong<N_klong_calc[i_delta_long];i_klong++){



			double delta_long=delta_long_list[i_delta_long]; //long-wavelength perturbation for CDM+b , for photons , and for massless neutrinos.


			double k_long=klong_list[i_klong];

			double ddelta_long_dt=delta_long*(-Hi*(1.+zi))*dTi_klong[i_klong]/Ti_klong[i_klong]; //d delta_long/dt, calculated as dT/dt*1/T*delta_long.

			double delta_short_min=delta_short_min_constant; //we restart the initial conditions, and redo the procedure.
			double delta_short_max=delta_short_max_constant;
			double delta_short=0.0,ddelta_short_dt; //private copies for each run
			double Ri, Ridot, Rpi, z_coll_iteration=0.0; //initial conditions and z of collapse of each iteration.




			if(debug_mode>1){
				printf("deltashort_min=%.2le, deltashort_max=%.2le \n", delta_short_min, delta_short_max);
			}

			for(;delta_short_max/delta_short_min-1.>tolerance_ini;){


				delta_short=(delta_short_max+delta_short_min)/2.0; //bisection, we update the min or max value.

				if(option_initial_velocity==0){
			///Option 1: through sigma(M)
					ddelta_short_dt=delta_short*(-Hi*(1.+zi))*d_sigmaM_i/sigmaM_i; //since: 		ddelta_dt/delta = dsigma_dt/sigma = -H*(1+z) dsigma_dz/sigma
				}
				else if(option_initial_velocity==1) {
			///Option 2: through transfer function (they agree)
					ddelta_short_dt=delta_short*(-Hi*(1.+zi))*dTi_kstar/Ti_kstar; //since:		ddelta_dt/delta = dT_dt/T = -H*(1+z) dT_dz/sigma
				}
				else {
					printf("ERROR: option_initial_velocity has to be either 0 or 1 \n");
					exit(0);
				}


			//these are the i.c. we will feed the function find_z_collapse_XXX.
			Ri= Ribar*(1.0-1./3.*(delta_short+delta_long)); //R_i in Mpc.
			Ridot= Ri*Hi*(1.0-1./3.*(ddelta_short_dt+ddelta_long_dt)/Hi); //dR/dt in Mpc/Mpc
			Rpi= -Ridot/((1.0+zi)*Hi); //dR/dz in Mpc.


				if(debug_mode > 0){
					printf("Ri=%.3le, dRi/dz=%.3le, delta_short=%.2le \n", Ri, Rpi, delta_short);
//					printf("i=%d, T_cdm_i=%.2le, T_g=%.2le, T_nu=%.2le \n",i_klong, Ti_klong[i_klong], transfer_gamma_klong[i_klong][33], transfer_nu_massless_klong[i_klong][33]);
				}


//we now call the collapse routine relevant for each case:

 				if(cosmo->Omega_extra + cosmo->Omeganu2 + cosmo->Omeganu1 == 0){//nothing extra
					if (debug_mode>0) printf("Collapse with nothing extra \n");
					z_coll_iteration=find_z_collapse_nothing(cosmo, Ri, Rpi, delta_long, zlist_transfer,
									Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong]);
				}
				else if(cosmo->Omega_extra + cosmo->Omeganu2 == 0){//only mnu1
					if (debug_mode>0) printf("Collapse with nu1 \n");
					z_coll_iteration=find_z_collapse_1nu(cosmo, Ri, Rpi, delta_long, zlist_transfer,
									Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
									transfer_nu1_klong[i_klong], zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS,
									Nz_solution, Rhalo_solution, Mnu1_solution);
				}
				else if (cosmo->Omega_extra == 0 && cosmo->Omeganu2 > 0){//mnu1 and mnu2
					if (debug_mode>0) printf("Collapse with nu1 and nu2 \n");
					z_coll_iteration=find_z_collapse_2nu(cosmo, Ri, Rpi, delta_long, zlist_transfer,
									Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong], transfer_nu1_klong[i_klong],transfer_nu2_klong[i_klong],
									zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_nu2_EoS, plist_nu2_EoS,
									Nz_solution, Rhalo_solution, Mnu1_solution, Mnu2_solution);

				}
				else if(cosmo->Omeganu1 + cosmo->Omeganu2 == 0){//only extra
					if (debug_mode>0) printf("Collapse with extra \n");
					z_coll_iteration=find_z_collapse_1nu(cosmo, Ri, Rpi, delta_long, zlist_transfer,
									Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
									transfer_extra_klong[i_klong], zmin_EoS, dz_EoS, Nz_EoS, rholist_extra_EoS, plist_extra_EoS,
									Nz_solution, Rhalo_solution, Mextra_solution);
				}
				else if (cosmo->Omega_extra > 0 && cosmo->Omeganu2 == 0){//mnu1 and extra, useful for CAMB
					if (debug_mode>0) printf("Collapse with nu1 and extra \n");
					z_coll_iteration=find_z_collapse_2nu(cosmo, Ri, Rpi, delta_long, zlist_transfer,
									Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong], transfer_nu1_klong[i_klong],
									transfer_extra_klong[i_klong],
									zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_extra_EoS, plist_extra_EoS,
									Nz_solution, Rhalo_solution, Mnu1_solution, Mextra_solution);

				}
				else{	//mnu1, mnu2 and extra
					if (debug_mode>0) printf("Collapse with nu1,nu2 & extra \n");
					z_coll_iteration=find_z_collapse_3nu(cosmo, Ri, Rpi, delta_long, zlist_transfer,
									Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
									transfer_nu1_klong[i_klong],transfer_nu2_klong[i_klong],transfer_extra_klong[i_klong],
									zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_nu2_EoS, plist_nu2_EoS, rholist_extra_EoS, plist_extra_EoS,
									Nz_solution, Rhalo_solution, Mnu1_solution, Mnu2_solution, Mextra_solution);

				}



				if(z_coll_iteration>cosmo->z_collapse){ //collapses too quickly
					delta_short_max=delta_short;
				}
				else{ 							//collapses too slowly
					delta_short_min=delta_short;
				}






			}	//end of delta_short loop

			if(debug_mode > 0){
					printf("delta_short=%le, delta_long=%le (with k_long=%le) \n",delta_short, delta_long, k_long);
					printf("z_coll_iteration= %le \n",z_coll_iteration);
				}


			if(abs(z_coll_iteration-cosmo->z_collapse)>tolerance_z){ //collapses too quickly
				printf("Not converged to z_collapse. Initial conditions too narrow, make boost_initial_conditions bigger. \n");
				printf("z_coll_it=%.3le and z_collapse=%.3le \n",z_coll_iteration,cosmo->z_collapse);
			}

			delta_short_collapse[i_delta_long][i_klong]=delta_short;


//			double percentage=1.*(i_klong+1.+cosmo->N_klong*i_delta_long)/N_delta_long/cosmo->N_klong;//not really a percentage, a one-centage.
//			printProgress(percentage);
//	we do not use percentage bar if parallel. It looks ugly.

		}	// end of k_long loop
	}	//end of delta_long loop



	// since the delta_long=0 does not depend on k we can just copy it for all ks.
	for(int i_delta_long=0;i_delta_long<1;i_delta_long++){
		for(int i_klong=1;i_klong<cosmo->N_klong;i_klong++){
			delta_short_collapse[i_delta_long][i_klong]=delta_short_collapse[i_delta_long][0];
		}
	}



	if (debug_mode > 1){
		for (i=0;i<Nz_solution;i++){
			printf("i=%d, z=%.1le, R(z)/Mpc=%.1le \n",i, logz_solution_min*exp(dlogz_solution*i),Rhalo_solution[i]);
		}
	}




if(do_clustering == 1 && (cosmo->Omega_extra + cosmo->Omeganu2 + cosmo->Omeganu1 > 0)){

	for(i_clustering_index=0;i_clustering_index<Nrecusions_nu_clustering;i_clustering_index++){


		printf("We add clustering of neutrinos/extra species: \n");

		//we compute the collapse of neutrinos/other species_counter.
		//we do it one at a time, since no more than one usually matters.

		if(cosmo->Omeganu1>0 && (cosmo->mnu1 > mnu_min_clustering*sqrt(Mhalo_min_clustering/cosmo->Mhalo)) ){
			nuclustering_check=findmass_1nu(cosmo, cosmo->mnu1, cosmo->T0_nu, zmin_EoS, dz_EoS, Nz_EoS,
					rholist_nu1_EoS, rholist_nu2_EoS, rholist_extra_EoS, etalist,
					logz_solution_min, dlogz_solution, Nz_solution, Rhalo_solution, Mnu1_solution);
		}
		if(cosmo->Omeganu2>0 && (cosmo->mnu2 > mnu_min_clustering*sqrt(Mhalo_min_clustering/cosmo->Mhalo)) ){
			nuclustering_check=findmass_1nu(cosmo, cosmo->mnu2, cosmo->T0_nu, zmin_EoS, dz_EoS, Nz_EoS,
					rholist_nu1_EoS, rholist_nu2_EoS, rholist_extra_EoS, etalist,
					logz_solution_min, dlogz_solution, Nz_solution, Rhalo_solution, Mnu2_solution);
		}
		if(cosmo->Omega_extra>0 && (cosmo->m_extra > cosmo->T0_extra/cosmo->T0_nu * mnu_min_clustering*sqrt(Mhalo_min_clustering/cosmo->Mhalo))){
			nuclustering_check=findmass_1nu(cosmo, cosmo->m_extra, cosmo->T0_extra, zmin_EoS, dz_EoS, Nz_EoS,
					rholist_nu1_EoS, rholist_nu2_EoS, rholist_extra_EoS, etalist,
					logz_solution_min, dlogz_solution, Nz_solution, Rhalo_solution, Mextra_solution);
		}
		do_check(nuclustering_check);


		///and we do the same thing now that we have calculated the nu_collapse, but with more precision.

		#pragma omp parallel for
		for(int i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
			#pragma omp parallel for
			for(int i_klong=0;i_klong<N_klong_calc[i_delta_long];i_klong++){

				double delta_long=delta_long_list[i_delta_long];


				double ddelta_long_dt=delta_long*(-Hi*(1.+zi))*dTi_klong[i_klong]/Ti_klong[i_klong];

				double delta_short_min=delta_short_min_constant; //we restart the initial conditions, and redo the procedure.
				double delta_short_max=delta_short_max_constant;
				double delta_short,ddelta_short_dt; //private copies for each run
				double Ri, Ridot, Rpi, z_coll_iteration; //initial conditions and z of collapse of each iteration.


				for(;delta_short_max/delta_short_min-1.>tolerance;){


					delta_short=(delta_short_max+delta_short_min)/2.; //bisection, we update the min or max value.

					if(option_initial_velocity==0){
				///Option 1: through sigma(M)
						ddelta_short_dt=delta_short*(-Hi*(1.+zi))*d_sigmaM_i/sigmaM_i; //since: 		ddelta_dt/delta = dsigma_dt/sigma = -H*(1+z) dsigma_dz/sigma
					}
					else if(option_initial_velocity==1) {
				///Option 2: through transfer function (they agree)
						ddelta_short_dt=delta_short*(-Hi*(1.+zi))*dTi_kstar/Ti_kstar; //since:		ddelta_dt/delta = dT_dt/T = -H*(1+z) dT_dz/sigma
					}
					else {
						printf("option_initial_velocity has to be either 0 or 1 \n");
						exit(0);
					}


				//these are the i.c. we will feed the function find_z_collapse_long.
				Ri= Ribar*(1.0-1./3.*(delta_short+delta_long)); //R_i in Mpc.
				Ridot= Ri*Hi*(1.0-1./3.*(ddelta_short_dt+ddelta_long_dt)/Hi); //dR/dt in Mpc/Mpc
				Rpi= -Ridot/((1.0+zi)*Hi); //dR/dz in Mpc.




		//we now call the collapse routine relevant for each case:

		 				if(cosmo->Omega_extra + cosmo->Omeganu2 + cosmo->Omeganu1 == 0){//nothing extra
							if (debug_mode>0) printf("Collapse with nothing extra \n");
							z_coll_iteration=find_z_collapse_nothing(cosmo, Ri, Rpi, delta_long, zlist_transfer,
											Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong]);
						}
						else if(cosmo->Omega_extra + cosmo->Omeganu2 == 0){//only mnu1
							if (debug_mode>0) printf("Collapse with nu1 \n");
							z_coll_iteration=find_z_collapse_1nu(cosmo, Ri, Rpi, delta_long, zlist_transfer,
											Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
											transfer_nu1_klong[i_klong], zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS,
											Nz_solution, Rhalo_solution, Mnu1_solution);
						}
						else if (cosmo->Omega_extra == 0 && cosmo->Omeganu2 > 0){//mnu1 and mnu2
							if (debug_mode>0) printf("Collapse with nu1 and nu2 \n");
							z_coll_iteration=find_z_collapse_2nu(cosmo, Ri, Rpi, delta_long, zlist_transfer,
											Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong], transfer_nu1_klong[i_klong],transfer_nu2_klong[i_klong],
											zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_nu2_EoS, plist_nu2_EoS,
											Nz_solution, Rhalo_solution, Mnu1_solution, Mnu2_solution);

						}
						else if(cosmo->Omeganu1 + cosmo->Omeganu2 == 0){//only extra
							if (debug_mode>0) printf("Collapse with extra \n");
							z_coll_iteration=find_z_collapse_1nu(cosmo, Ri, Rpi, delta_long, zlist_transfer,
											Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
											transfer_extra_klong[i_klong], zmin_EoS, dz_EoS, Nz_EoS, rholist_extra_EoS, plist_extra_EoS,
											Nz_solution, Rhalo_solution, Mextra_solution);
						}
						else if (cosmo->Omega_extra > 0 && cosmo->Omeganu2 == 0){//mnu1 and extra, useful for CAMB
							if (debug_mode>0) printf("Collapse with nu1 and extra \n");
							z_coll_iteration=find_z_collapse_2nu(cosmo, Ri, Rpi, delta_long, zlist_transfer,
											Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong], transfer_nu1_klong[i_klong],
											transfer_extra_klong[i_klong],
											zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_extra_EoS, plist_extra_EoS,
											Nz_solution, Rhalo_solution, Mnu1_solution, Mextra_solution);

						}
						else{	//mnu1, mnu2 and extra
							if (debug_mode>0) printf("Collapse with nu1,nu2 & extra \n");
							z_coll_iteration=find_z_collapse_3nu(cosmo, Ri, Rpi, delta_long, zlist_transfer,
											Ti_klong[i_klong], transfer_gamma_klong[i_klong], transfer_nu_massless_klong[i_klong],
											transfer_nu1_klong[i_klong],transfer_nu2_klong[i_klong],transfer_extra_klong[i_klong],
											zmin_EoS, dz_EoS, Nz_EoS, rholist_nu1_EoS, plist_nu1_EoS, rholist_nu2_EoS, plist_nu2_EoS, rholist_extra_EoS, plist_extra_EoS,
											Nz_solution, Rhalo_solution, Mnu1_solution, Mnu2_solution, Mextra_solution);

						}



						if(z_coll_iteration>cosmo->z_collapse){ //collapses too quickly
							delta_short_max=delta_short;
						}
						else{ 							//collapses too slowly
							delta_short_min=delta_short;
						}






					}	//end of delta_short loop



				if(abs(z_coll_iteration-cosmo->z_collapse)>tolerance_z){ //collapses too quickly
					printf("Not converged to z_collapse (after doing Mnu_collapse). Initial conditions too narrow, make delta_short_max and delta_short_min wider. \n");
					printf("z_coll_it=%.3le and z_collapse=%.3le \n",z_coll_iteration,cosmo->z_collapse);
				}

				delta_short_collapse[i_delta_long][i_klong]=delta_short;


			}	// end of k_long loop
		}	//end of delta_long loop

		printf("\n Iteration of clustering #%d \n",i_clustering_index+1);

	} //end of do_clustering loop

} //end of do_clustering if statement.



// since the delta_long=0 does not depend on k we can just copy it for all ks (again).
for(int i_delta_long=0;i_delta_long<1;i_delta_long++){
	for(int i_klong=1;i_klong<cosmo->N_klong;i_klong++){
		delta_short_collapse[i_delta_long][i_klong]=delta_short_collapse[i_delta_long][0];
	}
}






//////////////////////////////////////////////////////////////////////////////
//// we now extrapolate  to find\delta_crit to the redshift of collapse	/////
////////////////////////////////////////////////////////////////////////////
#pragma omp parallel for collapse(2)
for(int i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
		for(int i_klong=0;i_klong<cosmo->N_klong;i_klong++){
			double delta_long=delta_long_list[i_delta_long];
			double k_long=klong_list[i_klong];

			T_z_collapse_klong=interpol_2D(transfer_matter, zlist_transfer, Nz_transfer, k_transfer_array, length_transfer, cosmo->z_collapse, k_long);

			delta_long_collapse[i_delta_long][i_klong]=delta_long*T_z_collapse_klong/Ti_klong[i_klong];	//this is right as well //			printf("klong=%.1le, delta_long_collapse/delta_long= %.3le \n", k_long, delta_long_collapse[i_delta_long][i_klong]/delta_long);


			if(option_initial_velocity==0){
				delta_short_crit[i_delta_long][i_klong]=delta_short_collapse[i_delta_long][i_klong]*sigmaM_collapse/sigmaM_i;
			}
			else if(option_initial_velocity==1){
				delta_short_crit[i_delta_long][i_klong]=delta_short_collapse[i_delta_long][i_klong]*T_z_collapse_kstar/Ti_kstar;
			}
			else{
				printf("Error, select an option_initial_velocity \n");
			}
		}
	}


////////////////////////////////////////////
////and we save the results to files.	/////
//////////////////////////////////////////

		lengthname=sprintf(filename,"output/result-%d/delta_initial_z%.2f_M%.2f_Nk%d.dat",cosmo->file_tag,cosmo->z_collapse, log10(cosmo->Mhalo),cosmo->N_klong); //we save the delta_crit extrapolated to z=0.5
		fp=fopen(filename,"w");

		if(print_headers!=0){
			fprintf(fp,"delta_long_initial   k_long[1/Mpc]   delta_short_initial \n");
		}

		for(int i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
			double delta_long=delta_long_list[i_delta_long];
			for(int i_klong=0;i_klong<cosmo->N_klong;i_klong++){
				double k_long=klong_list[i_klong];
				fprintf(fp,"%le   %le   %le \n", delta_long, k_long, delta_short_collapse[i_delta_long][i_klong]);
			}
		}
		fclose(fp);


		lengthname=sprintf(filename,"output/result-%d/delta_crit_z%.2f_M%.2f_Nk%d.dat",cosmo->file_tag,cosmo->z_collapse, log10(cosmo->Mhalo),cosmo->N_klong); //we save the delta_crit extrapolated to z=0.5
		fp=fopen(filename,"w");

		if(print_headers!=0){
			fprintf(fp,"delta_long   k_long[1/Mpc]   delta_crit \n");
		}
		for(int i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
			for(int i_klong=0;i_klong<cosmo->N_klong;i_klong++){
				double k_long=klong_list[i_klong];
				fprintf(fp,"%le   %le   %le \n", delta_long_collapse[i_delta_long][i_klong], k_long, delta_short_crit[i_delta_long][i_klong]);
			}
		}
		fclose(fp);




} //end of Mhalo loop.










//////////////////////////////////////////
//// we free the allocated memory	/////
////////////////////////////////////////

	free(transfer_array_z0);
	free(k_transfer_array);

	free(transfer_array_gamma_zi);
	free(transfer_array_nu_massless_zi);
	free(transfer_array_nu1_zi);
	free(transfer_array_nu2_zi);
	free(transfer_array_extra_zi);


	free(transfer_array_zi);
	free(transfer_array_zip1);
	free(transfer_array_zip2);

	free(transfer_array_z_collapse);


	free(d_transfer_array_zi);

	free(filename);


	free_2D_array(transfer_matter, Nz_transfer);
	free_2D_array(transfer_gamma, Nz_transfer);
	free_2D_array(transfer_nu_massless, Nz_transfer);
	free_2D_array(transfer_nu1, Nz_transfer);
	free_2D_array(transfer_nu2, Nz_transfer);
	free_2D_array(transfer_extra, Nz_transfer);


	free_2D_array(transfer_gamma_klong, cosmo->N_klong);
	free_2D_array(transfer_nu_massless_klong, cosmo->N_klong);
	free_2D_array(transfer_nu1_klong, cosmo->N_klong);
	free_2D_array(transfer_nu2_klong, cosmo->N_klong);
	free_2D_array(transfer_extra_klong, cosmo->N_klong);



	free(klong_list);

	free(Ti_klong);
	free(dTi_klong);


	free(delta_long_list);

	free_2D_array(delta_short_crit,N_delta_long);
	free_2D_array(delta_long_collapse,N_delta_long);

	free_2D_array(delta_short_collapse,N_delta_long);



	free(zlist_EoS);
	free(plist_nu1_EoS);
	free(rholist_nu1_EoS);

	free(plist_nu2_EoS);
	free(rholist_nu2_EoS);

	free(plist_extra_EoS);
	free(rholist_extra_EoS);

	free(etalist);



	free(Rhalo_solution);
	free(Mnu1_solution);
	free(Mnu2_solution);
	free(Mextra_solution);



	return 1;
}






//now collapse functions:

double find_z_collapse_nothing
(Cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long, double * const zlist_transfer,
const double Tfm_klong, double * const transfer_gamma_klong, double * const transfer_nu_klong
){
// This function returns the z at which an overdensity collapses. global cosmological parameters assumed.
// In the presence of a long-wavelength PHOTON and MASSLESS NU perturbation.
//we read the instantaneous transfer function from CLASS/CAMB at each redshift.

	int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20


	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.



	double zf_code=cosmo->z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
	}

	const double zstep_log=(log(zf_code/zi))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood

	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;


	long i;
	double z=zi;
	double dE, OmL, OmR, H, OmM; //z-dependent.
	double H2; //H for previous redshift, for derivatives.

	double OmG, Omnu_massless; //photon, massless and massive neutrino energy density.
	double OmRbar, OmGbar, Omnu_masslessbar; //average for all of them



	double T_gamma, T_nu;//transfer functions at each k and z.

	const double T_matter=Tfm_klong;



//we set the initial H
	OmGbar= cosmo->OmegaG * pow(1.+zi,4.);
	Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+zi,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = cosmo->OmegaM * pow(1.+zi,3.);
	OmL = cosmo->OmegaL; //Omega_L(z), we take it as z-independent.
	H = cosmo->H0_Mpc* sqrt(OmL + OmM + OmRbar); //H(zi)


//and the initial OmR
	T_gamma = interpol(transfer_gamma_klong, zlist_transfer, Nz_transfer, zi);
	T_nu = interpol(transfer_nu_klong, zlist_transfer, Nz_transfer, zi);
	OmG= OmGbar * (1.+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	Omnu_massless= Omnu_masslessbar * (1.+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR2 = OmG + Omnu_massless;




	double z_next, dE2; //for the next values, to do Heun's method. We assume Omega_L does not change.

//symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	z_next = zi * exp(zstep_log);//first step
	OmGbar= cosmo->OmegaG * pow(1.+z_next,4.);
	Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = cosmo->OmegaM * pow(1.+z_next,3.);
	OmL = cosmo->OmegaL;
	H2  = cosmo->H0_Mpc * sqrt(OmL + OmM + OmRbar);// H(z_next)
	dE2 = (H2-H)/H2/(zi*zstep_log);

	double Rpp1, Rpp2; //d^2R(z)/dz^2


///////////////////////////////////////////////////////////
////    here we solve for the collapse								////
/////////////////////////////////////////////////////////
	for(i=1; i<npoints-1 && R2>0.; i++){
		//we update the values from the previous step
			R1 = R2;
			Rp1 = Rp2;
		//including all z-dependent quantities, for Heun's method. Remember 2 means next step.
			z=z_next;
			H=H2;
			OmR=OmR2;
			dE=dE2;

		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.

	//we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- cosmo->Mhalo_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR - 2*OmL)*cosmo->H0_Mpc*cosmo->H0_Mpc/pow(H*(1.+z),2.); //R''(R1,t1)



//we update all z-dependent stuff to the next z
		OmGbar= cosmo->OmegaG * pow(1.+z_next,4.);
		Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+z_next,4.);
		OmRbar = OmGbar + Omnu_masslessbar;
		OmM = cosmo->OmegaM * pow(1.+z_next,3.);


		T_gamma = interpol(transfer_gamma_klong, zlist_transfer, Nz_transfer, z_next);
		T_nu = interpol(transfer_nu_klong, zlist_transfer, Nz_transfer, z_next);

		OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
		OmR2 = OmG + Omnu_massless;


		H2  = cosmo->H0_Mpc * sqrt(OmL + OmM + OmRbar);// H(z_next)

		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value at the end.



//these are too annoying to keep even with debug_mode, activate manually if you want to debug:
//		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, cosmo->H0_Mpc*E_LCDM(cosmo, z), dE, dlogE_dz_LCDM(cosmo, z));
//		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);


//now we evolve R and R'.

		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;


//we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- cosmo->Mhalo_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR2 - 2*OmL)*cosmo->H0_Mpc*cosmo->H0_Mpc/pow(H2*(1.+z_next),2.); //R''(R2,t2)

//and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);





	}

	return z;

}







double find_z_collapse_1nu
(Cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long, double * const zlist_transfer,
const double Tfm_klong, double * const transfer_gamma_klong, double * const transfer_nu_klong, double * const transfer_nu_massive_klong,
const double zmin_EoS, const double dz_EoS, const long Nz_EoS, double * const rholist_EoS, double * const plist_EoS,
const long Nz_solution, double *R_solution, double *Mnu_solution
){
// This function returns the z at which an overdensity collapses. global cosmological parameters assumed.
// In the presence of a long-wavelength PHOTON, MASSLESS and 1 MASSIVE NU (OR OTHER EXTRA) perturbation.
//we read the instantaneous transfer function from CAMB at each redshift.
//and we save R(z) to perform BKT calculation of neutrino clustering, or we input Mnu(z) if we already calculated it.
//we only modify R_solution and Mnu_solution, the rest are read-only.


	const int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20

	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.



	double zf_code=cosmo->z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
	}

	const double zstep_log=(log(zf_code/zi))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood

	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;


	long i;
	double z=zi;



	const double T_matter=Tfm_klong;

	//the nu_massive can be any new species, just make sure that the Omegas are well defined here
	 	 if((cosmo->Omega_extra>0) && (cosmo->mnu1>0)){
	 	   printf("Using wrong collapse function, XXX_1nu only accepts one species. \n");
	 	 }
	double rhonu1_0=rholist_EoS[0];//rho at z=0. units are K^4




//we set the initial H. All these are average densities.
	double OmGbar= cosmo->OmegaG * pow(1.+zi,4.); //photon
	double Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+zi,4.); //massless nu
	double OmRbar = OmGbar + Omnu_masslessbar; //all radiation (massless nu + photon)
	double OmM = cosmo->OmegaM * pow(1.+zi,3.); //matter
	double OmL = cosmo->OmegaL; //Omega_Lambda(z), we take it as z-independent.
	double rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_EoS, Nz_EoS, zi); //rho of nu1 at z
	double rho_ratio_nu1=rhonu1_z/rhonu1_0; //ratio to rho of zero
	double Omnu1bar_z =  cosmo->Omeganu1 * rho_ratio_nu1 + cosmo->Omega_extra * rho_ratio_nu1; //omeganu1
	double H = cosmo->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z); //H(zi)



//and the initial OmR (we assme c_s^2=w=1/3 for photons and massless nu)
	double T_gamma = interpol(transfer_gamma_klong, zlist_transfer, Nz_transfer, zi);
	double T_nu = interpol(transfer_nu_klong, zlist_transfer, Nz_transfer, zi);
	double OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	double Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR, OmR2 = OmG + Omnu_massless;
//and initial perturbation in nu1	(massive) and pressure
	double pnu1_z = interpol_cubic(zmin_EoS, dz_EoS, plist_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu1_z = pnu1_z/rhonu1_z; //equation of state
	double T_nu1_z=interpol(transfer_nu_massive_klong, zlist_transfer, Nz_transfer, zi);
	double delta_nu1_z = delta_long*T_nu1_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next





//for the next values, to do Heun's method. We assume Omega_L does not change. Symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	double z_next = zi * exp(zstep_log);//first step
	OmGbar= cosmo->OmegaG * pow(1.+z_next,4.);
	Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = cosmo->OmegaM * pow(1.+z_next,3.);
	OmL = cosmo->OmegaL;
	rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_EoS, Nz_EoS, z_next);
	rho_ratio_nu1=rhonu1_z/rhonu1_0;
	Omnu1bar_z =  cosmo->Omeganu1 * rho_ratio_nu1 + cosmo->Omega_extra * rho_ratio_nu1;
	double H2 = cosmo->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z); // H(z_next)
	double dE, dE2 = (H2-H)/H2/(zi*zstep_log);

//and the w and c_s^2 of nu at first step we calculate.
	pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist_EoS, Nz_EoS, z_next);
	double wnu1_z2=pnu1_z/rhonu1_z;
	double d_wnu1_z = (wnu1_z2-wnu1_z)/(z_next-zi);
	double csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/(3.0*(1+wnu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).


	double Rpp1, Rpp2; //d^2R(z)/dz^2

	int i_solution=0; //we only check for the solution of neutrino clustering every few steps
	double Mhalo_nu1_Mpc=0; //M of the nu1 halo
	int collapse_steps= (int) npoints/Nz_solution;	//how many steps we take before saving R(z) and/or updating Mnu
	//we save manually the first value of R and mass to avoid running the loop at i=0.
  R_solution[i_solution]=R1;
  Mhalo_nu1_Mpc = Mnu_solution[i_solution];
  i_solution++;


///////////////////////////////////////////////////////////
////    here we solve for the collapse								////
/////////////////////////////////////////////////////////

	for(i=1; i<npoints-1 && R2>0.; i++){
		//we update the values from the previous step
			R1 = R2;
			Rp1 = Rp2;
		//including all z-dependent quantities, for Heun's method. Remember 2 means next step.
			z=z_next;
			H=H2;
			dE=dE2;
			OmR=OmR2;
			wnu1_z=wnu1_z2;

//we save R and read Mnu, only every few steps since it doesn't vary much.
			if (i % collapse_steps == 0){
				R_solution[i_solution]=R1;
				Mhalo_nu1_Mpc = Mnu_solution[i_solution];
				i_solution++;
		    if(debug_mode > 1){
		      printf("z=%.1le, R=%.1le, Mnu/Mhalo=%1le \n\n", z, R_solution[i_solution], Mhalo_nu1_Mpc/cosmo->Mhalo_Mpc);
		    }
		  }


		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.

	//we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- cosmo->Mhalo_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu1_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR + (1.+3.*wnu1_z + delta_nu1_z*(1+3.0*csq_ad_nu1_z))*Omnu1bar_z - 2*OmL)*cosmo->H0_Mpc*cosmo->H0_Mpc/pow(H*(1.+z),2.); //R''(R1,t1)



//we update all z-dependent stuff to the next z
		OmGbar= cosmo->OmegaG * pow(1.+z_next,4.);
		T_gamma = interpol(transfer_gamma_klong, zlist_transfer, Nz_transfer, z_next);
		OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)

		Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+z_next,4.);
		T_nu = interpol(transfer_nu_klong, zlist_transfer, Nz_transfer, z_next);
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)

		OmRbar = OmGbar + Omnu_masslessbar;
		OmR2 = OmG + Omnu_massless;

		OmM = cosmo->OmegaM * pow(1.+z_next,3.);



		rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_EoS, Nz_EoS, z_next);
		rho_ratio_nu1=rhonu1_z/rhonu1_0;
		Omnu1bar_z = cosmo->Omeganu1 * rho_ratio_nu1 + cosmo->Omega_extra * rho_ratio_nu1;
		T_nu1_z=interpol(transfer_nu_massive_klong, zlist_transfer, Nz_transfer, z_next);
		delta_nu1_z = delta_long*T_nu1_z/T_matter; //long-wavlength neutrino perturbation
		pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist_EoS, Nz_EoS, z_next);
		wnu1_z2=pnu1_z/rhonu1_z;
		d_wnu1_z = (wnu1_z2-wnu1_z)/zstep_lin;
		csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/3.0*(1.0+z_next)/(1.0+wnu1_z2); //adiabatic sound speed squared.


		H2  = cosmo->H0_Mpc * sqrt(OmL + OmM + OmRbar + Omnu1bar_z);// H(z_next)

		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value, although it's unnecessary.



		if(debug_mode > 1){ //check that derivatives are well defined, it should be for any reasonable value of npoints, but just in case:
			printf("w(z1)-w(z2)=%.3le , w(z2)=%.3le , dw/dz=%.3le \n", wnu1_z-wnu1_z2, wnu1_z2, d_wnu1_z);
			printf("H(z1)-H(z2)=%.3le , H(z2)=%.3le , dH/dz=%.3le \n", H-H2, H2, dE2*H2);
			printf("z=%.3le, cs^2/w(z2)=%.3le \n", z_next, csq_ad_nu1_z/wnu1_z2); //cs^2 and w are within 10% of each other.
	}


//these are too annoying to keep even with debug_mode, activate manually if you want:
//		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, cosmo->H0_Mpc*E_LCDM(cosmo, z), dE, dlogE_dz_LCDM(cosmo, z));
//		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);


//now we evolve R and R'.

		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;


//we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- cosmo->Mhalo_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu1_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR2 + (1.0+3.0*wnu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z - 2*OmL)*cosmo->H0_Mpc*cosmo->H0_Mpc/pow(H2*(1.+z_next),2.); //R''(R2,t2)


//and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);





	}



		return z;

}





double find_z_collapse_2nu
(Cosmology *const cosmo, const double Ri, const double Rpi, const double delta_long, double * const zlist_transfer,
const double Tfm_klong, double * const transfer_gamma_klong, double * const transfer_nu_klong, double * const transfer_nu1_klong, double * const transfer_nu2_klong,
const double zmin_EoS, const double dz_EoS, const long Nz_EoS, double * const rholist1_EoS, double * const plist1_EoS, double * const rholist2_EoS, double * const plist2_EoS,
const long Nz_solution, double *R_solution, double *Mnu1_solution, double *Mnu2_solution
){
// This function returns the z at which an overdensity collapses. global cosmological parameters assumed.
// In the presence of a long-wavelength PHOTON, MASSLESS and 2 MASSIVE NU (or 1 nu + 1 extra) perturbation.
//we read the instantaneous transfer function from CAMB at each redshift.
//and we save R(z) to perform BKT calculation of neutrino clustering, or we input Mnu(z) if we already calculated it.




	const int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20

	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.



	double zf_code=cosmo->z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
	}

	const double zstep_log=(log(zf_code/zi))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood

	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;


	long i;
	double z=zi;



	const double T_matter=Tfm_klong;

	//the nu1 is nu1, but nu2 can be any new species, just make sure that the Omegas are well defined here
	 	 if((cosmo->Omega_extra>0) && (cosmo->mnu2>0)){
	 	   printf("Using wrong collapse function, XXX_2nu only accepts nu1 + one other species. \n");
	 	 }
	double rhonu1_0=rholist1_EoS[0];//rho1 at z=0. units are K^4
  double rhonu2_0=rholist2_EoS[0];//rho2 at z=0. units are K^4




//we set the initial H. All these are average densities.
	double OmGbar= cosmo->OmegaG * pow(1.+zi,4.); //photon
	double Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+zi,4.); //massless nu
	double OmRbar = OmGbar + Omnu_masslessbar; //all radiation (massless nu + photon)
	double OmM = cosmo->OmegaM * pow(1.+zi,3.); //matter
	double OmL = cosmo->OmegaL; //Omega_Lambda(z), we take it as z-independent.

	double rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, zi); //rho of nu1 at z
	double rho_ratio_nu1=rhonu1_z/rhonu1_0; //ratio to rho of zero
	double Omnu1bar_z =  cosmo->Omeganu1 * rho_ratio_nu1; //omeganu1

	double rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, zi); //rho of nu1 at z
	double rho_ratio_nu2=rhonu2_z/rhonu2_0; //ratio to rho of zero
	double Omnu2bar_z =  cosmo->Omeganu2 * rho_ratio_nu2 + cosmo->Omega_extra * rho_ratio_nu2; //omeganu2 or omega extra, only one!

	double H = cosmo->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z); //H(zi)


//and the initial OmR (we assme c_s^2=w=1/3 for photons and massless nu)
	double T_gamma = interpol(transfer_gamma_klong, zlist_transfer, Nz_transfer, zi);
	double T_nu = interpol(transfer_nu_klong, zlist_transfer, Nz_transfer, zi);
	double OmG= OmGbar * (1.+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	double Omnu_massless= Omnu_masslessbar * (1.+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR, OmR2 = OmG + Omnu_massless;
//and initial perturbation in nu1	(massive) and pressure
	double pnu1_z = interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu1_z = pnu1_z/rhonu1_z; //equation of state
	double T_nu1_z=interpol(transfer_nu1_klong, zlist_transfer, Nz_transfer, zi);
	double delta_nu1_z = delta_long*T_nu1_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//and initial perturbation in nu2	(massive) and pressure
	double pnu2_z = interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu2_z = pnu2_z/rhonu2_z; //equation of state
	double T_nu2_z=interpol(transfer_nu2_klong, zlist_transfer, Nz_transfer, zi);
	double delta_nu2_z = delta_long*T_nu2_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next



//for the next values, to do Heun's method. We assume Omega_L does not change. Symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	double z_next = zi * exp(zstep_log);//first step
	OmGbar= cosmo->OmegaG * pow(1.+z_next,4.);
	Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = cosmo->OmegaM * pow(1.+z_next,3.);
	OmL = cosmo->OmegaL;

	rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, z_next);
	rho_ratio_nu1=rhonu1_z/rhonu1_0;
	Omnu1bar_z =  cosmo->Omeganu1 * rho_ratio_nu1;

	rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, z_next);
	rho_ratio_nu2=rhonu2_z/rhonu2_0;
	Omnu2bar_z =  cosmo->Omeganu2 * rho_ratio_nu2 + cosmo->Omega_extra * rho_ratio_nu2;

	double H2 = cosmo->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z); // H(z_next)
	double dE, dE2 = (H2-H)/H2/(zi*zstep_log);

//and the w and c_s^2 of nus at first step we calculate.
//for nu1
		pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, z_next);
		double wnu1_z2=pnu1_z/rhonu1_z;
		double d_wnu1_z = (wnu1_z2-wnu1_z)/(z_next-zi);
		double csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/(3.0*(1+wnu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//for nu2
		pnu2_z=interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, z_next);
		double wnu2_z2=pnu2_z/rhonu2_z;
		double d_wnu2_z = (wnu2_z2-wnu2_z)/(z_next-zi);
		double csq_ad_nu2_z = wnu2_z2 + d_wnu2_z/(3.0*(1+wnu2_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).


	double Rpp1, Rpp2; //d^2R(z)/dz^2

	int i_solution=0; //we only check for the solution of neutrino clustering every few steps
	double Mhalo_nu1_Mpc=0, Mhalo_nu2_Mpc=0; //M of the nu1 and nu2 haloes
	int collapse_steps= (int) npoints/Nz_solution;	//how many steps we take before saving R(z) and/or updating Mnu
	//we save manually the first value of R and mass to avoid running the loop at i=0.
  R_solution[i_solution]=R1;
	Mhalo_nu1_Mpc = Mnu1_solution[i_solution];
	Mhalo_nu2_Mpc = Mnu2_solution[i_solution];
  i_solution++;


///////////////////////////////////////////////////////////
////    here we solve for the collapse								////
/////////////////////////////////////////////////////////
	for(i=1; i<npoints-1 && R2>0.; i++){
		//we update the values from the previous step
			R1 = R2;
			Rp1 = Rp2;
		//including all z-dependent quantities, for Heun's method. Remember 2 means next step.
			z=z_next;
			H=H2;
			dE=dE2;
			OmR=OmR2;
			wnu1_z=wnu1_z2;
			wnu2_z=wnu2_z2;


//we save R and read Mnu, only every few steps since it doesn't vary much.
			if (i % collapse_steps == 0){
				R_solution[i_solution]=R1;
				Mhalo_nu1_Mpc = Mnu1_solution[i_solution];
				Mhalo_nu2_Mpc = Mnu2_solution[i_solution];
				i_solution++;
		    if(debug_mode > 1){
		      printf("z=%.1le, R=%.1le, Mnu1=%1le , Mnu2=%1le \n\n", z, R_solution[i_solution], Mhalo_nu1_Mpc, Mhalo_nu2_Mpc);
		    }
		  }


		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.


//these are the terms that go inside the sum, define them outside for clarity.
		double Onu1 = (1.0+3.0*wnu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
		double Onu2 = (1.0+3.0*wnu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;
	//we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- cosmo->Mhalo_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu1_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu2_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR + Onu1 + Onu2 - 2*OmL)*cosmo->H0_Mpc*cosmo->H0_Mpc/pow(H*(1.+z),2.); //R''(R1,t1)



//we update all z-dependent stuff to the next z
		OmGbar= cosmo->OmegaG * pow(1.+z_next,4.);
		T_gamma = interpol(transfer_gamma_klong, zlist_transfer, Nz_transfer, z_next);
		OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)

		Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+z_next,4.);
		T_nu = interpol(transfer_nu_klong, zlist_transfer, Nz_transfer, z_next);
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)

		OmRbar = OmGbar + Omnu_masslessbar;
		OmR2 = OmG + Omnu_massless;

		OmM = cosmo->OmegaM * pow(1.+z_next,3.);


		rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, z_next);
		rho_ratio_nu1=rhonu1_z/rhonu1_0;
		Omnu1bar_z = cosmo->Omeganu1 * rho_ratio_nu1;
		T_nu1_z=interpol(transfer_nu1_klong, zlist_transfer, Nz_transfer, z_next);
		delta_nu1_z = delta_long*T_nu1_z/T_matter; //long-wavlength neutrino perturbation
		pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, z_next);
		wnu1_z2=pnu1_z/rhonu1_z;
		d_wnu1_z = (wnu1_z2-wnu1_z)/zstep_lin;
		csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/3.0*(1.0+z_next)/(1.0+wnu1_z2); //adiabatic sound speed squared.

		rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, z_next);
		rho_ratio_nu2=rhonu2_z/rhonu2_0;
		Omnu2bar_z = cosmo->Omeganu2 * rho_ratio_nu2 + cosmo->Omega_extra * rho_ratio_nu2;
		T_nu2_z=interpol(transfer_nu2_klong, zlist_transfer, Nz_transfer, z_next);
		delta_nu2_z = delta_long*T_nu2_z/T_matter; //long-wavlength neutrino perturbation
		pnu2_z=interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, z_next);
		wnu2_z2=pnu2_z/rhonu2_z;
		d_wnu2_z = (wnu2_z2-wnu2_z)/zstep_lin;
		csq_ad_nu2_z = wnu2_z2 + d_wnu2_z/3.0*(1.0+z_next)/(1.0+wnu2_z2); //adiabatic sound speed squared.

//we update the Onu terms that go inside the integral.
		Onu1 = (1.0+3.0*wnu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
		Onu2 = (1.0+3.0*wnu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;

		H2  = cosmo->H0_Mpc * sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z);// H(z_next)

		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value at the end.



//these are too annoying to keep even with debug_mode, activate manually if you want:
//		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, cosmo->H0_Mpc*E_LCDM(cosmo, z), dE, dlogE_dz_LCDM(cosmo, z));
//		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);


//now we evolve R and R'.

		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;


//we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- cosmo->Mhalo_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu1_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu2_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR2 + Onu1 + Onu2 - 2*OmL)*cosmo->H0_Mpc*cosmo->H0_Mpc/pow(H2*(1.+z_next),2.); //R''(R2,t2)


//and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);





	}



		return z;

}



double find_z_collapse_3nu
(Cosmology *const cosmo, const double Ri, const double Rpi, const double delta_long, double * const  zlist_transfer,
const double Tfm_klong, double * const transfer_gamma_klong, double * const transfer_nu_klong, double * const transfer_nu1_klong, double * const transfer_nu2_klong, double * const transfer_nu3_klong,
const double zmin_EoS, const double dz_EoS, const long Nz_EoS, double * const rholist1_EoS, double * const plist1_EoS, double * const rholist2_EoS, double * const plist2_EoS, double * const rholist3_EoS, double * const plist3_EoS,
const long Nz_solution, double *R_solution, double *Mnu1_solution, double *Mnu2_solution, double *Mnu3_solution
){
// This function returns the z at which an overdensity collapses. global cosmological parameters assumed.
// In the presence of a long-wavelength PHOTON, MASSLESS and 2 MASSIVE NU perturbation, and an EXTRA SPECIES.
//we read the instantaneous transfer function from CAMB at each redshift.
//and we save R(z) to perform BKT calculation of neutrino clustering, or we input Mnu(z) if we already calculated it.





	const int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20

	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.



	double zf_code=cosmo->z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
	}

	const double zstep_log=(log(zf_code/zi))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood

	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;


	long i;
	double z=zi;



	const double T_matter=Tfm_klong;

	//the nu1 is nu1, but nu2 can be any new species, just make sure that the Omegas are well defined here
	 	 if((cosmo->Omega_extra==0) || (cosmo->mnu2==0) || (cosmo->mnu1==0)){
	 	   printf("Using wrong collapse function, XXX_3nu works with 3 nus (or 2+extra). \n");
	 	 }
	double rhonu1_0=rholist1_EoS[0];//rho1 at z=0. units are K^4
	double rhonu2_0=rholist2_EoS[0];//rho2 at z=0. units are K^4
	double rhonu3_0=rholist3_EoS[0];//rho3 at z=0. units are K^4





//we set the initial H. All these are average densities.
	double OmGbar= cosmo->OmegaG * pow(1.+zi,4.); //photon
	double Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+zi,4.); //massless nu
	double OmRbar = OmGbar + Omnu_masslessbar; //all radiation (massless nu + photon)
	double OmM = cosmo->OmegaM * pow(1.+zi,3.); //matter
	double OmL = cosmo->OmegaL; //Omega_Lambda(z), we take it as z-independent.

	double rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, zi); //rho of nu1 at z
	double rho_ratio_nu1=rhonu1_z/rhonu1_0; //ratio to rho of zero
	double Omnu1bar_z =  cosmo->Omeganu1 * rho_ratio_nu1; //omeganu1

	double rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, zi); //rho of nu2 at z
	double rho_ratio_nu2=rhonu2_z/rhonu2_0; //ratio to rho of zero
	double Omnu2bar_z =  cosmo->Omeganu2 * rho_ratio_nu2; //omeganu2

	double rhonu3_z=interpol_cubic(zmin_EoS, dz_EoS, rholist3_EoS, Nz_EoS, zi); //rho of nu3 at z
	double rho_ratio_nu3=rhonu3_z/rhonu3_0; //ratio to rho of zero
	double Omnu3bar_z = cosmo->Omega_extra * rho_ratio_nu3; //omeganu3

	double H = cosmo->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z + Omnu3bar_z); //H(zi)


//and the initial OmR
	double T_gamma = interpol(transfer_gamma_klong, zlist_transfer, Nz_transfer, zi);
	double T_nu = interpol(transfer_nu_klong, zlist_transfer, Nz_transfer, zi);
	double OmG= OmGbar * (1.+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	double Omnu_massless= Omnu_masslessbar * (1.+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR, OmR2 = OmG + Omnu_massless;
//and initial perturbation in nu1	(massive) and pressure
	double pnu1_z = interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu1_z = pnu1_z/rhonu1_z; //equation of state
	double T_nu1_z=interpol(transfer_nu1_klong, zlist_transfer, Nz_transfer, zi);
	double delta_nu1_z = delta_long*T_nu1_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//and initial perturbation in nu2	(massive) and pressure
	double pnu2_z = interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu2_z = pnu2_z/rhonu2_z; //equation of state
	double T_nu2_z=interpol(transfer_nu2_klong, zlist_transfer, Nz_transfer, zi);
	double delta_nu2_z = delta_long*T_nu2_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//and initial perturbation in nu3	(massive) and pressure
	double pnu3_z = interpol_cubic(zmin_EoS, dz_EoS, plist3_EoS, Nz_EoS, zi); //pressure of nu1 at z
	double wnu3_z = pnu3_z/rhonu3_z; //equation of state
	double T_nu3_z=interpol(transfer_nu3_klong, zlist_transfer, Nz_transfer, zi);
	double delta_nu3_z = delta_long*T_nu3_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next


//for the next values, to do Heun's method. We assume Omega_L does not change. Symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	double z_next = zi * exp(zstep_log);//first step
	OmGbar= cosmo->OmegaG * pow(1.+z_next,4.);
	Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = cosmo->OmegaM * pow(1.+z_next,3.);
	OmL = cosmo->OmegaL;

	rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, z_next);
	rho_ratio_nu1=rhonu1_z/rhonu1_0;
	Omnu1bar_z =  cosmo->Omeganu1 * rho_ratio_nu1;

	rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, z_next);
	rho_ratio_nu2=rhonu2_z/rhonu2_0;
	Omnu2bar_z =  cosmo->Omeganu2 * rho_ratio_nu2;

	rhonu3_z=interpol_cubic(zmin_EoS, dz_EoS, rholist3_EoS, Nz_EoS, z_next);
	rho_ratio_nu3=rhonu3_z/rhonu3_0;
	Omnu3bar_z = cosmo->Omega_extra * rho_ratio_nu3;

	double H2 = cosmo->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z + Omnu3bar_z); // H(z_next)
	double dE, dE2 = (H2-H)/H2/(zi*zstep_log);

//and the w and c_s^2 of nus at first step we calculate.
//for nu1
	pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, z_next);
	double wnu1_z2=pnu1_z/rhonu1_z;
	double d_wnu1_z = (wnu1_z2-wnu1_z)/(z_next-zi);
	double csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/(3.0*(1+wnu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//for nu2
	pnu2_z=interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, z_next);
	double wnu2_z2=pnu2_z/rhonu2_z;
	double d_wnu2_z = (wnu2_z2-wnu2_z)/(z_next-zi);
	double csq_ad_nu2_z = wnu2_z2 + d_wnu2_z/(3.0*(1+wnu2_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//for nu3
	pnu3_z=interpol_cubic(zmin_EoS, dz_EoS, plist3_EoS, Nz_EoS, z_next);
	double wnu3_z2=pnu3_z/rhonu3_z;
	double d_wnu3_z = (wnu3_z2-wnu3_z)/(z_next-zi);
	double csq_ad_nu3_z = wnu3_z2 + d_wnu3_z/(3.0*(1+wnu3_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).






	double Rpp1, Rpp2; //d^2R(z)/dz^2

	int i_solution=0; //we only check for the solution of neutrino clustering every few steps
	double Mhalo_nu1_Mpc=0, Mhalo_nu2_Mpc=0, Mhalo_nu3_Mpc=0; //M of the nu1, nu2, and nu3 haloes
	int collapse_steps= (int) npoints/Nz_solution;	//how many steps we take before saving R(z) and/or updating Mnu
	//we save manually the first value of R and mass to avoid running the loop at i=0.
  R_solution[i_solution]=R1;
	Mhalo_nu1_Mpc = Mnu1_solution[i_solution];
	Mhalo_nu2_Mpc = Mnu2_solution[i_solution];
	Mhalo_nu3_Mpc = Mnu3_solution[i_solution];
  i_solution++;




///////////////////////////////////////////////////////////
////    here we solve for the collapse								////
/////////////////////////////////////////////////////////
	for(i=1; i<npoints-1 && R2>0.; i++){
		//we update the values from the previous step
			R1 = R2;
			Rp1 = Rp2;
		//including all z-dependent quantities, for Heun's method. Remember 2 means next step.
			z=z_next;
			H=H2;
			dE=dE2;
			OmR=OmR2;
			wnu1_z=wnu1_z2;
			wnu2_z=wnu2_z2;
			wnu3_z=wnu3_z2;

//we save R and read Mnu, only every few steps since it doesn't vary much.
			if (i % collapse_steps == 0){
				R_solution[i_solution]=R1;
				Mhalo_nu1_Mpc = Mnu1_solution[i_solution];
				Mhalo_nu2_Mpc = Mnu2_solution[i_solution];
				Mhalo_nu3_Mpc = Mnu3_solution[i_solution];
				i_solution++;
		    if(debug_mode > 1){
		      printf("z=%.1le, R=%.1le, Mnu1=%1le , Mnu2=%1le , Mnu3=%1le \n\n", z, R_solution[i_solution], Mhalo_nu1_Mpc, Mhalo_nu2_Mpc, Mhalo_nu3_Mpc);
		    }
		  }


		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.

//these are the terms that go inside the sum, define them outside for clarity.
		double Onu1 = (1.0+3.0*wnu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
		double Onu2 = (1.0+3.0*wnu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;
		double Onu3 = (1.0+3.0*wnu3_z + delta_nu3_z*(1.0+3.0*csq_ad_nu3_z))*Omnu3bar_z;

	//we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- cosmo->Mhalo_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu1_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu2_Mpc/R1/R1/pow(H*(1.+z),2.)
		- Mhalo_nu3_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR + Onu1 + Onu2 + Onu3 - 2*OmL)*cosmo->H0_Mpc*cosmo->H0_Mpc/pow(H*(1.+z),2.); //R''(R1,t1)



//we update all z-dependent stuff to the next z
		OmGbar= cosmo->OmegaG * pow(1.+z_next,4.);
		T_gamma = interpol(transfer_gamma_klong, zlist_transfer, Nz_transfer, z_next);
		OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)

		Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+z_next,4.);
		T_nu = interpol(transfer_nu_klong, zlist_transfer, Nz_transfer, z_next);
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)

		OmRbar = OmGbar + Omnu_masslessbar;
		OmR2 = OmG + Omnu_massless;

		OmM = cosmo->OmegaM * pow(1.+z_next,3.);


		rhonu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist1_EoS, Nz_EoS, z_next);
		rho_ratio_nu1=rhonu1_z/rhonu1_0;
		Omnu1bar_z = cosmo->Omeganu1 * rho_ratio_nu1;
		T_nu1_z=interpol(transfer_nu1_klong, zlist_transfer, Nz_transfer, z_next);
		delta_nu1_z = delta_long*T_nu1_z/T_matter; //long-wavlength neutrino perturbation
		pnu1_z=interpol_cubic(zmin_EoS, dz_EoS, plist1_EoS, Nz_EoS, z_next);
		wnu1_z2=pnu1_z/rhonu1_z;
		d_wnu1_z = (wnu1_z2-wnu1_z)/zstep_lin;
		csq_ad_nu1_z = wnu1_z2 + d_wnu1_z/3.0*(1.0+z_next)/(1.0+wnu1_z2); //adiabatic sound speed squared.

		rhonu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist2_EoS, Nz_EoS, z_next);
		rho_ratio_nu2=rhonu2_z/rhonu2_0;
		Omnu2bar_z = cosmo->Omeganu2 * rho_ratio_nu2;
		T_nu2_z=interpol(transfer_nu2_klong, zlist_transfer, Nz_transfer, z_next);
		delta_nu2_z = delta_long*T_nu2_z/T_matter; //long-wavlength neutrino perturbation
		pnu2_z=interpol_cubic(zmin_EoS, dz_EoS, plist2_EoS, Nz_EoS, z_next);
		wnu2_z2=pnu2_z/rhonu2_z;
		d_wnu2_z = (wnu2_z2-wnu2_z)/zstep_lin;
		csq_ad_nu2_z = wnu2_z2 + d_wnu2_z/3.0*(1.0+z_next)/(1.0+wnu2_z2); //adiabatic sound speed squared.

		rhonu3_z=interpol_cubic(zmin_EoS, dz_EoS, rholist3_EoS, Nz_EoS, z_next);
		rho_ratio_nu3=rhonu3_z/rhonu3_0;
		Omnu3bar_z = cosmo->Omega_extra * rho_ratio_nu3;
		T_nu3_z=interpol(transfer_nu3_klong, zlist_transfer, Nz_transfer, z_next);
		delta_nu3_z = delta_long*T_nu3_z/T_matter; //long-wavlength neutrino perturbation
		pnu3_z=interpol_cubic(zmin_EoS, dz_EoS, plist3_EoS, Nz_EoS, z_next);
		wnu3_z2=pnu3_z/rhonu3_z;
		d_wnu3_z = (wnu3_z2-wnu3_z)/zstep_lin;
		csq_ad_nu3_z = wnu3_z2 + d_wnu3_z/3.0*(1.0+z_next)/(1.0+wnu3_z2); //adiabatic sound speed squared.

//we update the Onu terms that go inside the integral.
		Onu1 = (1.0+3.0*wnu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
		Onu2 = (1.0+3.0*wnu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;
		Onu3 = (1.0+3.0*wnu3_z + delta_nu3_z*(1.0+3.0*csq_ad_nu3_z))*Omnu3bar_z;


		H2  = cosmo->H0_Mpc * sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z + Omnu3bar_z);// H(z_next)

		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value at the end.




//these are too annoying to keep even with debug_mode, activate manually if you want:
//		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, cosmo->H0_Mpc*E_LCDM(cosmo, z), dE, dlogE_dz_LCDM(cosmo, z));
//		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);


//now we evolve R and R'.

		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;


//we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- cosmo->Mhalo_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu1_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu2_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- Mhalo_nu3_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR2 + Onu1 + Onu2 + Onu3 - 2*OmL)*cosmo->H0_Mpc*cosmo->H0_Mpc/pow(H2*(1.+z_next),2.); //R''(R2,t2)


//and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);





	}


		return z;

}






int findmass_1nu (Cosmology *cosmo, double mass_0, double temp_0, double zmin_EoS, double dz_EoS, long Nz_EoS,
	double* rho1list_EoS, double* rho2list_EoS, double* rho3list_EoS, double* etalist,
	double logz_solution_min, double dlogz_solution, long Nz_solution, double* R_solution, double* Mnu_solution){
//this function finds Mnu(t) given the R(t) input. We use the BKT approximation for neutrino clustering.
//mass_0 and temp_0 are the mass and temperature (at z=0) for the particle (neutrino or extra), we specify it to distinguish which case we do.


  double precision_findmass = fmin(fmax(precision_clustering,0.5),2.0); //precision parameter to regulate speed vs. precision.

	long nzpoints=400*precision_findmass;//how many points we integrate over in z (or t)
	double zstep;//= (z-zi)/(nzpoints-1.); //linear binning in z, since we're interested in early times

  int Nx = (int) 12*precision_findmass; //how many positions
	int Nq = 10; //how many momenta, 10 is good.
	int Nmu = (int) 8*precision_findmass; //how many x \cdot p cosines.

	double ***f1;	//this is f1(x,p,mu) at different times
	f1=allocate_3D_array(Nx,Nq,Nmu);

	double OmRbar, OmM, OmGbar, Omnu_masslessbar, Omnu1, Omnu2, Omnu3; //average for all of them, nu3=extra.
  double OmL=cosmo->OmegaL;



	double rho_ratio_nu1; //ratio of rho(z)/rho(z=0)
	double rho1_0=rho1list_EoS[0] + (cosmo->Omeganu1 == 0) * 1;//rho at z=0. units are K^4. add 1 if not present to avoid NaNs
	double rho1_z;

	double rho_ratio_nu2; //ratio of rho(z)/rho(z=0)
	double rho2_0=rho2list_EoS[0] + (cosmo->Omeganu2 == 0) * 1;//rho at z=0. units are K^4
	double rho2_z;

	double rho_ratio_nu3; //ratio of rho(z)/rho(z=0)
	double rho3_0=rho3list_EoS[0] + (cosmo->Omega_extra == 0) * 1;//rho at z=0. units are K^4
	double rho3_z;


  //all temperatures, mases, and momenta in eV unless otherwise specified (we'd say: _Mpc):
	double temp = temp_0*KtoeV;//at z=0 in eV, we will use comoving momenta
  double mass = mass_0;//
	//we do each species clustering every time.




	long i, iq, ix, imu, iz;



  int lengthname=200;
 	char *filename; //To open files.
 	filename=(char *)malloc((lengthname+1)*sizeof(char));
	FILE *fp; //to save Mnu(z) in case
	if(debug_mode>0){
	  lengthname=sprintf(filename,"tests/Mnu_z-%d.dat",cosmo->file_tag); //we save escape fraction in a file
	  fp=fopen(filename,"w");
	}


//we first find \eta(z), the superconformal time, \eta(z) = \int dz/(H(z)*(1+z))

	double H;


  double Mnusmooth;

	double eta_in, eta_out; //eta(t) and eta(t'), we integrate over t'.
	double z_in, z_out; //z insider integral and outside


  double factor = 2.0*mass/temp; //x, q, and z-independent factor, 2 because 2 degrees of freedom in f0. unitless.
  double factorx; //1/x^2
  double factorq; //df0/dq = e^q/T/(1+e^q/T)^2

  double x, q, mu;
  double parenthesis; //(x-xc)^2 (vectors).
  double R_halo, x_halo; //Rhalo(t) and xhalo=Rhalo(t)/(a(t));

  double integrand, xc;
  double deltaM, Mhalo_smooth;

  double xmin,xmax,dx;
  double qmin,qmax,dq_log; //q is logspaced around temp.
  double mumin,mumax,dmu;

  xmin=0.; //xmax will be z-dependent, since we go up to R_halo.

  mumin=-1.;//cosine always goes between -1 and 1, we do linear spacing.
  mumax=1.;
  dmu=(mumax-mumin)/(Nmu-1.);

  qmin=temp/20.;
  qmax=temp*10.;
  dq_log = log(qmax/qmin)/(Nq-1.); //logspaced is better.

  z_out=zi; //we set it first at z_initial~200 at the begining, for the first test.

  long Nz_solution_max = (log(cosmo->z_collapse)-logz_solution_min)/dlogz_solution; //to make sure we don't continue after collapse.



	for (i=0;(i<Nz_solution) && (i<Nz_solution_max);i++){ //we stop when we cross z_collapse, which is before zf_code=z_collapse/2. We keep the same array than in collapse code for convenience.
		Mnu_solution[i] = 0.0;
		z_out = exp(logz_solution_min+dlogz_solution*i); //that's where R_solution is saved, and where we save Mnu_solution.
		eta_out=interpol_cubic(zmin_EoS, dz_EoS, etalist, Nz_EoS, z_out);

//    xmax = (1+z_out)*interpol_cubic(logz_solution_min, dlogz_solution, R_solution, Nz_solution, log(z_out)); //comoving radius of halo at z_out.
    xmax = (1+z_out)*R_solution[i];
    dx=(xmax-xmin)/(Nx-1.);

		if (debug_mode>1) 	printf("i=%ld, z=%.1le, eta=%.1le \n",i,z_out,eta_out);


    for(ix=1;ix<Nx;ix++){//comoving position, we avoid x=0.
      x=xmin+dx*ix;
      factorx=1.0/x/x;
      for(iq=0;iq<Nq;iq++){//comoving momentum
        q=qmin*exp(dq_log*iq);
        factorq = exp(q/temp)/pow(1.+exp(q/temp),2.); //assumming FERMIONS.
        for(imu=0;imu<Nmu;imu++){//cosine between them
          mu=mumin + dmu*imu;
          f1[ix][iq][imu] = 0.;//just in case.

    		    for (iz=0;iz<nzpoints-1;iz++){
                zstep = (zi-z_out)/(nzpoints-1.); //we integrate z_inside (z_in) from z_initial (z_i) to z_outside (z_out).
    		        z_in = z_out + zstep * iz;
    		        eta_in=interpol_cubic(zmin_EoS, dz_EoS, etalist, Nz_EoS, z_in);

                OmGbar= cosmo->OmegaG * pow(1.+z_in,4.);
                Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+z_in,4.);
                OmRbar = OmGbar + Omnu_masslessbar;
                OmM = cosmo->OmegaM * pow(1.+z_in,3.);

								rho1_z=interpol_cubic(zmin_EoS, dz_EoS, rho1list_EoS, Nz_EoS, z_in);
                rho_ratio_nu1=rho1_z/rho1_0;
								Omnu1 =  cosmo->Omeganu1 * rho_ratio_nu1;

								rho2_z=interpol_cubic(zmin_EoS, dz_EoS, rho2list_EoS, Nz_EoS, z_in);
                rho_ratio_nu2=rho2_z/rho2_0;
								Omnu2 =  cosmo->Omeganu2 * rho_ratio_nu2;

								rho3_z=interpol_cubic(zmin_EoS, dz_EoS, rho3list_EoS, Nz_EoS, z_in);
                rho_ratio_nu3=rho3_z/rho3_0;
								Omnu3 =  cosmo->Omega_extra * rho_ratio_nu3;

								H = cosmo->H0_Mpc * sqrt(OmL + OmM + OmRbar + Omnu1 + Omnu2 + Omnu3);


    		        xc=(eta_out-eta_in)*q/mass;
                R_halo = interpol_cubic(logz_solution_min, dlogz_solution, R_solution, Nz_solution, log(z_in));
                x_halo = (1+z_in) * R_halo;

                Mhalo_smooth = cosmo->H0_Mpc*cosmo->H0_Mpc*pow(R_halo,3.)*OmM/2.;
                deltaM = cosmo->Mhalo_Mpc - Mhalo_smooth;//Mhalo-Mhalo_smooth. In Mpc.
                integrand = zstep/H * deltaM * (xc/x-mu);

                parenthesis = x*x + xc*xc - 2*x*xc*mu;

                if(parenthesis<x_halo*x_halo){
                  f1[ix][iq][imu]+= integrand*pow(x/x_halo,3.);
                }
                else{
                  f1[ix][iq][imu]+= integrand/pow(parenthesis/x/x,3/2.);
                }

//               printf("z_in=%.1le, deltaM=%.1le, xc=%.1le, x_halo=%.1le, parenthesis=%.1le \n",z_in,deltaM,xc,x_halo,parenthesis);

            }//z-loop
            f1[ix][iq][imu] *= factorq * factorx * factor;

            Mnu_solution[i] += (2*PI*x*x*dx) * (dmu) * (dq_log * q * q * q /(2.*PI*PI) ) * f1[ix][iq][imu]; //\int d^3x d^3q = dx x^2 (2PI) dmu dp p^2/(2PI^2).
        }//mu-loop
      }//q-loop
    }//x-loop

    Mnu_solution[i]*=mass*eVtoMpc/pow(hbareVMpc,3.); //to get Mnu in Mpc we need hbar^3 (that is why there is a 1/(2pi)^3 in p integral). There are 3 q in eV to convert to Mpc-1 (since x is in Mpc),
    //we have to use G to convert mass to Mpc -> eVtoMpc.



//smooth mass of ALL non-CDM components. Just for comparison.
  Mnusmooth = cosmo->H0_Mpc*cosmo->H0_Mpc*pow(R_solution[i],3.)*(cosmo->Omeganu1+cosmo->Omeganu2+cosmo->Omega_extra)/2.*pow(1+z_out,3.);



	if(debug_mode>0){
		printf("z=%.1le, Mnu/Mhalo=%1le \n",z_out,Mnu_solution[i]/cosmo->Mhalo_Mpc);
		fprintf(fp, "%le  %le  %le \n",z_out, Mnu_solution[i]/cosmo->Mhalo_Mpc, (Mnu_solution[i]+Mnusmooth)/cosmo->Mhalo_Mpc);
	}



}


	if(debug_mode>0){
  	fclose(fp);
	}

  	free(filename);






	free_3D_array(f1,Nx,Nq);

  return 1;

}
