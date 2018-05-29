////////////////////////////////////////////////////////////////////////////////////////
//
//	Auxiliar functions for RelicFast.
//	By Julian B Munoz (05/2018 @Harvard)
//
////////////////////////////////////////////////////////////////////////////////////////

#include "auxiliar.h"

//////////////////////////////////
/////	ARRAY UTILITIES:	/////
////////////////////////////////


double *allocate_1D_array(int n1){
   double *array = (double *) calloc(n1, sizeof(double));
   if (array == NULL) {
      printf("Error in allocate_1D_array with length %d\n", n1);
      exit(1);
    }
   return array;
}

int *allocate_1D_array_int(int n1){
   int *array = (int *) calloc(n1, sizeof(int));
   if (array == NULL) {
      printf("Error in allocate_1D_array with length %d\n", n1);
      exit(1);
    }
   return array;
}



double **allocate_2D_array(int n1, int n2){
   double **array = (double **) calloc(n1, sizeof(double *));
   if (array == NULL){
      printf("Error in allocate_2D_array with lengths %d, %d \n", n1, n2);
      exit(1);
    }
   for (int i = 0; i < n1; i++){
     array[i] = allocate_1D_array(n2);
   }

   return array;
}

void free_2D_array(double **array, int n1){
    for (int i = 0; i < n1; i++){
      free(array[i]);
    }
    free(array);
}


double ***allocate_3D_array(int n1, int n2, int n3){
   double ***array = (double ***) calloc(n1, sizeof(double **));
   if (array == NULL){
      printf("Error in allocate_2D_array with lengths %d, %d, %d \n", n1, n2, n3);
      exit(1);
    }
   for (int i = 0; i < n1; i++){
     array[i] = allocate_2D_array(n2,n3);
   }

   return array;
}


void free_3D_array(double ***array, int n1, int n2){

    for (int i = 0; i < n1; i++){
    	free_2D_array(array[i],n2);
    }

    free(array);
}






void free_cosmology(Cosmology *cosmo){
//we free a cosmology structure pointer, and all the pointers inside it.

  free(cosmo->Mhalo_array);
  free(cosmo->j_collapse_array);
  free(cosmo->z_collapse_array);
  free(cosmo->klong_list_input);

  free(cosmo);

}



void do_check(int var_to_check){
//this just checks if variables have the expected value or no.
  if(var_to_check <= 0){
    printf("do_check() found unexpected value. Revise check_X variables. \n");
  }
}

















//////////////////////////////////////
/////	MATHEMATICAL FUNCTIONS:	/////
////////////////////////////////////


double interpol(double data[], double xtab[],int length, double x){
  //Linear interpolation algorithm

    double d1,res,v1,v2,x1,xmin,xmax;
    xmin=xtab[0];
    xmax=xtab[length-1];


    if((x<xmin)) {
        printf("interpol(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
		x=xmin;
    }

    if((x>xmax)) {
        printf("interpol(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
		x=xmax;
    }



	int i=floor(length/2),count;

    for(count=2; i<length-1 && (x>=xtab[i+1] || x<xtab[i]);++count){
    	if (x>=xtab[i+1])
    		i+=1;
//    		i=i+length/(count);
    	else if (x<xtab[i])
			i-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol \n");
    }

    x1=xtab[i];
    v1=data[i];


    double step ;

//we avoid the edge by interpolating downwards if it's there.
    if(i==length-1){
		i--;
    }


    step = xtab[i+1]-xtab[i];
    v2=data[i+1];

    d1=(v2-v1)/step;

    res= v1 + d1*(x-x1);


//     printf("%e,%e, %f, %d \n",v1,v2,step,length);
//     printf("x=%e,v1=%e, res=%e \n",x1,v1,res);

    return res;

}




double interpol_2D(double **data, double xtab[], int lengthx, double ytab[], int lengthy, double x, double y){
  //Linear interpolation algorithm in 2D
  //assumes both lists are ordered (upwards).

    double d1x,d1y,res,v1,v2x,v2y,x1,y1,xmin,xmax, ymin,ymax;


    xmin=xtab[0];
    xmax=xtab[lengthx-1];
    if((x<xmin)) {
        printf("interpol_2D(x,y) is out of range in x. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
		x=xmin;
    }

    if((x>xmax)) {
        printf("interpol_2D(x,y) is out of range in x. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
		x=xmax;
    }


    ymin=ytab[0];
    ymax=ytab[lengthy-1];
    if((y<ymin)) {
        printf("interpol_2D(x,y) is out of range in y. Range is y=(%f %f) and y=%le \n",ymin,ymax,y);
//        return 0;
		y=ymin;
    }

    if((y>ymax)) {
        printf("interpol_2D(x,y) is out of range in y. Range is y=(%f %f) and y=%le \n",ymin,ymax,y);
//        return 0;
		y=ymax;
    }

    int i, j, count;


	i=floor(lengthx/2);
    for(count=2; i<lengthx-1 && (x>=xtab[i+1] || x<xtab[i]);++count){
    	if (x>=xtab[i+1])
    		i+=1;
//    		i=i+length/(count);
    	else if (x<xtab[i])
			i-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol_2D in x \n");
    }


	j=floor(lengthy/2);
    for(count=2; j<lengthy-1 && (y>=ytab[j+1] || y<ytab[j]);++count){
    	if (y>=ytab[j+1])
    		j+=1;
//    		i=i+length/(count);
    	else if (y<ytab[j])
			j-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on interpol_2D in y \n");
    }

    v1=data[i][j];
    x1=xtab[i];
    y1=ytab[j];


    double stepx, stepy;

//we avoid the edge by interpolating downwards if it's there.
    if(i==lengthx-1){
		i--;
    }

    if(j==lengthy-1){
		j--;
	}


    stepx = xtab[i+1]-xtab[i];
    v2x=data[i+1][j];
	  stepy = ytab[j+1]-ytab[j];
    v2y=data[i][j+1];


    d1x=(v2x-v1)/stepx;
    d1y=(v2y-v1)/stepy;

    res= v1 + d1x*(x-x1) +  d1y*(y-y1);

//    printf(" ix=%d (x-x1)=%.3le, iy=%d (y-y1)=%.3le \n", i ,(x-x1) ,j , (y-y1));


//	printf("%e,%e, %f, %d \n",v1,v2,step,length);
//	printf("x=%e, y=%e, v1=%e, res=%e \n",x1,y1,v1,res);

    return res;

}






long find_value(long length, double xtab[], double x){
  //finds the position in the list xtab closest to x

    double xmin,xmax;
    xmin=xtab[0];
    xmax=xtab[length-1];


    if((x<xmin)) {
        printf("find_value(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
		      return 0;
    }

    if((x>xmax)) {
        printf("find_value(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
		      return length -1;
    }

    if(x==xmin){//same but with not warnings.
      return 0;
    }
    if(x==xmax){
      return length-1;
    }



	long i=floor(length/2),count;

    for(count=2; x>=xtab[i+1] || x<xtab[i];++count){
    	if (x>=xtab[i+1])
    		i+=1;
//    		i=i+length/(count);
    	else if (x<xtab[i])
			i-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on find_value \n");
    }

    if (xtab[i+1]-x < x - xtab[i]){
     i++; //if it's closer to the one above add one.
    }



    return i;

}

long find_value_reverse(long length, double xtab[], double x){
  //finds the position in the list xtab closest to x
  //here we assume the xlist is reversed in order (goes down)

    double xmin,xmax;
    xmin=xtab[length-1];
    xmax=xtab[0];


    if((x<xmin)) {
        printf("find_value_reverse(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
		x=xmin;
    }

    if((x>xmax)) {
        printf("find_value_reverse(x) is out of range. Range is x=(%f %f) and x=%le \n",xmin,xmax,x);
//        return 0;
		x=xmax;
    }



	long i=floor(length/2),count;

    for(count=2; x<=xtab[i+1] || x>xtab[i];++count){
    	if (x<=xtab[i+1])
    		i+=1;
//    		i=i+length/(count);
    	else if (x>xtab[i])
			i-=1;
//    		i=i-length/(count);
    	else
    		printf("ERROR on find_value_reverse \n");
    }

    if (xtab[i+1]-x > x - xtab[i]){
     i++; //if it's closer to the one above add one.
    }



    return i;

}




double nintegrate(double data[],double xtab[], int length){

	double sum = 0;
	long double step; //To avoid possible overflows when two values are too close.

	long i;

    for (i = 0; i < length-2; ++i) {
    	step = xtab[i+1]-xtab[i];
        sum +=  0.5*(data[i] + data[i+1]) * step;
//        printf("i= %d, sum=%lf \n",i,sum);
    }

return sum;

}


int triangle(double l1,double l2,double l3){
	//Checking for triangle identities. Returns 0 if not valid and 1 if valid.
	int res=1;

	if((l1+l2-l3<0) || (l1-l2+l3<0) || (-l1+l2+l3<0)){
		res=0;
		}

	return res;
}


double interpol_cubic(double x0, double dx, double *ytab, int Nx, double x) {
  //cubic interpolation routine, assumes x equally spaced.
  long ix;
  double frac;

  if (Nx < 4) {
    fprintf(stderr, "Error: interpol_cubic: Table needs to be of dimension 4 at least\n");
    exit(1);
  }

  // Check if in range
  if (  (dx > 0 && (x<x0 || x>x0+dx*(Nx-1)))
      ||(dx < 0 && (x>x0 || x<x0+dx*(Nx-1))) ) {
	fprintf(stderr, "Error: interpol_cubic: value out of range. Range is (%le,%le) and x=%le. \n",x0,x0+dx*(Nx-1),x);
 //   exit(1);
  }

  // Identify location to interpolate
  ix = (long)floor((x-x0)/dx);
  if (ix<1) ix=1;
  if (ix>Nx-3) ix=Nx-3;
  frac = (x-x0)/dx-ix;
  ytab += ix-1;

  // Return value
  return(
    -ytab[0]*frac*(1.-frac)*(2.-frac)/6.
    +ytab[1]*(1.+frac)*(1.-frac)*(2.-frac)/2.
    +ytab[2]*(1.+frac)*frac*(2.-frac)/2.
    -ytab[3]*(1.+frac)*frac*(1.-frac)/6.
  );
}







void reverse(double *list,int N){
//reverse a list. Pretty straightforward.
	double aux;
	int j;
	for(j=0;j<=floor(N/2);j++){
		aux=list[j];
		list[j]=list[N-1-j];
		list[N-1-j]=aux;
	}

}




//for progressbar:
//PBWIDTH and PBSTR defined in .h
void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}











//////////////////////////////////////
/////	COSMOLOGICAL FUNCTIONS:	/////
////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//unused, although they are useful for checks://
double E_LCDM(Cosmology *cosmo, double z){
//computes E(z) = H(z)/H0 for given cosmological parameters.

	double E1 = sqrt(cosmo->OmegaL + cosmo->OmegaM * pow(1.+z,3.)+ cosmo->OmegaR * pow(1.+z,4.));

	return E1;

}
double dlogE_dz_LCDM(Cosmology *cosmo, double z){
//computes dE(z)/dz * 1/E(z)  for given cosmological parameters.
// for H(z)= Sqrt[OmL + OmM (1 + z)^3 + OmR (1 + z)^4].
//H'(z)/H(z) = (3 OmM (1 + z)^2 + 4 OmR (1 + z)^3)/(2 (OmL + OmM (1 + z)^3 + OmR (1 + z)^4))
//saves the trouble of doing it numerically

	double dE1 = (3.*cosmo->OmegaM*(1+z)*(1+z) + 4.*cosmo->OmegaR*pow(1.+z,3.))/(2.*(cosmo->OmegaL + cosmo->OmegaM *pow(1.+z,3.) + cosmo->OmegaR *pow(1.+z,4.)));

	return dE1;

}
////////////////////////////////////////////////////////////////////////





void findeta(Cosmology *cosmo, double *etalist, double zmin_EoS, double dz_EoS, long Nz_EoS,
            double* rholist_nu1_EoS, double* rholist_nu2_EoS, double* rholist_extra_EoS){
  //finds superconformal time etalist, at z_EoS list. Given LCDM + neutrinos and extra species.
  //used for clustering of neutrinos in the BKT approximation.

  double z, H, rho_ratio_nu1, rho_ratio_nu2, rho_ratio_extra;
  double rho_nu1_z, rho_nu1_z0;
  double rho_nu2_z, rho_nu2_z0;
  double rho_extra_z, rho_extra_z0;

  double Omnuextra, OmM, OmL, OmGbar, OmRbar, Omnu_masslessbar;

  long i;

  rho_nu1_z0 = (cosmo->Omeganu1>0) * rholist_nu1_EoS[0] + (cosmo->Omeganu1==0) * 1; //we make sure that dividing by it does not blow anything up.
  rho_nu2_z0 = (cosmo->Omeganu2>0) * rholist_nu2_EoS[0] + (cosmo->Omeganu2==0) * 1;
  rho_extra_z0 = (cosmo->Omega_extra>0) * rholist_extra_EoS[0] + (cosmo->Omega_extra==0) * 1;

		etalist[0]=0.; //everything will be a function of eta-eta'. so initial value does not matter.

  OmL=cosmo->OmegaL; //z independent.

		for (i=1;i<Nz_EoS;i++){

			z=zmin_EoS + dz_EoS*i;

			OmGbar= cosmo->OmegaG * pow(1.+z,4.);
			Omnu_masslessbar= cosmo->Omeganu_massless * pow(1.+z,4.);
			OmRbar = OmGbar + Omnu_masslessbar;
			OmM = cosmo->OmegaM * pow(1.+z,3.);

			rho_nu1_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_nu1_EoS, Nz_EoS, z);
			rho_nu2_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_nu2_EoS, Nz_EoS, z);
			rho_extra_z=interpol_cubic(zmin_EoS, dz_EoS, rholist_extra_EoS, Nz_EoS, z);
			rho_ratio_nu1=rho_nu1_z/rho_nu1_z0;
			rho_ratio_nu2=rho_nu2_z/rho_nu2_z0;
			rho_ratio_extra=rho_extra_z/rho_extra_z0;

			Omnuextra =  cosmo->Omeganu1 * rho_ratio_nu1 + cosmo->Omeganu2 * rho_ratio_nu2 + cosmo->Omega_extra * rho_ratio_extra;


			H = cosmo->H0_Mpc* sqrt(OmL + OmM + OmRbar + Omnuextra); //H(zi) in Mpc

			etalist[i]=etalist[i-1] - dz_EoS / H * (1.+z);//d eta = d t/a^2, - sign because eta grows at low z (it's super-conformal time)
//      printf("etalist[i]=%.1le , ratioextra=%.1le \n",etalist[i],rho_ratio_extra);
    } //cubic interpolator for eta.


}



///////////////////////////////////////////////////////////////////////////////////////
//		Calculate sigma(M) and sigma_8 at some redshift
//
//	input Mhalo_Mpc and transfer-function data
///////////////////////////////////////////////////////////////////////////////////////
double getsigma_M(Cosmology *cosmo, double *k_transfer, double *transfer_function){

	int j, Nj=500; //number of kpoints. Log-spaced

	double kmin_log=log(k_transfer[0]);
	double kmax_log=log(k_transfer[length_transfer-1]);
	double kstep_log=(kmax_log-kmin_log)/(Nj-1.);

	double sigmaM_2=0.; //sigma^2 = \int dlog(k) k^3  W(k,M)^2 Pcc(k) / (2pi^2)

	double R = pow(cosmo->OmegaM/2.*cosmo->H0_Mpc*cosmo->H0_Mpc/cosmo->Mhalo_Mpc,-1./3); //OmegaM=Omegac+Omegab. in Mpc.


	double k, window, transfer, Powcc;// k, W(kR), T(k), Pmm(k)

	for(j=1;j<Nj-1;j++){
		k=exp(kmin_log+j*kstep_log);//in Mpc
		window=3.*(sin(k*R)/(k*R)-cos(k*R))/(k*R)/(k*R);
		transfer=interpol(transfer_function,k_transfer,length_transfer,k);
		Powcc=(2.*PI*PI)*cosmo->As*pow(k/kpivot,cosmo->ns-1.)/k/k/k*transfer*transfer;
		sigmaM_2+=kstep_log * k*k*k * window*window * Powcc / (2.*PI*PI);	//		printf("k=%le, Powcc=%le, transfer= %le \n",k,Powcc, transfer);
	}

	return sqrt(sigmaM_2);
}



double getsigma_8(Cosmology *cosmo, double *k_transfer, double *transfer_function){
	//for R=8 Mpc/h. It's sigma_8=0.83 for Planck2015 best-fit
	int j, Nj=500; //number of kpoints. Log-spaced

	double kmin_log=log(k_transfer[0]);
	double kmax_log=log(k_transfer[length_transfer-1]);
	double kstep_log=(kmax_log-kmin_log)/(Nj-1.);

	double sigmaM_2=0.; //sigma^2 = \int dlog(k) k^3  W(k,M)^2 Pcc(k) / (2pi^2)

	double R = 8./cosmo->h; //In Mpc/h.



	double k, window, transfer, Powcc;// k, W(kR), T(k), Pmm(k)

	for(j=1;j<Nj-1;j++){
		k=exp(kmin_log+j*kstep_log);//in Mpc
		window=3.*(sin(k*R)/(k*R)-cos(k*R))/(k*R)/(k*R);
		transfer=interpol(transfer_function,k_transfer,length_transfer,k);
		Powcc=(2.*PI*PI)*cosmo->As*pow(k/kpivot,cosmo->ns-1.)/k/k/k*transfer*transfer;
		sigmaM_2+=kstep_log * k*k*k * window*window * Powcc / (2.*PI*PI);		//		printf("k=%le, Powcc=%le, integrand= %le \n",k,Powcc, k*k * window*window * Powcc / (2.*PI*PI));
	}

	return sqrt(sigmaM_2);
}




double getsigma_M_z(Cosmology *cosmo, double z, double *z_transfer, double *k_transfer, double **transfer_function){
//calculates sigma(M) @z.
	int j, Nj=1000; //number of kpoints. Log-spaced

	double kmin_log=log(k_transfer[0]);
	double kmax_log=log(k_transfer[length_transfer-1]);
	double kstep_log=(kmax_log-kmin_log)/(Nj-1.);

	double sigmaM_2=0.; //sigma^2 = \int dlog(k) k^3  W(k,M)^2 Pcc(k) / (2pi^2)

	double R = pow(cosmo->OmegaM/2.*cosmo->H0_Mpc*cosmo->H0_Mpc/cosmo->Mhalo_Mpc,-1./3); //OmegaM=Omegac+Omegab. in Mpc.


	double k, window, transfer, Powcc;// k, W(kR), T(k), Pmm(k)

	for(j=1;j<Nj-1;j++){
		k=exp(kmin_log+j*kstep_log);//in Mpc
		window=3.*(sin(k*R)/(k*R)-cos(k*R))/(k*R)/(k*R);
		transfer=interpol_2D(transfer_function,z_transfer, Nz_transfer, k_transfer,length_transfer,z,k);
		Powcc=(2.*PI*PI)*cosmo->As*pow(k/kpivot,cosmo->ns-1.)/k/k/k*transfer*transfer;
		sigmaM_2+=kstep_log * k*k*k * window*window * Powcc / (2.*PI*PI);	//		printf("k=%le, Powcc=%le, transfer= %le \n",k,Powcc, transfer);
	}

	return sqrt(sigmaM_2);
}


double getsigma_8_z(Cosmology *cosmo, double z, double *z_transfer, double *k_transfer, double **transfer_function){
//calculates sigma(M) @z.
	int j, Nj=1000; //number of kpoints. Log-spaced

	double kmin_log=log(k_transfer[0]);
	double kmax_log=log(k_transfer[length_transfer-1]);
	double kstep_log=(kmax_log-kmin_log)/(Nj-1.);

	double sigmaM_2=0.; //sigma^2 = \int dlog(k) k^3  W(k,M)^2 Pcc(k) / (2pi^2)

	double R = 8./cosmo->h; //In Mpc/h.


	double k, window, transfer, Powcc;// k, W(kR), T(k), Pmm(k)

	for(j=1;j<Nj-1;j++){
		k=exp(kmin_log+j*kstep_log);//in Mpc
		window=3.*(sin(k*R)/(k*R)-cos(k*R))/(k*R)/(k*R);
		transfer=interpol_2D(transfer_function,z_transfer, Nz_transfer, k_transfer,length_transfer,z ,k );
		Powcc=(2.*PI*PI)*cosmo->As*pow(k/kpivot,cosmo->ns-1.)/k/k/k*transfer*transfer;
		sigmaM_2+=kstep_log * k*k*k * window*window * Powcc / (2.*PI*PI);	//		printf("k=%le, Powcc=%le, transfer= %le \n",k,Powcc, transfer);
	}

	return sqrt(sigmaM_2);
}
