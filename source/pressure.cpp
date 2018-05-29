////////////////////////////////////////////////////////////////////////////////////
//
//	Code to compute the pressure and energy density of an additional non-CDM component
//  given its mass and temperature, assuming it's a thermal relic, you can get rho and p.
//	By Julian B Munoz (05/2018 @Harvard)
//
//////////////////////////////////////////////////////////////////////////////////


#include "pressure.h"

#define Nkpoints 30 //number of p (momentum) points. Log-spaced, 30 is enough.

double pressure_WDM(double mass, double Temp){
//calculates Pbar(mass,Temp) in Kelvin^4 as in Eq. (19) of Loverde14.
//for a WDM particle of mass mass (in eV) and temperature Temp (in K)

	double moverT=mass/(Temp*KtoeV); //mass in units of Temperature

	int j, Nj=Nkpoints;

	double deltaPmax=30.;
	double deltaPmin=10.;
	double pmin_log=log(1./deltaPmin);
	double pmax_log=log(1.*deltaPmax);//in units of Temp
	double pstep_log=(pmax_log-pmin_log)/(Nj-1.);

	double Pbar=0.; // P_nu = 2 (\int d^3p/(2pi^3) p^2/(3 sqrt(p^2+m^2)) * 1/(exp(p/Tnu)+1)).   //   	printf("%le %le %le \n",Toverm,pmin_log,pmin_log);

	double p, fermi, integrand;// p (in units of Temp), fermi factor f, integrand


	for(j=0;j<Nj-1;j++){
		p=exp(pmin_log+j*pstep_log);//in units of mnu
		fermi=1./(1.+exp(p)); //p in units of Tnu
		integrand= p*p/(3.*sqrt(p*p+moverT*moverT));
		Pbar+= pstep_log * p * p * p / (2.*PI*PI) * fermi * 2 * integrand; //in units of Temp^4, 2 because 2 degrees of freedom
//		printf("p/mnu=%le, fermi=%le, integrand= %le \n",p, fermi, integrand);
	}

	Pbar*=pow(Temp,4.); //in K^4

	return Pbar;
}


double density_WDM(double mass, double Temp){
//calculates rhobar(mass,Temp) in Kelvin^4 as in Eq. (19) of Loverde14.
//for a WDM particle of mass mass (in eV) and temperature Temp (in K)

	double moverT=mass/(Temp*KtoeV); //mnu in units of Tnu.


	int j, Nj=Nkpoints; //number of kpoints. Log-spaced, 30 seems enough.

	double deltaPmax=30.;
	double deltaPmin=10.;
	double pmin_log=log(1./deltaPmin);
	double pmax_log=log(1.*deltaPmax);//in units of Temp
	double pstep_log=(pmax_log-pmin_log)/(Nj-1.); //p will be logspaced in units of Tnu.

	double rhobar=0.; // rho_nu = 2 (\int d^3p/(2pi^3) sqrt(p^2+m^2) * 1/(exp(p/Tnu)+1)).

	double p, fermi, integrand;// p (in units of mnu), fermi factor f, integrand


	for(j=0;j<Nj-1;j++){
		p=exp(pmin_log+j*pstep_log);//in units of mnu
		fermi=1./(1.+exp(p)); //p in units of Tnu
		integrand= sqrt(p*p+moverT*moverT);
		rhobar+= pstep_log * p * p * p / (2.*PI*PI) * fermi * 2 * integrand; //in units of Temp^4, 2 because 2 degrees of freedom	//		printf("p/mnu=%le, fermi=%le, integrand= %le \n",p, fermi, integrand);
	}

	rhobar*=pow(Temp,4.); //in K^4

	return rhobar;
}




double EoS_WDM(double mass, double Temp){
//calculates w, Equation of state, of WDM with mass and Temp, (in eV and K, respectively)

	double w=pressure_WDM(mass,Temp)/density_WDM(mass,Temp);



	return w;
}
