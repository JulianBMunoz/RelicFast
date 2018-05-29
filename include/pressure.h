#ifndef __PRESS_H_INCLUDED__
#define __PRESS_H_INCLUDED__


#include <stdio.h>
#include <math.h>

#include "common.h"





double pressure_WDM(double mass, double Temp);

double density_WDM(double mass, double Temp);

double EoS_WDM(double mass, double Temp);






#endif
