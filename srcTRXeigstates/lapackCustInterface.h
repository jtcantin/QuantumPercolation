//
//  lapackCustInterface.h
//  
//
//  Created by Joshua Tyler Cantin on 2015-06-03.
//
//

#ifndef ____lapackCustInterface__
#define ____lapackCustInterface__

#include <cstdlib>
#include <iostream>
#include <new>

#include <iomanip>
#include <ctime>
#include <cmath>

#ifdef LINUX_MKL
#include "mkl.h"
//#include "mkl_lapack.h"
#include "mkl_lapacke.h"
#endif

#ifdef OSX_LAPACKE
#include <lapacke.h>
#endif

void dsyevrEigvalsInterface(double* Hami, int N, int* numEigvals, double* eigvals);
void dsyevrEigstatesInterface(double* Hami, int N, int* numEigvals, double* eigvals, double* eigStateMatrix);


#endif /* defined(____lapackCustInterface__) */
