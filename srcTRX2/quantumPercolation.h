//
//  quantumPercolation.h
//  
//
//  Created by Joshua Tyler Cantin on 2015-04-24.
//
//

#ifndef ____quantumPercolation__
#define ____quantumPercolation__

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

#include "lapackCustInterface.h"

#endif /* defined(____quantumPercolation__) */
