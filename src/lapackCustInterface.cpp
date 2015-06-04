//
//  lapackCustInterface.cpp
//  
//
//  Created by Joshua Tyler Cantin on 2015-06-03.
//
//

#include "lapackCustInterface.h"

using namespace std;

//This function takes a Hamiltonian of size N by N and returns
//  the number of eigenvalues in addition to placing the
//  eigenvalues into eigvals, which must be an array of doubles
//  of size N. The eigenvalues are determined using the
//  dsyevr routine from Lapack.
//The Hamiltonian must be an array of doubles of size (N*N) that
//  is stored using the column major format.
//  ie. H(i,j) = Hami[i + N*j]
void dsyevrEigvalsInterface(double* Hami, int N, int* numEigvals, double* eigvals)
{

    // print out Hamiltonian
    //            for (i=0; i<N; i++) {
    //                for (j=0; j<N; j++) {
    //                    cout << Hami[i][j] ;
    //                }
    //                cout << endl;
    //            }


    //Prepare Lapack dsyevr call
    lapack_int N_l;
    lapack_int dummy_int = 1;

    //Output parameters
    lapack_int numEval,isuppz;
    double** dummy_matrix = new double* [1];
    dummy_matrix[0] = new double [1];

    lapack_int ldz = N;
    lapack_int lda = N;
    int info;

    double abstol;

    #ifdef OSX_LAPACKE
    //Get the machine epsilon
    abstol = LAPACKE_dlamch('E');
    #endif

    #ifdef LINUX_MKL
    //Arbitrary value
    abstol = 1.0E-10;
    #endif

    N_l = N;

    //I compared the below to the results from dstevr, and they are identical, which were identical (relDiff < 1E-14) to the python values
    info = LAPACKE_dsyevr(LAPACK_COL_MAJOR,'N','A','U',N_l,Hami,lda,dummy_int,dummy_int,dummy_int,dummy_int,abstol, &numEval,eigvals,*dummy_matrix,ldz,&isuppz);

    if (info == 0) {
        cout << "LAPACKE_dsyevr completed successfully." << endl;
        }
    else {
        cout << "!!!LAPACKE_dsyevr ERROR!!!" << endl;
        cout << "dstevr info: " << info << endl;
        exit (EXIT_FAILURE);
    }
    
    *numEigvals = numEval;
    
//    cout << "NumEvals: " << *numEigvals << endl;
//    cout << "Eigenvalues: " << endl;
//    cout << setprecision(presInt);
//    int i;
//    for (i=0; i<*numEigvals; i++) {
//        cout << eigvals[i] << endl;
//    }
}