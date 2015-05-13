//
//  quantumPercolation.cpp
//  
//
//  Created by Joshua Tyler Cantin on 2015-04-24.
//
//

#include "quantumPercolation.h"

using namespace std;

int main(int argc, char **argv) {
    cout << "Hello World." << endl;
    
    //Make Hamiltonian
    int N;
    double E;
    double t;
    int randTest;
    
    
    N = atoi(argv[1]);
    E = atof(argv[2]);
    t = atof(argv[3]);
    randTest = atoi(argv[4]);
    
    double** Hami;
    
    Hami = new double* [N];
    
    int i,j;
    
    for (i=0; i<N; i++) {
        Hami[i] = new double [N];
    }
    
    //Add in t and on-site energy
    
    double LO, HI;
    if (randTest == 1) {
        srand (static_cast <unsigned> (time(0)));
        LO = -1*E;
        HI = E;
    }
    
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (i==j) {
                Hami[i][j] = E;
                if (randTest == 1) {
                    Hami[i][j] = LO + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI-LO)));
                }
            }
            else if (i-j == 1) {
                Hami[i][j] = t;
            }
            else if (i-j == -1) {
                Hami[i][j] = t;
            }
        }
    }
    
    
    
    // print out Hamiltonian
//    for (i=0; i<N; i++) {
//        for (j=0; j<N; j++) {
//            cout << Hami[i][j] ;
//        }
//        cout << endl;
//    }
    
    //Prepare Lapack dsterv call
    lapack_int N_l;
    double dummy = 0.0;
    lapack_int dummy_int = 1;
    
    //Output parameters
    lapack_int numEval,isuppz;
    double* evals = new double [N];
    double** dummy_matrix = new double* [1];
    dummy_matrix[0] = new double [1];
    
    lapack_int ldz = N;
    int info;
    
//    double reltol = 1E-7;
//    double minMag = ?;
    double abstol;
    
    //Get the machine epsilon
    abstol = LAPACKE_dlamch('E');
    
    double* diag;
    double* subDiag;
    
    diag = new double [N];
    subDiag = new double [N-1];
    
    //Assign diagonal elements
//    cout << "diag: ";
    for (i=0; i<N; i++) {
        diag[i] = Hami[i][i];
//        cout << diag[i];
    }
//    cout << endl;
    
    //Assign subdiagonal elements
//    cout << "subDiag: ";
    for (i=0; i<(N-1); i++) {
        subDiag[i] = Hami[i][i+1];
//        cout << subDiag[i];
    }
//    cout << endl;
    
    N_l = N;
    
    //I compared the below to the results from Numpy, and they agree to less than a relative error of 1E-14 for (N,E,t) = (10,1,2)
    info = LAPACKE_dstevr(LAPACK_ROW_MAJOR,'N','A',N_l,diag, subDiag,dummy,dummy,dummy_int,dummy_int,abstol,&numEval,evals,*dummy_matrix,ldz,&isuppz);
    
    if (info == 0) {
        cout << "LAPACKE_dstevr completed successfully." << endl;
    }
    else {
        cout << "!!!LAPACKE_dstevr ERROR!!!" << endl;
        cout << "dstevr info: " << info << endl;
        exit (EXIT_FAILURE);
    }
    
    
    cout << "NumEvals: " << numEval << endl;
    cout << "Eigenvalues: " << endl;
    cout << setprecision(10);
    for (i=0; i<numEval; i++) {
        cout << evals[i] << endl;
    }
    
    for (i=0; i<N; i++) {
        delete[] Hami[i];
    }
    delete[] Hami;
    delete[] diag;
    delete[] subDiag;
}