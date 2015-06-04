//
//  quantumPercolation.cpp
//  
//
//  Created by Joshua Tyler Cantin on 2015-04-24.
//
// To run:
// ./exec N E t randE_i[0(no) or 1(yes)] solver[0(dstevr), 1(dsyevr)]

#include "quantumPercolation.h"

using namespace std;

int main(int argc, char **argv) {
    cout << "Hello World." << endl;
    
    int presInt = 32;
    
    //Make Hamiltonian
    int N;
    double E;
    double t;
    int randTest;
    int solver;
    
    int alpha = 3;
    int longRange = 0;
    
    
    if (argc != 6) {
        cout << "Not enough arguments!" << endl;
        exit (EXIT_FAILURE);
    }
    
    N = atoi(argv[1]);
    E = atof(argv[2]);
    t = atof(argv[3]);
    randTest = atoi(argv[4]);
    solver = atoi(argv[5]);
    
    double** Hami;
    double* Hami2;
    
    Hami = new double* [N];
    Hami2 = new double [N*N];
    
    int i,j;
    
    for (i=0; i<N; i++) {
        Hami[i] = new double [N];
        for (j=0; j<N; j++) {
            Hami[i][j] = 0.0;
//            Hami2[N*i+j] = 0.0; //Row major order
            Hami2[i+N*j] = 0.0; //Column major order
        }
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
    
    if (longRange == 0)
    {
    
        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                if (i==j) {
                    Hami2[i+N*j] = E;
                    if (randTest == 1) {
                        Hami2[i+N*j] = LO + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI-LO)));
                    }
                }
                else if (i-j == 1) {
                    Hami2[i+N*j] = t;
                }
                else if (i-j == -1) {
                    Hami2[i+N*j] = t;
                }
            }
        }
    }
    else if (longRange==1)
    {
        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                if (i==j) {
                    Hami2[i+N*j] = E;
                    if (randTest == 1) {
                        Hami2[i+N*j] = LO + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(HI-LO)));
                    }
                }
                else {
                    Hami2[i+N*j] = t / pow( fabs(i - j), alpha);
                }
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
//    
//    cout << setprecision(3) << fixed;
//    for (i=0; i<N; i++) {
//        for (j=0; j<N; j++) {
//            cout << Hami2[i+N*j] << " ";
//        }
//        cout << endl;
//    }
    
    switch (solver) {
        case 0:
        {
            
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
           
#ifdef OSX_LAPACKE 
            //Get the machine epsilon
            abstol = LAPACKE_dlamch('E');
#endif

#ifdef LINUX_MKL
	    //Arbitrary value
	    abstol = 1.0E-10;
#endif            

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
            cout << setprecision(presInt);
            for (i=0; i<numEval; i++) {
                cout << evals[i] << endl;
            }
            
            delete[] diag;
            delete[] subDiag;
            delete [] evals;
        }
            break;
            
        case 1:
        {
            int numEigvals;
            double* eigvals;
            eigvals = new double [N];
            
            dsyevrEigvalsInterface(Hami2, N, &numEigvals, eigvals);
            
            
            cout << "NumEvals: " << numEigvals << endl;
            cout << "Eigenvalues: " << endl;
            cout << setprecision(presInt);
            int i;
            for (i=0; i<numEigvals; i++) {
                cout << eigvals[i] << endl;
            }
            
            delete [] eigvals;
            
//            // print out Hamiltonian
////            for (i=0; i<N; i++) {
////                for (j=0; j<N; j++) {
////                    cout << Hami[i][j] ;
////                }
////                cout << endl;
////            }
//            
//            
//            //Prepare Lapack dsyevr call
//            lapack_int N_l;
////            double dummy = 0.0;
//            lapack_int dummy_int = 1;
//            
//            //Output parameters
//            lapack_int numEval,isuppz;
//            double* evals = new double [N];
//            double** dummy_matrix = new double* [1];
//            dummy_matrix[0] = new double [1];
//            
//            lapack_int ldz = N;
//            lapack_int lda = N;
//            int info;
//            
//            //    double reltol = 1E-7;
//            //    double minMag = ?;
//            double abstol;
//
//#ifdef OSX_LAPACKE 
//            //Get the machine epsilon
//            abstol = LAPACKE_dlamch('E');
//#endif
//
//#ifdef LINUX_MKL
//	    //Arbitrary value
//	    abstol = 1.0E-10;
//#endif            
//            
//            N_l = N;
//            
//            //I compared the below to the results from dstevr, and they are identical
//            info = LAPACKE_dsyevr(LAPACK_COL_MAJOR,'N','A','U',N_l,Hami2,lda,dummy_int,dummy_int,dummy_int,dummy_int,abstol, &numEval,evals,*dummy_matrix,ldz,&isuppz);
//            
//            if (info == 0) {
//                cout << "LAPACKE_dsyevr completed successfully." << endl;
//            }
//            else {
//                cout << "!!!LAPACKE_dsyevr ERROR!!!" << endl;
//                cout << "dstevr info: " << info << endl;
//                exit (EXIT_FAILURE);
//            }
//            
//            cout << "NumEvals: " << numEval << endl;
//            cout << "Eigenvalues: " << endl;
//            cout << setprecision(presInt);
//            for (i=0; i<numEval; i++) {
//                cout << evals[i] << endl;
//            }
//            delete [] evals;
        }

            
            break;
        default:
            cout << "Solver not defined! No work done." << endl;
            break;
        }
    
    
    
    
    for (i=0; i<N; i++) {
        delete[] Hami[i];
    }
    delete[] Hami;
    delete[] Hami2;
    
}
