//
//  main.cpp
//  disorders
//
//  Created by Tianrui Xu on 2015-05-03.
//  Copyright (c) 2015 TianruiXu. All rights reserved.
//


#include <iostream>
#include <vector>
#include <stdlib.h>
#include <Eigen/Dense>
#include <time.h>
#include <Eigen/Eigenvalues>
#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "matrix_set.h"
#include "quantumPercolation.h"

#define I complex<double>(0.0,1.0);
const double pi=3.141592653589793238462643383279502884197169399375105820974944;
const double eps=8.85419*pow(10.0,-12);
const double h=6.62606957*pow(10.0,-34);
const double hbar=1.054571726*pow(10.0,-34);

using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;

int zhop, size, occnum;
double Gamma;
double range,pvac;
string occfile="rd3do",rndfile="rd3dw",outfile="rd3d";

int main(int argc, const char * argv[])
{
    /** %%%%%%%%%%%%%%%%%%%% INITIALIZING %%%%%%%%%%%%%%%%%%%% */
    const char *dis=argv[argc-1];
    //string c="dis_in.txt";
    const char* inputFile = argv[1];
    string line;
    int dum=0;
    VectorXd a;
    ifstream input (inputFile);
    if (input.is_open())
    {
        while (getline(input,line))
        {
            if (line[0]=='!'){dum=dum+1;}
            else
            {
                if (dum==1)
                {
                    size=atoi(line.c_str());
                    occfile=occfile.append(line.c_str());
                    rndfile=rndfile.append(line.c_str());
                    outfile=outfile.append(line.c_str());
                    outfile=outfile.append("_");                    
                }
                else if (dum==2)
                {
                    if (line.compare("Long")==0)
			{
			    range=size*2.0;
			    outfile=outfile.append("LR");
			    outfile=outfile.append("_"); 
			}
                    else 
			{
			    range=atoi(line.c_str());
			    outfile=outfile.append(line.c_str());
			    outfile=outfile.append("NN_"); 
			}
                }
                else if (dum==3){Gamma=atof(line.c_str());}
                else if (dum==4){zhop=atoi(line.c_str());}
                else if (dum==5)
                {
                    pvac=atoi(line.c_str());
                    occfile=occfile.append(line.c_str());
                    rndfile=rndfile.append(line.c_str());
                    outfile=outfile.append(line.c_str());
                    outfile=outfile.append("_");
                }

/*                else if (dum==6)
                {
                    cout<<'\t'<<dum<<endl;
                    occfile=occfile.append(line.c_str());
                    rndfile=rndfile.append(line.c_str());
                }
*/
            }
        }
        input.close();
    }
    else {cout<<"Unable to open file: "<<inputFile<< " \n";};
    
    occfile=occfile.append(dis);
    rndfile=rndfile.append(dis);
    outfile=outfile.append(dis);
    
    occfile=occfile.append(".txt");
    rndfile=rndfile.append(".txt");
    
    cout<<"lattice size = "<<size<<endl;
    cout<<"range = "<<range<<endl;
    cout<<"decay = "<<Gamma<<endl;
    cout<<"symmetrical? "<<zhop<<endl;
    cout<<"percent of vacancy: "<<pvac<<endl;
    cout<<"occupation file: "<<occfile<<endl;
    cout<<"on-site energy file: "<<rndfile<<endl;
    /** %%%%%%%%%%%%%%%%%%%% END OF INITIALIZATION %%%%%%%%%%%%%%%%%%%% */

    Hamiltonian C;
    double* H;
    C.Hami(size, occfile, rndfile, range, zhop, Gamma, H, occnum);

//Lapack
    int numEigvals;
    double* eigvals;
    int presInt = 16;
    eigvals = new double [occnum];
            
    dsyevrEigvalsInterface(H, occnum, &numEigvals, eigvals);
       
    cout << "NumEvals: " << numEigvals << endl;
    cout << "Eigenvalues: " << endl;
    cout << setprecision(presInt);
    ofstream myfile2;
    myfile2.open((outfile+"_Eval.txt").c_str());
    myfile2 << setprecision(presInt);
    for (int i=0; i<numEigvals; i++) {
        cout << eigvals[i] << endl;
	myfile2<< eigvals[i] <<endl;
    } 
    myfile2.close();    
    delete [] eigvals;
    delete [] H;

    
}
