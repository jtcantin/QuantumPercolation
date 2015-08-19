//
//  matrix_set.h
//  disorders
//
//  Created by Tianrui Xu on 2015-05-08.
//  Copyright (c) 2015 TianruiXu. All rights reserved.
//

#ifndef __disorders__matrix_set__
#define __disorders__matrix_set__

#include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;

class Hamiltonian {
    public:
        /** set size and lattice index of the Hamiltonian */
        //void set_Hamiltonian(int size);
        //int lsz() const;
    
        /** returns basis of the Hamiltonian */
        void H_Basis(int lattice_sz, string occfile);
    
        /** returns diagonal part of the Hamiltonian */
        VectorXd H_Diagonal(string onsite);
    
        /** returns the interaction of the system in between of molecules on two sites */
        double Vdd(Vector3d v1, Vector3d v2, int t_z, double gamma);
    
        /** returns off-diagonal part of the Hamiltonian w random vacancy */
        MatrixXd Ham_offdiag(int lattice_sz, string vacfile, double hop_range, int zhop, double gamma);
    
        /** returns off-diagonal part of the Hamiltonian w random OCCUPATION */
        MatrixXd Ham_offdiag2(int lattice_sz, string occfile, double hop_range, int zhop, double gamma);
    
        /** returns the Hamiltonian in double vector form*/
        void Hami(int lattice_sz, string occfile, string onsite, double hop_range, int zhop, double gamma, double*& H, int &occ_num);

        /** load a .bin matrix */
//        void loadMatrixbin(string filename, MatrixXd& m);
    
        /** save a .bin matrix */
        void saveMatrixbin(string filename, const MatrixXd& m);
    
        /** load a .txt matrix */
//        void loadMatrixtxt(string filename, MatrixXd& m);
    
        /** load a .txt vector */
        void loadVectortxt(string filename, VectorXd& m);
    
};


#endif /* defined(__disorders__matrix_set__) */
