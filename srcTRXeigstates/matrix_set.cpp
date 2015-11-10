//
//  matrix_set.cpp
//  disorders
//
//  Created by Tianrui Xu on 2015-05-08.
//  Copyright (c) 2015 TianruiXu. All rights reserved.
//

#include "matrix_set.h"
#include <stdio.h>
#include <string>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;

double Hamiltonian::Vdd(Vector3d v1, Vector3d v2, int t_z, double gamma)
{
    double Vdd=0.0;
    const double J3_2h=1.0;//-52.1258738608;
    if (t_z==1)
    {
        Vector3d DC(0,0,1.0);
        double prod, r, costheta;
        Vector3d vd;
        vd=v2-v1;
        prod=vd.dot(DC);
        r=vd.norm()*DC.norm();
        costheta=pow(prod/r,2);
        Vdd=J3_2h/pow(vd.norm(),gamma)*(1-3*costheta);
    }
    else if (t_z==0)
    {
        Vector3d vd;
        vd=v2-v1;
/*        cout<<v1.transpose()<<endl;
        cout<<v2.transpose()<<endl;
        cout<<vd.norm()<<endl;*/
        Vdd=J3_2h/pow(vd.norm(),gamma);
    }
    else
    {cout<<"please check input file and resubmit the job!!! \n";};
    return Vdd;
}

void Hamiltonian::H_Basis(int lattice_sz, string occfile)
{
    const int sites=lattice_sz*lattice_sz*lattice_sz;
    double ind_x[sites];
    double ind_y[sites];
    double ind_z[sites];
    
    //basis
    int cs=0;
    for(int u=0; u<lattice_sz; u++)
    {
        for(int v=0; v<lattice_sz; v++)
        {
            for(int w=0; w<lattice_sz; w++)
            {
                ind_z[cs]=u;
                ind_y[cs]=v;
                ind_x[cs]=w;
                //cout<<ind_x[cs]<<'\t'<<ind_y[cs]<<'\t'<<ind_z[cs]<<endl;
                cs++;
            }
        }
    }
    
    Hamiltonian C;
    VectorXd occ;
    C.loadVectortxt(occfile,occ);
    const int nocc=int(occ.size());
    for (int i=0; i<nocc; i++)
    {
        int ind=occ[i];
        cout<<ind<<'\t'<<ind_x[ind]<<'\t'<<ind_y[ind]<<'\t'<<ind_z[ind]<<endl;
    }
}



MatrixXd Hamiltonian::Ham_offdiag(int lattice_sz, string vacfile, double hop_range, int zhop, double gamma)
{
    MatrixXd H_off;
    const int sites=lattice_sz*lattice_sz*lattice_sz;
    double ind_x[sites];
    double ind_y[sites];
    double ind_z[sites];
    
    //basis
    int cs=0;
    for(int u=0; u<lattice_sz; u++)
    {
        for(int v=0; v<lattice_sz; v++)
        {
            for(int w=0; w<lattice_sz; w++)
            {
                ind_z[cs]=u;
                ind_y[cs]=v;
                ind_x[cs]=w;
                //cout<<ind_x[cs]<<'\t'<<ind_y[cs]<<'\t'<<ind_z[cs]<<endl;
                cs++;
            }
        }
    }
    
    Hamiltonian C;
    VectorXd vac;
    C.loadVectortxt(vacfile,vac);
    const int nvac=int(vac.size());
    int fil[sites-nvac];
    
    for (int c=0; c<sites-nvac;c++)
    {
        fil[c]=0;
    }
    
    H_off.resize(sites-nvac,sites-nvac);
    
    cout<<"vac size "<<vac.size()<<endl;
    cout<<"fil size "<<sites-nvac<<endl;
    
    //Fill random
    cs=0;
    if (nvac==0)
    {
        for (int l=0; l<sites; l++)
        {
            fil[cs]=l;
            cs++;
        }
    }
    else
    {
        int vcn=0;//count number of vacancy
        for(int l=0; l<sites; l++)
        {
            //cout<<"chkpt0"<<endl;
            if(l!=vac[vcn])
            {
                fil[cs]=l;
                //cout<<"fil["<<cs<<"]= "<<fil[cs]<<'\n';
                cs++;
            }
            else
            {
                //cout<<"chkpt2"<<endl;
                //cout<<"vac["<<vcn<<"]= "<<vac[vcn]<<'\n';
                if (vcn==vac.size()-1)
                {
                    cout<<"exhaust all vacancy sites \n";
                }
                else
                {
                    vcn++;
                }
            }
        }
    }

    //fill H
    int sm, sn;
    for(int j=0; j<sites-nvac; j++)
    {
        sm=fil[j];
        for(int l=0; l<sites-nvac; l++)
        {
            sn=fil[l];
            Vector3d Hj(ind_x[sm],ind_y[sm],ind_z[sm]);
            Vector3d Hl(ind_x[sn],ind_y[sn],ind_z[sn]);
            Vector3d Hd=Hj-Hl;
            if (Hd.norm()<=hop_range)
            {
                H_off(j,l)=C.Vdd(Hj,Hl,zhop,gamma);//H(row,col)
            }
            else {H_off(j,l)=0.0;}
        }
        H_off(j,j)=0;
    }
    
    return H_off;
}

MatrixXd Hamiltonian::Ham_offdiag2(int lattice_sz, string occfile, double hop_range, int zhop, double gamma)
{
    MatrixXd H_off;
    const int sites=lattice_sz*lattice_sz*lattice_sz;
    double ind_x[sites];
    double ind_y[sites];
    double ind_z[sites];
    
    //basis
    int cs=0;
    for(int u=0; u<lattice_sz; u++)
    {
        for(int v=0; v<lattice_sz; v++)
        {
            for(int w=0; w<lattice_sz; w++)
            {
                ind_z[cs]=u;
                ind_y[cs]=v;
                ind_x[cs]=w;
                //cout<<ind_x[cs]<<'\t'<<ind_y[cs]<<'\t'<<ind_z[cs]<<endl;
                cs++;
            }
        }
    }
    
    Hamiltonian C;
    VectorXd occ;
    C.loadVectortxt(occfile,occ);
    const int nocc=int(occ.size());
    
    H_off.resize(nocc,nocc);
    
    cout<<"occ size "<<occ.size()<<endl;
    
    //fill H
    int sm, sn;
    for(int j=0; j<nocc; j++)
    {
        sm=occ[j];
        for(int l=0; l<nocc; l++)
        {
            sn=occ[l];
            Vector3d Hj(ind_x[sm],ind_y[sm],ind_z[sm]);
            Vector3d Hl(ind_x[sn],ind_y[sn],ind_z[sn]);
            Vector3d Hd=Hj-Hl;
            if (Hd.norm()<=hop_range)
            {
                H_off(j,l)=C.Vdd(Hj,Hl,zhop,gamma);//H(row,col)
            }
            else {H_off(j,l)=0.0;}
        }
        H_off(j,j)=0;
    }
    
    return H_off;
}

void Hamiltonian::loadVectortxt(string filename, VectorXd& m)
{
    int length;
    double tmp;
    ifstream vecfile (filename.c_str());
    if (vecfile.is_open())
    {
        vecfile>>length;
        m.resize(length);
        for(int j=0; j<length; j++)
        {
            vecfile>>tmp;
            m(j)=tmp;
        }
        vecfile.close();
    }
    else {cout << "Unable to open file "+filename+". \n";};
}



 void Hamiltonian::saveMatrixbin(string filename, const MatrixXd& m)
{
    ofstream f(filename.c_str(), ios::binary);
    // write the row_count and col_count into the file
    int rows = m.rows();
    int cols = m.cols();
    f.write((char *)&rows, sizeof(rows));
    f.write((char *)&cols, sizeof(cols));
    // write the matrix elements into the file
    f.write((char *)m.data(), sizeof(double)*rows*cols);
    f.close();
}

VectorXd Hamiltonian::H_Diagonal(string onsite)
{
    Hamiltonian Hdiag;
    VectorXd onsiteE;
    Hdiag.loadVectortxt(onsite,onsiteE);
    return onsiteE;
}

void Hamiltonian::Hami(int lattice_sz, string occfile, string onsite, double hop_range, int zhop, double gamma, double*& H, int &occ_num)
{
    //basis
    int N = lattice_sz*lattice_sz*lattice_sz;
    double ind_x[N];
    double ind_y[N];
    double ind_z[N];
    
    int cs=0;
    for(int u=0; u<lattice_sz; u++)
    {
        for(int v=0; v<lattice_sz; v++)
        {
            for(int w=0; w<lattice_sz; w++)
            {
                ind_z[cs]=u;
                ind_y[cs]=v;
                ind_x[cs]=w;
                //cout<<ind_x[cs]<<'\t'<<ind_y[cs]<<'\t'<<ind_z[cs]<<endl;
                cs++;
            }
        }
    }
    
    //fill H
    Hamiltonian C;
    VectorXd onsiteE,occ;
    C.loadVectortxt(onsite,onsiteE);
    C.loadVectortxt(occfile,occ);
    const int nocc=int(occ.size());
    occ_num=nocc;
    
    long long occVal = (long long) nocc;
    
    long long numElements = occVal*occVal;

    H = new double[numElements];
//    H = new double[nocc*nocc];
    for (int i=0; i<nocc; i++)
    {
        for (int j=0; j<nocc; j++)
        {
            H[i+nocc*j]=0.0;
        }
    }
    
    int sm, sn;
    for(int i=0; i<nocc; i++)
    {
        sm=occ[i];
        for(int j=0; j<nocc; j++)
        {
            sn=occ[j];
            Vector3d Hj(ind_x[sm],ind_y[sm],ind_z[sm]);
            Vector3d Hl(ind_x[sn],ind_y[sn],ind_z[sn]);
            Vector3d Hd=Hj-Hl;
            if (Hd.norm()<=hop_range)
            {
                H[i+nocc*j]=C.Vdd(Hj,Hl,zhop,gamma);
            }
            else {H[i+nocc*j]=0.0;}
        }
        H[i+nocc*i]=onsiteE[i];
    }
}

