#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <iostream>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <vector>
#include "input.h"
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include "dco.hpp"

using namespace Eigen;

/**
 * This header file contains all the needed
 * helper functions for the 2d heat diffusion
 * problem to be solved
 **/


template<class active>
void initialize_T(std::vector<active> &T, input my_input) {
        for (int i = 0; i < my_input.nx * my_input.ny; i++) {
                if (i % my_input.nx == (my_input.nx - 1))
                        T[i] = my_input.T_right;
                else if (i % my_input.nx == 0)
                        T[i] = my_input.T_left;
                else
                        T[i] = my_input.T_init;
        }
}

template<class active>
void printArray(std::vector<active> &A, input my_input) {
        for (int i = 0; i < my_input.nx*my_input.ny; i++){
                std::cout << A[i] << "\t";
                if(i%my_input.nx == (my_input.nx-1))
                    std::cout << std::endl;
        }

}

template<class active>
void printMatrix(std::vector<std::vector<active>> &A, int nx, int ny) {
        for (int i = 0; i < nx*ny; i++) {
                for (int j = 0; j < nx*ny; j++)
                        std::cout << A[i][j] << "\t";
                std::cout << "\n";
        }
        std::cout << "\n";
}


template<class active>
void identityMatrix(std::vector<std::vector<active>> &A, int nx, int ny) {
        for (int i = 0; i < nx * ny; i++) {
                for (int j = 0; j < nx * ny; j++) {
                        if (i == j)
                                A[i][i] = 1;
                        else
                                A[i][j] = 0;
                }
        }
}


template<class active>
void makeR(std::vector<std::vector<active>> &R, std::vector<active> &c, input my_input) {
        for (int i = 0; i < my_input.nx * my_input.ny; i++) {
                if (i % my_input.nx == 0)
                        R[i][i] = 0;
                else if (i % my_input.nx == (my_input.nx-1))
                        R[i][i] = 0;
                else {
                        if ((i - my_input.nx) >= 0)
                                R[i][i - my_input.nx] = (c[i] + c[i - my_input.nx]) / 2 / pow(my_input.delta_y, 2);

                        R[i][i - 1] = (c[i] + c[i - 1]) / 2 / pow(my_input.delta_x, 2);

                        if ( i-my_input.nx >=0 && i+my_input.nx <= my_input.nx*my_input.ny)
                                R[i][i] = -(c[i - 1] + 2 * c[i] + c[i + 1]) / 2 / pow(my_input.delta_x, 2) -
                                        (c[i + my_input.nx] + 2 * c[i] + c[i - my_input.nx]) / 2 / pow(my_input.delta_y, 2);
                        else if (i-my_input.nx >= 0)
                                R[i][i] = -(c[i - 1] + 2 * c[i] + c[i + 1]) / 2 / pow(my_input.delta_x, 2) -
                                        (c[i] + c[i - my_input.nx]) / 2 / pow(my_input.delta_y, 2);
                        else
                                R[i][i] = -(c[i - 1] + 2 * c[i] + c[i + 1]) / 2 / pow(my_input.delta_x, 2) -
                                        (c[i + my_input.nx] + c[i]) / 2 / pow(my_input.delta_y, 2);

                        R[i][i + 1] = (c[i] + c[i + 1]) / 2 / pow(my_input.delta_x, 2);

                        if ( (i + my_input.nx) <= (my_input.nx * my_input.ny))
                                R[i][i + my_input.nx] = (c[i] + c[i + my_input.nx]) / 2 / pow(my_input.delta_y, 2);
                }
        }

}


template<class active>
void scalarMultiplication(std::vector<std::vector<active>> &A, double k, int nx, int ny) {
        for (int i = 0; i < nx*ny; i++)
                for (int j = 0; j < nx*ny; j++)
                        A[i][j] = k * A[i][j];
}


template<class active>
void subtractMatrix(std::vector<std::vector<active>> &A, std::vector<std::vector<active>> &B, int nx, int ny) {
        // A - B and store in B
        for (int i = 0; i < nx*ny; i++)
                for (int j = 0; j < nx*ny; j++)
                        B[i][j] = A[i][j] - B[i][j];
}


// Solves Linear System by Gauss-elimination
template<class active>
void solve_ls(std::vector<std::vector<active>> &A, std::vector<active> &b, std::vector<active> &x, int nx, int ny) {
        // Ax = b
        //double C[13*13][13*13 + 1] = { 0 }; // here C[n*n][n*n+1]
        std::vector<std::vector<active>> C(13 * 13,std::vector<active> ((13 * 13)+1,0));
        for (int i = 0; i < nx*ny; i++)
                for (int j = 0; j < nx*ny + 1; j++) {
                        if (j < nx*ny)
                                C[i][j] = A[i][j];
                        else
                                C[i][j] = b[i];
                }


        active e;
        for (int j = 0; j < nx*ny; j++) // loop for the generation of upper triangular matrix
        {
                for (int i = 0; i < nx*ny; i++)
                {
                        if (i > j)
                        {
                                e = C[i][j] / C[j][j];
                                for (int k = 0; k < nx*ny + 1; k++)
                                {
                                        C[i][k] = C[i][k] - e * C[j][k];
                                }
                        }
                }
        }

        active sum;
        x[nx*ny - 1] = C[nx*ny - 1][nx*ny] / C[nx*ny - 1][nx*ny - 1];
        // this loop is for backward substitution
        for (int i = nx*ny - 2; i >= 0; i--)
        {
                sum = 0;
                for (int j = i + 1; j < nx*ny; j++)
                {
                        sum = sum + C[i][j] * x[j];
                }
                x[i] = (C[i][nx*ny] - sum) / C[i][i];
        }

}

template<class active>
void coordinates(std::vector<active> &x, std::vector<active> &y, input my_input){
    for(int i = 0; i < my_input.nx*my_input.ny; i++){
         x[i] = (i%my_input.nx)/(my_input.nx-1.0);
    }
    for(int i=0;i<my_input.nx*my_input.ny;i++){
        y[i] = floor(i/my_input.ny)/(my_input.ny-1);
    }
}

template<class active>
void initialize_c(std::vector<active> &C, std::vector<active> &x, std::vector<active> &y, input my_input){

    // function to initialize c for c1 and c2 using x,y coordinates
    // c1,c2 - c, x_start, x_end, y_start, y_end
    if (my_input.nc == 1){
        for (int i = 0; i< my_input.nx*my_input.ny; i++)
            C[i] = my_input.c1[0];
    }
    if (my_input.nc == 2){
        for(int i=0; i<my_input.nx*my_input.ny; i++){
            if( (x[i]>=my_input.c1[1] && x[i]<=my_input.c1[2]) && (y[i]>=my_input.c1[3] && y[i]<=my_input.c1[4])){
                C[i] = my_input.c1[0];
            }
            if( (x[i]>=my_input.c2[1] && x[i]<=my_input.c2[2]) && (y[i]>=my_input.c2[3] && y[i]<=my_input.c2[4])){
                C[i] = my_input.c2[0];
            }
        }
    }
    if (my_input.nc == 3){
        for(int i=0; i<my_input.nx*my_input.ny; i++){
            if( (x[i]>=my_input.c1[1] && x[i]<=my_input.c1[2]) && (y[i]>=my_input.c1[3] && y[i]<=my_input.c1[4])){
                C[i] = my_input.c1[0];
            }
            if( (x[i]>=my_input.c2[1] && x[i]<=my_input.c2[2]) && (y[i]>=my_input.c2[3] && y[i]<=my_input.c2[4])){
                C[i] = my_input.c2[0];
            }
            if( (x[i]>=my_input.c3[1] && x[i]<=my_input.c3[2]) && (y[i]>=my_input.c3[3] && y[i]<=my_input.c3[4])){
                C[i] = my_input.c3[0];
            }
        }
    }
}


template<class active>
void initialize_q(std::vector<active> &Q, std::vector<active> &x, std::vector<active> &y, input my_input){

    for(int i=0; i<my_input.nx * my_input.ny; i++){
        if( (x[i]>= my_input.q[1] && x[i]<=my_input.q[2]) && (y[i]>=my_input.q[3] && y[i]<=my_input.q[4])){
            Q[i] = my_input.q[0];
        }
        else
            Q[i] = 0;
    }

    // boundary conditions
    for (int i = 0; i < my_input.nx * my_input.ny; i++)
            if ( (i % my_input.nx == (my_input.nx - 1)) || (i % my_input.nx == 0) )
                    Q[i] = 0;


}

template<class active>
void scalarMult_array(std::vector<active> &A, double k, int nx, int ny) {
        for (int i = 0; i < nx*ny; i++)
                        A[i] = k * A[i];
}

template<class active>
void add_array(std::vector<active> &C, std::vector<active> &A, std::vector<active> &B, int nx, int ny){
    // C = A + B
    for (int i = 0; i < nx*ny; i++)
                    C[i] = A[i] + B[i];
}

template<class active, class passive>
void copy_vector(std::vector<passive> &T_init, std::vector<active> &T, input my_input){
    for (int i=0; i<my_input.nx*my_input.ny;i++)
        T[i] = T_init[i];
}


// Solves linear system by Eigen Dense solver
template<class active>
void solve_ls_eigen(std::vector<std::vector<active>> &A, std::vector<active> &b, std::vector<active> &x, int nx, int ny) {

    typedef Matrix<active,Dynamic,Dynamic> matrix_dco;
    typedef Matrix<active,Dynamic,1> vector_dco;

    matrix_dco C(nx*ny,nx*ny);
    vector_dco d(nx*ny);
    vector_dco y(nx*ny);

    for(int i = 0; i < nx * ny; i++){
        d(i) = b[i];
        for(int j = 0; j < nx * ny; j++){
            C(i,j) = A[i][j];
        }
    }

    y = C.partialPivLu().solve(d);
    for(int i = 0; i < nx * ny; i++)
        x[i] = y(i);
}



// Solves linear system by Eigen Sparse solver
template<class active>
void solve_ls_sparse(std::vector<std::vector<active>> &A, std::vector<active> &b, std::vector<active> &x, int nx, int ny) {

    typedef Eigen::SparseMatrix<active> SpMat;
    typedef Matrix<active,Dynamic,1> vector_dco;

    SpMat C(nx*ny,nx*ny);
    vector_dco d(nx*ny);

    C.reserve(VectorXi::Constant(nx*ny,10));

    for(int i = 0; i < nx * ny; i++){
        d(i) = b[i];
    }

    for(int i = 0; i < nx*ny; i++){
        for(int j = 0; j < nx*ny; j++){
            if(A[i][j] != 0){
                C.insert(i,j) = A[i][j];
            }
        }
    }

    C.makeCompressed();


    Eigen::SparseLU<SpMat> LUD(C);
    vector_dco y = LUD.solve(d);

    for(int i = 0; i < nx * ny; i++)
        x[i] = y(i);
}


#endif // DIFFUSION_H

