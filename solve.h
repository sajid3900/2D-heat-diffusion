#ifndef SOLVE_H
#define SOLVE_H

#include "diffusion.h"
#include <cmath>
#include "visualize.h"
#include <vector>
#include "input.h"
#include "dco.hpp"

using namespace dco;

// Initializing DCO variables
typedef double DCO_BASE_TYPE;
typedef ga1s<DCO_BASE_TYPE> DCO_MODE;
typedef DCO_MODE::type DCO_TYPE;
typedef DCO_MODE::tape_t DCO_TAPE_TYPE;



/** Solve for temperature
 * INPUT:
 * T_init   ---> Initial temperature values
 * Q_init   ---> Heat source values
 * C        ---> Heat diffusitivity values
 * my_input ---> The input class
 * OUTPUT:
 * T_new    --> The calculated temperature values
 **/
template<class active, class passive>
std::vector<active> solve(std::vector<passive> &T_init, std::vector<active> &C, std::vector<passive> &Q_init,input my_input){

    vector<active> T_new(my_input.nx * my_input.ny,0);
    vector<vector<active>> R(my_input.nx * my_input.ny,vector<active> (my_input.nx * my_input.ny,0));
    vector<vector<active>> I(my_input.nx * my_input.ny,vector<active> (my_input.nx * my_input.ny,0));
    identityMatrix(I,my_input.nx,my_input.ny);
    vector<active> T(my_input.nx * my_input.ny,0);
    vector<active> Q(my_input.nx * my_input.ny,0);
    vector<active> TQ(my_input.nx * my_input.ny,0);
    copy_vector(T_init,T,my_input);
    copy_vector(Q_init,Q,my_input);
    makeR(R, C, my_input);
    scalarMultiplication(R, my_input.delta_t,my_input.nx,my_input.ny);
    subtractMatrix(I,R,my_input.nx,my_input.ny);
    scalarMult_array(Q,my_input.delta_t,my_input.nx,my_input.ny);

    int count = 0;
    // Time loop
    while(count < my_input.m){

        for (int j = 0; j < my_input.nx*my_input.ny; j++)
            T_new[j] = 0;
        add_array(TQ,T,Q,my_input.nx,my_input.ny);


        // By Gauss Elimination
        //solve_ls(R, TQ, T_new,my_input.nx,my_input.ny);
        // Using Dense Eigen solver
        //solve_ls_eigen(R,TQ,T_new,my_input.nx,my_input.ny);
        // Using Sparse Eigen solver
        solve_ls_sparse(R,TQ,T_new,my_input.nx,my_input.ny);
        for (int j = 0; j < my_input.nx*my_input.ny; j++)
            T[j] = T_new[j];
        count++;
        }
        return T_new;

}


#endif // SOLVE_H
