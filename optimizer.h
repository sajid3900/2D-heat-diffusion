#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <cmath>
#include "solve.h"
#include "dco.hpp"
#include <vector>
#include "diffusion.h"
#include "input.h"
#include "visualize.h"

using namespace dco;

/** The cost function
 * INPUT:
 * T_obs    ---> Observed temperature values
 * T_comp   ---> Computed Temperature values
 * nx, ny   ---> Mesh size parameters
 * OUTPUT:
 * J        ---> The cost function value (scalar)
 **/
template<class active>
active cost_function(std::vector<active> &T_obs, std::vector<active> &T_comp, int nx, int ny){

    active n = (nx-2)*ny; // to remove boundary values
    active sum = 0;
    active J;

    for(int i=0; i<nx*ny;i++){
        sum = sum + (((T_obs[i]))-(T_comp[i]))*((T_obs[i])-(T_comp[i]));
    }

    J = sum / n;
    return J;
}




/** The computation of gradients by adjoint
 * INPUT:
 * T        ---> Initial temperature values
 * Q        ---> Heat source values
 * C        ---> Heat diffusitivity values
 * T_obs    ---> Observed temperature values
 * dJdc     ---> To store the gradients
 * my_input ---> The input class
 * count    ---> required for vtk plots
 * OUTPUT:
 * J        ---> The cost function value (scalar)
 **/
template<class active>
DCO_TYPE compute_grad(std::vector<active> &T, std::vector<active> &Q, std::vector<active> &C,
                  std::vector<active> &T_obs, std::vector<active> &dJdc, input my_input, int count){


    int nx= my_input.nx, ny = my_input.ny;
    typedef active DCO_BASE_TYPE;
    typedef ga1s<DCO_BASE_TYPE> DCO_MODE;
    typedef DCO_MODE::type DCO_TYPE;
    typedef DCO_MODE::tape_t DCO_TAPE_TYPE;
    DCO_MODE::global_tape=DCO_TAPE_TYPE::create();
    std::vector<DCO_TYPE> CC(nx*ny), TT_obs(nx*ny);
    typename dco::ga1s<active>::type J_adjoint;
    std::vector<DCO_TYPE> Tnew(my_input.nx * my_input.ny);

    // Initialization of the variables
    for(int i = 0; i < my_input.nx * my_input.ny; i++){
        CC[i] = C[i];
        TT_obs[i] = T_obs[i];
        dJdc[i] = 0.0;
        DCO_MODE::global_tape->register_variable(CC[i]);
    }


    // Solve for Tnew
    Tnew = solve(T,CC,Q,my_input);

    // Write vtk file for each iteration
    char *vtk_filename;
    vtk_filename = "./vtk/visualize_";
    createVTK(Tnew,CC,my_input.nx,my_input.ny,count,vtk_filename);

    // Seeding
    DCO_MODE::global_tape->register_output_variable(J_adjoint);
    J_adjoint = cost_function(TT_obs, Tnew, my_input.nx, my_input.ny);
    derivative(J_adjoint) = 1.0;

    // Harvest the derivatives
    DCO_MODE::global_tape->interpret_adjoint();

    // Store the derivatives
    for(int i = 0; i < my_input.nx * my_input.ny; i++)
        dJdc[i] = derivative(CC[i]);
    DCO_TAPE_TYPE::remove(DCO_MODE::global_tape);

    return J_adjoint;

}



/** The computation of gradients
 * INPUT:
 * T        ---> Initial temperature values
 * Q        ---> Heat source values
 * C        ---> Heat diffusitivity values
 * T_obs    ---> Observed temperature values
 * my_input ---> The input class
 **/
template<class active>
void optimizer(std::vector<active> &T,std::vector<active> &C, std::vector<active> &Q, std::vector<active> &T_obs, input my_input){


    int nx= my_input.nx, ny = my_input.ny;
    std::vector<active> dJdc(nx*ny);
    int count = 0;
    double grad_norm = 10;
    DCO_TYPE J;
    while(count < my_input.opt_steps && grad_norm > my_input.opt_err){
        J = compute_grad(T,Q,C,T_obs,dJdc,my_input,count); // Adjoint mode
        // compute_grad_tangent(T,Q,C,T_obs,dJdc,my_input,count); // Tangent mode
        grad_norm = 0;
        for(int i = 0; i < nx*ny; i++){
            C[i] = C[i] - my_input.alpha* dJdc[i];
            grad_norm += dJdc[i]*dJdc[i];

        }
        grad_norm = sqrt(grad_norm/(nx*ny));
        count++;
        if (count % 10 == 0)
           std::cout << "Iteration = " << count << ", J = " << value(J) << ", Grad_norm = " << grad_norm << std::endl;
         //std::cout << "Iteration = " << count << ", Grad_norm = " << grad_norm << std::endl; // use for tangent mode
   }

}

/** The computation of gradients by Tangent mode
 * INPUT:
 * T        ---> Initial temperature values
 * Q        ---> Heat source values
 * C        ---> Heat diffusitivity values
 * T_obs    ---> Observed temperature values
 * dJdc     ---> To store the gradients
 * my_input ---> The input class
 * count    ---> required for vtk plots
 **/
template<class active>
void compute_grad_tangent(std::vector<active> &T, std::vector<active> &Q, std::vector<active> &C,
                  std::vector<active> &T_obs, std::vector<active> &dJdc, input my_input, int count){


    int nx= my_input.nx, ny = my_input.ny;
    typedef active DCO_BASE_TYPE;
    typedef gt1s<DCO_BASE_TYPE> DCO_MODE;
    typedef DCO_MODE::type DCO_TYPE;
    std::vector<DCO_TYPE> CC(nx*ny), TT_obs(nx*ny);
    typename dco::gt1s<active>::type J;
    std::vector<DCO_TYPE> Tnew(my_input.nx * my_input.ny);

    // Initialization of the variables
    for(int i = 0; i < my_input.nx * my_input.ny; i++){
        CC[i] = C[i];
        TT_obs[i] = T_obs[i];
        dJdc[i] = 0.0;
    }

    for(int i = 0; i<my_input.nx*my_input.ny;i++){
        derivative(CC[i]) = 1.0;
        Tnew = solve(T,CC,Q,my_input);
        J = cost_function(TT_obs, Tnew, my_input.nx, my_input.ny);
        dJdc[i] = derivative(J);
        derivative(CC[i]) = 0.0;
    }

}



#endif // OPTIMIZER_H
