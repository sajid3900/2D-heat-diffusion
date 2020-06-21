/*------------------------------------------------------------------------------------
SiSc Lab Project:   Parameter Estimation of 2D Heat Distribution in Heterogeneous Media

Authors:            Muhammad Sajid Ali
                    Shubhaditya Burela
                    Aneesh Futane
                    Pourya Pilva

--------------------------------------------------------------------------------------*/


#include <iostream>
#include <cmath>
#include <cstdio>
#include "diffusion.h"
#include "visualize.h"
#include "solve.h"
#include "optimizer.h"
#include <vector>
#include "input.h"
#include "dco.hpp"

int main()
{
    // Reading the input file. Change the filename(or destination) to the required case.
    char *filename;
    filename = "../SiSc_lab2d/case_1.json"; // Change to case 2 or 3 as required
    input my_input;
    my_input.read_input(filename); // Reading inputs
    my_input.print_input(); // Print inputs read
    int nx = my_input.nx;
    int ny = my_input.ny;

    // Setting the mesh coordinates
    std::vector<double> x(nx*ny,0), y(nx*ny,0);
    coordinates(x,y,my_input);

    // Initializing vectors
    std::vector<double> T(nx*ny,0); // Temperature
    std::vector<double> C(nx*ny,0); // Heat diffuitivity
    std::vector<double> Q(nx*ny,0); // Heat source
    std::vector<double> T_new(nx*ny,0); // Solved temperature distribution stored in this
    initialize_T(T,my_input);
    std::cout << "Initial Temperature:" << std::endl;
    printArray(T,my_input);
    initialize_c(C,x,y,my_input);
    std::cout << "\nDiffusivity, c: " << std::endl;
    printArray(C,my_input);
    initialize_q(Q,x,y,my_input);
    std::cout << "\nHeat source, q: " << std::endl;
    printArray(Q,my_input);

    // Solve
    T_new = solve(T,C,Q,my_input);

    // Print result
    std::cout << "Result: " << std::endl;
    printArray(T_new,my_input);

    // Writre vtk file
    char *vtk_filename;
    vtk_filename = "./vtk/visualize_final_";
    createVTK(T_new,C,my_input.nx,my_input.ny,0,vtk_filename);

    // -------------------------------------OPTIMIZATION-----------------------------------------------------------------------

    // Initiializing
    vector<double> T_obs(nx*ny,0); // Observed temperature
    vector<double> C_pert(nx*ny,0); // Perturbed c
   // Intitializing Observed values and perturbed c
    for(int i=0;i<nx*ny;i++){
        T_obs[i] = T_new[i];
        C_pert[i] = my_input.c_init; // Initial guess
    }

    // Optimizing by gradient descent
    optimizer(T,C_pert,Q,T_obs,my_input);

   // Write optimized c values
   std::vector<double> C_new(nx*ny,0);
   for(int i=0;i<nx*ny;i++)
       C_new[i]= C_pert[i];

   // Solve for T_new
   T_new = solve(T,C_new,Q,my_input);

   // Print results
   std::cout << "Result C_new: " << std::endl;
   printArray(C_new,my_input);
   std::cout << "Result T_new: " << std::endl;
   printArray(T_new,my_input);


   //Write vtk
   createVTK(T_new,C_new,my_input.nx,my_input.ny,1,vtk_filename);

}


