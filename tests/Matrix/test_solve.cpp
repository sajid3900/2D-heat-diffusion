#include <gtest/gtest.h>
#include "solve.h"
#include "diffusion.h"
#include "input.h"
#include "optimizer.h"
#include "dco.hpp"


TEST(DIFFUSION,initialize_T)
{
    vector<double> T(3*3);
    input my_input;
    my_input.nx = 3;
    my_input.ny = 3;
    my_input.T_left = 300;
    my_input.T_right = 330;
    my_input.T_init = 300;
    initialize_T(T,my_input);
    vector<double> TT{300,300,330,300,300,330,300,300,330};
    EXPECT_TRUE(T==TT);
}

TEST(DIFFUSION,initialize_c)
{
    vector<double> C(3*3);
    vector<double> x(3*3);
    vector<double> y(3*3);
    x = {0,0.5,1,0,0.5,1,0,0.5,1};
    y = {0,0,0,0.5,0.5,0.5,1,1,1};
    input my_input;
    my_input.nx = 3;
    my_input.ny = 3;
    my_input.nc = 2;
    my_input.c1 = {1.0,0,1,0,0.5};
    my_input.c2 = {2.0,0,1,0.5,1};
    initialize_c(C,x,y,my_input);
    vector<double> CC{1.0,1.0,1.0,2.0,2.0,2.0,2.0,2.0,2.0};
    EXPECT_TRUE(C==CC);
}

TEST(DIFFUSION,initialize_q)
{
    vector<double> Q(3*3);
    vector<double> x(3*3);
    vector<double> y(3*3);
    x = {0,0.5,1,0,0.5,1,0,0.5,1};
    y = {0,0,0,0.5,0.5,0.5,1,1,1};
    input my_input;
    my_input.nx = 3;
    my_input.ny = 3;
    my_input.q = {1.0,0.45,0.55,0.45,0.55};
    initialize_q(Q,x,y,my_input);
    vector<double> QQ{0,0,0,0,1,0,0,0,0};
    EXPECT_TRUE(Q==QQ);
}

TEST(DIFFUSION,make_R)
{
    input my_input;
    my_input.nx = 3;
    my_input.ny = 3;
    my_input.delta_x = 0.5,
    my_input.delta_y = 0.5;
    vector<double> C = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};;
    vector<vector<double>> R(3*3,vector<double> (3*3,0));
    makeR(R,C,my_input);
    vector<double> QQ{0,0,0,0,1,0,0,0,0};
    vector<vector<double>> RR{ {0,0,0,0,0,0,0,0,0},
                {0.0004,-0.0012,0.0004,0,0.0004,0,0,0,0},
                {0,0,0,0,0,0,0,0,0},
                {0,0,0,0,0,0,0,0,0},
                {0,0.0004,0,0.0004,-0.0016,0.0004,0,0.0004,0},
                {0,0,0,0,0,0,0,0,0},
                {0,0,0,0,0,0,0,0,0},
                {0,0,0,0,0.0004,0,0.0004,-0.0012,0.0004},
                {0,0,0,0,0,0,0,0,0} };
    for (int i= 0;i<9;i++)
        for (int j=0;j<9;j++)
            EXPECT_NEAR(R[i][j], RR[i][j], 1e-1);

}

TEST(SOLVE,solve)
{
    vector<double> Q = {0,0,0,0,1,0,0,0,0};
    vector<double> C = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001};
    vector<double> T = {300,300,330,300,300,330,300,300,330};
    input my_input;
    my_input.nx = 3;
    my_input.ny = 3;
    my_input.delta_x = 0.5,
    my_input.delta_y = 0.5;
    my_input.delta_t = 0.1;
    my_input.m = 400;
    vector<double> T_new(3*3);
    T_new = solve(T,C,Q,my_input);
    vector<double> T_expected = {300,300.781,330,300,339.22,330,300,300.781,330};
    for (int i= 0;i<9;i++)
        EXPECT_NEAR(T_new[i], T_expected[i], 1e-1);
}


TEST(OPTIMIZER,cost_function)
{
    vector<double> T = {300,300,330,300,300,330,300,300,330};
    vector<double> T_obs = {300,300.781,330,300,339.22,330,300,300.781,330};
    int nx = 3;
    int ny = 3;
    double J = cost_function(T_obs,T,nx,ny);
    EXPECT_NEAR(J,513.143, 1e-1);
}

/*
TEST(OPTIMIZER,optimizer)
{
    char *filename;
    filename = "../SiSc_lab2d/case_4.json"; // Change to case 2 or 3 as required
    input my_input;
    my_input.read_input(filename);
    int nx = my_input.nx, ny = my_input.ny;
    std::vector<double> x(nx*ny,0), y(nx*ny,0);
    coordinates(x,y,my_input);

    // Initializing vectors
    std::vector<double> T(nx*ny,0); // Temperature
    std::vector<double> C(nx*ny,0); // Heat diffuitivity
    std::vector<double> Q(nx*ny,0); // Heat source
    std::vector<double> T_new(nx*ny,0); // Solved temperature distribution stored in this
    initialize_T(T,my_input);
    initialize_c(C,x,y,my_input);
    initialize_q(Q,x,y,my_input);

    // Solve
    T_new = solve(T,C,Q,my_input);

    vector<double> T_obs(nx*ny,0); // Observed temperature
    vector<double> C_pert(nx*ny,0); // Perturbed c
   // Intitializing Observed values and perturbed c
    for(int i=0;i<nx*ny;i++){
        T_obs[i] = T_new[i];
        C_pert[i] = my_input.c_init; // Initial guess
    }

    // Optimizing by gradient descent
    optimizer(T,C_pert,Q,T_obs,my_input);

    std::vector<double> CC = {0.000201157,8.28209e-05,0.000132067,
                              0.000194345,9.58933e-05,0.000202355,
                              0.000201157,8.28209e-05,0.000132067};

    for(int i = 0; i < 9; i++)
        EXPECT_NEAR(C_pert[i],CC[i],1e-3);
}
*/

int main(int argc, char* argv[]) {


    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

}
