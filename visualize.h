#ifndef VISUALIZE_H
#define VISUALIZE_H

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace std;


/** Writes vtk file
 * INPUT:
 * T_new    ---> Temperature values
 * C        ---> Heat diffusitivity values
 * nx,ny    ---> Mesh size parameters
 * count    ---> Time iteration number
 * filename ---> location to write the vtk file
 **/
template<class active>
void createVTK(std::vector<active> &T_new, std::vector<active> &C, int nx, int ny, int count, char *filename)
{

        // get mesh information from myinput

        //myinput = argInput;
        double Lx,Ly;
        Lx = 1.0;
        Ly = 1.0;

        // calculate node distances
        double dx = Lx/(nx-1);
        double dy = Ly/(ny-1);

        // Writing a vtk file named visualize.vtk
        ofstream vtkFile;
        //vtkFile.open("./vtk/visualize_" + std::to_string(count) + ".vtk", ios::out);
        vtkFile.open(filename + std::to_string(count) + ".vtk", ios::out);
        //cout << endl << "*************** VISUALIZATION ****************" << endl;
        if(vtkFile.is_open())
        {
                // vtk file for version 4.0 as an ASCII file
                vtkFile << "# vtk DataFile Version 4.0\n";
                vtkFile << "vtk output\n";
                vtkFile << "ASCII\n";

                //  A Rectilinear grid with dimensions
                vtkFile << "DATASET RECTILINEAR_GRID\n";
                vtkFile << "DIMENSIONS " << nx << " " << ny << " 1\n";

                // Write coordinates in x direction
                vtkFile << "X_COORDINATES " << 	nx << " float\n";
                vtkFile << "0 ";
                for(int i=1; i<nx-1; i++)
                {
                        vtkFile << dx*i << " ";
                }
                vtkFile << Lx << "\n";

                // Write coordinated in y direction
                vtkFile << "Y_COORDINATES " << 	ny << " float\n";
                vtkFile << "0 ";
                for(int i=1; i<ny-1; i++)
                {
                        vtkFile << dy*i << " ";
                }
                vtkFile << Ly << "\n";

                // Write z-coordinates (zero, we have a 2D mesh)
                vtkFile << "Z_COORDINATES 1 float\n";
                vtkFile << "0\n";

                /* Number of nodes and elements
                 *  	CELL_DATA: number of elements
                 *  	POINT_DATA: number of nodes
                 */
                vtkFile << "CELL_DATA " << (nx-1)*(ny-1) << "\n";
                vtkFile << "POINT_DATA " << nx*ny << "\n";

                vtkFile << "FIELD FieldData 2\n";

                // Field1: Temperature
                vtkFile << "Temperature 1 " << nx*ny << " float\n";
                vtkFile << T_new[0];
                for(int i=1; i<nx*ny; i++)
                {
                        vtkFile << " " << T_new[i];
                }
                vtkFile << "\n";

                // Field2: Heat conductivity
                vtkFile << "Heat_conductivity 1 " << nx*ny << " float\n";
                vtkFile << C[0];
                for(int i=1; i<nx*ny; i++)
                {
                        vtkFile << " " << C[i];
                }
                vtkFile << "\n";


                // close vtk file
                vtkFile.close();
                //cout << "Visualization file successfully saved in visualize.vtk file! Check it with Paraview!" << endl;
        }
        else
        {
                // Show error if not able to open the file
                std::cout << "Error: Unable to open vtk file for visualization" << endl;
        }

        return;
}

#endif // VISUALIZE_H
