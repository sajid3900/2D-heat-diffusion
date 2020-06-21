#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <cmath>
#include <cstdio>
#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"
#include <vector>
#include <fstream>
#include <string>


using namespace rapidjson;

/** Input class defined to read json file**/
class input{

public:

    int nx, ny, m, nc, opt_steps;
    double Lx, Ly, tf, T_init, T_left, T_right, delta_x, delta_y, delta_t, c_init, alpha, opt_err;
    std::string case_name;
    std::vector<double> c1,c2,c3,q;


    void read_input(char* filename){

        FILE *fp;
        fp = fopen(filename, "rb"); // non-Windows use "r"

        char readBuffer[65536];
        FileReadStream is(fp, readBuffer, sizeof(readBuffer));

        Document d;
        d.ParseStream(is);

        fclose(fp);

        case_name = d["name"].GetString();
        Lx = d["Lx"].GetDouble();
        Ly = d["Ly"].GetDouble();
        nx = d["nx"].GetInt();
        ny = d["ny"].GetInt();
        tf = d["tf"].GetDouble();
        m = d["m"].GetInt();
        nc = d["nc"].GetInt();
        T_init = d["T_init"].GetDouble();
        T_left = d["T_left"].GetDouble();
        T_right = d["T_right"].GetDouble();
        c_init = d["c_init"].GetDouble();
        alpha = d["alpha"].GetDouble();
        opt_steps = d["opt_steps"].GetInt();
        opt_err = d["opt_err"].GetDouble();

        delta_x = Lx / (nx-1);
        delta_y = Ly / (ny-1);
        delta_t = tf / m;



        // Reading heat source and heat diffusitivity values
        const Value& z = d["q"]; // heat source
        assert(z.IsArray());
        const Value& a = d["c1"];  // first c
        assert(a.IsArray());
        for(SizeType i = 0; i < a.Size(); i++){
            c1.push_back(a[i].GetDouble()); // c,x_start,x_end,y_start,y_end
            q.push_back(z[i].GetDouble());
        }

        if(nc>1){
            const Value& b = d["c2"]; // second c
            assert(b.IsArray());
            for(SizeType i = 0; i < a.Size(); i++)
                    c2.push_back(b[i].GetDouble());
        }
        if(nc>2){
            const Value& e = d["c3"]; // third c
            assert(e.IsArray());
            for(SizeType i = 0; i < a.Size(); i++)
                    c3.push_back(e[i].GetDouble());
        }


    }

    void print_input(){

        std::cout << "Case: " << case_name << std::endl;
        std::cout << "Lx = "<< Lx << std::endl;
        std::cout << "Ly = "<< Ly << std::endl;
        std::cout << "nx = "<< nx << std::endl;
        std::cout << "ny = "<< ny << std::endl;
        std::cout << "Final time = "<< tf << std::endl;
        std::cout << "Number of time steps = "<< m << std::endl;
        std::cout << "T_init = "<< T_init << std::endl;
        std::cout << "T_left = "<< T_left << std::endl;
        std::cout << "T_right = "<< T_right << std::endl;
        std::cout << "delta_x = " << delta_x << std::endl;
        std::cout << "delta_y = " << delta_y << std::endl;
        std::cout << "delta_t = " << delta_t << std::endl;
    }


};




#endif // INPUT_H
