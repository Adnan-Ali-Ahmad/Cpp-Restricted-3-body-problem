#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "functions.h"
#include <object.h>
#include <ctime>
#include <distances.h>

using namespace std;

double pi = acos(-1.);
double G = 6.67408e-11;

int main(){

    srand(time(NULL)); //for random generation purposes

    object sun(1.989e30);
    object jup(1.89813e27,5.20336301,0.04839266,1.30530,14.75385,100.55615,34.40438); //Jupiter
    object COM(0.); //center of mass

    object M1 = sun; //copies
    object M2 = jup;

    double GM = G*(sun.m+jup.m);
    double T = 25.*orbital_period(GM,jup.a); //simulation time
    int it = 0;
    double dt; //time step
    int snapshot = 50;
    double t = 0.; //secs
    double t0 = 0.; //J2000 epoch

    get_cartesian(GM,t,t0,sun);
    get_cartesian(GM,t,t0,jup);
    update_COM(COM,jup,sun);

    ofstream file1("../Outputs_Kirkwood3D/a0.txt");
    ofstream file2("../Outputs_Kirkwood3D/a.txt");
    ofstream file3("../Outputs_Kirkwood3D/initial_distribution.txt");
    ofstream file4("../Outputs_Kirkwood3D/final_distribution.txt");
    ofstream file5("../Outputs_Kirkwood3D/jup.txt");


    int N_ast = 5000; //number of asteroids
    object ast[N_ast];
    object M3[N_ast]; //asteroid copies
    distances SIM[N_ast];

    for(int i=0;i<N_ast;i++){
        ast[i].a = fRand(2.,jup.a);
        ast[i].M = fRand(0.,360.);
        ast[i].i = fRand(0.,5.);
        get_cartesian(GM,t,t0,ast[i]);
        ast[i].write_out_keplerian(file1);
        //write_in_COM(ast[i],COM,file3);
        ofstream file;
        file.open("../Outputs_Kirkwood3D/traj/ast"+to_string(i)+".txt");
    }

    //write_in_COM(jup,COM,file5);
    dt = adapt_dt(GM,ast,N_ast,1000.);
    cout<<"Timestep = "<<dt<<" seconds."<<endl;

    while(t < T){

        if((it%snapshot) == 0){
            cout<<"Simulation Progress: "<<floor(100.*t/T)+1<<"%\r";
            update_COM(COM,sun,jup);
            write_in_COM(jup,sun,file5);
            for(int i=0;i<N_ast;i++){
                ofstream file;
                file.open("../Outputs_Kirkwood3D/traj/ast"+to_string(i)+".txt",fstream::app);
                write_in_COM(ast[i],sun,file);
            }
        }
        RK4(sun,jup,M1,M2,M3,ast,SIM,N_ast,dt);
        t+=dt;
        it+=1;
    }
    update_COM(COM,jup,sun);
    //write_in_COM(jup,COM,file5);
    for(int i =0;i<N_ast;i++){
        get_keplerian(GM,ast[i],COM);
        ast[i].write_out_keplerian(file2);
        //write_in_COM(ast[i],COM,file4);
    }

    return 0;
}
