#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#include <object.h>
#include <distances.h>

void RK4(object &m1, object &m2,object &M1,object &M2,object M3[],object ast[],distances SIM[],int N_ast, double dt);
void write_in_COM(object m,object COM,ofstream& file);
void get_cartesian(double GM,double t,double t0,object &m);
void get_keplerian(double GM,object &m,object COM);
double orbital_period(double GM, double a);
void write_hamiltonian(object m1,object m2,object m3,double t,ofstream& file);
double fRand(double fMin, double fMax);
void update_COM(object &COM,object m1, object m2);
double adapt_dt(double GM,object ast[],int N_ast,double precision);

#endif // FUNCTIONS_H_INCLUDED
