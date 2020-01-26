#ifndef OBJECT_H
#define OBJECT_H
#include <iostream>
#include <stdio.h>
#include <fstream>

using namespace std;

class object
{
    public:
        object(double m=0.,double a=0.,double e=0.,double i=0.,double w=0.,double Om=0.,double M=0.);
        virtual ~object();

        double m;
        double a,e,i,w,Om,M;
        double x,y,z;
        double vx,vy,vz;

        void write_out(ofstream& file){
            file<<x<<"\t"<<y<<"\t"<<z<<"\t"<<vx<<"\t"<<vy<<"\t"<<vz<<endl;
        }

        void print_out(){
        cout <<"x= "<<x<<"\t"<<"y= "<<y<<"\t"<<"z= "<<z<<"\t"<<"vx= "<<vx<<"\t"<<"vy= "<<vy<<"\t"<<"vz= "<<vz<<endl;
        }

        void write_out_keplerian(ofstream& file){
            file<<a<<"\t"<<e<<"\t"<<i<<"\t"<<w<<"\t"<<Om<<"\t"<<M<<endl;
        }

        void print_out_keplerian(){
            cout<<a<<"\t"<<e<<"\t"<<i<<"\t"<<w<<"\t"<<Om<<"\t"<<M<<endl;
        }

    protected:

    private:
};

#endif // OBJECT_H
