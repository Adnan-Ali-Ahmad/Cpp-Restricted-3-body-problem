#ifndef DISTANCES_H
#define DISTANCES_H
#include "object.h"
#include <math.h>

extern double G;

class distances
{
    public:
        distances(double r12=0.,double r13=0.,double r23=0.);
        virtual ~distances();

        double r12, r13, r23;
        double r12_3, r13_3, r23_3;

        double rij(object &mi, object &mj){
            double x2 = (mj.x-mi.x)*(mj.x-mi.x);
            double y2 = (mj.y-mi.y)*(mj.y-mi.y);
            double z2 = (mj.z-mi.z)*(mj.z-mi.z);
            return sqrt(x2+y2+z2);
        }

        void get_distance_12(object &m1, object &m2){
            r12 = rij(m1,m2);
            r12_3 = r12*r12*r12;
        }

        void get_distances_13_23(object &m1, object &m2, object &m3){
            r13 = rij(m1,m3);
            r23 = rij(m2,m3);
            r13_3 = r13*r13*r13;
            r23_3 = r23*r23*r23;
        }

        void print_out(){
            cout <<"r12= "<<r12<<"\t"<<"r13= "<<r13<<"\t"<<"r23= "<<r23<<endl;
        }

    protected:

    private:
};

#endif // DISTANCES_H
