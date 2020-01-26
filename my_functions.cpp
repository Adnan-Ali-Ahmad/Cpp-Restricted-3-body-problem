#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include "functions.h"
#include <object.h>
#include <distances.h>

using namespace std;
extern double G;
extern double pi;

double rij(object &mi, object &mj){
    double x2 = (mj.x-mi.x)*(mj.x-mi.x);
    double y2 = (mj.y-mi.y)*(mj.y-mi.y);
    return sqrt(x2+y2);
}

double acc1(double m2,double x2, double x1, double r12_3){
    return G*m2*(x2-x1)/r12_3;
}

double acc2(double acc1, double m1, double m2){
    return -acc1*m1/m2;
}

double acc3x(object &M1,object &M2,object &M3,distances &SIM){
    return -G*M1.m*(M3.x-M1.x)/SIM.r13_3 -G*M2.m*(M3.x-M2.x)/SIM.r23_3;
}

double acc3y(object &M1,object &M2,object &M3,distances &SIM){
    return -G*M1.m*(M3.y-M1.y)/SIM.r13_3 -G*M2.m*(M3.y-M2.y)/SIM.r23_3;
}

double acc3z(object &M1,object &M2,object &M3,distances &SIM){
    return -G*M1.m*(M3.z-M1.z)/SIM.r13_3 -G*M2.m*(M3.z-M2.z)/SIM.r23_3;
}

void RK4(object &m1, object &m2,object &M1,object &M2,object M3[], object ast[],distances SIM[],int N_ast, double dt){

    M1 = m1; //temporary objects for RK purposes
    M2 = m2;

    double M1x,M2x,M1y,M2y,M1z,M2z;
    double M1vx,M1vy,M2vx,M2vy,M2vz,M1vz;
    double M3x[N_ast],M3y[N_ast],M3z[N_ast];
    double M3vx[N_ast],M3vy[N_ast],M3vz[N_ast];
    M1x = M1.x; M1y = M1.y; M1z = M1.z; M1vx = M1.vx; M1vy = M1.vy; M1vz = M1.vz;
    M2x = M2.x; M2y = M2.y; M2z = M2.z; M2vx = M2.vx; M2vy = M2.vy; M2vz = M2.vz;

    M3[0] = ast[0];
    M3x[0] = ast[0].x;
    M3y[0] = ast[0].y;
    M3z[0] = ast[0].z;
    M3vx[0] = ast[0].vx;
    M3vy[0] = ast[0].vy;
    M3vz[0] = ast[0].vz;
    SIM[0].get_distance_12(M1,M2);
    SIM[0].get_distances_13_23(M1,M2,M3[0]);

    for(int i=1;i<N_ast;i++){
        M3[i] = ast[i];
        M3x[i] = ast[i].x;
        M3y[i] = ast[i].y;
        M3z[i] = ast[i].z;
        M3vx[i] = ast[i].vx;
        M3vy[i] = ast[i].vy;
        M3vz[i] = ast[i].vz;
        SIM[i].r12 = SIM[0].r12;
        SIM[i].get_distances_13_23(M1,M2,M3[i]);
    }

    double k1x1, k2x1, k3x1, k4x1;
    double l1x1, l2x1, l3x1, l4x1;
    double k1y1, k2y1, k3y1, k4y1;
    double l1y1, l2y1, l3y1, l4y1;
    double k1z1, k2z1, k3z1, k4z1;
    double l1z1, l2z1, l3z1, l4z1;

    double k1x2, k2x2, k3x2, k4x2;
    double l1x2, l2x2, l3x2, l4x2;
    double k1y2, k2y2, k3y2, k4y2;
    double l1y2, l2y2, l3y2, l4y2;
    double k1z2, k2z2, k3z2, k4z2;
    double l1z2, l2z2, l3z2, l4z2;

    double k1x3[N_ast], k2x3[N_ast], k3x3[N_ast], k4x3[N_ast];
    double l1x3[N_ast], l2x3[N_ast], l3x3[N_ast], l4x3[N_ast];
    double k1y3[N_ast], k2y3[N_ast], k3y3[N_ast], k4y3[N_ast];
    double l1y3[N_ast], l2y3[N_ast], l3y3[N_ast], l4y3[N_ast];
    double k1z3[N_ast], k2z3[N_ast], k3z3[N_ast], k4z3[N_ast];
    double l1z3[N_ast], l2z3[N_ast], l3z3[N_ast], l4z3[N_ast];

    double ax1,ay1,az1,ax2,ay2,az2,ax3,ay3,az3;

    ax1 = acc1(M2.m,M2.x,M1.x,SIM[0].r12_3);
    ay1 = acc1(M2.m,M2.y,M1.y,SIM[0].r12_3);
    az1 = acc1(M2.m,M2.z,M1.z,SIM[0].r12_3);
    k1x1 = dt*ax1;
    k1y1 = dt*ay1;
    k1z1 = dt*az1;
    l1x1 = dt*M1.vx;
    l1y1 = dt*M1.vy;
    l1z1 = dt*M1.vz;

    ax2 = acc2(ax1,M1.m,M2.m);
    ay2 = acc2(ay1,M1.m,M2.m);
    az2 = acc2(az1,M1.m,M2.m);
    k1x2 = dt*ax2;
    k1y2 = dt*ay2;
    k1z2 = dt*az2;
    l1x2 = dt*M2.vx;
    l1y2 = dt*M2.vy;
    l1z2 = dt*M2.vz;

    for(int i=0;i<N_ast;i++){
        ax3 = acc3x(M1,M2,M3[i],SIM[i]);
        ay3 = acc3y(M1,M2,M3[i],SIM[i]);
        az3 = acc3z(M1,M2,M3[i],SIM[i]);
        k1x3[i] = dt*ax3;
        k1y3[i] = dt*ay3;
        k1z3[i] = dt*az3;
        l1x3[i] = dt*M3[i].vx;
        l1y3[i] = dt*M3[i].vy;
        l1z3[i] = dt*M3[i].vz;
        M3[i].x = M3x[i] + 0.5*l1x3[i];
        M3[i].y = M3y[i] + 0.5*l1y3[i];
        M3[i].z = M3z[i] + 0.5*l1z3[i];
        M3[i].vx = M3vx[i] + 0.5*k1x3[i];
        M3[i].vy = M3vy[i] + 0.5*k1y3[i];
        M3[i].vz = M3vz[i] + 0.5*k1z3[i];
        SIM[i].get_distances_13_23(M1,M2,M3[i]);
    }

    M1.x = M1x + 0.5*l1x1; M1.y = M1y + 0.5*l1y1; M1.z = M1z + 0.5*l1z1;
    M1.vx = M1vx + 0.5*k1x1; M1.vy = M1vy + 0.5*k1y1; M1.vz = M1vz + 0.5*k1z1;
    M2.x = M2x + 0.5*l1x2; M2.y = M2y + 0.5*l1y2; M2.z = M2z + 0.5*l1z2;
    M2.vx = M2vx + 0.5*k1x2; M2.vy = M2vy + 0.5*k1y2; M2.vz = M2vz + 0.5*k1z2;
    SIM[0].get_distance_12(M1,M2);

    ax1 = acc1(M2.m,M2.x,M1.x,SIM[0].r12_3);
    ay1 = acc1(M2.m,M2.y,M1.y,SIM[0].r12_3);
    az1 = acc1(M2.m,M2.z,M1.z,SIM[0].r12_3);
    k2x1 = dt*ax1;
    k2y1 = dt*ay1;
    k2z1 = dt*az1;
    l2x1 = dt*M1.vx;
    l2y1 = dt*M1.vy;
    l2z1 = dt*M1.vz;

    ax2 = acc2(ax1,M1.m,M2.m);
    ay2 = acc2(ay1,M1.m,M2.m);
    az2 = acc2(az1,M1.m,M2.m);
    k2x2 = dt*ax2;
    k2y2 = dt*ay2;
    k2z2 = dt*az2;
    l2x2 = dt*M2.vx;
    l2y2 = dt*M2.vy;
    l2z2 = dt*M2.vz;

    for(int i=0;i<N_ast;i++){
        ax3 = acc3x(M1,M2,M3[i],SIM[i]);
        ay3 = acc3y(M1,M2,M3[i],SIM[i]);
        az3 = acc3z(M1,M2,M3[i],SIM[i]);
        k2x3[i] = dt*ax3;
        k2y3[i] = dt*ay3;
        k2z3[i] = dt*az3;
        l2x3[i] = dt*M3[i].vx;
        l2y3[i] = dt*M3[i].vy;
        l2z3[i] = dt*M3[i].vz;
        M3[i].x = M3x[i] + 0.5*l2x3[i];
        M3[i].y = M3y[i] + 0.5*l2y3[i];
        M3[i].z = M3z[i] + 0.5*l2z3[i];
        M3[i].vx = M3vx[i] + 0.5*k2x3[i];
        M3[i].vy = M3vy[i] + 0.5*k2y3[i];
        M3[i].vz = M3vz[i] + 0.5*k2z3[i];
        SIM[i].get_distances_13_23(M1,M2,M3[i]);
    }

    M1.x = M1x + 0.5*l2x1; M1.y = M1y + 0.5*l2y1; M1.z = M1z + 0.5*l2z1;
    M1.vx = M1vx + 0.5*k2x1; M1.vy = M1vy + 0.5*k2y1; M1.vz = M1vz + 0.5*k2z1;
    M2.x = M2x + 0.5*l2x2; M2.y = M2y + 0.5*l2y2; M2.z = M2z + 0.5*l2z2;
    M2.vx = M2vx + 0.5*k2x2; M2.vy = M2vy + 0.5*k2y2; M2.vz = M2vz + 0.5*k2z2;
    SIM[0].get_distance_12(M1,M2);

    ax1 = acc1(M2.m,M2.x,M1.x,SIM[0].r12_3);
    ay1 = acc1(M2.m,M2.y,M1.y,SIM[0].r12_3);
    az1 = acc1(M2.m,M2.z,M1.z,SIM[0].r12_3);
    k3x1 = dt*ax1;
    k3y1 = dt*ay1;
    k3z1 = dt*az1;
    l3x1 = dt*M1.vx;
    l3y1 = dt*M1.vy;
    l3z1 = dt*M1.vz;

    ax2 = acc2(ax1,M1.m,M2.m);
    ay2 = acc2(ay1,M1.m,M2.m);
    az2 = acc2(az1,M1.m,M2.m);
    k3x2 = dt*ax2;
    k3y2 = dt*ay2;
    k3z2 = dt*az2;
    l3x2 = dt*M2.vx;
    l3y2 = dt*M2.vy;
    l3z2 = dt*M2.vz;

    for(int i=0;i<N_ast;i++){
        ax3 = acc3x(M1,M2,M3[i],SIM[i]);
        ay3 = acc3y(M1,M2,M3[i],SIM[i]);
        az3 = acc3z(M1,M2,M3[i],SIM[i]);
        k3x3[i] = dt*ax3;
        k3y3[i] = dt*ay3;
        k3z3[i] = dt*az3;
        l3x3[i] = dt*M3[i].vx;
        l3y3[i] = dt*M3[i].vy;
        l3z3[i] = dt*M3[i].vz;
        M3[i].x = M3x[i] + l3x3[i];
        M3[i].y = M3y[i] + l3y3[i];
        M3[i].z = M3z[i] + l3z3[i];
        M3[i].vx = M3vx[i] + k3x3[i];
        M3[i].vy = M3vy[i] + k3y3[i];
        M3[i].vz = M3vz[i] + k3z3[i];
        SIM[i].get_distances_13_23(M1,M2,M3[i]);
    }

    M1.x = M1x + l3x1; M1.y = M1y + l3y1; M1.z = M1z + l3z1;
    M1.vx = M1vx + k3x1; M1.vy = M1vy + k3y1; M1.vz = M1vz + k3z1;
    M2.x = M2x + l3x2; M2.y = M2y + l3y2; M2.z = M2z + l3z2;
    M2.vx = M2vx + k3x2; M2.vy = M2vy + k3y2; M2.vz = M2vz + k3z2;
    SIM[0].get_distance_12(M1,M2);

    ax1 = acc1(M2.m,M2.x,M1.x,SIM[0].r12_3);
    ay1 = acc1(M2.m,M2.y,M1.y,SIM[0].r12_3);
    az1 = acc1(M2.m,M2.z,M1.z,SIM[0].r12_3);
    k4x1 = dt*ax1;
    k4y1 = dt*ay1;
    k4z1 = dt*az1;
    l4x1 = dt*M1.vx;
    l4y1 = dt*M1.vy;
    l4z1 = dt*M1.vz;

    ax2 = acc2(ax1,M1.m,M2.m);
    ay2 = acc2(ay1,M1.m,M2.m);
    az2 = acc2(az1,M1.m,M2.m);
    k4x2 = dt*ax2;
    k4y2 = dt*ay2;
    k4z2 = dt*az2;
    l4x2 = dt*M2.vx;
    l4y2 = dt*M2.vy;
    l4z2 = dt*M2.vz;

    for(int i=0;i<N_ast;i++){
        ax3 = acc3x(M1,M2,M3[i],SIM[i]);
        ay3 = acc3y(M1,M2,M3[i],SIM[i]);
        az3 = acc3z(M1,M2,M3[i],SIM[i]);
        k4x3[i] = dt*ax3;
        k4y3[i] = dt*ay3;
        k4z3[i] = dt*az3;
        l4x3[i] = dt*M3[i].vx;
        l4y3[i] = dt*M3[i].vy;
        l4z3[i] = dt*M3[i].vz;
    }

    m1.vx += 1./6. * (k1x1+2*k2x1+2*k3x1+k4x1);
    m1.vy += 1./6. * (k1y1+2*k2y1+2*k3y1+k4y1);
    m1.vz += 1./6. * (k1z1+2*k2z1+2*k3z1+k4z1);
    m1.x += 1./6. * (l1x1+2*l2x1+2*l3x1+l4x1);
    m1.y += 1./6. * (l1y1+2*l2y1+2*l3y1+l4y1);
    m1.z += 1./6. * (l1z1+2*l2z1+2*l3z1+l4z1);

    m2.vx += 1./6. * (k1x2+2*k2x2+2*k3x2+k4x2);
    m2.vy += 1./6. * (k1y2+2*k2y2+2*k3y2+k4y2);
    m2.vz += 1./6. * (k1z2+2*k2z2+2*k3z2+k4z2);
    m2.x += 1./6. * (l1x2+2*l2x2+2*l3x2+l4x2);
    m2.y += 1./6. * (l1y2+2*l2y2+2*l3y2+l4y2);
    m2.z += 1./6. * (l1z2+2*l2z2+2*l3z2+l4z2);

    for(int i=0;i<N_ast;i++){
        ast[i].vx += 1./6. * (k1x3[i]+2.*k2x3[i]+2.*k3x3[i]+k4x3[i]);
        ast[i].vy += 1./6. * (k1y3[i]+2.*k2y3[i]+2.*k3y3[i]+k4y3[i]);
        ast[i].vz += 1./6. * (k1z3[i]+2.*k2z3[i]+2.*k3z3[i]+k4z3[i]);
        ast[i].x += 1./6. * (l1x3[i]+2.*l2x3[i]+2.*l3x3[i]+l4x3[i]);
        ast[i].y += 1./6. * (l1y3[i]+2.*l2y3[i]+2.*l3y3[i]+l4y3[i]);
        ast[i].z += 1./6. * (l1z3[i]+2.*l2z3[i]+2.*l3z3[i]+l4z3[i]);
    }

}

void get_cartesian(double GM,double t,double t0,object &m){
//get cartesian elements from keplerian elements
//MAKE SURE THAT t AND t0 ARE IN SECONDS
    double a = m.a*1.496e11;
    double e = m.e;
    //1: setting M(t)
    double M;
    if(t==t0){
        M = m.M*pi/180.;
    }
    else {
        double dt = t-t0;
        if (a == 0.){
            M = 0.;
        } else {
            M = m.M*pi/180.+dt*sqrt(GM/(a*a*a));
        }
    }
    //2: Computing eccentric anomaly with Newton-Raphson
    double E = M; double E0 = E;
    for(int i=0;i<500;i++){
        E0 = E;
        E = E - (E-e*sin(E)-M)/(1-e*cos(E));
        if(fabs(E-E0) < 1e-8){
            break;
        }
    }
    //3: Obtaining true anomaly
    double mu = 2.*atan2(sqrt(1.+e)*sin(0.5*E),sqrt(1.-e)*cos(0.5*E));
    //4: Obtaining distance to central body
    double rc = a*(1.-e*cos(E));
    //5: Obtaining x,y,z,vx,vy,vz in the orbital frame
    m.x = rc*(cos(mu)); m.vx = sqrt(GM*a)/rc *-sin(E);
    m.y = rc*(sin(mu)); m.vy = sqrt(GM*a)/rc *sqrt(1-e*e)*cos(E);
    m.z = 0.; m.vz = 0.;
    if(rc == 0.){
        m.vx = 0.;
        m.vy = 0.;
        m.vz = 0.;
    }
    //6: Obtaining x,y,z,vx,vy,vz in the inertial frame
    double w = m.w*pi/180.;
    double Om = m.Om*pi/180.;
    double i = m.i*pi/180.;
    double x = m.x; double y = m.y;
    double vx = m.vx; double vy = m.vy;

    m.x = x*(cos(w)*cos(Om)-sin(w)*cos(i)*sin(Om)) - y*(sin(w)*cos(Om)+cos(w)*cos(i)*sin(Om));
    m.y = x*(cos(w)*sin(Om)+sin(w)*cos(i)*cos(Om)) + y*(cos(w)*cos(i)*cos(Om)-sin(w)*sin(Om));
    m.z = x*(sin(w)*sin(i)) + y*(cos(w)*sin(i));
    m.vx = vx*(cos(w)*cos(Om)-sin(w)*cos(i)*sin(Om)) - vy*(sin(w)*cos(Om)+cos(w)*cos(i)*sin(Om));
    m.vy = vx*(cos(w)*sin(Om)+sin(w)*cos(i)*cos(Om)) + vy*(cos(w)*cos(i)*cos(Om)-sin(w)*sin(Om));
    m.vz = vx*(sin(w)*sin(i)) + vy*(cos(w)*sin(i));

}

void get_keplerian(double GM,object &m,object COM){
//get keplerian elements from cartesian elements
    //1: Preparations
    double x = m.x-COM.x;
    double y = m.y-COM.y;
    double z = m.z-COM.z;
    double vx = m.vx-COM.vx;
    double vy = m.vy-COM.vy;
    double vz = m.vz-COM.vz;
    double r = sqrt(x*x+y*y+z*z);
    double v = sqrt(vx*vx+vy*vy+vz*vz);
    //angular momentum:
    double hx = y*vz-z*vy;
	double hy = z*vx-x*vz;
	double hz = x*vy-y*vx;
	double h = sqrt(hx*hx+hy*hy+hz*hz);
	//obtaining eccentricity:
	double ex = (vy*hz-vz*hy)/GM-x/r;
	double ey = (vz*hx-vx*hz)/GM-y/r;
	double ez = (vx*hy-vy*hx)/GM-z/r;
	m.e = sqrt(ex*ex+ey*ey+ez*ez);
	//obtaining node vector:
	double nx = hy;
	double ny = -hx;
	double n=sqrt(nx*nx+ny*ny);
	//obtaining semimajor axis:
	m.a = -GM/2./(v*v/2. - GM/r);
	//2: Obtaining inclination
	m.i = acos(hz/h);
	//3: Obtaining longitude of the ascending node and argument of periapsis
	if(fabs(sin(m.i)) < 1e-8){
        m.Om = 0.;
	} else {
        m.Om = acos(-nx/n);
        if(m.Om && ny>0.){
            m.Om = 2.*pi-m.Om;
        }
	}
    if(m.e == 0.){
        m.w = 0.;
    }
	else {
		m.w = acos((cos(m.Om)*ex+sin(m.Om)*ey)/m.e);
		if(m.w && ez<0.){
            m.w = 2.*pi-m.w;
		}
	}
	//4: Obtaining true anomaly and eccentric anomaly
	double cosi = cos(hz/h);
	double comega = cos(m.Om);
	double somega = sin(m.Om);
	double cosR = x*comega+y*somega;
	double sinR = (y*comega-x*somega)*cosi+z*sin(m.i);
	double norm = sqrt(cosR*cosR+sinR*sinR);
    cosR /= norm;
    sinR /= norm;
	double T = acos(cosR);
	if(sinR < 0){
        T = 2*pi-T;
	}
	T -= m.w;
	double E = acos((cos(T)+m.e)/(1.+m.e*cos(T)));
	if(sin(T) < 0.){
        E = 2.*pi-E;
	}
	//5: Obtaining Mean anomaly
	m.M = E-m.e*sin(E);

	//6: Converting
	m.a /= 1.496e11;
	m.i *= 180./pi;
	m.w *= 180./pi;
	m.Om *= 180./pi;
	m.M *= 180./pi;
}

void write_in_COM(object m,object COM,ofstream& file){
    //writes coordinates of m1 m2 and m3 in the center of mass
    m.x -= COM.x; m.y -= COM.y; m.z -= COM.z;
    m.vx -= COM.vx; m.vy -= COM.vy; m.vz -= COM.vz;
    m.write_out(file);
}

double orbital_period(double GM, double a){
    a *= 1.496e11;
    double T = sqrt(4*pi*pi*a*a*a/GM);
    return T;
}

void update_COM(object &COM,object m1, object m2){
    COM.x = (m1.m*m1.x+m2.m*m2.x)/(m1.m+m2.m);
    COM.y = (m1.m*m1.y+m2.m*m2.y)/(m1.m+m2.m);
    COM.z = (m1.m*m1.z+m2.m*m2.z)/(m1.m+m2.m);
    COM.vx = (m1.m*m1.vx+m2.m*m2.vx)/(m1.m+m2.m);
    COM.vy = (m1.m*m1.vy+m2.m*m2.vy)/(m1.m+m2.m);
    COM.vz = (m1.m*m1.vz+m2.m*m2.vz)/(m1.m+m2.m);
}

double fRand(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double adapt_dt(double GM,object ast[],int N_ast,double precision){
    double amin = ast[0].a;
    for(int i=1;i<N_ast;i++){
        if(ast[i].a<amin){
            amin = ast[i].a;
        }
    }
    return orbital_period(GM,amin)/precision;
}
