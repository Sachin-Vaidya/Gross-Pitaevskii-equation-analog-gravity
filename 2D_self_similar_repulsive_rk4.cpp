#include<iostream>
#include<cmath>
#include<complex>
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include <fstream>
#include <omp.h>
#include<random>
#include<chrono>
//#include<sys/types.h>
//#include<sys/stat.h>
#include</home/vaidya2/build_dir/source_dir/Eigen/Dense>
//#include "home/vaidya2/ondemand.halstead/data/sys/myjobs/projects/default/1/build_dir/Eigen/Core"

#include<iomanip>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef complex<double> dcomp;

#define Zmax 100

double B, r0;

double ode1(double v, double rho, double drho);
double ode2(double v, double rho, double drho);
void RK4(double v0, double rho0, double drho0, double h);



// -----------------------------CITATION-------------------------------------//
//https://stackoverflow.com/questions/49041945/runge-kutta-4th-order-to-solve-2nd-order-ode-using-c
// Autor    : Carlos Eduardo da Silva Lima
// Tema     : Runge-Kutta de quarta ordem
// Linguagem: C++ (ANSI)
// IDE      : Online GBD
// Data     : 19/11/2022
// --------------------------End of CITATION---------------------------------//


    
int main()
{
    
    RK4(Zmax, 2./Zmax, -2.*pow(Zmax,-2.), 0.001);
    
    return 0;
}


double ode1(double v, double rho, double drho)
{
    return drho;
}

double ode2(double v, double rho, double drho)
{
    return (-drho / v - rho * v * v / 0.16e2 + B * B * pow(rho, -0.3e1) * pow(v, -0.2e1) + pow(rho, 0.3e1));
}

void RK4(double v0, double rho0, double drho0, double h)
{
    cout<<"h = "<<h<<endl;
    
    int N = abs(round(Zmax/h));
    
    cout<<N<<endl;
    
    std::string B_val("B_"), h_val("_h_"), folder("2D_self_similar_new_sol/"), rk4("_self_similar_repulsive_rk4_new_sol_1_-2"), file_for_B_r0_val("r0_vs_B2"), extension(".txt");
    
    string r0_B_file = file_for_B_r0_val+extension;
    
    std::ifstream r0_B;
    r0_B.open(r0_B_file.c_str());
    
    char num[1000], num1[1000], num2[1000];
    
    int num_files = 0;
    do
    {
        r0_B>>B>>r0;
        
        cout<<B<<"\t"<<r0<<endl;
    
        double v[N+1], rho[N+1], drho[N+1];
        
        v[N] = v0;
        rho[N] = 1.;//rho0*sqrt(B);
        drho[N] = -2.;//drho0*sqrt(B);
        
        cout<<v[N]<<"\t"<<rho[N]<<"\t"<<drho[N]<<endl;
        
        int i = N-1;
        do
        {
            
            double k11 = ode1(v[i+1],rho[i+1],drho[i+1]);
            double k12 = ode2(v[i+1],rho[i+1],drho[i+1]);
            
            double k21 = ode1(v[i+1]-h*(1./2),rho[i+1]-h*(k11/2),drho[i+1]-h*(k12/2));
            double k22 = ode2(v[i+1]-h*(1./2),rho[i+1]-h*(k11/2),drho[i+1]-h*(k12/2));
            
            double k31 = ode1(v[i+1]-h*(1./2),rho[i+1]-h*(k21/2),drho[i+1]-h*(k22/2));
            double k32 = ode2(v[i+1]-h*(1./2),rho[i+1]-h*(k21/2),drho[i+1]-h*(k22/2));
            
            double k41 = ode1(v[i+1]-h,rho[i+1]-h*k31,drho[i+1]-h*k32);
            double k42 = ode2(v[i+1]-h,rho[i+1]-h*k31,drho[i+1]-h*k32);
            
            rho[i] = rho[i+1] - h*((k11+2.*(k21+k31)+k41)/6);
            drho[i] = drho[i+1] - h*((k12+2.*(k22+k32)+k42)/6);
            v[i] = v[i+1] - h;
            
            i--;
        }
        while(i>=0);
        
        sprintf(num1,"%.3f",B);
        sprintf(num2,"%.5f",abs(h));
        string rk4_file = folder+B_val+num1+h_val+num2+rk4+extension;
    
        std::ofstream rk4_output;
        rk4_output.open(rk4_file.c_str());
        
        int j = 0;
        do
        {
            rk4_output <<std::setprecision(16)<<v[j]<<"\t"<<rho[j]<<"\t"<<drho[j]<<endl;
            j++;
        }
        while(j<=N);
        rk4_output.close();
        
        num_files+=1;
    }
    while(num_files<2);
    r0_B.close();
}