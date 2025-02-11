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
//#include "home/vaidya2/build_dir/Eigen/Core"

#include<iomanip>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef complex<double> dcomp;

//----------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------//
               //-------------------solving in x----------------------//
//----------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------//

// Laplace Conformal Borel Pade integral

int main()
{
    double dummy;
    double Rmin = 0.;
    double Rmax = 100.;
    double h;
    double h2 = 0.001;
    int N;
    int n = round(10./h2);
    double pi = 3.141592653589793238;
    double epsilon;
    
    double Res0, Res1, Res2, Res3, Res4, Res5, Res6, Res7, Res8, Res9, Res10; 
    
    int Order[4] = {96,100,192,200};
    
    int p0, p1, p2;
    double R_0to1_p0, R_0to1_p1, R_0to1_p2, R_0to1_factor_p0, R_0to1_factor_p1, R_0to1_factor_p2;
    double R_1toInf_p0, R_1toInf_p1, R_1toInf_p2, R_1toInf_factor_p0, R_1toInf_factor_p1, R_1toInf_factor_p2;
    
    //B = 1
    
    /*//B = 2.5
    Res0 = 11459.173690383833;
    Res1 = 11459.173690383833;
    Res2 = 8594.380267787874;
    Res3 = 5729.586845191916;
    Res4 = 3580.991778244948;
    Res5 = 2148.595066946969;*/
    
    double f0_Re[n], f0_Im[n], f1_Re[n], f1_Im[n], f2_Re[n], f2_Im[n], f3_Re[n], f3_Im[n], f4_Re[n], f4_Im[n], f5_Re[n], f5_Im[n], f6_Re[n], f6_Im[n], f7_Re[n], f7_Im[n], f8_Re[n], f8_Im[n], f9_Re[n], f9_Im[n], f10_Re[n], f10_Im[n];
    double f0_Re_LCBP[n], f0_Im_LCBP[n], f1_Re_LCBP[n], f1_Im_LCBP[n], f2_Re_LCBP[n], f2_Im_LCBP[n], f3_Re_LCBP[n], f3_Im_LCBP[n], f4_Re_LCBP[n], f4_Im_LCBP[n], f5_Re_LCBP[n], f5_Im_LCBP[n], f6_Re_LCBP[n], f6_Im_LCBP[n], f7_Re_LCBP[n], f7_Im_LCBP[n], f8_Re_LCBP[n], f8_Im_LCBP[n], f9_Re_LCBP[n], f9_Im_LCBP[n], f10_Re_LCBP[n], f10_Im_LCBP[n];
    double H_sing_CBP_Re[n], H_sing_CBP_Im[n];
    
    std::string file_input("B1_CBP_singularity_subtraction_full_range_epsilon_0_utoz_0to1_O"), file_input2("B1_CBP_singularity_subtraction_full_range_epsilon_0_u_1toInf_x_0to1_O"), file_sing_input("B1_CBP_singularity_subtraction_epsilon_0"), file_output("B1_LCBP_singularity_subtraction_full_range_epsilon_0_Simpsons_O"), extension(".txt");
    
    char num[1000];
    
    char order[1][1000] = { "" };
    char NEW[1][10] = { "1" };
    
    for(int p=0;p<4;p++)
    {
        sprintf(num,"%d",Order[p]);
        for(int q=0;q<1;q++)
        {
            if(p==0)
            {
                cout<<"O96"<<endl;
                //O(u^96)
                Res0 = -8.8071673034896043222;
                Res1 = -8.8071673034896043213;
                Res2 = -6.6053754776172032424;
                Res3 = -4.4035836517448021705;
                Res4 = -2.7522397823405012216;
                Res5 = -1.6513438694043008087;
                Res6 = -0.96328392381917564859;
                Res7 = -0.55044795646810029732;
                Res8 = -0.30962697551330622269;
                Res9 = -0.17201498639628141773;
                Res10 = -0.094608242517954546170;
            }
            else if(p==1)
            {
                cout<<"O100"<<endl;
                //O(u^100)
                Res0 = -8.8071673034896043222; 
                Res1 = -8.8071673034896043222; 
                Res2 = -6.6053754776172032414; 
                Res3 = -4.4035836517448021628; 
                Res4 = -2.7522397823405013499; 
                Res5 = -1.6513438694043008157;
                Res6 = -0.96328392381917547545; 
                Res7 = -0.55044795646810028597; 
                Res8 = -0.30962697551330639240; 
                Res9 = -0.17201498639628135286;
                Res10 = -0.094608242517954746890;
            }
            else if(p==2)
            {
                cout<<"O192"<<endl;
                //O(u^192)
                Res0 = -8.8071673034896043222;
                Res1 = -8.8071673034896043222; 
                Res2 = -6.6053754776172032417; 
                Res3 = -4.4035836517448021611; 
                Res4 = -2.7522397823405013507; 
                Res5 = -1.6513438694043008104; 
                Res6 = -0.96328392381917547274; 
                Res7 = -0.55044795646810027014; 
                Res8 = -0.30962697551330640195; 
                Res9 = -0.17201498639628133442; 
                Res10 = -0.094608242517954733930;
            }
            else
            {
                cout<<"O200"<<endl;
                //O(u^200)
                Res0 = -8.8071673034896043222;
                Res1 = -8.8071673034896043222;
                Res2 = -6.6053754776172032417;
                Res3 = -4.4035836517448021611; 
                Res4 = -2.7522397823405013507; 
                Res5 = -1.6513438694043008104; 
                Res6 = -0.96328392381917547274; 
                Res7 = -0.55044795646810027014; 
                Res8 = -0.30962697551330640195; 
                Res9 = -0.17201498639628133442; 
                Res10 = -0.094608242517954733930;
            }
            
            if(q==0)
            {
                epsilon = 0.;
                N = 100000;
                h = 1./N;
            }
            /*else if(q==1)
            {
                epsilon = 0.;
                N = round((1./6)*(180 + 999 + 1000 + 17800 + 19980));
            }*/
            
            double R_0to1[N+1], H0_CBP_Re_0to1[N+1], H0_CBP_Im_0to1[N+1], H1_CBP_Re_0to1[N+1], H1_CBP_Im_0to1[N+1], H2_CBP_Re_0to1[N+1], H2_CBP_Im_0to1[N+1], H3_CBP_Re_0to1[N+1], H3_CBP_Im_0to1[N+1], H4_CBP_Re_0to1[N+1], H4_CBP_Im_0to1[N+1], H5_CBP_Re_0to1[N+1], H5_CBP_Im_0to1[N+1], H6_CBP_Re_0to1[N+1], H6_CBP_Im_0to1[N+1], H7_CBP_Re_0to1[N+1], H7_CBP_Im_0to1[N+1], H8_CBP_Re_0to1[N+1], H8_CBP_Im_0to1[N+1], H9_CBP_Re_0to1[N+1], H9_CBP_Im_0to1[N+1], H10_CBP_Re_0to1[N+1], H10_CBP_Im_0to1[N+1];
            double R_1toInf[N+1], H0_CBP_Re_1toInf[N+1], H0_CBP_Im_1toInf[N+1], H1_CBP_Re_1toInf[N+1], H1_CBP_Im_1toInf[N+1], H2_CBP_Re_1toInf[N+1], H2_CBP_Im_1toInf[N+1], H3_CBP_Re_1toInf[N+1], H3_CBP_Im_1toInf[N+1], H4_CBP_Re_1toInf[N+1], H4_CBP_Im_1toInf[N+1], H5_CBP_Re_1toInf[N+1], H5_CBP_Im_1toInf[N+1], H6_CBP_Re_1toInf[N+1], H6_CBP_Im_1toInf[N+1], H7_CBP_Re_1toInf[N+1], H7_CBP_Im_1toInf[N+1], H8_CBP_Re_1toInf[N+1], H8_CBP_Im_1toInf[N+1], H9_CBP_Re_1toInf[N+1], H9_CBP_Im_1toInf[N+1], H10_CBP_Re_1toInf[N+1], H10_CBP_Im_1toInf[N+1];
            
            if(q==0)
            {
                string input_file = file_input+num+extension;
                std::ifstream input(input_file.c_str());
                cout<<input_file<<endl;
                if(input.is_open())
                {
                    cout<<"opened the file"<<endl;
                    for(int number=0;number<N+1;number++)
                    {
                        input>>R_0to1[number]>>H0_CBP_Re_0to1[number]>>H0_CBP_Im_0to1[number]>>H1_CBP_Re_0to1[number]>>H1_CBP_Im_0to1[number]>>H2_CBP_Re_0to1[number]>>H2_CBP_Im_0to1[number]>>H3_CBP_Re_0to1[number]>>H3_CBP_Im_0to1[number]>>H4_CBP_Re_0to1[number]>>H4_CBP_Im_0to1[number]>>H5_CBP_Re_0to1[number]>>H5_CBP_Im_0to1[number]>>H6_CBP_Re_0to1[number]>>H6_CBP_Im_0to1[number]>>H7_CBP_Re_0to1[number]>>H7_CBP_Im_0to1[number]>>H8_CBP_Re_0to1[number]>>H8_CBP_Im_0to1[number]>>H9_CBP_Re_0to1[number]>>H9_CBP_Im_0to1[number]>>H10_CBP_Re_0to1[number]>>H10_CBP_Im_0to1[number];
                    }
                    input.close();
                }
                else
                {
                    cout<<"can't open the file"<<endl;
                }
                
                string input2_file = file_input2+num+extension;
                std::ifstream input2(input2_file.c_str());
                cout<<input2_file<<endl;
                if(input2.is_open())
                {
                    cout<<"opened the file"<<endl;
                    for(int number=0;number<N+1;number++)
                    {
                        input2>>R_1toInf[number]>>H0_CBP_Re_1toInf[number]>>H0_CBP_Im_1toInf[number]>>H1_CBP_Re_1toInf[number]>>H1_CBP_Im_1toInf[number]>>H2_CBP_Re_1toInf[number]>>H2_CBP_Im_1toInf[number]>>H3_CBP_Re_1toInf[number]>>H3_CBP_Im_1toInf[number]>>H4_CBP_Re_1toInf[number]>>H4_CBP_Im_1toInf[number]>>H5_CBP_Re_1toInf[number]>>H5_CBP_Im_1toInf[number]>>H6_CBP_Re_1toInf[number]>>H6_CBP_Im_1toInf[number]>>H7_CBP_Re_1toInf[number]>>H7_CBP_Im_1toInf[number]>>H8_CBP_Re_1toInf[number]>>H8_CBP_Im_1toInf[number]>>H9_CBP_Re_1toInf[number]>>H9_CBP_Im_1toInf[number]>>H10_CBP_Re_1toInf[number]>>H10_CBP_Im_1toInf[number];
                    }
                    input2.close();
                }
                else
                {
                    cout<<"can't open the file"<<endl;
                }
                
                string sing_input_file = file_sing_input+extension;
                std::ifstream sing_input(sing_input_file.c_str());
                cout<<sing_input_file<<endl;
                if(sing_input.is_open())
                {
                    cout<<"opened the file"<<endl;
                    for(int number=0;number<=n-1;number++)
                    {
                        sing_input>>dummy>>H_sing_CBP_Re[number]>>H_sing_CBP_Im[number];
                    }
                    sing_input.close();
                }
                else
                {
                    cout<<"can't open the file"<<endl;
                }
            }
            
            for(int k=1;k<n+1;k++)
            {
                f0_Re[k-1] = 0.;
                f0_Im[k-1] = 0.;
                f1_Re[k-1] = 0.;
                f1_Im[k-1] = 0.;
                f2_Re[k-1] = 0.;
                f2_Im[k-1] = 0.;
                f3_Re[k-1] = 0.;
                f3_Im[k-1] = 0.;
                f4_Re[k-1] = 0.;
                f4_Im[k-1] = 0.;
                f5_Re[k-1] = 0.;
                f5_Im[k-1] = 0.;
                f6_Re[k-1] = 0.;
                f6_Im[k-1] = 0.;
                f7_Re[k-1] = 0.;
                f7_Im[k-1] = 0.;
                f8_Re[k-1] = 0.;
                f8_Im[k-1] = 0.;
                f9_Re[k-1] = 0.;
                f9_Im[k-1] = 0.;
                f10_Re[k-1] = 0.;
                f10_Im[k-1] = 0.;
                
                for(int j=0;j<=round(N/2)-1;j++)
                {
                    p0 = 2*j;
                    p1 = 2*j+1;
                    p2 = 2*j+2;
                    
                    R_0to1_p0 = 2.*R_0to1[p0]/(1+pow(R_0to1[p0],2.));
                    R_0to1_p1 = 2.*R_0to1[p1]/(1+pow(R_0to1[p1],2.));
                    R_0to1_p2 = 2.*R_0to1[p2]/(1+pow(R_0to1[p2],2.));
                    
                    R_0to1_factor_p0 = 2.*(1-pow(R_0to1[p0],2.))*pow(1+pow(R_0to1[p0],2.),-2.);
                    R_0to1_factor_p1 = 2.*(1-pow(R_0to1[p1],2.))*pow(1+pow(R_0to1[p1],2.),-2.);
                    R_0to1_factor_p2 = 2.*(1-pow(R_0to1[p2],2.))*pow(1+pow(R_0to1[p2],2.),-2.);
                    
                    f0_Re[k-1] += h*(1./3)*(H0_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H0_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H0_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f0_Im[k-1] += h*(1./3)*(H0_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H0_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H0_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f1_Re[k-1] += h*(1./3)*(H1_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H1_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H1_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f1_Im[k-1] += h*(1./3)*(H1_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H1_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H1_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f2_Re[k-1] += h*(1./3)*(H2_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H2_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H2_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f2_Im[k-1] += h*(1./3)*(H2_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H2_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H2_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f3_Re[k-1] += h*(1./3)*(H3_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H3_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H3_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f3_Im[k-1] += h*(1./3)*(H3_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H3_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H3_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f4_Re[k-1] += h*(1./3)*(H4_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H4_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H4_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f4_Im[k-1] += h*(1./3)*(H4_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H4_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H4_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f5_Re[k-1] += h*(1./3)*(H5_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H5_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H5_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f5_Im[k-1] += h*(1./3)*(H5_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H5_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H5_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f6_Re[k-1] += h*(1./3)*(H6_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H6_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H6_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f6_Im[k-1] += h*(1./3)*(H6_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H6_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H6_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f7_Re[k-1] += h*(1./3)*(H7_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H7_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H7_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f7_Im[k-1] += h*(1./3)*(H7_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H7_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H7_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f8_Re[k-1] += h*(1./3)*(H8_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H8_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H8_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f8_Im[k-1] += h*(1./3)*(H8_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H8_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H8_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f9_Re[k-1] += h*(1./3)*(H9_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H9_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H9_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f9_Im[k-1] += h*(1./3)*(H9_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H9_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H9_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f10_Re[k-1] += h*(1./3)*(H10_CBP_Re_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H10_CBP_Re_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H10_CBP_Re_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                    f10_Im[k-1] += h*(1./3)*(H10_CBP_Im_0to1[p0]*R_0to1_factor_p0*exp(-R_0to1_p0/(k*h2)) + 4.*H10_CBP_Im_0to1[p1]*R_0to1_factor_p1*exp(-R_0to1_p1/(k*h2)) + H10_CBP_Im_0to1[p2]*R_0to1_factor_p2*exp(-R_0to1_p2/(k*h2)));
                }
                
                for(int j=0;j<=round(N/2)-1;j++)
                {
                    p0 = 2*j;
                    p1 = 2*j+1;
                    p2 = 2*j+2;
                    
                    if(j==0)
                    {
                        R_1toInf_p0 = 0.;
                    }
                    else
                    {
                        R_1toInf_p0 = 1./R_1toInf[p0];
                    }
                    R_1toInf_p1 = 1./R_1toInf[p1];
                    R_1toInf_p2 = 1./R_1toInf[p2];
                    
                    if(j==0)
                    {
                        R_1toInf_factor_p0 = 0.;
                    }
                    else
                    {
                        R_1toInf_factor_p0 = pow(R_1toInf[p0],-2.);
                    }
                    R_1toInf_factor_p1 = pow(R_1toInf[p1],-2.);
                    R_1toInf_factor_p2 = pow(R_1toInf[p2],-2.);
                    
                    f0_Re[k-1] += h*(1./3)*(H0_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H0_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H0_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f0_Im[k-1] += h*(1./3)*(H0_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H0_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H0_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f1_Re[k-1] += h*(1./3)*(H1_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H1_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H1_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f1_Im[k-1] += h*(1./3)*(H1_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H1_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H1_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f2_Re[k-1] += h*(1./3)*(H2_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H2_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H2_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f2_Im[k-1] += h*(1./3)*(H2_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H2_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H2_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f3_Re[k-1] += h*(1./3)*(H3_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H3_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H3_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f3_Im[k-1] += h*(1./3)*(H3_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H3_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H3_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f4_Re[k-1] += h*(1./3)*(H4_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H4_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H4_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f4_Im[k-1] += h*(1./3)*(H4_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H4_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H4_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f5_Re[k-1] += h*(1./3)*(H5_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H5_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H5_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f5_Im[k-1] += h*(1./3)*(H5_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H5_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H5_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f6_Re[k-1] += h*(1./3)*(H6_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H6_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H6_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f6_Im[k-1] += h*(1./3)*(H6_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H6_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H6_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f7_Re[k-1] += h*(1./3)*(H7_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H7_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H7_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f7_Im[k-1] += h*(1./3)*(H7_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H7_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H7_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f8_Re[k-1] += h*(1./3)*(H8_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H8_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H8_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f8_Im[k-1] += h*(1./3)*(H8_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H8_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H8_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f9_Re[k-1] += h*(1./3)*(H9_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H9_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H9_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f9_Im[k-1] += h*(1./3)*(H9_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H9_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H9_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f10_Re[k-1] += h*(1./3)*(H10_CBP_Re_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H10_CBP_Re_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H10_CBP_Re_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));
                    f10_Im[k-1] += h*(1./3)*(H10_CBP_Im_1toInf[p0]*R_1toInf_factor_p0*exp(-R_1toInf_p0/(k*h2)) + 4.*H10_CBP_Im_1toInf[p1]*R_1toInf_factor_p1*exp(-R_1toInf_p1/(k*h2)) + H10_CBP_Im_1toInf[p2]*R_1toInf_factor_p2*exp(-R_1toInf_p2/(k*h2)));  
                }
                
                f0_Re_LCBP[k-1] = (1+f0_Re[k-1]+Res0*H_sing_CBP_Re[k-1])*exp(-0./(k*h2))*pow(k*h2,0.);
                f0_Im_LCBP[k-1] = (f0_Im[k-1]+Res0*H_sing_CBP_Im[k-1])*exp(-0./(k*h2))*pow(k*h2,0.);
                f1_Re_LCBP[k-1] = (1+f1_Re[k-1]+Res1*H_sing_CBP_Re[k-1])*exp(-1./(k*h2))*pow(k*h2,1./2);
                f1_Im_LCBP[k-1] = (f1_Im[k-1]+Res1*H_sing_CBP_Im[k-1])*exp(-1./(k*h2))*pow(k*h2,1./2);
                f2_Re_LCBP[k-1] = ((1./2)+f2_Re[k-1]+Res2*H_sing_CBP_Re[k-1])*exp(-2./(k*h2))*pow(k*h2,2./2);
                f2_Im_LCBP[k-1] = (f2_Im[k-1]+Res2*H_sing_CBP_Im[k-1])*exp(-2./(k*h2))*pow(k*h2,2./2);
                f3_Re_LCBP[k-1] = ((1./4)+f3_Re[k-1]+Res3*H_sing_CBP_Re[k-1])*exp(-3./(k*h2))*pow(k*h2,3./2);
                f3_Im_LCBP[k-1] = (f3_Im[k-1]+Res3*H_sing_CBP_Im[k-1])*exp(-3./(k*h2))*pow(k*h2,3./2);
                f4_Re_LCBP[k-1] = ((1./8)+f4_Re[k-1]+Res4*H_sing_CBP_Re[k-1])*exp(-4./(k*h2))*pow(k*h2,4./2);
                f4_Im_LCBP[k-1] = (f4_Im[k-1]+Res4*H_sing_CBP_Im[k-1])*exp(-4./(k*h2))*pow(k*h2,4./2);
                f5_Re_LCBP[k-1] = ((1./16)+f5_Re[k-1]+Res5*H_sing_CBP_Re[k-1])*exp(-5./(k*h2))*pow(k*h2,5./2);
                f5_Im_LCBP[k-1] = (f5_Im[k-1]+Res5*H_sing_CBP_Im[k-1])*exp(-5./(k*h2))*pow(k*h2,5./2);
                f6_Re_LCBP[k-1] = ((1./32)+f6_Re[k-1]+Res6*H_sing_CBP_Re[k-1])*exp(-6./(k*h2))*pow(k*h2,6./2);
                f6_Im_LCBP[k-1] = (f6_Im[k-1]+Res6*H_sing_CBP_Im[k-1])*exp(-6./(k*h2))*pow(k*h2,6./2);
                f7_Re_LCBP[k-1] = ((1./64)+f7_Re[k-1]+Res7*H_sing_CBP_Re[k-1])*exp(-7./(k*h2))*pow(k*h2,7./2);
                f7_Im_LCBP[k-1] = (f7_Im[k-1]+Res7*H_sing_CBP_Im[k-1])*exp(-7./(k*h2))*pow(k*h2,7./2);
                f8_Re_LCBP[k-1] = ((1./128)+f8_Re[k-1]+Res8*H_sing_CBP_Re[k-1])*exp(-8./(k*h2))*pow(k*h2,8./2);
                f8_Im_LCBP[k-1] = (f8_Im[k-1]+Res8*H_sing_CBP_Im[k-1])*exp(-8./(k*h2))*pow(k*h2,8./2);
                f9_Re_LCBP[k-1] = ((1./256)+f9_Re[k-1]+Res9*H_sing_CBP_Re[k-1])*exp(-9./(k*h2))*pow(k*h2,9./2);
                f9_Im_LCBP[k-1] = (f9_Im[k-1]+Res9*H_sing_CBP_Im[k-1])*exp(-9./(k*h2))*pow(k*h2,9./2);
                f10_Re_LCBP[k-1] = ((1./512)+f10_Re[k-1]+Res10*H_sing_CBP_Re[k-1])*exp(-10./(k*h2))*pow(k*h2,10./2);
                f10_Im_LCBP[k-1] = (f10_Im[k-1]+Res10*H_sing_CBP_Im[k-1])*exp(-10./(k*h2))*pow(k*h2,10./2);
            }
            
            cout<<"f0_Re_LCBP[0] = "<<f0_Re_LCBP[0]<<endl;
            
            string output_file = file_output+num+extension;
            std::ofstream output(output_file.c_str());
            
            if(output.is_open())
            {
                for(int count = 0; count < n; count ++)
                {
                    output <<std::setprecision(32)<<h2*(count+1)<<"\t"<<f0_Re_LCBP[count]<<"\t"<<f0_Im_LCBP[count]<<"\t"<<f1_Re_LCBP[count]<<"\t"<<f1_Im_LCBP[count]<<"\t"<<f2_Re_LCBP[count]<<"\t"<<f2_Im_LCBP[count]<<"\t"<<f3_Re_LCBP[count]<<"\t"<<f3_Im_LCBP[count]<<"\t"<<f4_Re_LCBP[count]<<"\t"<<f4_Im_LCBP[count]<<"\t"<<f5_Re_LCBP[count]<<"\t"<<f5_Im_LCBP[count]<<"\t"<<f6_Re_LCBP[count]<<"\t"<<f6_Im_LCBP[count]<<"\t"<<f7_Re_LCBP[count]<<"\t"<<f7_Im_LCBP[count]<<"\t"<<f8_Re_LCBP[count]<<"\t"<<f8_Im_LCBP[count]<<"\t"<<f9_Re_LCBP[count]<<"\t"<<f9_Im_LCBP[count]<<"\t"<<f10_Re_LCBP[count]<<"\t"<<f10_Im_LCBP[count]<<endl;
                }
                output.close();
            }
            else
            {
                cout<<"can't write into a file"<<endl;
            }
        }
    }
}
