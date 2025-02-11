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

//----------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------//
               //-------------------solving in x----------------------//
//----------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------//

// Non-singular GPE solution from background funnel or potentials: Black/White Hole 

int main()
{
    double bc_l,bc_r,d_bc_r,F_diff,F_diff2,series_left,series_right,flag,ode,ode2;
    double omega = 1.;
    double n = 10.;
    double B, R0;
    double Zmin = -1.;
    double Zmax = 1.;
    double r0;
    int N = 200;
    int N2 = N;
    double d = 2;
    double pi = 3.141592653589793238;
    
    double full_D1[N2+1][N2+1], full_D2[N2+1][N2+1], full_R[N2+1], c[N2+1], del_f[N2+1], del2_f[N2+1];
    
    for(int k=0;k<=N2;k++)
    {
        full_R[k] = cos(pi*k/N2);
    }
    
    for(int k=0;k<=N2;k++)
    {
        if(k==0 || k==N2)
        {
            c[k] = 2.;
        }
        else
        {
            c[k] = 1.;
        }
    }
    
    for(int k=0;k<=N2;k++)
    {
        for(int j=0;j<=N2;j++)
        {
            if(k==j)
            {
                if(k==0)
                {
                    full_D1[k][j] = (2.*pow(N2,2)+1)*pow(6.,-1);
                }
                else if(k==N2)
                {
                    full_D1[k][j] = -(2.*pow(N2,2)+1)*pow(6.,-1);
                }
                else
                {
                    full_D1[k][j] = -0.5*full_R[k]*pow(1-pow(full_R[k],2),-1);
                }
            }
            else
            {
                if((k+j+2) % 2 == 0)
                {
                    flag=1.;
                }
                else
                {
                    flag=-1.;
                }
                full_D1[k][j] = c[k]*flag*pow(c[j]*(full_R[k]-full_R[j]),-1.);
            }
            
            //cout<<"D["<<k<<"]["<<j<<"] = "<<D[k][j]<<"\t";
        }
        //cout<<endl;
    }
    
    for(int k=0;k<N2+1;k++)
    {
        for(int j=0;j<N2+1;j++)
        {
            full_D2[k][j]=0.;
            for(int i=0;i<N2+1;i++)
            {
                full_D2[k][j] += full_D1[k][i]*full_D1[i][j];
            }
        }
    }
    
    MatrixXd Diag_function_of_F(N,N);
    MatrixXd M(N,N);
    MatrixXd MInv(N,N);
    MatrixXd D1(N,N);
    MatrixXd D2(N,N);
    
    VectorXd R(N);
    VectorXd F(N);
    VectorXd F_p(N);
    VectorXd F_p2(N);
    VectorXd F_update(N);
    VectorXd function_of_F(N);
    VectorXd ODE(N);
    VectorXd function_of_F_for_theta(N);
    VectorXd bv_for_fd(N); //boundary_val_inclusion_for_finite_diff
    VectorXd F1(N);
    
    VectorXd f(N);
    VectorXd h(N);
    VectorXd f_prime(N);
    VectorXd h_prime(N);
    
    VectorXd g(N);
    VectorXd V(N);
    VectorXd MetricRR(N);
    VectorXd DetMetric(N);
    VectorXd del_DetMetric(N);
    
    VectorXd theta(N);
    
    
    bc_r = 1.;
    
    cout<<"N2 = "<<N2<<endl;
    cout<<"N = "<<N<<endl;
    cout<<"Zmax = "<<Zmax<<endl;
    
    //--------------Time loop for evolution of fields----------------//
    for(int k=0;k<N;k++)
    {
        R(k) = full_R[N-k];
    }
    
    
    for(int k=0;k<N;k++)
    {
        for(int j=0;j<N;j++)
        {
            D2(k,j) = full_D2[N-k][N-j];
            D1(k,j) = full_D1[N-k][N-j];
        }
    }
    
    std::string B_val("B_"), bcl("_bcl_"), N_val("_N_"), iterative_cheb_shifted_scaled_coordinate_sym_R0("_2D_iterative_full_cheb_shifted_scaled_coordinate_sym_R0_"), iterative_cheb_shifted_scaled_coordinate_asym_R0("_2D_iterative_full_cheb_shifted_scaled_coordinate_asym_R0_"), iterative_cheb_shifted_scaled_coordinate_R0("_2D_iterative_full_cheb_shifted_scaled_coordinate_R0_"), file_for_B_r0_val("r0_vs_B"), extension(".txt");
    
    string r0_B_file = file_for_B_r0_val+extension;
    
    std::ifstream r0_B;
    r0_B.open(r0_B_file.c_str());
    
    char num[1000], num1[1000], num2[1000], num3[1000];
    char oldname[1000] = "file.txt";
    char newname[1000];
    
    int result;
    
    int num_files = 0;
    do
    {
        r0_B>>B>>r0;
        
        cout<<"B = "<<B<<", r0 = "<<r0<<endl;
        
        if(B==1.)
        {   
            B = 1.;
            R0 = 1.4827;
            /*if(r0>=0.0002)
            {
                r0 = r0 - 0.0002;
            }*/
            
            bc_l = 1.;
            r0 = 0.;
            
            do
            {
                n = R0;
                
                for(int k=0;k<N;k++)
                {
                    F(k) = 1.;//(bc_r * (1 + R(k))) / 0.2e1 + (bc_l * (1 - R(k))) / 0.2e1;//-0.19e2 / 0.2e1 * R(k) * R(k) + R(k) / 0.2e1 + 0.10e2;//0.1e1 - pow(R(k) / 0.2e1 - 0.1e1 / 0.2e1, 0.10e2) + 0.375e3 / 0.64e2 * pow(R(k) * R(k) - 0.1e1, 0.10e2) * pow(R(k) - 0.6e0, 0.20e2);// + 3.*pow(1-pow(R(k),2.),20.);
                    F_update(k) = F(k);
                    
                    g(k) = 1.;
                    V(k) = 0.;
                    
                    h(k) = 0.1e1 + pow(R0, 0.4e1) * pow((r0 + n * (1 + R(k)) / (1 - R(k))), (-4));
                    h_prime(k) = 0.8e1 * n * pow((-1 + R(k)), 3) * pow(R0, 0.4e1) * pow((n * R(k) - r0 * R(k) + n + r0), (-5));
                    
                    f(k) = 1.;
                    f_prime(k) = 0.;
                    
                    /*f(k) = h(k);
                    f_prime(k) = h_prime(k);*/
                }
            
                for(int k=0;k<N;k++)
                {
                    F_p(k) = bc_r*full_D1[N-k][0];
                    F_p2(k) = bc_r*full_D2[N-k][0];
                    for(int j=0;j<N;j++)
                    {
                        F_p(k) += D1(k,j)*F(j);
                        F_p2(k) += D2(k,j)*F(j);
                    }
                }
                
                function_of_F(0) = F_p(0);
                for(int k=1;k<N;k++)
                {
                    function_of_F(k) = ((-0.8e1 * (pow((n - r0) * R(k) + n + r0, 0.2e1) * (g(k) * F(k) * F(k) + V(k) - 1) * pow(F(k), 4) * h(k) + B * B * pow(-0.1e1 + R(k), 0.2e1)) * n * n * pow(f(k), 0.2e1) + 0.2e1 * ((n - r0) * R(k) + n + r0) * ((F_p2(k) * (n - r0) * R(k) * R(k) + (0.2e1 * n * F_p(k) - 0.2e1 * r0 * (F_p(k) - F_p2(k))) * R(k) - n * F_p2(k) + 0.2e1 * r0 * (F_p(k) - F_p2(k) / 0.2e1)) * h(k) + F_p(k) * h_prime(k) * ((n - r0) * R(k) + n + r0) * (-0.1e1 + R(k)) / 0.2e1) * pow(-0.1e1 + R(k), 0.3e1) * pow(F(k), 3) * f(k) - f_prime(k) * F_p(k) * pow((n - r0) * R(k) + n + r0, 0.2e1) * h(k) * pow(-0.1e1 + R(k), 0.4e1) * pow(F(k), 3)) * pow(n, -0.2e1) * pow(f(k), -0.2e1) * pow((n - r0) * R(k) + n + r0, -0.2e1) / h(k) * pow(F(k), (-3)) / 0.8e1);
                }
            
                int t=0;
                
                do
                {
                    /*if(t<3)
                    {
                        F_diff2 = F_diff;
                        ode2 = ode;
                        
                        for(int k=0;k<N;k++)
                        {
                            F(k)=F_update(k);
                            ODE(k) = function_of_F(k);
                        }
                    }
                    else if(t>=3 && ode<ode2)
                    {
                        F_diff2 = F_diff;
                        ode2 = ode;
                        
                        for(int k=0;k<N;k++)
                        {
                            F(k)=F_update(k);
                            ODE(k) = function_of_F(k);
                        }
                    }
                    else if(t>=3)
                    {
                        omega = omega/2;
                    }*/
                    
                    F_diff2 = F_diff;
                    ode2 = ode;
                    
                    for(int k=0;k<N;k++)
                    {
                        F(k)=F_update(k);
                        ODE(k) = function_of_F(k);
                    }
                    
                    for(int k=1;k<N;k++)
                    {
                        for(int j=0;j<N;j++)
                        {
                            M(k,j) = (pow((-1 + R(k)), 4) / f(k) * pow(n, -0.2e1) / 0.4e1)*D2(k,j) + ((pow((-1 + R(k)), 3) * ((((4 * n - 4 * r0) * R(k) + 4 * r0) * h(k) + ((n - r0) * R(k) + n + r0) * (-1 + R(k)) * h_prime(k)) * f(k) - ((n - r0) * R(k) + n + r0) * h(k) * (-1 + R(k)) * f_prime(k)) / h(k) / ((n - r0) * R(k) + n + r0) * pow(n, (-2)) * pow(f(k), (-2))) / 0.8e1)*D1(k,j);
                            if(k==j)
                            {
                               M(k,j) += (-3 * g(k) * F(k) * F(k) - V(k) + 1 + 3 * B * B * pow((-1 + R(k)), 2) * pow(((n - r0) * R(k) + n + r0), (-2)) / h(k) * pow(F(k), (-4)));
                            }
                        }
                    }
                    
                    for(int j=0;j<N;j++)
                    {
                        M(0,j) = D1(0,j);
                    }
                    
                    MInv=M.inverse();
                    
                    F_update = F - (omega*MInv*function_of_F);
                    
                    for(int k=0;k<N;k++)
                    {
                        F_p(k) = bc_r*full_D1[N-k][0];
                        F_p2(k) = bc_r*full_D2[N-k][0];
                        for(int j=0;j<N;j++)
                        {
                            F_p(k) += D1(k,j)*F_update(j);
                            F_p2(k) += D2(k,j)*F_update(j);
                        }
                    }
                    
                    function_of_F(0) = F_p(0);
                    for(int k=1;k<N;k++)
                    {
                        function_of_F(k) = ((-0.8e1 * (pow((n - r0) * R(k) + n + r0, 0.2e1) * (g(k) * F_update(k) * F_update(k) + V(k) - 1) * pow(F_update(k), 4) * h(k) + B * B * pow(-0.1e1 + R(k), 0.2e1)) * n * n * pow(f(k), 0.2e1) + 0.2e1 * ((n - r0) * R(k) + n + r0) * ((F_p2(k) * (n - r0) * R(k) * R(k) + (0.2e1 * n * F_p(k) - 0.2e1 * r0 * (F_p(k) - F_p2(k))) * R(k) - n * F_p2(k) + 0.2e1 * r0 * (F_p(k) - F_p2(k) / 0.2e1)) * h(k) + F_p(k) * h_prime(k) * ((n - r0) * R(k) + n + r0) * (-0.1e1 + R(k)) / 0.2e1) * pow(-0.1e1 + R(k), 0.3e1) * pow(F_update(k), 3) * f(k) - f_prime(k) * F_p(k) * pow((n - r0) * R(k) + n + r0, 0.2e1) * h(k) * pow(-0.1e1 + R(k), 0.4e1) * pow(F_update(k), 3)) * pow(n, -0.2e1) * pow(f(k), -0.2e1) * pow((n - r0) * R(k) + n + r0, -0.2e1) / h(k) * pow(F_update(k), (-3)) / 0.8e1);
                    }
                    
                    F_diff=pow(F(0)-F_update(0),2.);
                    for(int k=1;k<N;k++)
                    {
                        F_diff+=pow(F(k)-F_update(k),2.);
                    }
                    F_diff=sqrt(F_diff);
                    
                    ode=pow(function_of_F(0),2.);
                    for(int k=1;k<N;k++)
                    {
                        ode+=pow(function_of_F(k),2.);
                    }
                    ode=sqrt(ode);
                    
                    cout<<t<<" ==> F_diff = "<<F_diff<<", ode = "<<ode<<", omega = "<<omega<<endl;
                    
                    t=t+1;
                }
                //while(t<10 || ode2>1.e-9 || F_diff2>1.e-9);
                while(t<10 || (ode2>1.e-6 && ode<ode2) || ode<ode2 || (F_diff2>1.e-6 && F_diff<F_diff2) || F_diff<F_diff2);
                //while(t<1000);
                
                cout<<t<<" ==> F_diff = "<<F_diff<<", ode = "<<ode<<endl;
                
                cout<<std::setprecision(7)<<"bc_l = "<<bc_l<<", F(0) = "<<F(0)<<", F'(0) = "<<F_p(0)<<", R0 = "<<R0<<endl;
                
                R0=R0+0.0001;
            }
            while(ode2>1.e-7 || F_diff2>1.e-7);
            
            R0=R0-0.0001;
            
            for(int k=0;k<=N2;k++)
            {
                del_f[k] = bc_r*full_D1[N-k][0];
                del2_f[k] = bc_r*full_D2[N-k][0];
                for(int j=0;j<N;j++)
                {
                    del_f[k] += full_D1[N-k][N-j]*F(j);
                    del2_f[k] += full_D2[N-k][N-j]*F(j);
                }
            }
            
            std::ofstream infile("file.txt");
            
            if(infile.is_open())
            {
                for(int count = 0; count < N; count ++)
                {
                    infile <<std::setprecision(32)<<R(count)<<"\t"<<F(count)<<"\t"<<ODE(count)<<"\t"<<del_f[count]<<"\t"<<del2_f[count]<<endl;
                }
                infile <<std::setprecision(32)<<Zmax<<"\t"<<bc_r<<"\t"<<0<<"\t"<<del_f[N]<<"\t"<<del2_f[N]<<endl;
                
                infile.close();
            }
            else
            {
                cout<<"can't write into a file"<<endl;
            }
            
            sprintf(num,"%.4f",R0);
            sprintf(num1,"%.1f",B);
            sprintf(num2,"%d",N);
            sprintf(num3,"%.4f",F(0));
            string r0_val = B_val+num1+N_val+num2+iterative_cheb_shifted_scaled_coordinate_asym_R0+num+bcl+num3+extension;
            //string r0_val = B_val+num1+N_val+num2+iterative_cheb_shifted_scaled_coordinate_sym_R0+num+bcl+num3+extension;
            
            strcpy(newname, r0_val.c_str());
            
            result= rename( oldname , newname );
        }
        num_files += 1;
    }
    while(num_files<22);
    r0_B.close();
}
