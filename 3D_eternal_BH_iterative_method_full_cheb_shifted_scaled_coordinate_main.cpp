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

// 3d singular GPE solution: Black/White hole

int main()
{
    double bc_l,bc_r,d_bc_r,F_diff,F_diff2,series_left,series_right,flag,ode,ode2;
    double omega = 1.;
    double n = 10.;
    double B;
    double Zmin = -1.;
    double Zmax = 1.;
    double r0;
    int N = 200;
    int N2 = N;
    double d = 3;
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
    
    MatrixXd Diag_function_of_F(N-1,N-1);
    MatrixXd M(N-1,N-1);
    MatrixXd MInv(N-1,N-1);
    MatrixXd D1(N-1,N-1);
    MatrixXd D2(N-1,N-1);
    
    VectorXd R(N-1);
    VectorXd F(N-1);
    VectorXd F_p(N-1);
    VectorXd F_p2(N-1);
    VectorXd F_update(N-1);
    VectorXd function_of_F(N-1);
    VectorXd ODE(N-1);
    VectorXd function_of_F_for_theta(N);
    VectorXd bv_for_fd(N-1); //boundary_val_inclusion_for_finite_diff
    VectorXd F1(N-1);
    
    VectorXd theta(N);
    
    bc_l = 0.;
    bc_r = 1.;
    
    cout<<"N2 = "<<N2<<endl;
    cout<<"N = "<<N<<endl;
    cout<<"Zmax = "<<Zmax<<endl;
    
    //--------------Time loop for evolution of fields----------------//
    for(int k=0;k<N-1;k++)
    {
        R(k) = full_R[N-1-k];
    }
    
    
    for(int k=0;k<N-1;k++)
    {
        for(int j=0;j<N-1;j++)
        {
            D2(k,j) = full_D2[N-1-k][N-1-j];
            D1(k,j) = full_D1[N-1-k][N-1-j];
        }
    }
    
    std::string B_val("B_"), N_val("_N_"), iterative_cheb_shifted_scaled_coordinate_r0("_3D_iterative_full_cheb_shifted_scaled_coordinate_r0_"), file_for_B_r0_val("r0_vs_B"), extension(".txt");
    
    string r0_B_file = file_for_B_r0_val+extension;
    
    std::ifstream r0_B;
    r0_B.open(r0_B_file.c_str());
    
    char num[1000], num1[1000], num2[1000];
    char oldname[1000] = "file.txt";
    char newname[1000];
    
    int result;
    
    int num_files = 0;
    do
    {
        r0_B>>B>>r0;
        
        cout<<"B = "<<B<<", r0 = "<<r0<<endl;
        
        if(B==500)
        {   
            /*if(r0>=0.0002)
            {
                r0 = r0 - 0.0002;
            }*/;
            
            do
            {
                for(int k=0;k<N-1;k++)
                {
                    F(k) = 1.-(0.25*pow(R(k)-1,2));// + 3.*pow(1-pow(R(k),2.),20.);
                    F_update(k) = F(k);
                }
            
                for(int k=0;k<N-1;k++)
                {
                    F_p(k) = bc_l*full_D1[N-k-1][N] + bc_r*full_D1[N-k-1][0];
                    F_p2(k) = bc_l*full_D2[N-k-1][N] + bc_r*full_D2[N-k-1][0];
                    for(int j=0;j<N-1;j++)
                    {
                        F_p(k) += D1(k,j)*F(j);
                        F_p2(k) += D2(k,j)*F(j);
                    }
                }
                
                for(int k=0;k<N-1;k++)
                {
                    function_of_F(k) = (B * B * pow((n - r0) * R(k) + n + r0, 0.3e1) * pow(((-n + r0) * R(k) - n - r0) / (-0.1e1 + R(k)), -(2 * d)) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k))) * pow(-0.1e1 + R(k), -0.2e1) * pow(F(k), 0.5e1) - ((n - r0) * R(k) + n + r0) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k))) * F(k) + ((n - r0) * R(k) + (-d + 2) * n + r0) * pow(-0.1e1 + R(k), 0.3e1) * pow(n, -0.2e1) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k))) * F_p(k) / 0.2e1 + ((n - r0) * R(k) + n + r0) * pow(-0.1e1 + R(k), 0.4e1) * pow(n, -0.2e1) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k))) * F_p2(k) / 0.4e1 + (-((n - r0) * R(k) + n + r0) * pow(-0.1e1 + R(k), 0.4e1) * pow(n, -0.2e1) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k))) * F_p(k) * F_p(k) / 0.2e1 + ((n - r0) * R(k) + n + r0) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k)))) / F(k));
                }
            
                int t=0;
                
                do
                {
                    F_diff2 = F_diff;
                    ode2 = ode;
                    
                    for(int k=0;k<N-1;k++)
                    {
                        F(k)=F_update(k);
                        ODE(k) = function_of_F(k);
                    }
                    
                    for(int k=0;k<N-1;k++)
                    {
                        for(int j=0;j<N-1;j++)
                        {
                            M(k,j) = (pow((-1 + R(k)), 4) * pow(n, -0.2e1) / 0.4e1)*D2(k,j) + (-(((2 * R(k) * R(k) * F_p(k) - R(k) * F(k) + (d - 2) * F(k) - 2 * F_p(k)) * n) - 0.2e1 * ((R(k) * F_p(k)) - F(k) / 0.2e1 - F_p(k)) * (-1 + R(k)) * r0) * pow((-1 + R(k)), 3) / F(k) * pow(n, (-2)) / (((1 + R(k)) * n) - r0 * (-1 + R(k))) / 0.2e1)*D1(k,j);
                            if(k==j)
                            {
                               M(k,j) += (0.5e1 * B * B * pow((n - r0) * R(k) + n + r0, 0.2e1) * pow(-(n * R(k) - r0 * R(k) + n + r0) / (-0.1e1 + R(k)), -(2 * d)) * pow(F(k), 0.4e1) * pow(-0.1e1 + R(k), -0.2e1) - 0.1e1 + (pow(-0.1e1 + R(k), 0.4e1) * F_p(k) * F_p(k) - 0.2e1 * n * n) * pow(F(k), -0.2e1) * pow(n, -0.2e1) / 0.2e1);
                            }
                        }
                    }
                    
                    MInv=M.inverse();
                    
                    F_update = F - (omega*MInv*function_of_F);
                    
                    for(int k=0;k<N-1;k++)
                    {
                        F_p(k) = bc_l*full_D1[N-k-1][N] + bc_r*full_D1[N-k-1][0];
                        F_p2(k) = bc_l*full_D2[N-k-1][N] + bc_r*full_D2[N-k-1][0];
                        for(int j=0;j<N-1;j++)
                        {
                            F_p(k) += D1(k,j)*F_update(j);
                            F_p2(k) += D2(k,j)*F_update(j);
                        }
                    }
                    
                    for(int k=0;k<N-1;k++)
                    {
                        function_of_F(k) = (B * B * pow((n - r0) * R(k) + n + r0, 0.3e1) * pow(((-n + r0) * R(k) - n - r0) / (-0.1e1 + R(k)), -(2 * d)) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k))) * pow(-0.1e1 + R(k), -0.2e1) * pow(F_update(k), 0.5e1) - ((n - r0) * R(k) + n + r0) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k))) * F_update(k) + ((n - r0) * R(k) + (-d + 2) * n + r0) * pow(-0.1e1 + R(k), 0.3e1) * pow(n, -0.2e1) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k))) * F_p(k) / 0.2e1 + ((n - r0) * R(k) + n + r0) * pow(-0.1e1 + R(k), 0.4e1) * pow(n, -0.2e1) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k))) * F_p2(k) / 0.4e1 + (-((n - r0) * R(k) + n + r0) * pow(-0.1e1 + R(k), 0.4e1) * pow(n, -0.2e1) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k))) * F_p(k) * F_p(k) / 0.2e1 + ((n - r0) * R(k) + n + r0) / ((0.1e1 + R(k)) * n - r0 * (-0.1e1 + R(k)))) / F_update(k));
                    }
                    
                    F_diff=pow(F(0)-F_update(0),2.);
                    for(int k=1;k<N-1;k++)
                    {
                        F_diff+=pow(F(k)-F_update(k),2.);
                    }
                    F_diff=sqrt(F_diff);
                    
                    ode=pow(function_of_F(0),2.);
                    for(int k=1;k<N-1;k++)
                    {
                        ode+=pow(function_of_F(k),2.);
                    }
                    ode=sqrt(ode);
                    
                    cout<<t<<" ==> F_diff = "<<F_diff<<", ode = "<<ode<<endl;
                    
                    t=t+1;
                }
                while(t<10 || (ode2>1.e-5 && ode<=ode2) || ode<=ode2 || (F_diff2>1.e-5 && F_diff<=F_diff2) || F_diff<=F_diff2);
                //while(t<1000);
                
                cout<<t<<" ==> F_diff = "<<F_diff<<", ode = "<<ode<<endl;
                
                cout<<std::setprecision(7)<<"r0 = "<<r0<<endl;
                
                r0=r0+0.0001;
            }
            while(ode2>1.e-5 || F_diff2>1.e-5);
            
            r0=r0-0.0001;
            
            for(int k=0;k<=N2;k++)
            {
                del_f[k] = bc_l*full_D1[N-k][N] + bc_r*full_D1[N-k][0];
                del2_f[k] = bc_l*full_D2[N-k][N] + bc_r*full_D2[N-k][0];
                for(int j=0;j<N2-1;j++)
                {
                    del_f[k] += full_D1[N-k][N-j-1]*F(j);
                    del2_f[k] += full_D2[N-k][N-j-1]*F(j);
                }
            }
            
            std::ofstream infile("file.txt");
            
            if(infile.is_open())
            {
                infile <<std::setprecision(32)<<Zmin<<"\t"<<bc_l<<"\t"<<0<<"\t"<<del_f[0]<<"\t"<<del2_f[0]<<endl;
                for(int count = 0; count < N-1; count ++)
                {
                    infile <<std::setprecision(32)<<R(count)<<"\t"<<F(count)<<"\t"<<ODE(count)<<"\t"<<del_f[count+1]<<"\t"<<del2_f[count+1]<<endl;
                }
                infile <<std::setprecision(32)<<Zmax<<"\t"<<bc_r<<"\t"<<0<<"\t"<<del_f[N]<<"\t"<<del2_f[N]<<endl;
                
                infile.close();
            }
            else
            {
                cout<<"can't write into a file"<<endl;
            }
            
            sprintf(num,"%.4f",r0);
            sprintf(num1,"%.2f",B);
            sprintf(num2,"%d",N);
            string r0_val = B_val+num1+N_val+num2+iterative_cheb_shifted_scaled_coordinate_r0+num+extension;
            
            strcpy(newname, r0_val.c_str());
            
            result= rename( oldname , newname );
        }
        num_files += 1;
    }
    while(num_files<=20);
    r0_B.close();
}
