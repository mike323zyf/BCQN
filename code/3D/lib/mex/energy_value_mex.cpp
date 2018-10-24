#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <iostream>


using namespace Eigen;
using namespace std;


double amips_param;
double c1_param;
double c2_param;
double d1_param;


double energy_value(int type, Vector3d S)
{
    
    double value = 0;
    double J;
    double l1;
    double l2;
    
    switch(type)
    {
        case 0: //arap
            value = (S[0] - 1) * (S[0] - 1) + (S[1] - 1) * (S[1] - 1) + (S[2] - 1) * (S[2] - 1);
            break;
            
        case 1: // mips
            value = S[0] / S[1] + S[1] / S[0];
            break;
            
        case 2: // iso
            value = S[0] * S[0] + 1.0 / (S[0] * S[0]) + S[1] * S[1] + 1.0 / (S[1] * S[1]) + S[2] * S[2] + 1.0 / (S[2] * S[2]);
            break;
            
        case 3: // amips
            value = exp(amips_param * (S[0] / S[1] + S[1] / S[0]));
            break;
            
        case 4: // conf
            value = S[0] / S[1];
            value *= value;
            break;
            
        case 5:
            J = S[0] * S[1] * S[2];
            l1 = S[0] * S[0] + S[1] * S[1] + S[2] * S[2];
            l2 = S[0] * S[0] * S[1] * S[1] + S[1] * S[1] * S[2] * S[2] + S[2] * S[2] * S[0] * S[0];
            //value = (J - 1) * (J - 1);
            value = c1_param * (pow(J, -2.0 / 3.0) * l1 - 3) + c2_param * (pow(J, -4.0 / 3.0) * l2 - 3) + d1_param * (J - 1) * (J - 1);
            break;
    }
    
    return value;

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *output_mex;

    double *tet_num, *X_g_inv, *tet_vols, *obj_tet, *q_target, *type, *amips_s, *c1, *c2, *d1;
    double *output;
    
    tet_num = mxGetPr(prhs[0]);
    X_g_inv = mxGetPr(prhs[1]);
    tet_vols = mxGetPr(prhs[2]);
    obj_tet = mxGetPr(prhs[3]);
    q_target = mxGetPr(prhs[4]);
    type = mxGetPr(prhs[5]);
    amips_s = mxGetPr(prhs[6]);
    c1 = mxGetPr(prhs[7]);
    c2 = mxGetPr(prhs[8]);
    d1 = mxGetPr(prhs[9]);
    
    output_mex = plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    output = mxGetPr(output_mex);
    
    output[0] = 0;
    
    int tet_n = tet_num[0];
    int energy_type = type[0];
    amips_param = amips_s[0];
    c1_param = c1[0];
    c2_param = c2[0];
    d1_param = d1[0];
    
    int tet[4];
    double tet_vol;
    
    Matrix3d B, X_f, A;
    Vector3d S;
    Matrix3d U, V;
 
    for(int i = 0; i < tet_n; i++)
    {
        tet[0] = obj_tet[i] - 1; 
        tet[1] = obj_tet[i + tet_n] - 1; 
        tet[2] = obj_tet[i + 2 * tet_n] - 1; 
        tet[3] = obj_tet[i + 3 * tet_n] - 1; 
        
        tet_vol = tet_vols[i];
        
        B(0, 0) = X_g_inv[i];
        B(1, 0) = X_g_inv[i + tet_n];
        B(2, 0) = X_g_inv[i + tet_n * 2]; 
        B(0, 1) = X_g_inv[i + tet_n * 3];
        B(1, 1) = X_g_inv[i + tet_n * 4];
        B(2, 1) = X_g_inv[i + tet_n * 5];
        B(0, 2) = X_g_inv[i + tet_n * 6];
        B(1, 2) = X_g_inv[i + tet_n * 7];
        B(2, 2) = X_g_inv[i + tet_n * 8];
        
        X_f(0, 0) = q_target[3 * tet[1]] - q_target[3 * tet[0]];
        X_f(1, 0) = q_target[3 * tet[1] + 1] - q_target[3 * tet[0] + 1]; 
        X_f(2, 0) = q_target[3 * tet[1] + 2] - q_target[3 * tet[0] + 2]; 
        X_f(0, 1) = q_target[3 * tet[2]] - q_target[3 * tet[0]];
        X_f(1, 1) = q_target[3 * tet[2] + 1] - q_target[3 * tet[0] + 1]; 
        X_f(2, 1) = q_target[3 * tet[2] + 2] - q_target[3 * tet[0] + 2];         
        X_f(0, 2) = q_target[3 * tet[3]] - q_target[3 * tet[0]];
        X_f(1, 2) = q_target[3 * tet[3] + 1] - q_target[3 * tet[0] + 1]; 
        X_f(2, 2) = q_target[3 * tet[3] + 2] - q_target[3 * tet[0] + 2]; 
        
        A = X_f * B;
        
        JacobiSVD<Matrix3d> svd(A, ComputeFullU | ComputeFullV);
        
        S = svd.singularValues();
        U = svd.matrixU();
        V = svd.matrixV();
        
        /*if(U.determinant() <= 0)
        {
            U(0, 2) *= -1;
            U(1, 2) *= -1;
            U(2, 2) *= -1;
        }
        
        if(V.determinant() <= 0)
        {
            V(0, 2) *= -1;
            V(1, 2) *= -1;
            V(2, 2) *= -1;
        }
        
        if(A.determinant() <= 0)
        {
            //mexPrintf("inverted\n");
            S[2] *= -1;
        }*/
        
        
        double deviation = energy_value(energy_type, S);
        
        if(A.determinant() <= 0 && energy_type > 0)
        {
            output[0] = 1e20;
            break;
        }
        
        output[0] += tet_vol * deviation;
        
    }
    
    return;
    
}