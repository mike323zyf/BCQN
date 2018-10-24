#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>


using namespace Eigen;
using namespace std;

double amips_param;
double c1_param;
double c2_param;
double d1_param;


Vector3d energy_derivative(int type, Vector3d S)
{
    
    Vector3d es;
    double J;
    double l1;
    double l2;
    
    switch(type)
    {
        case 0: //arap
            es[0] = 2 * (S[0] - 1);
            es[1] = 2 * (S[1] - 1);
            es[2] = 2 * (S[2] - 1);
            break;
            
        case 1: // mips
            es[0] = 1.0 / S[1] - S[1] / (S[0] * S[0]);
            es[1] = 1.0 / S[0] - S[0] / (S[1] * S[1]);
            break;
            
        case 2: // iso
            es[0] = 2 * S[0] - 2.0 / (S[0] * S[0] * S[0]);
            es[1] = 2 * S[1] - 2.0 / (S[1] * S[1] * S[1]);
            es[2] = 2 * S[2] - 2.0 / (S[2] * S[2] * S[2]);
            break;
            
        case 3: // amips
            es[0] = amips_param * exp(amips_param * (S[0] / S[1] + S[1] / S[0])) * (1.0 / S[1] - S[1] / (S[0] * S[0]));
            es[1] = amips_param * exp(amips_param * (S[0] / S[1] + S[1] / S[0])) * (1.0 / S[0] - S[0] / (S[1] * S[1]));
            break;
            
        case 4: // conf
            es[0] = 2 * S[0] / (S[1] * S[1]);
            es[1] = -2 * S[0] * S[0] / (S[1] * S[1] * S[1]);
            break;
            
        case 5:
            J = S[0] * S[1] * S[2];
            
            //J = max(J, 1e-8);
            
            l1 = S[0] * S[0] + S[1] * S[1] + S[2] * S[2];
            l2 = S[0] * S[0] * S[1] * S[1] + S[1] * S[1] * S[2] * S[2] + S[2] * S[2] * S[0] * S[0];
            
            //es[0] = 2 * (J - 1) * S[1];
            //es[1] = 2 * (J - 1) * S[0];
            
            es[0] = c1_param * (-2.0 / 3.0 * pow(J, -5.0 / 3.0) * S[1] * S[2] * l1 + pow(J, -2.0 / 3.0) * 2 * S[0]) 
                    + c2_param * (-4.0 / 3.0 * pow(J, -7.0 / 3.0) * S[1] * S[2] * l2 + pow(J, -4.0 / 3.0) * 2 * (S[1] * S[1] + S[2] * S[2]) * S[0]) 
                    + d1_param * 2 * (J - 1) * S[1] * S[2];
            
            es[1] = c1_param * (-2.0 / 3.0 * pow(J, -5.0 / 3.0) * S[2] * S[0] * l1 + pow(J, -2.0 / 3.0) * 2 * S[1]) 
                    + c2_param * (-4.0 / 3.0 * pow(J, -7.0 / 3.0) * S[2] * S[0] * l2 + pow(J, -4.0 / 3.0) * 2 * (S[0] * S[0] + S[2] * S[2]) * S[1]) 
                    + d1_param * 2 * (J - 1) * S[2] * S[0];
            
            es[2] = c1_param * (-2.0 / 3.0 * pow(J, -5.0 / 3.0) * S[0] * S[1] * l1 + pow(J, -2.0 / 3.0) * 2 * S[2]) 
                    + c2_param * (-4.0 / 3.0 * pow(J, -7.0 / 3.0) * S[0] * S[1] * l2 + pow(J, -4.0 / 3.0) * 2 * (S[0] * S[0] + S[1] * S[1]) * S[2]) 
                    + d1_param * 2 * (J - 1) * S[0] * S[1];

            break;
    }
    
    return es;

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *output_mex, *J_value_mex, *JT_value_mex, *wu_mex, *bu_mex;
    
    double *tet_num, *X_g_inv, *tet_vols, *obj_tet, *q_target, *type, *amips_s, *F_dot, *ver_num, *c1_g, *c2_g, *d1_g, *J_index, *JT_index, *JTJ_info, *x2u;
    double *output, *J_value, *JT_value, *wu, *bu;
    
    tet_num = mxGetPr(prhs[0]);
    X_g_inv = mxGetPr(prhs[1]);
    tet_vols = mxGetPr(prhs[2]);
    obj_tet = mxGetPr(prhs[3]);
    q_target = mxGetPr(prhs[4]);
    type = mxGetPr(prhs[5]);
    amips_s = mxGetPr(prhs[6]);
    F_dot = mxGetPr(prhs[7]);
    ver_num = mxGetPr(prhs[8]);
    c1_g = mxGetPr(prhs[9]);
    c2_g = mxGetPr(prhs[10]);
    d1_g = mxGetPr(prhs[11]);
    J_index = mxGetPr(prhs[12]);
    JT_index = mxGetPr(prhs[13]);
    JTJ_info = mxGetPr(prhs[14]);
    x2u = mxGetPr(prhs[15]);
    
    int tet_n = tet_num[0];
    int ver_n = ver_num[0];
    int energy_type = type[0];
    amips_param = amips_s[0];
    c1_param = c1_g[0];
    c2_param = c2_g[0];
    d1_param = d1_g[0];
    int J_rn = JTJ_info[0];
    int J_cn = JTJ_info[1];
    int JT_rn = JTJ_info[2];
    int JT_cn = JTJ_info[3];
    
    output_mex = plhs[0] = mxCreateDoubleMatrix(3 * ver_n, 1, mxREAL);
    J_value_mex = plhs[1] = mxCreateDoubleMatrix(J_rn * J_cn, 1, mxREAL);
    JT_value_mex = plhs[2] = mxCreateDoubleMatrix(JT_rn * JT_cn, 1, mxREAL);
    wu_mex = plhs[3] = mxCreateDoubleMatrix(tet_n, 1, mxREAL);
    bu_mex = plhs[4] = mxCreateDoubleMatrix(tet_n, 1, mxREAL);
    
    output = mxGetPr(output_mex);
    J_value = mxGetPr(J_value_mex);
    JT_value = mxGetPr(JT_value_mex);
    wu = mxGetPr(wu_mex);
    bu = mxGetPr(bu_mex);
    
    memset(J_value, 0, J_rn * J_cn * sizeof(double));
    memset(JT_value, 0, JT_rn * JT_cn * sizeof(double));
    memset(output, 0, 3 * ver_n * sizeof(double));
    
    vector<int> offset;
    
    offset.resize(3 * ver_n, 0);
    
    int tet[4];
    double tet_vol;
    
    Matrix3d B, X_f, A, F_d;
    Vector3d S, u0, u1, u2, v0, v1, v2;
    Matrix3d U, V;
    
    Vector3d es;
    Vector3d cs;
    
    double es1, es2, es3;
    double d1, d2, d3;
    double dphi;
 
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
        }*/
        
        if(A.determinant() <= 0)
        {
            mexPrintf("element inverted\n");
        }
        
        bu[i] = S[0] * S[1] * S[2];
        
        es = energy_derivative(energy_type, S);
        
        cs = S;
        
        es1 = es[0];
        es2 = es[1];
        es3 = es[2];
        
        ////////////////////////////////////////////
        
        wu[i] = 0;
        
        for(int j = 0; j < 12; j++)
        {
            
            F_d(0, 0) = F_dot[i + j * tet_n];
            F_d(1, 0) = F_dot[i + j * tet_n + 12 * tet_n];
            F_d(2, 0) = F_dot[i + j * tet_n + 12 * 2 * tet_n];
            
            F_d(0, 1) = F_dot[i + j * tet_n + 12 * 3 * tet_n];
            F_d(1, 1) = F_dot[i + j * tet_n + 12 * 4 * tet_n];
            F_d(2, 1) = F_dot[i + j * tet_n + 12 * 5 * tet_n];
            
            F_d(0, 2) = F_dot[i + j * tet_n + 12 * 6 * tet_n];
            F_d(1, 2) = F_dot[i + j * tet_n + 12 * 7 * tet_n];
            F_d(2, 2) = F_dot[i + j * tet_n + 12 * 8 * tet_n];
            
            u0 = U.col(0);
            u1 = U.col(1);
            u2 = U.col(2);
            v0 = V.col(0);
            v1 = V.col(1);
            v2 = V.col(2);
            
            d1 = u0.transpose() * F_d * v0;
            d2 = u1.transpose() * F_d * v1;
            d3 = u2.transpose() * F_d * v2;
            
            dphi = d1 * es1 + d2 * es2 + d3 * es3;
            
            output[3 * tet[j / 3] + (j % 3)] += tet_vol * dphi;
            
            if(x2u[3 * tet[j / 3] + (j % 3)] > 0)
            {
                int idx = x2u[3 * tet[j / 3] + (j % 3)] - 1;
                
                J_value[idx * J_cn + offset[3 * tet[j / 3] + (j % 3)]] = d1 * cs[1] * cs[2] + d2 * cs[0] * cs[2] + d3 * cs[0] * cs[1];//d1 * cs[1] + d2 * cs[0];
                
                offset[3 * tet[j / 3] + (j % 3)]++;
                
                JT_value[i * JT_cn + j] = d1 * cs[1] * cs[2] + d2 * cs[0] * cs[2] + d3 * cs[0] * cs[1];//d1 * cs[1] + d2 * cs[0];
                
                wu[i] += JT_value[i * JT_cn + j] * JT_value[i * JT_cn + j];
            }
            
        }
        
        
        
    }
    
    return;
    
}