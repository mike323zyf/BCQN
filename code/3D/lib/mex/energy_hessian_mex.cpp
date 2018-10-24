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
            
            J = max(J, 1e-8);
            
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


double energy_hessian(double s1, double s2, double s3, double s1j, double s1k, double s2j, double s2k, double s3j, double s3k, double s1jk, double s2jk, double s3jk, int type)
{
    
    double hessian;
    double v1, v2, v3;
     
    switch(type)
    {
        case 0: //arap
            v1 = 2 * ((s1 - 1) * s1jk + (s2 - 1) * s2jk + (s3 - 1) * s3jk);
            v2 = 2 * (s1j * s1k + s2j * s2k + s3j * s3k);
            hessian = v1 + v2;
            break;
            
        case 1: // mips
            
            break;
            
        case 2: // iso
            
            v1 = 2 * s1k * s1j + 2 * s1 * s1jk;
            v2 = 2 * s2k * s2j + 2 * s2 * s2jk;
            v3 = 2 * s3k * s3j + 2 * s3 * s3jk;
            
            v1 += 6 * s1k * s1j / pow(s1, 4) - 2 * s1jk / pow(s1, 3);
            v2 += 6 * s2k * s2j / pow(s2, 4) - 2 * s2jk / pow(s2, 3);
            v3 += 6 * s3k * s3j / pow(s3, 4) - 2 * s3jk / pow(s3, 3);
            
            hessian = v1 + v2 + v3;
            break;
            
        case 3: // amips
            
            break;
            
        case 4: // conf
            
            break;
            
        case 5: // gmr
            
            double J = s1 * s2 * s3;
            double l1 = s1 * s1 + s2 * s2 + s3 * s3;
            double l2 = s1 * s1 * s2 * s2 + s2 * s2 * s3 * s3 + s3 * s3 * s1 * s1;
            
            double Jj = s1j * s2 * s3 + s1 * s2j * s3 + s1 * s2 * s3j;
            double Jk = s1k * s2 * s3 + s1 * s2k * s3 + s1 * s2 * s3k;
            
            double l1j = 2 * (s1 * s1j + s2 * s2j + s3 * s3j);
            double l1k = 2 * (s1 * s1k + s2 * s2k + s3 * s3k);
            
            double l2j = 2 * (s1 * s1j * s2 * s2 + s2 * s2j * s1 * s1 + s2 * s2j * s3 * s3 + s3 * s3j * s2 * s2 + s3 * s3j * s1 * s1 + s1 * s1j * s3 * s3);
            double l2k = 2 * (s1 * s1k * s2 * s2 + s2 * s2k * s1 * s1 + s2 * s2k * s3 * s3 + s3 * s3k * s2 * s2 + s3 * s3k * s1 * s1 + s1 * s1k * s3 * s3);
            
            double Jjk = (s1jk * s2 * s3 + s1j * s2k * s3 + s1j * s2 * s3k) + (s1k * s2j * s3 + s1 * s2jk * s3 + s1 * s2j * s3k) + (s1k * s2 * s3j + s1 * s2k * s3j + s1 * s2 * s3jk);
            
            double l1jk = 2 * (s1k * s1j + s1 * s1jk + s2k * s2j + s2 * s2jk + s3k * s3j + s3 * s3jk);
            
            double l2jk = s1k * s1j * s2 * s2 + s1 * s1jk * s2 * s2 + 2 * s1 * s1j * s2 * s2k;
            l2jk += s2k * s2j * s1 * s1 + s2 * s2jk * s1 * s1 + 2 * s2 * s2j * s1 * s1k;
            l2jk += s2k * s2j * s3 * s3 + s2 * s2jk * s3 * s3 + 2 * s2 * s2j * s3 * s3k;
            l2jk += s3k * s3j * s2 * s2 + s3 * s3jk * s2 * s2 + 2 * s3 * s3j * s2 * s2k;
            l2jk += s3k * s3j * s1 * s1 + s3 * s3jk * s1 * s1 + 2 * s3 * s3j * s1 * s1k;
            l2jk += s1k * s1j * s3 * s3 + s1 * s1jk * s3 * s3 + 2 * s1 * s1j * s3 * s3k;
            l2jk *= 2;
                         
            hessian = c1_param * (-2.0 / 3.0 * pow(J, -5.0 / 3.0) * Jk * l1j + pow(J, -2.0 / 3.0) * l1jk + 10.0 / 9.0 * pow(J, -8.0 / 3.0) * Jk * Jj * l1 - 2.0 / 3.0 * pow(J, -5.0 / 3.0) * Jjk * l1 - 2.0 / 3.0 * pow(J, -5.0 / 3.0) * Jj * l1k);              
            hessian += c2_param * (-4.0 / 3.0 * pow(J, -7.0 / 3.0) * Jk * l2j + pow(J, -4.0 / 3.0) * l2jk + 28.0 / 9.0 * pow(J, -10.0 / 3.0) * Jk * Jj * l2 - 4.0 / 3.0 * pow(J, -7.0 / 3.0) * Jjk * l2 - 4.0 / 3.0 * pow(J, -7.0 / 3.0) * Jj * l2k);              
            hessian += 2 * d1_param * (Jk * Jj + J * Jjk - Jjk);              
            
            break;
    }
    
    return hessian;

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *row_mex, *col_mex, *val_mex, *grad_mex, *J_value_mex, *JT_value_mex, *wu_mex, *bu_mex;
    
    double *tet_num, *X_g_inv, *tet_vols, *obj_tet, *q_target, *type, *amips_s, *F_dot, *ver_num, *c1_g, *c2_g, *d1_g, *J_index, *JT_index, *JTJ_info, *x2u, *clamp;
    double *row, *col, *val, *grad, *J_value, *JT_value, *wu, *bu;
    
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
    clamp = mxGetPr(prhs[16]);
    
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
    double proj_threshold = clamp[0];
    
    grad_mex = plhs[0] = mxCreateDoubleMatrix(3 * ver_n, 1, mxREAL);
    J_value_mex = plhs[1] = mxCreateDoubleMatrix(J_rn * J_cn, 1, mxREAL);
    JT_value_mex = plhs[2] = mxCreateDoubleMatrix(JT_rn * JT_cn, 1, mxREAL);
    wu_mex = plhs[3] = mxCreateDoubleMatrix(tet_n, 1, mxREAL);
    bu_mex = plhs[4] = mxCreateDoubleMatrix(tet_n, 1, mxREAL);
    row_mex = plhs[5] = mxCreateDoubleMatrix(144 * tet_n, 1, mxREAL);
    col_mex = plhs[6] = mxCreateDoubleMatrix(144 * tet_n, 1, mxREAL);
    val_mex = plhs[7] = mxCreateDoubleMatrix(144 * tet_n, 1, mxREAL);
    
    grad = mxGetPr(grad_mex);
    J_value = mxGetPr(J_value_mex);
    JT_value = mxGetPr(JT_value_mex);
    wu = mxGetPr(wu_mex);
    bu = mxGetPr(bu_mex);
    row = mxGetPr(row_mex);
    col = mxGetPr(col_mex);
    val = mxGetPr(val_mex);
    
    memset(J_value, 0, J_rn * J_cn * sizeof(double));
    memset(JT_value, 0, JT_rn * JT_cn * sizeof(double));
    memset(grad, 0, 3 * ver_n * sizeof(double));
    
    vector<int>* offset = new vector<int>;
    
    offset->resize(3 * ver_n, 0);
    
    int tet[4];
    double tet_vol;
    int index[12];
    
    Matrix3d B, X_f, A, F_d, T, T_ddot;
    Vector3d S, u0, u1, u2, v0, v1, v2;
    Matrix3d U, V;
    Vector2d temp_w;
    
    MatrixXd ELS(12, 12);
    VectorXd ELS_S;
    MatrixXd ELS_X;
    DiagonalMatrix<double, 12> ELS_D;
    
    
    Matrix3d T_dot[12];
    Matrix3d w_U_dot[12];
    Matrix3d w_V_dot[12];
    
    
    Matrix2d temp_A, inv_A, temp_B, inv_B, temp_C, inv_C;
    
    Vector3d es;
    Vector3d cs;
    
    double es1, es2, es3;
    double d1, d2, d3;
    double dphi;
    
    Matrix2d II;
    II(0, 0) = 1;
    II(1, 1) = 1;
    II(0, 1) = 0;
    II(1, 0) = 0;
    
    double tao = 0.01;
    
    double value;
 
    for(int i = 0; i < tet_n; i++)
    {
        tet[0] = obj_tet[i] - 1; 
        tet[1] = obj_tet[i + tet_n] - 1; 
        tet[2] = obj_tet[i + 2 * tet_n] - 1; 
        tet[3] = obj_tet[i + 3 * tet_n] - 1; 
        
        index[0] = 3 * tet[0];
        index[1] = index[0] + 1;
        index[2] = index[0] + 2;
        
        index[3] = 3 * tet[1];
        index[4] = index[3] + 1;
        index[5] = index[3] + 2;
        
        index[6] = 3 * tet[2];
        index[7] = index[6] + 1;
        index[8] = index[6] + 2;
        
        index[9] = 3 * tet[3];
        index[10] = index[9] + 1;
        index[11] = index[9] + 2;
        
        
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
        
        bu[i] = S[0] * S[1] * S[2];
        
        es = energy_derivative(energy_type, S);
        
        cs = S;
        
        es1 = es[0];
        es2 = es[1];
        es3 = es[2];
        
        ////////////////////////////////////////////
        
        T(0, 0) = S[0];
        T(0, 1) = 0;
        T(0, 2) = 0;
        T(1, 0) = 0;
        T(1, 1) = S[1];
        T(1, 2) = 0;
        T(2, 0) = 0;
        T(2, 1) = 0;
        T(2, 2) = S[2];
        
        temp_A(0, 0) = S[1];
        temp_A(0, 1) = S[0];
        temp_A(1, 0) = S[0];
        temp_A(1, 1) = S[1];
        
        temp_B(0, 0) = S[2];
        temp_B(0, 1) = S[1];
        temp_B(1, 0) = S[1];
        temp_B(1, 1) = S[2];
        
        temp_C(0, 0) = S[2];
        temp_C(0, 1) = S[0];
        temp_C(1, 0) = S[0];
        temp_C(1, 1) = S[2];
        
        
        u0 = U.col(0);
        u1 = U.col(1);
        u2 = U.col(2);
        v0 = V.col(0);
        v1 = V.col(1);
        v2 = V.col(2);
        
        if(abs(S[0] - S[1]) < 1e-5)
        {
            inv_A = temp_A.transpose() * temp_A + tao * II;
            inv_A = inv_A.inverse() * temp_A.transpose();
            //mexPrintf("ill conditioned\n");
        }else{
            inv_A = temp_A.inverse();
        }
        
        
        if(abs(S[1] - S[2]) < 1e-5)
        {
            inv_B = temp_B.transpose() * temp_B + tao * II;
            inv_B = inv_B.inverse() * temp_B.transpose();
            //mexPrintf("ill conditioned\n");
        }else{
            inv_B = temp_B.inverse();
        }
        
        
        if(abs(S[2] - S[0]) < 1e-5)
        {
            inv_C = temp_C.transpose() * temp_C + tao * II;
            inv_C = inv_C.inverse() * temp_C.transpose();
            //mexPrintf("ill conditioned\n");
        }else{
            inv_C = temp_C.inverse();
        }
        
        
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
            
            
            d1 = u0.transpose() * F_d * v0;
            d2 = u1.transpose() * F_d * v1;
            d3 = u2.transpose() * F_d * v2;
            
            dphi = d1 * es1 + d2 * es2 + d3 * es3;
            
            grad[3 * tet[j / 3] + (j % 3)] += tet_vol * dphi;
            
            if(x2u[3 * tet[j / 3] + (j % 3)] > 0)
            {
                int idx = x2u[3 * tet[j / 3] + (j % 3)] - 1;
                
                J_value[idx * J_cn + (*offset)[3 * tet[j / 3] + (j % 3)]] = d1 * cs[1] * cs[2] + d2 * cs[0] * cs[2] + d3 * cs[0] * cs[1];//d1 * cs[1] + d2 * cs[0];
                
                (*offset)[3 * tet[j / 3] + (j % 3)]++;
                
                JT_value[i * JT_cn + j] = d1 * cs[1] * cs[2] + d2 * cs[0] * cs[2] + d3 * cs[0] * cs[1];//d1 * cs[1] + d2 * cs[0];
                
                wu[i] += JT_value[i * JT_cn + j] * JT_value[i * JT_cn + j];
            }
            
            
            T_dot[j](0, 0) = d1;
            T_dot[j](0, 1) = 0;
            T_dot[j](0, 2) = 0;
            T_dot[j](1, 0) = 0;
            T_dot[j](1, 1) = d2;
            T_dot[j](1, 2) = 0;
            T_dot[j](2, 0) = 0;
            T_dot[j](2, 1) = 0;
            T_dot[j](2, 2) = d3;
            
            
            w_U_dot[j](0, 0) = 0;
            w_U_dot[j](1, 1) = 0;
            w_U_dot[j](2, 2) = 0;
            
            w_V_dot[j](0, 0) = 0;
            w_V_dot[j](1, 1) = 0; 
            w_V_dot[j](2, 2) = 0; 
            
            
            temp_w[0] = u0.transpose() * F_d * v1;
            temp_w[1] = -1.0 * u1.transpose() * F_d * v0;
            temp_w = inv_A * temp_w;
          
            w_U_dot[j](0, 1) = temp_w[0];
            w_U_dot[j](1, 0) = -1 * temp_w[0];
      
            w_V_dot[j](0, 1) = temp_w[1];
            w_V_dot[j](1, 0) = -1 * temp_w[1];
            
            
            temp_w[0] = u1.transpose() * F_d * v2;
            temp_w[1] = -1.0 * u2.transpose() * F_d * v1;
            temp_w = inv_B * temp_w;
            
            w_U_dot[j](1, 2) = temp_w[0];
            w_U_dot[j](2, 1) = -1 * temp_w[0];
            
            w_V_dot[j](1, 2) = temp_w[1];
            w_V_dot[j](2, 1) = -1 * temp_w[1];
            
            
            temp_w[0] = u0.transpose() * F_d * v2;
            temp_w[1] = -1.0 * u2.transpose() * F_d * v0;
            temp_w = inv_C * temp_w;
            
            w_U_dot[j](0, 2) = temp_w[0];
            w_U_dot[j](2, 0) = -1 * temp_w[0];
            
            w_V_dot[j](0, 2) = temp_w[1];
            w_V_dot[j](2, 0) = -1 * temp_w[1];      
        }
        
        
        
        for(int j = 0; j < 12; j++)
        {
            for(int k = j; k < 12; k++)
            {
                T_ddot = w_U_dot[k].transpose() * w_U_dot[j] * T + T * w_V_dot[j] * w_V_dot[k].transpose() - w_U_dot[j] * T * w_V_dot[k] - w_U_dot[k] * T * w_V_dot[j];
                
                value = energy_hessian(S[0], S[1], S[2], T_dot[j](0, 0), T_dot[k](0, 0), T_dot[j](1, 1), T_dot[k](1, 1), T_dot[j](2, 2), T_dot[k](2, 2), T_ddot(0, 0), T_ddot(1, 1), T_ddot(2, 2), energy_type);
                
                //val[144 * i + j * 12 + k] = tet_vol * value;
                //val[144 * i + k * 12 + j] = tet_vol * value;
                
                ELS(j, k) = tet_vol * value;
                ELS(k, j) = tet_vol * value;
                
                row[144 * i + j * 12 + k] = index[j];
                col[144 * i + j * 12 + k] = index[k];
                
                row[144 * i + k * 12 + j] = index[k];
                col[144 * i + k * 12 + j] = index[j];
            }
        }
        
        
        SelfAdjointEigenSolver<MatrixXd> eig(ELS);
        
        ELS_S = eig.eigenvalues();
        
        ELS_X = eig.eigenvectors();
        
        for(int j = 0; j < 12; j++)
        {
            ELS_S[j] = fmax(ELS_S[j], proj_threshold);
        }
        
        ELS_D.diagonal() << ELS_S[0], ELS_S[1], ELS_S[2], ELS_S[3], ELS_S[4], ELS_S[5], ELS_S[6], ELS_S[7], ELS_S[8], ELS_S[9], ELS_S[10], ELS_S[11];
        
        ELS = ELS_X * ELS_D * ELS_X.transpose();
        
        for(int j = 0; j < 12; j++)
        {
            for(int k = j; k < 12; k++)
            {         
                val[144 * i + j * 12 + k] = ELS(j, k);
                val[144 * i + k * 12 + j] = ELS(k, j);
            }
        }
        
    }
    
    offset->clear();
    delete offset;
    
    return;
    
}