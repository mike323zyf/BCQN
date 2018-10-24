#include <mex.h>
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <string>


using namespace Eigen;
using namespace std;

double amips_param;
double c1_param;
double c2_param;
double d1_param;


Vector2d energy_derivative(int type, Vector2d S)
{
    
    Vector2d es;
    double J;
    double l1;
    double l2;
    
    switch(type)
    {
        case 0: //arap
            es[0] = 2 * (S[0] - 1);
            es[1] = 2 * (S[1] - 1);
            break;
            
        case 1: // mips
            es[0] = 1.0 / S[1] - S[1] / (S[0] * S[0]);
            es[1] = 1.0 / S[0] - S[0] / (S[1] * S[1]);
            //es[0] *= 100;
            //es[1] *= 100;
            break;
            
        case 2: // iso
            es[0] = 2 * S[0] - 2.0 / (S[0] * S[0] * S[0]);
            es[1] = 2 * S[1] - 2.0 / (S[1] * S[1] * S[1]);
            break;
            
        case 3: // amips
            es[0] = amips_param * exp(amips_param * (S[0] / S[1] + S[1] / S[0])) * (1.0 / S[1] - S[1] / (S[0] * S[0]));
            es[1] = amips_param * exp(amips_param * (S[0] / S[1] + S[1] / S[0])) * (1.0 / S[0] - S[0] / (S[1] * S[1]));
            break;
            
        case 4: // conf
            es[0] = 2 * S[0] / (S[1] * S[1]);
            es[1] = -2 * S[0] * S[0] / (S[1] * S[1] * S[1]);
            break;
            
        case 5: // gmr
						//            J = S[0] * S[1];
						//            l1 = S[0] * S[0] + S[1] * S[1];
						//            l2 = S[0] * S[0] * S[1] * S[1];
						//            
						//            es[0] = c1_param * (-2.0 / 3.0 * pow(J, -5.0 / 3.0) * S[1] * l1 + pow(J, -2.0 / 3.0) * 2 * S[0]) 
						//                    + c2_param * (-4.0 / 3.0 * pow(J, -7.0 / 3.0) * S[1] * l2 + pow(J, -4.0 / 3.0) * 2 * S[1] * S[1] * S[0]) 
						//                    + d1_param * 2 * (J - 1) * S[1];
						//            
						//            es[1] = c1_param * (-2.0 / 3.0 * pow(J, -5.0 / 3.0) * S[0] * l1 + pow(J, -2.0 / 3.0) * 2 * S[1]) 
						//                    + c2_param * (-4.0 / 3.0 * pow(J, -7.0 / 3.0) * S[0] * l2 + pow(J, -4.0 / 3.0) * 2 * S[0] * S[0] * S[1]) 
						//                    + d1_param * 2 * (J - 1) * S[0];

						es[0] = c1_param * (1.0 / S[1] - S[1] / (S[0] * S[0])) - d1_param *2.*S[1]*(1.0 - S[0]*S[1]);
						es[1] = c1_param * (1.0 / S[0] - S[0] / (S[1] * S[1])) - d1_param *2.*S[0]*(1.0 - S[0]*S[1]);

						break;
				
				
        case 6: // olg
            
            if(S[0] > 1.0 / S[1])
            {
                es[0] = 1;
                es[1] = 0;
            }else{
                es[0] = 0;
                es[1] = -1.0 / (S[1] * S[1]);
            }
            
            break;
    }
    
    return es;

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *output_mex, *J_value_mex, *JT_value_mex, *wu_mex, *bu_mex;
    const int *dims;
    double *tri_num, *X_g_inv, *tri_areas, *obj_tri, *q_target, *type, *amips_s, *F_dot, *ver_num, *c1_g, *c2_g, *d1_g, *J_index, *JT_index, *JTJ_info, *x2u;
    double *output, *J_value, *JT_value, *wu, *bu;
    
    tri_num = mxGetPr(prhs[0]);
    X_g_inv = mxGetPr(prhs[1]);
    tri_areas = mxGetPr(prhs[2]);
    obj_tri = mxGetPr(prhs[3]);
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
    
    int tri_n = tri_num[0];
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
    
    output_mex = plhs[0] = mxCreateDoubleMatrix(2 * ver_n, 1, mxREAL);
    J_value_mex = plhs[1] = mxCreateDoubleMatrix(J_rn * J_cn, 1, mxREAL);
    JT_value_mex = plhs[2] = mxCreateDoubleMatrix(JT_rn * JT_cn, 1, mxREAL);
    wu_mex = plhs[3] = mxCreateDoubleMatrix(tri_n, 1, mxREAL);
    bu_mex = plhs[4] = mxCreateDoubleMatrix(tri_n, 1, mxREAL);
    
    output = mxGetPr(output_mex);
    J_value = mxGetPr(J_value_mex);
    JT_value = mxGetPr(JT_value_mex);
    wu = mxGetPr(wu_mex);
    bu = mxGetPr(bu_mex);
    
    memset(J_value, 0, J_rn * J_cn * sizeof(double));
    memset(JT_value, 0, JT_rn * JT_cn * sizeof(double));
    memset(output, 0, 2 * ver_n * sizeof(double));

    vector<int> offset;
    
    offset.resize(2 * ver_n, 0);
    
    int tri[3];
    double tri_area;
    
    Matrix2d B, X_f, A, F_d;
    Vector2d S, u0, u1, v0, v1;
    Matrix2d U, V;
    
    Vector2d es;
    Vector2d cs;
    
    double es1, es2;
    double d1, d2;
    double dphi;
  
 
    for(int i = 0; i < tri_n; i++)
    {
        tri[0] = obj_tri[i] - 1; 
        tri[1] = obj_tri[i + tri_n] - 1; 
        tri[2] = obj_tri[i + 2 * tri_n] - 1; 
        tri_area = tri_areas[i];
        
        B(0, 0) = X_g_inv[i];
        B(0, 1) = X_g_inv[i + tri_n * 2];
        B(1, 0) = X_g_inv[i + tri_n];
        B(1, 1) = X_g_inv[i + tri_n * 3];
        
        X_f(0, 0) = q_target[2 * tri[1]] - q_target[2 * tri[0]];
        X_f(1, 0) = q_target[2 * tri[1] + 1] - q_target[2 * tri[0] + 1];  
        X_f(0, 1) = q_target[2 * tri[2]] - q_target[2 * tri[0]];
        X_f(1, 1) = q_target[2 * tri[2] + 1] - q_target[2 * tri[0] + 1];
        
        A = X_f * B;
        
        JacobiSVD<Matrix2d> svd(A, ComputeFullU | ComputeFullV);
        
        S = svd.singularValues();
        U = svd.matrixU();
        V = svd.matrixV();
        
        if(A.determinant() <= 0)
        {
            //mexPrintf("inverted\n");
        }
        
        bu[i] = S[0] * S[1];
        
        es = energy_derivative(energy_type, S);
        
        cs = S;
        
        es1 = es[0];
        es2 = es[1];
        
        ////////////////////////////////////////////
        
        wu[i] = 0;
        
        u0 = U.col(0);
        u1 = U.col(1);
        v0 = V.col(0);
        v1 = V.col(1);
        
        for(int j = 0; j < 6; j++)
        {
            
            F_d(0, 0) = F_dot[i + j * tri_n];
            F_d(0, 1) = F_dot[i + j * tri_n + 6 * 2 * tri_n];
            F_d(1, 0) = F_dot[i + j * tri_n + 6 * tri_n];
            F_d(1, 1) = F_dot[i + j * tri_n + 6 * 3 * tri_n];
            
            d1 = u0.transpose() * F_d * v0;
            d2 = u1.transpose() * F_d * v1;
            
            dphi = d1 * es1 + d2 * es2;
            
            int cache_1 = 2 * tri[j / 2] + (j % 2);
            int cache_2 = x2u[cache_1];
            
            output[cache_1] += tri_area * dphi;
            
            if(cache_2 > 0)
            {
                int idx = cache_2 - 1;
                
                double value = d1 * cs[1] + d2 * cs[0];
                
                J_value[idx * J_cn + offset[cache_1]] = value;
                
                offset[cache_1]++;
                
                JT_value[i * JT_cn + j] = value;
                
                wu[i] += value * value;
            }
        }  
    }
    
    return;
    
}