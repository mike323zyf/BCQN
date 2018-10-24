#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <iostream>


using namespace Eigen;
using namespace std;


double amips_param;
double energy_min;
double energy_max = 10;
double c1_param;
double c2_param;
double d1_param;


Vector3d color_map(double v)
{
	Vector3d c;
    
    c[0] = c[1] = c[2] = 1;

	double dv;

	if (v < energy_min)
		v = energy_min;
	if (v > energy_max)
		v = energy_max;
	dv = energy_max - energy_min;

	if (v < (energy_min + 0.25 * dv)) {
		c[0] = 0;
		c[1] = 4 * (v - energy_min) / dv;
	}
	else if (v < (energy_min + 0.5 * dv)) {
		c[0] = 0;
		c[2] = 1 + 4 * (energy_min + 0.25 * dv - v) / dv;
	}
	else if (v < (energy_min + 0.75 * dv)) {
		c[0] = 4 * (v - energy_min - 0.5 * dv) / dv;
		c[2] = 0;
	}
	else {
		c[1] = 1 + 4 * (energy_min + 0.75 * dv - v) / dv;
		c[2] = 0;
	}

	return c;
}


void set_energy_min(int type)
{

    switch(type)
    {
        case 0: //arap
            energy_min = 0;
            break;
            
        case 1: // mips
            energy_min = 2;
            break;
            
        case 2: // iso
            energy_min = 4;
            break;
            
        case 3: // amips
            energy_min = exp(amips_param * 2);
            break;
            
        case 4: // conf
            energy_min = 1;
            break;
            
        case 5: 
            energy_min = -c1_param - 2 * c2_param;
            break;
    }
    
}


double energy_value(int type, Vector2d S)
{
    
    double value = 0;
    double J;
    double l1;
    double l2;
    
    switch(type)
    {
        case 0: //arap
            value = (S[0] - 1) * (S[0] - 1) + (S[1] - 1) * (S[1] - 1);
            break;
            
        case 1: // mips
            value = S[0] / S[1] + S[1] / S[0];
            break;
            
        case 2: // iso
            value = S[0] * S[0] + 1.0 / (S[0] * S[0]) + S[1] * S[1] + 1.0 / (S[1] * S[1]);
            break;
            
        case 3: // amips
            value = exp(amips_param * (S[0] / S[1] + S[1] / S[0]));
            break;
            
        case 4: // conf
            value = S[0] / S[1];
            value *= value;
            break;
            
        case 5: // gmr
            J = S[0] * S[1];
            l1 = S[0] * S[0] + S[1] * S[1];
            l2 = S[0] * S[0] * S[1] * S[1];
            
            value = c1_param * (pow(J, -2.0 / 3.0) * l1 - 3) + c2_param * (pow(J, -4.0 / 3.0) * l2 - 3) + d1_param * (J - 1) * (J - 1);
            break;
            
        case 6: // olg
            
            value = max(S[0], 1.0 / S[1]);
            break;
    }
    
    return value;

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *output_mex;
    const int *dims;
    double *tri_num, *X_g_inv, *tri_areas, *obj_tri, *q_target, *type, *amips_s, *c1, *c2, *d1;
    double *output;
    
    tri_num = mxGetPr(prhs[0]);
    X_g_inv = mxGetPr(prhs[1]);
    tri_areas = mxGetPr(prhs[2]);
    obj_tri = mxGetPr(prhs[3]);
    q_target = mxGetPr(prhs[4]);
    type = mxGetPr(prhs[5]);
    amips_s = mxGetPr(prhs[6]);
    c1 = mxGetPr(prhs[7]);
    c2 = mxGetPr(prhs[8]);
    d1 = mxGetPr(prhs[9]);   
    
    int tri_n = tri_num[0];
    int energy_type = type[0];
    amips_param = amips_s[0];
    c1_param = c1[0];
    c2_param = c2[0];
    d1_param = d1[0];
    
    output_mex = plhs[0] = mxCreateDoubleMatrix(tri_n, 1, mxREAL);
    
    output = mxGetPr(output_mex);
    
    
    int tri[3];
    double tri_area;
    
    
    Matrix2d B, X_f, A;
    Vector2d S;
    Matrix2d U, V;
 
    for(int i = 0; i < tri_n; i++)
    {
        //mexPrintf("%d %d %d\n", (int)obj_tri[i], (int)obj_tri[i + tri_n], (int)obj_tri[i + 2 * tri_n]);
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
        
        if(A.determinant() <= 0)
        {
            //mexPrintf("element inverted");
        }
        
        JacobiSVD<Matrix2d> svd(A, ComputeFullU | ComputeFullV);
        
        S = svd.singularValues();
        U = svd.matrixU();
        V = svd.matrixV();
        
        output[i] = energy_value(energy_type, S);
        
    }
    
    return;
    
}