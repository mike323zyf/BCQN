#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <iostream>


using namespace Eigen;
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *output_mex;
    
    double *tri_num, *obj_tri, *q_target;
    double *output;
     
    tri_num = mxGetPr(prhs[0]);
    obj_tri = mxGetPr(prhs[1]);
    q_target = mxGetPr(prhs[2]);
    
    int tri_n = tri_num[0];
    
    output_mex = plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    output = mxGetPr(output_mex);
    
    int tri[3];
    
    output[0] = 1;
     
    tri[0] = obj_tri[0] - 1; 
    tri[1] = obj_tri[0 + tri_n] - 1; 
    tri[2] = obj_tri[0 + 2 * tri_n] - 1; 
    
    Vector3d va(q_target[2 * tri[0]], q_target[2 * tri[0] + 1], 0);
    Vector3d vb(q_target[2 * tri[1]], q_target[2 * tri[1] + 1], 0);
    Vector3d vc(q_target[2 * tri[2]], q_target[2 * tri[2] + 1], 0);
    
    Vector3d e1 = va - vb;
    Vector3d e2 = vc - vb;
    Vector3d flag = e1.cross(e2);
     
    for(int i = 0; i < tri_n; i++)
    {   
        tri[0] = obj_tri[i] - 1; 
        tri[1] = obj_tri[i + tri_n] - 1; 
        tri[2] = obj_tri[i + 2 * tri_n] - 1; 
        
        Vector3d va(q_target[2 * tri[0]], q_target[2 * tri[0] + 1], 0);
        Vector3d vb(q_target[2 * tri[1]], q_target[2 * tri[1] + 1], 0);
        Vector3d vc(q_target[2 * tri[2]], q_target[2 * tri[2] + 1], 0);
        
        Vector3d e1 = va - vb;
        Vector3d e2 = vc - vb;
        Vector3d f1 = e1.cross(e2);
        
        if(f1.dot(flag) < 0)
        {
            output[0] = 0;
        }
    }
    
    return;
    
}