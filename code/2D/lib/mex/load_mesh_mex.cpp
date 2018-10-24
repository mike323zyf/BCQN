#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <igl/readOBJ.h>
#include <iostream>


using namespace Eigen;
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    mxArray *vers_mex, *tris_mex, *uvs_mex;
   
    char *name;
    double *vers, *tris, *uvs;
    
    name = mxArrayToString(prhs[0]);
    
    mexPrintf("%s\n", name);
    
    MatrixXd v, tc, n;
    MatrixXi f, ftc, fn;
    
    igl::readOBJ(name, v, tc, n, f, ftc, fn);
    
    mexPrintf("ver #: %d; tri #: %d; uv #: %d.\n", v.rows(), f.rows(), tc.rows());
    
    vers_mex = plhs[0] = mxCreateDoubleMatrix(v.rows() * 3, 1, mxREAL);
    tris_mex = plhs[1] = mxCreateDoubleMatrix(f.rows() * 3, 1, mxREAL);
    uvs_mex = plhs[2] = mxCreateDoubleMatrix(tc.rows() * 2, 1, mxREAL);
    
    vers = mxGetPr(vers_mex);
    tris = mxGetPr(tris_mex);
    uvs = mxGetPr(uvs_mex);
    
    for(int i = 0; i < v.rows(); i++)
    {
        vers[i] = v(i, 0);
        vers[v.rows() + i] = v(i, 1);
        vers[2 * v.rows() + i] = v(i, 2);
    }
    
    for(int i = 0; i < f.rows(); i++)
    {
        tris[i] = f(i, 0);
        tris[f.rows() + i] = f(i, 1);
        tris[2 * f.rows() + i] = f(i, 2);
    }
    
    for(int i = 0; i < tc.rows(); i++)
    {
        uvs[i] = tc(i, 0);
        uvs[tc.rows() + i] = tc(i, 1);
    }
    
    return;
    
}