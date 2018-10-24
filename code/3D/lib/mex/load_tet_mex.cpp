#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <stdio.h>


using namespace Eigen;
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    mxArray *vers_mex, *tets_mex;
   
    char *name;
    double *vers, *tets;
    
    name = mxArrayToString(prhs[0]);
    
    mexPrintf("%s\n", name);
    
    int ver_num, tet_num;
    float x, y, z;
    int i1, i2, i3, i4;
    
    FILE* tet_file = fopen(name, "rt");
    
    fscanf(tet_file, "#num of vertex: %d\n", &ver_num);
    fscanf(tet_file, "#num of tetrahedra: %d\n\n", &tet_num);
    
    MatrixXd ver = MatrixXd::Zero(ver_num, 3);
    MatrixXi tet = MatrixXi::Zero(tet_num, 4);
    
    for(int i = 0; i < ver_num; i++)
    {
        fscanf(tet_file, "%f %f %f\n", &x, &y, &z);
        
        ver(i, 0) = x;
        ver(i, 1) = y;
        ver(i, 2) = z;
    }
    
    fscanf(tet_file, "\n");
    
    for(int i = 0; i < tet_num; i++)
    {
        fscanf(tet_file, "%d %d %d %d\n", &i1, &i2, &i3, &i4);
        
        tet(i, 0) = i1;
        tet(i, 1) = i2;
        tet(i, 2) = i3;
        tet(i, 3) = i4;
    }
    
    fclose(tet_file);
    
    vers_mex = plhs[0] = mxCreateDoubleMatrix(ver.rows() * 3, 1, mxREAL);
    tets_mex = plhs[1] = mxCreateDoubleMatrix(tet.rows() * 4, 1, mxREAL);
    
    vers = mxGetPr(vers_mex);
    tets = mxGetPr(tets_mex);
    
    for(int i = 0; i < ver.rows(); i++)
    {
        vers[i] = ver(i, 0);
        vers[ver.rows() + i] = ver(i, 1);
        vers[2 * ver.rows() + i] = ver(i, 2);
    }
    
    for(int i = 0; i < tet.rows(); i++)
    {
        tets[i] = tet(i, 0);
        tets[tet.rows() + i] = tet(i, 1);
        tets[2 * tet.rows() + i] = tet(i, 2);
        tets[3 * tet.rows() + i] = tet(i, 3);
    }
    
    return;
    
}