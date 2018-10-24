#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <stdio.h>


using namespace Eigen;
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    char *file_name;
    
    double *tet_num, *ver_num, *obj_tet, *q_target;
    
    tet_num = mxGetPr(prhs[0]);
    ver_num = mxGetPr(prhs[1]);
    obj_tet = mxGetPr(prhs[2]);
    q_target = mxGetPr(prhs[3]);
    file_name = mxArrayToString(prhs[4]);
    
    int tet_n = tet_num[0];
    int ver_n = ver_num[0];
    
    FILE* tet_file = fopen(file_name, "wt");
    
    fprintf(tet_file, "#num of vertex: %d\n", ver_n);
    fprintf(tet_file, "#num of tetrahedra: %d\n\n", tet_n);
    
    for(int i = 0; i < ver_n; i++)
    {
        fprintf(tet_file, "%f %f %f\n", q_target[3 * i], q_target[3 * i + 1], q_target[3 * i + 2]);
    }
    
    fprintf(tet_file, "\n");
    
    int tet[4];
    
    for(int i = 0; i < tet_n; i++)
    {
        tet[0] = obj_tet[i] - 1; 
        tet[1] = obj_tet[i + tet_n] - 1; 
        tet[2] = obj_tet[i + 2 * tet_n] - 1; 
        tet[3] = obj_tet[i + 3 * tet_n] - 1; 
        
        fprintf(tet_file, "%d %d %d %d\n", tet[0], tet[1], tet[2], tet[3]);
    }

    fclose(tet_file);
    
    return;
    
}