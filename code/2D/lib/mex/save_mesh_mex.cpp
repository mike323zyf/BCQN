#include <mex.h> 
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <vector>


using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
      
    char *file_name;
    
    double *simplex_num, *ver_num, *simplex, *ver, *simplex_dim, *ver_dim;
    
    file_name = mxArrayToString(prhs[6]);
    
    simplex_num = mxGetPr(prhs[0]);
    ver_num = mxGetPr(prhs[1]);
    simplex = mxGetPr(prhs[2]);
    ver = mxGetPr(prhs[3]);
    simplex_dim = mxGetPr(prhs[4]);
    ver_dim = mxGetPr(prhs[5]);
    
    int simplex_n = simplex_num[0];
    int ver_n = ver_num[0];
    int simplex_d = simplex_dim[0];
    int ver_d = ver_dim[0];
    
    vector<int> sim;
    
    for(int i = 0; i < simplex_d; i++)
    {
        sim.push_back(0);
    }
 
    FILE* obj_file = fopen(file_name, "wt");
    
    for(int i = 0; i < ver_n; i++)
    {
        if(ver_d == 2)
        {
            fprintf(obj_file, "v %f %f 0.0\n", ver[2 * i + 0], ver[2 * i + 1]);
        }else{
            fprintf(obj_file, "v %f %f %f\n", ver[3 * i + 0], ver[3 * i + 1], ver[3 * i + 2]);
        }   
    }
    
    fprintf(obj_file, "\n");
    
    for(int i = 0; i < simplex_n; i++)
    {
        sim[0] = simplex[i]; 
        sim[1] = simplex[i + simplex_n]; 
        sim[2] = simplex[i + 2 * simplex_n]; 
        
        fprintf(obj_file, "f %d %d %d\n", sim[0], sim[1], sim[2]);
    }
    
    fclose(obj_file);
    
    return;
    
}