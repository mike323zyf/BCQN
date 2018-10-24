#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "robert/lapamg.h"


using namespace Eigen;
using namespace std;


LaplacianAMGStuff* S;
int expect_cycle;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *output_mex;
    
    double *nnz, *dof, *rows, *cols, *vals, *msg, *rhs;
    
    double *output;
    
    nnz = mxGetPr(prhs[0]);
    dof = mxGetPr(prhs[1]);
    rows = mxGetPr(prhs[2]);
    cols = mxGetPr(prhs[3]);
    vals = mxGetPr(prhs[4]);
    msg = mxGetPr(prhs[5]);
    rhs = mxGetPr(prhs[6]);
    
    int nnz_num = nnz[0];
    int dof_num = dof[0];
    int status = msg[0];
    
    output_mex = plhs[0] = mxCreateDoubleMatrix(dof_num, 1, mxREAL);
    
    output = mxGetPr(output_mex);
    
    if(status == 0)
    {
        int* rowbegin = (int*)malloc((dof_num + 1) * sizeof(int));
        int* column = (int*)malloc(nnz_num * sizeof(int));
        double* value = vals;
        int nmin = 1000;
        
        for(int i = 0; i < nnz_num; i++)
        {
            column[i] = cols[i] - 1;
            rows[i]--;
        }
        
        
        ///////////////////////////////////////////////////////////////
        
        FILE* file_a = fopen("matlab.txt", "wt");
        
        for(int i = 0; i < nnz_num; i++)
        {
            int tr = rows[i];
            int tc = cols[i] - 1;
            fprintf(file_a, "%d %d %f\n", tr, tc, vals[i]);
            
        }
        
        fclose(file_a);
        
        ///////////////////////////////////////////////////////////////
        
            
        rowbegin[0] = 0;
        
        int index = 0;
        int count = 0;
        
        for(int i = 0; i < nnz_num; i++)
        {
            if(rows[i] - index > 0.5)
            {
                index++;
                rowbegin[index] = rowbegin[index - 1] + count;
                count = 0;
            }
    
            count++;        
        }
        
        rowbegin[dof_num] = nnz_num;
        
        
        ///////////////////////////////////////////////////////////////
        
        FILE* file_b = fopen("csr.txt", "wt");
        
        for(int i = 0; i < dof_num + 1; i++)
        {
            
            fprintf(file_b, "%d ", rowbegin[i]);
            
        }
        
        fprintf(file_b, "\n\n");
        
        for(int i = 0; i < nnz_num; i++)
        {
            
            fprintf(file_b, "%d ", column[i]);
            
        }
        
        fprintf(file_b, "\n\n");
        
        for(int i = 0; i < nnz_num; i++)
        {
            
            fprintf(file_b, "%f ", value[i]);
            
        }
        
        fclose(file_b);
        
        ///////////////////////////////////////////////////////////////
        
         
        S = prepareLaplacianAMG(dof_num, rowbegin, column, value, nmin);
        
        expect_cycle = solveLaplacianAMG(S, rhs, output, 0, 100, 1e-8);
        
        printf("iteration %d\n", expect_cycle);

        free(rowbegin);
        free(column);
        
    }else{
        
        solveLaplacianAMG(S, rhs, output, expect_cycle, expect_cycle, 1e-8);
    }
     
    return;
    
}