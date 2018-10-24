#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <time.h>


using namespace Eigen;
using namespace std;


typedef Triplet<double> eigen_entry;
typedef vector<eigen_entry> entry_list;


double get_Fischer_Burmeister(VectorXd &x, VectorXd &pj, int num)
{
    
    double fb = 0;
    
    for (int i = 0; i < num; i++)
    {
        double ent = x[i] + pj[i] - sqrt(x[i] * x[i] + pj[i] * pj[i]);
        fb += ent * ent;
    }
    
    return sqrt(fb);
    
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *output_mex;

    double *J_index, *JT_index, *J_value, *JT_value, *JTJ_info, *pu, *wu, *bu, *max_iter, *step_size;
    double *output;
    
    J_index = mxGetPr(prhs[0]);
    JT_index = mxGetPr(prhs[1]);
    J_value = mxGetPr(prhs[2]);
    JT_value = mxGetPr(prhs[3]);
    JTJ_info = mxGetPr(prhs[4]); 
    pu = mxGetPr(prhs[5]);
    wu = mxGetPr(prhs[6]);
    bu = mxGetPr(prhs[7]);
    max_iter = mxGetPr(prhs[8]);
    step_size = mxGetPr(prhs[9]);
    
    // JTJ_info : dof_n, max_valence, tri_n, 3 * 2
    
    int J_row = JTJ_info[0];
    int J_col = JTJ_info[2];
    int max_it = max_iter[0];
    double step = step_size[0];
    
    output_mex = plhs[0] = mxCreateDoubleMatrix(J_row, 1, mxREAL);
    
    output = mxGetPr(output_mex);

    SparseMatrix<double> J_u = SparseMatrix<double>(J_row, J_col);
    SparseMatrix<double> JT_u = SparseMatrix<double>(J_col, J_row);
    VectorXd b_u = VectorXd::Zero(J_col);
     
    entry_list en_list;
    
    for(int i = 0; i < JTJ_info[0]; i++)
    {
        for(int j = 0; j < JTJ_info[1]; j++)
        {
            int index = i * JTJ_info[1] + j;
            
            if(J_index[index] >= 0)
            {
                en_list.push_back(eigen_entry(i, J_index[index], J_value[index]));
            }
        }
    }
    
    J_u.setFromTriplets(en_list.begin(), en_list.end());
    
    
    en_list.clear();
    
    for(int i = 0; i < JTJ_info[2]; i++)
    {
        for(int j = 0; j < JTJ_info[3]; j++)
        {
            int index = i * JTJ_info[3] + j;
            
            if(JT_index[index] >= 0)
            {
                en_list.push_back(eigen_entry(i, JT_index[index], JT_value[index]));
            }
        }
    }
    
    JT_u.setFromTriplets(en_list.begin(), en_list.end());
    
    
    VectorXd p_u = VectorXd::Zero(J_row);
    
    for(int i = 0; i < J_row; i++)
    {
        p_u[i] = pu[i];
    }
    
    for(int i = 0; i < J_col; i++)
    {
        b_u[i] = bu[i];
    }
   
    VectorXd b = b_u + JT_u * p_u;
    
    
    if(b.minCoeff() < 0)
    {
    
        VectorXd lambda = VectorXd::Zero(J_col);
    
        VectorXd pi;
    
        double pre_FB, post_FB;
    
        int iteration = -1;
    
        ///////////////////////////////////////////////////////////////////////////////////////////////
    
        pi = JT_u * (J_u * lambda) + b;
        
        post_FB = get_Fischer_Burmeister(lambda, pi, J_col);
    
        for(int outer = 0; outer < max_it; outer++)
        {
            iteration = outer;
        
        
            pre_FB = post_FB;
        
            for(int iter = 0; iter < J_col; iter++)
            {
                lambda[iter] = fmax(lambda[iter] - step * pi[iter] / wu[iter], 0.0);
            }
        
            pi = JT_u * (J_u * lambda) + b;
        
            post_FB = get_Fischer_Burmeister(lambda, pi, J_col);
            
            if(abs(post_FB) < 1e-3)
            {
                break;
            }
        
            if((pre_FB - post_FB) / pre_FB < 1e-3)
            {
                break;
            }
        
        }
    
        //mexPrintf("%d\n", iteration);
    
        VectorXd p_mod = J_u * lambda;
    
        for(int i = 0; i < J_row; i++)
        {
            output[i] = p_mod[i];
        }
        
    }else{
        for(int i = 0; i < J_row; i++)
        {
            output[i] = 0;
        }
    }
   
    return;
    
}















