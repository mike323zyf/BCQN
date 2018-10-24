#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <time.h>
#include <tbb/tbb.h>


using namespace Eigen;
using namespace std;
using namespace tbb;


double *J_index, *JT_index, *J_value, *JT_value, *JTJ_info, *pu, *wu, *bu, *max_iter, *step_size;
double *b, *lambda, *pi, *Jx_in, *p_mod, *dummy;

int J_row;
int J_col;
int row_spy;
int col_spy;
int max_it;
double step;

class spmv_mul_add
{

public:
	spmv_mul_add(double* A_index, double* A_value, double* input, double* addition, double* output, int spy)
	{

		this->A_index = A_index;	
		this->A_value = A_value;	
		this->input = input;
        this->addition = addition;
		this->output = output;
		this->spy = spy;

	};

	~spmv_mul_add() {};

	void operator() (const blocked_range<int>& r) const
	{
		for (int i = r.begin(); i != r.end(); i++)
		{
			output[i] = 0;

			for (int j = 0; j < spy; j++)
			{
				int index = i * spy + j;
                
                if(A_index[index] >= 0)
                {
                    int idx = A_index[index];
                    output[i] += A_value[index] * input[idx];
                }	
			} 
            
            output[i] += addition[i];
		}
	};

	double* A_index;
	double* A_value;
	double* input;
    double* addition;
	double* output;
	int spy;

};


class jacobi_update
{

public:
	jacobi_update(double* dof, double* update, double* divisor, double alpha)
	{
		this->dof = dof;	
		this->update = update;	
		this->divisor = divisor;
        this->alpha = alpha;
	};

	~jacobi_update() {};

	void operator() (const blocked_range<int>& r) const
	{
		for (int i = r.begin(); i != r.end(); i++)
		{
            dof[i] = fmax(dof[i] - alpha * update[i] / divisor[i], 0.0);
		}
	};

	double* dof;
	double* update;
	double* divisor;
    double alpha;

};


class set_zero
{

public:
	set_zero(double* data)
	{
		this->data = data;	
	};

	~set_zero() {};

	void operator() (const blocked_range<int>& r) const
	{
		for (int i = r.begin(); i != r.end(); i++)
		{
			data[i] = 0;
		}
	};

	double* data;

};


double get_Fischer_Burmeister(double *x, double *pj, int num)
{
    
    double fb = 0;
    
    for (int i = 0; i < num; i++)
    {
        double ent = x[i] + pj[i] - sqrt(x[i] * x[i] + pj[i] * pj[i]);
        fb += ent * ent;
    }
    
    return sqrt(fb);
    
}


void damped_Jacobi()
{
   
    task_scheduler_init init(8);
    
    parallel_for(blocked_range<int>(0, J_col), set_zero(lambda));
    
    parallel_for(blocked_range<int>(0, J_row), set_zero(dummy));
    
    parallel_for(blocked_range<int>(0, J_col), spmv_mul_add(JT_index, JT_value, pu, bu, b, col_spy));
    
    double pre_FB, post_FB;
    
    parallel_for(blocked_range<int>(0, J_row), spmv_mul_add(J_index, J_value, lambda, dummy, Jx_in, row_spy));
    
    parallel_for(blocked_range<int>(0, J_col), spmv_mul_add(JT_index, JT_value, Jx_in, b, pi, col_spy));
    
    post_FB = get_Fischer_Burmeister(lambda, pi, J_col);
    
    int iteration = -1;
    
    for(int outer = 0; outer < max_it; outer++)
    {
        iteration = outer;
        pre_FB = post_FB;
        
        parallel_for(blocked_range<int>(0, J_col), jacobi_update(lambda, pi, wu, step));
        
        ///////////////////////////////////////////////////////////////////////////////
        
        parallel_for(blocked_range<int>(0, J_row), spmv_mul_add(J_index, J_value, lambda, dummy, Jx_in, row_spy));
        
        parallel_for(blocked_range<int>(0, J_col), spmv_mul_add(JT_index, JT_value, Jx_in, b, pi, col_spy));
        
        ///////////////////////////////////////////////////////////////////////////////
        
        post_FB = get_Fischer_Burmeister(lambda, pi, J_col);
        
        if((pre_FB - post_FB) / pre_FB < 1e-3)
        {
            break;
        }
        
    }
    
    //mexPrintf("%d\n", iteration);
    
    parallel_for(blocked_range<int>(0, J_row), spmv_mul_add(J_index, J_value, lambda, dummy, p_mod, row_spy));
  
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mxArray *output_mex;

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
    
    J_row = JTJ_info[0];
    J_col = JTJ_info[2];
    row_spy = JTJ_info[1];
    col_spy = JTJ_info[3];
    max_it = max_iter[0];
    step = step_size[0];
    
    b = (double*)malloc(J_col * sizeof(double));
    lambda = (double*)malloc(J_col * sizeof(double));
    pi = (double*)malloc(J_col * sizeof(double));
    Jx_in = (double*)malloc(J_row * sizeof(double));
    dummy = (double*)malloc(J_row * sizeof(double));
    
    output_mex = plhs[0] = mxCreateDoubleMatrix(J_row, 1, mxREAL);   
    p_mod = mxGetPr(output_mex);
    
    damped_Jacobi();
       
    free(b);
    free(lambda);
    free(pi);
    free(Jx_in);
    free(dummy);
   
    return;
    
}















