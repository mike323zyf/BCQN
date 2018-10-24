#include "mex.h"
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "mexHelpers.cpp"
#include "panozo/Param_State.h"
#include "panozo/GlobalLocalParametrization.h"
#include "panozo/StateManager.h"
#include <igl/components.h>
#include <iostream>

using namespace std;
using namespace Eigen;


void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray*prhs[])
{
	// assign input
	int n_tri = mxGetM(prhs[0]); // # rows of F
	int d_simplex = mxGetN(prhs[0]); // # cols of F
	int n_vert = mxGetM(prhs[1]); // # rows of V
	int dim = mxGetN(prhs[1]); // # cols of V
	const Map<MatrixXd, Aligned> Fmatlab(mxGetPr(prhs[0]), n_tri, d_simplex);
	const Map<MatrixXd, Aligned> V(mxGetPr(prhs[1]), n_vert, dim);
	
	// update index numbers to 0-base
	MatrixXd Fd (Fmatlab);
	Fd = Fd.array() - 1;
    MatrixXi F = MatrixXi::Zero(Fd.rows(), Fd.cols());
    
    for(int i = 0; i < Fd.rows(); i++)
    {
        for(int j = 0; j < Fd.cols(); j++)
        {
            F(i, j) = Fd(i, j);
        }
    }

	// compute
    
    Param_State state;
    
    state.method = Param_State::GLOBAL_ARAP_IRLS;
    state.flips_linesearch = true;
    state.update_all_energies = false;
    state.proximal_p = 0.0001;
    
    state.V = V;
    state.F = F;
    state.v_num = state.V.rows();
    state.f_num = state.F.rows();
    
    igl::doublearea(state.V,state.F, state.M); state.M /= 2.;

    state.global_local_energy = Param_State::SYMMETRIC_DIRICHLET;
    state.cnt_flips = false;
    
    state.mesh_area = state.M.sum();
    //state.V /= sqrt(state.mesh_area);
    //state.mesh_area = 1;
    
    StateManager state_manager;
    GlobalLocalParametrization param(state_manager, &state);
    
    param.init_parametrization();
    
    MatrixXd UV = state.uv;

	// assign outputs
	mapDenseMatrixToMex(UV, &(plhs[0]));
}