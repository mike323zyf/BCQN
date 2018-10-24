#include <mex.h> 
#include <math.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>


using namespace Eigen;
using namespace std;


double get_tet_vol(Vector3d p0, Vector3d p1, Vector3d p2, Vector3d p3)
{
    
    Vector3d u = p1 - p0;
    Vector3d v = p2 - p0;
    Vector3d w = p3 - p0;
    
    return abs((u.cross(v)).dot(w)) / 6.0;
    
}


double get_tri_area(Vector3d p0, Vector3d p1, Vector3d p2)
{
    
    Vector3d u = p1 - p0;
    Vector3d v = p2 - p0;
    
    return 0.5 * sqrt(u.dot(u) * v.dot(v) - u.dot(v) * u.dot(v));
    
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    mxArray *X_g_inv_mex, *tet_vols_mex, *F_dot_mex, *row_mex, *col_mex, *val_mex, *x2u_mex, *J_mex, *J_info_mex, *JT_mex, *JT_info_mex, *perimeter_mex;
    
    double *tet_num, *tet_list, *ver_num, *ver_list, *dirichlet;
    double *X_g_inv, *tet_vols, *F_dot, *row, *col, *val, *x2u, *J, *J_info, *JT, *JT_info, *perimeter;

    tet_num = mxGetPr(prhs[0]);
    tet_list = mxGetPr(prhs[1]);
    ver_num = mxGetPr(prhs[2]);
    ver_list = mxGetPr(prhs[3]);
    dirichlet = mxGetPr(prhs[4]);
    
    int tet_n = tet_num[0];
    int ver_n = ver_num[0];
    
    X_g_inv_mex = plhs[0] = mxCreateDoubleMatrix(9 * tet_n, 1, mxREAL);
    tet_vols_mex = plhs[1] = mxCreateDoubleMatrix(tet_n, 1, mxREAL);
    F_dot_mex = plhs[2] = mxCreateDoubleMatrix(108 * tet_n, 1, mxREAL);
    
    
    row_mex = plhs[3] = mxCreateDoubleMatrix(48 * tet_n, 1, mxREAL);
    col_mex = plhs[4] = mxCreateDoubleMatrix(48 * tet_n, 1, mxREAL);
    val_mex = plhs[5] = mxCreateDoubleMatrix(48 * tet_n, 1, mxREAL);
    
    x2u_mex = plhs[6] = mxCreateDoubleMatrix(3 * ver_n, 1, mxREAL);
    
    X_g_inv = mxGetPr(X_g_inv_mex);
    tet_vols = mxGetPr(tet_vols_mex);
    F_dot = mxGetPr(F_dot_mex);
    row = mxGetPr(row_mex);
    col = mxGetPr(col_mex);
    val = mxGetPr(val_mex);
    
    x2u = mxGetPr(x2u_mex);
    
    int count = 1;
    int fixed_num = 0;
    int tet[4];
    
    for(int i = 0; i < 3 * ver_n; i++)
    {
        if(dirichlet[i] == 0)
        {
            x2u[i] = count;
            count++;
        }else
        {
            fixed_num++;
            x2u[i] = 0;
        }
    }
    
    vector<vector<int>*> ver_tet_map;
    
    for(int i = 0; i < ver_n; i++)
    {
        vector<int>* tets = new vector<int>();
        ver_tet_map.push_back(tets);
    }
    
    for(int i = 0; i < tet_n; i++)
    {
        tet[0] = tet_list[i] - 1; 
        tet[1] = tet_list[i + tet_n] - 1; 
        tet[2] = tet_list[i + 2 * tet_n] - 1;
        tet[3] = tet_list[i + 3 * tet_n] - 1;
        
        ver_tet_map[tet[0]]->push_back(i);
        ver_tet_map[tet[1]]->push_back(i);
        ver_tet_map[tet[2]]->push_back(i);
        ver_tet_map[tet[3]]->push_back(i);
    }
    
    int max_valence = 0;
    
    for(int i = 0; i < ver_n; i++)
    {
        if(max_valence < ver_tet_map[i]->size())
        {
            max_valence = ver_tet_map[i]->size();
        }
    }
    
    J_mex = plhs[7] = mxCreateDoubleMatrix((3 * ver_n - fixed_num) * max_valence, 1, mxREAL);
    J_info_mex = plhs[8] = mxCreateDoubleMatrix(2, 1, mxREAL);
    JT_mex = plhs[9] = mxCreateDoubleMatrix(tet_n * 4 * 3, 1, mxREAL);
    JT_info_mex = plhs[10] = mxCreateDoubleMatrix(2, 1, mxREAL);
    perimeter_mex = plhs[11] = mxCreateDoubleMatrix(ver_n, 1, mxREAL);

    J = mxGetPr(J_mex);
    J_info = mxGetPr(J_info_mex);
    JT = mxGetPr(JT_mex);
    JT_info = mxGetPr(JT_info_mex);
    perimeter = mxGetPr(perimeter_mex);
    
    count = 0;
    
    J_info[0] = 3 * ver_n - fixed_num;
    J_info[1] = max_valence;
    
    for(int i = 0; i < ver_n; i++)
    {
        perimeter[i] = 0;
        
        if(dirichlet[3 * i] == 0)
        {
            for(int j = 0; j < max_valence; j++)
            {
                if(j < ver_tet_map[i]->size())
                {
                    J[count * max_valence + j] = ver_tet_map[i]->at(j);
                }else{
                    J[count * max_valence + j] = -1;
                }
            }
            
            count++;
            
            for(int j = 0; j < max_valence; j++)
            {
                if(j < ver_tet_map[i]->size())
                {
                    J[count * max_valence + j] = ver_tet_map[i]->at(j);
                }else{
                    J[count * max_valence + j] = -1;
                }
            }
            
            count++;
            
            for(int j = 0; j < max_valence; j++)
            {
                if(j < ver_tet_map[i]->size())
                {
                    J[count * max_valence + j] = ver_tet_map[i]->at(j);
                }else{
                    J[count * max_valence + j] = -1;
                }
            }
            
            count++;
        }
    }
    
    JT_info[0] = tet_n;
    JT_info[1] = 4 * 3;

    for(int i = 0; i < tet_n; i++)
    { 
        tet[0] = tet_list[i] - 1; 
        tet[1] = tet_list[i + tet_n] - 1; 
        tet[2] = tet_list[i + 2 * tet_n] - 1;
        tet[3] = tet_list[i + 3 * tet_n] - 1;
        
        for(int j = 0; j < 4; j++)
        {
            if(x2u[3 * tet[j]] > 0)
            {
                JT[i * 3 * 4 + 3 * j + 0] = x2u[3 * tet[j] + 0] - 1;
                JT[i * 3 * 4 + 3 * j + 1] = x2u[3 * tet[j] + 1] - 1;
                JT[i * 3 * 4 + 3 * j + 2] = x2u[3 * tet[j] + 2] - 1;
            }else{
                JT[i * 3 * 4 + 3 * j + 0] = -1;
                JT[i * 3 * 4 + 3 * j + 1] = -1;
                JT[i * 3 * 4 + 3 * j + 2] = -1;
            }     
        }
    }

    for(int i = 0; i < ver_n; i++)
    {
        ver_tet_map[i]->clear();
        delete ver_tet_map[i];
    }
    
    ver_tet_map.clear();
 
    ///////////////////////////////////////////////////////////////////////
    
    for(int i = 0; i < 108 * tet_n; i++)
    {
        F_dot[i] = 0;
    }
  
    Matrix3d X_g, B;
    
    double tet_vol;
    
    for(int i = 0; i < tet_n; i++)
    {
        tet[0] = tet_list[i] - 1; 
        tet[1] = tet_list[i + tet_n] - 1; 
        tet[2] = tet_list[i + 2 * tet_n] - 1;
        tet[3] = tet_list[i + 3 * tet_n] - 1;
        
        Vector3d p0(ver_list[3 * tet[0]], ver_list[3 * tet[0] + 1], ver_list[3 * tet[0] + 2]);
        Vector3d p1(ver_list[3 * tet[1]], ver_list[3 * tet[1] + 1], ver_list[3 * tet[1] + 2]);
        Vector3d p2(ver_list[3 * tet[2]], ver_list[3 * tet[2] + 1], ver_list[3 * tet[2] + 2]);
        Vector3d p3(ver_list[3 * tet[3]], ver_list[3 * tet[3] + 1], ver_list[3 * tet[3] + 2]);
        
        perimeter[tet[0]] += get_tri_area(p1, p2, p3);
        perimeter[tet[1]] += get_tri_area(p0, p2, p3);
        perimeter[tet[2]] += get_tri_area(p1, p0, p3);
        perimeter[tet[3]] += get_tri_area(p1, p2, p0);
        
        tet_vol = tet_vols[i] = get_tet_vol(p0, p1, p2, p3);
        
        Vector3d e0 = p1 - p0;
        Vector3d e1 = p2 - p0;
        Vector3d e2 = p3 - p0;
        
        X_g << e0[0], e1[0], e2[0], e0[1], e1[1], e2[1], e0[2], e1[2], e2[2];
        
        //cout << X_g.determinant() << endl;
        
        B = X_g.inverse();
        
        X_g_inv[i] = B(0, 0);
        X_g_inv[tet_n + i] = B(1, 0);
        X_g_inv[2 * tet_n + i] = B(2, 0);
        X_g_inv[3 * tet_n + i] = B(0, 1);
        X_g_inv[4 * tet_n + i] = B(1, 1);
        X_g_inv[5 * tet_n + i] = B(2, 1);      
        X_g_inv[6 * tet_n + i] = B(0, 2);
        X_g_inv[7 * tet_n + i] = B(1, 2);
        X_g_inv[8 * tet_n + i] = B(2, 2);
    
        
        F_dot[0 * tet_n * 12 + 0 * tet_n + i] = -1.0 * (B(0, 0) + B(1, 0) + B(2, 0));// 0 0
        F_dot[3 * tet_n * 12 + 0 * tet_n + i] = -1.0 * (B(0, 1) + B(1, 1) + B(2, 1));// 0 1
        F_dot[6 * tet_n * 12 + 0 * tet_n + i] = -1.0 * (B(0, 2) + B(1, 2) + B(2, 2));// 0 2
        
        F_dot[1 * tet_n * 12 + 1 * tet_n + i] = -1.0 * (B(0, 0) + B(1, 0) + B(2, 0));// 1 0
        F_dot[4 * tet_n * 12 + 1 * tet_n + i] = -1.0 * (B(0, 1) + B(1, 1) + B(2, 1));// 1 1
        F_dot[7 * tet_n * 12 + 1 * tet_n + i] = -1.0 * (B(0, 2) + B(1, 2) + B(2, 2));// 1 2
        
        F_dot[2 * tet_n * 12 + 2 * tet_n + i] = -1.0 * (B(0, 0) + B(1, 0) + B(2, 0));// 2 0
        F_dot[5 * tet_n * 12 + 2 * tet_n + i] = -1.0 * (B(0, 1) + B(1, 1) + B(2, 1));// 2 1
        F_dot[8 * tet_n * 12 + 2 * tet_n + i] = -1.0 * (B(0, 2) + B(1, 2) + B(2, 2));// 2 2
        
        
        F_dot[0 * tet_n * 12 + 3 * tet_n + i] = B(0, 0);// 0 0
        F_dot[3 * tet_n * 12 + 3 * tet_n + i] = B(0, 1);// 0 1
        F_dot[6 * tet_n * 12 + 3 * tet_n + i] = B(0, 2);// 0 2
        
        F_dot[1 * tet_n * 12 + 4 * tet_n + i] = B(0, 0);// 1 0
        F_dot[4 * tet_n * 12 + 4 * tet_n + i] = B(0, 1);// 1 1
        F_dot[7 * tet_n * 12 + 4 * tet_n + i] = B(0, 2);// 1 2
        
        F_dot[2 * tet_n * 12 + 5 * tet_n + i] = B(0, 0);// 2 0
        F_dot[5 * tet_n * 12 + 5 * tet_n + i] = B(0, 1);// 2 1
        F_dot[8 * tet_n * 12 + 5 * tet_n + i] = B(0, 2);// 2 2
        
        
        F_dot[0 * tet_n * 12 + 6 * tet_n + i] = B(1, 0);// 0 0
        F_dot[3 * tet_n * 12 + 6 * tet_n + i] = B(1, 1);// 0 1
        F_dot[6 * tet_n * 12 + 6 * tet_n + i] = B(1, 2);// 0 2
        
        F_dot[1 * tet_n * 12 + 7 * tet_n + i] = B(1, 0);// 1 0
        F_dot[4 * tet_n * 12 + 7 * tet_n + i] = B(1, 1);// 1 1
        F_dot[7 * tet_n * 12 + 7 * tet_n + i] = B(1, 2);// 1 2
        
        F_dot[2 * tet_n * 12 + 8 * tet_n + i] = B(1, 0);// 2 0
        F_dot[5 * tet_n * 12 + 8 * tet_n + i] = B(1, 1);// 2 1
        F_dot[8 * tet_n * 12 + 8 * tet_n + i] = B(1, 2);// 2 2
        
        
        F_dot[0 * tet_n * 12 + 9 * tet_n + i] = B(2, 0);// 0 0
        F_dot[3 * tet_n * 12 + 9 * tet_n + i] = B(2, 1);// 0 1
        F_dot[6 * tet_n * 12 + 9 * tet_n + i] = B(2, 2);// 0 2
        
        F_dot[1 * tet_n * 12 + 10 * tet_n + i] = B(2, 0);// 1 0
        F_dot[4 * tet_n * 12 + 10 * tet_n + i] = B(2, 1);// 1 1
        F_dot[7 * tet_n * 12 + 10 * tet_n + i] = B(2, 2);// 1 2
        
        F_dot[2 * tet_n * 12 + 11 * tet_n + i] = B(2, 0);// 2 0
        F_dot[5 * tet_n * 12 + 11 * tet_n + i] = B(2, 1);// 2 1
        F_dot[8 * tet_n * 12 + 11 * tet_n + i] = B(2, 2);// 2 2
        
        
        // row 1
        double a1 = -1.0 * (B(0, 0) + B(1, 0) + B(2, 0));
        double a2 = a1;
        double b1 = -1.0 * (B(0, 1) + B(1, 1) + B(2, 1));
        double b2 = b1;
        double c1 = -1.0 * (B(0, 2) + B(1, 2) + B(2, 2));
        double c2 = c1;
        row[48 * i + 0] = 3 * tet[0] + 0;
        col[48 * i + 0] = 3 * tet[0] + 0;
        val[48 * i + 0] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a1 = B(0, 0);
        b1 = B(0, 1);
        c1 = B(0, 2);
        row[48 * i + 1] = 3 * tet[0] + 0;
        col[48 * i + 1] = 3 * tet[1] + 0;
        val[48 * i + 1] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a1 = B(1, 0);
        b1 = B(1, 1);
        c1 = B(1, 2);
        row[48 * i + 2] = 3 * tet[0] + 0;
        col[48 * i + 2] = 3 * tet[2] + 0;
        val[48 * i + 2] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a1 = B(2, 0);
        b1 = B(2, 1);
        c1 = B(2, 2);
        row[48 * i + 3] = 3 * tet[0] + 0;
        col[48 * i + 3] = 3 * tet[3] + 0;
        val[48 * i + 3] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        
        // row 2
        row[48 * i + 4] = 3 * tet[0] + 1;
        col[48 * i + 4] = 3 * tet[0] + 1;
        val[48 * i + 4] = (a2 * a2 + b2 * b2 + c2 * c2) * tet_vol;  
        
        a1 = B(0, 0);
        b1 = B(0, 1);
        c1 = B(0, 2);
        row[48 * i + 5] = 3 * tet[0] + 1;
        col[48 * i + 5] = 3 * tet[1] + 1;
        val[48 * i + 5] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a1 = B(1, 0);
        b1 = B(1, 1);
        c1 = B(1, 2);
        row[48 * i + 6] = 3 * tet[0] + 1;
        col[48 * i + 6] = 3 * tet[2] + 1;
        val[48 * i + 6] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a1 = B(2, 0);
        b1 = B(2, 1);
        c1 = B(2, 2);
        row[48 * i + 7] = 3 * tet[0] + 1;
        col[48 * i + 7] = 3 * tet[3] + 1;
        val[48 * i + 7] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        
        // row 3
        row[48 * i + 8] = 3 * tet[0] + 2;
        col[48 * i + 8] = 3 * tet[0] + 2;
        val[48 * i + 8] = (a2 * a2 + b2 * b2 + c2 * c2) * tet_vol;  
        
        a1 = B(0, 0);
        b1 = B(0, 1);
        c1 = B(0, 2);
        row[48 * i + 9] = 3 * tet[0] + 2;
        col[48 * i + 9] = 3 * tet[1] + 2;
        val[48 * i + 9] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a1 = B(1, 0);
        b1 = B(1, 1);
        c1 = B(1, 2);
        row[48 * i + 10] = 3 * tet[0] + 2;
        col[48 * i + 10] = 3 * tet[2] + 2;
        val[48 * i + 10] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a1 = B(2, 0);
        b1 = B(2, 1);
        c1 = B(2, 2);
        row[48 * i + 11] = 3 * tet[0] + 2;
        col[48 * i + 11] = 3 * tet[3] + 2;
        val[48 * i + 11] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        
        
        
        // row 4
        a1 = B(0, 0);
        b1 = B(0, 1);
        c1 = B(0, 2);
        row[48 * i + 12] = 3 * tet[1] + 0;
        col[48 * i + 12] = 3 * tet[0] + 0;
        val[48 * i + 12] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        row[48 * i + 13] = 3 * tet[1] + 0;
        col[48 * i + 13] = 3 * tet[1] + 0;
        val[48 * i + 13] = (a1 * a1 + b1 * b1 + c1 * c1) * tet_vol;  
        
        double a3 = B(1, 0);
        double b3 = B(1, 1);
        double c3 = B(1, 2);
        row[48 * i + 14] = 3 * tet[1] + 0;
        col[48 * i + 14] = 3 * tet[2] + 0;
        val[48 * i + 14] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        a3 = B(2, 0);
        b3 = B(2, 1);
        c3 = B(2, 2);
        row[48 * i + 15] = 3 * tet[1] + 0;
        col[48 * i + 15] = 3 * tet[3] + 0;
        val[48 * i + 15] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        
        
        
        // row 5
        row[48 * i + 16] = 3 * tet[1] + 1;
        col[48 * i + 16] = 3 * tet[0] + 1;
        val[48 * i + 16] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        row[48 * i + 17] = 3 * tet[1] + 1;
        col[48 * i + 17] = 3 * tet[1] + 1;
        val[48 * i + 17] = (a1 * a1 + b1 * b1 + c1 * c1) * tet_vol;  
        
        a3 = B(1, 0);
        b3 = B(1, 1);
        c3 = B(1, 2);
        row[48 * i + 18] = 3 * tet[1] + 1;
        col[48 * i + 18] = 3 * tet[2] + 1;
        val[48 * i + 18] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        a3 = B(2, 0);
        b3 = B(2, 1);
        c3 = B(2, 2);
        row[48 * i + 19] = 3 * tet[1] + 1;
        col[48 * i + 19] = 3 * tet[3] + 1;
        val[48 * i + 19] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        
        
        // row 6
        row[48 * i + 20] = 3 * tet[1] + 2;
        col[48 * i + 20] = 3 * tet[0] + 2;
        val[48 * i + 20] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        row[48 * i + 21] = 3 * tet[1] + 2;
        col[48 * i + 21] = 3 * tet[1] + 2;
        val[48 * i + 21] = (a1 * a1 + b1 * b1 + c1 * c1) * tet_vol;  
        
        a3 = B(1, 0);
        b3 = B(1, 1);
        c3 = B(1, 2);
        row[48 * i + 22] = 3 * tet[1] + 2;
        col[48 * i + 22] = 3 * tet[2] + 2;
        val[48 * i + 22] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        a3 = B(2, 0);
        b3 = B(2, 1);
        c3 = B(2, 2);
        row[48 * i + 23] = 3 * tet[1] + 2;
        col[48 * i + 23] = 3 * tet[3] + 2;
        val[48 * i + 23] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        
        
        // row 7
        a1 = B(1, 0);
        b1 = B(1, 1);
        c1 = B(1, 2);
        row[48 * i + 24] = 3 * tet[2] + 0;
        col[48 * i + 24] = 3 * tet[0] + 0;
        val[48 * i + 24] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a3 = B(0, 0);
        b3 = B(0, 1);
        c3 = B(0, 2);
        row[48 * i + 25] = 3 * tet[2] + 0;
        col[48 * i + 25] = 3 * tet[1] + 0;
        val[48 * i + 25] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        row[48 * i + 26] = 3 * tet[2] + 0;
        col[48 * i + 26] = 3 * tet[2] + 0;
        val[48 * i + 26] = (a1 * a1 + b1 * b1 + c1 * c1) * tet_vol;  
        
        a3 = B(2, 0);
        b3 = B(2, 1);
        c3 = B(2, 2);
        row[48 * i + 27] = 3 * tet[2] + 0;
        col[48 * i + 27] = 3 * tet[3] + 0;
        val[48 * i + 27] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        
        
        
        // row 8
        row[48 * i + 28] = 3 * tet[2] + 1;
        col[48 * i + 28] = 3 * tet[0] + 1;
        val[48 * i + 28] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a3 = B(0, 0);
        b3 = B(0, 1);
        c3 = B(0, 2);
        row[48 * i + 29] = 3 * tet[2] + 1;
        col[48 * i + 29] = 3 * tet[1] + 1;
        val[48 * i + 29] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        row[48 * i + 30] = 3 * tet[2] + 1;
        col[48 * i + 30] = 3 * tet[2] + 1;
        val[48 * i + 30] = (a1 * a1 + b1 * b1 + c1 * c1) * tet_vol;  
        
        a3 = B(2, 0);
        b3 = B(2, 1);
        c3 = B(2, 2);
        row[48 * i + 31] = 3 * tet[2] + 1;
        col[48 * i + 31] = 3 * tet[3] + 1;
        val[48 * i + 31] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        
        
        // row 9
        row[48 * i + 32] = 3 * tet[2] + 2;
        col[48 * i + 32] = 3 * tet[0] + 2;
        val[48 * i + 32] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a3 = B(0, 0);
        b3 = B(0, 1);
        c3 = B(0, 2);
        row[48 * i + 33] = 3 * tet[2] + 2;
        col[48 * i + 33] = 3 * tet[1] + 2;
        val[48 * i + 33] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        row[48 * i + 34] = 3 * tet[2] + 2;
        col[48 * i + 34] = 3 * tet[2] + 2;
        val[48 * i + 34] = (a1 * a1 + b1 * b1 + c1 * c1) * tet_vol;  
        
        a3 = B(2, 0);
        b3 = B(2, 1);
        c3 = B(2, 2);
        row[48 * i + 35] = 3 * tet[2] + 2;
        col[48 * i + 35] = 3 * tet[3] + 2;
        val[48 * i + 35] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        
        
        // row 10
        a1 = B(2, 0);
        b1 = B(2, 1);
        c1 = B(2, 2);
        row[48 * i + 36] = 3 * tet[3] + 0;
        col[48 * i + 36] = 3 * tet[0] + 0;
        val[48 * i + 36] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a3 = B(0, 0);
        b3 = B(0, 1);
        c3 = B(0, 2);
        row[48 * i + 37] = 3 * tet[3] + 0;
        col[48 * i + 37] = 3 * tet[1] + 0;
        val[48 * i + 37] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        a3 = B(1, 0);
        b3 = B(1, 1);
        c3 = B(1, 2);
        row[48 * i + 38] = 3 * tet[3] + 0;
        col[48 * i + 38] = 3 * tet[2] + 0;
        val[48 * i + 38] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        row[48 * i + 39] = 3 * tet[3] + 0;
        col[48 * i + 39] = 3 * tet[3] + 0;
        val[48 * i + 39] = (a1 * a1 + b1 * b1 + c1 * c1) * tet_vol;  
        
        
        
        // row 11
        row[48 * i + 40] = 3 * tet[3] + 1;
        col[48 * i + 40] = 3 * tet[0] + 1;
        val[48 * i + 40] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a3 = B(0, 0);
        b3 = B(0, 1);
        c3 = B(0, 2);
        row[48 * i + 41] = 3 * tet[3] + 1;
        col[48 * i + 41] = 3 * tet[1] + 1;
        val[48 * i + 41] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        a3 = B(1, 0);
        b3 = B(1, 1);
        c3 = B(1, 2);
        row[48 * i + 42] = 3 * tet[3] + 1;
        col[48 * i + 42] = 3 * tet[2] + 1;
        val[48 * i + 42] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        row[48 * i + 43] = 3 * tet[3] + 1;
        col[48 * i + 43] = 3 * tet[3] + 1;
        val[48 * i + 43] = (a1 * a1 + b1 * b1 + c1 * c1) * tet_vol; 
        
        
        
        // row 12
        row[48 * i + 44] = 3 * tet[3] + 2;
        col[48 * i + 44] = 3 * tet[0] + 2;
        val[48 * i + 44] = (a1 * a2 + b1 * b2 + c1 * c2) * tet_vol;  
        
        a3 = B(0, 0);
        b3 = B(0, 1);
        c3 = B(0, 2);
        row[48 * i + 45] = 3 * tet[3] + 2;
        col[48 * i + 45] = 3 * tet[1] + 2;
        val[48 * i + 45] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        a3 = B(1, 0);
        b3 = B(1, 1);
        c3 = B(1, 2);
        row[48 * i + 46] = 3 * tet[3] + 2;
        col[48 * i + 46] = 3 * tet[2] + 2;
        val[48 * i + 46] = (a1 * a3 + b1 * b3 + c1 * c3) * tet_vol;  
        
        row[48 * i + 47] = 3 * tet[3] + 2;
        col[48 * i + 47] = 3 * tet[3] + 2;
        val[48 * i + 47] = (a1 * a1 + b1 * b1 + c1 * c1) * tet_vol; 
     
    }
     
    return;
    
}