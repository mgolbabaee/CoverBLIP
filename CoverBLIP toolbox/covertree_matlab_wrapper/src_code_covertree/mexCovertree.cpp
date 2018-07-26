#include <math.h>
#include <mex.h> 

// define input indices for creating the covertree
#define INPUT_Command_String          0
#define INPUT_Points                           1

/* include classes and other sub-routines */
#include "../src_hdr/class_handle.hpp"
#include "../src_hdr/MATLAB_Matrix_ReadOnly.cc"
#include "../src_hdr/MATLAB_Matrix_ReadWrite.cc"
#include "covertree.cc"

//#include <omp.h>

// note: 'prhs' represents the Right-Hand-Side arguments from MATLAB (inputs)
//       'plhs' represents the  Left-Hand-Side arguments from MATLAB (outputs)

/***************************************************************************************/
// define the "gateway" function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // get the command string
    char cmd[64];
    if (nrhs < 1 || mxGetString(prhs[INPUT_Command_String], cmd, sizeof(cmd)))
        {
        mexPrintf("mexCovertree:\n");
        mexPrintf("\n");
        mexPrintf("ERROR: first input should be a command string less than 64 characters long.\n");
        mexPrintf("\n");
        mexPrintf("       List of Acceptable Command Strings:\n");
        mexPrintf("       -----------------------------------\n");
        mexPrintf("       'new'\n");
        mexPrintf("       'delete'\n");
        mexPrintf("       'enn_search_batch'\n");
        mexPrintf("       'enn_search_loop\n");
        mexPrintf("\n");
        mexPrintf("       The number of additional inputs depends on the command string.\n");
        mexErrMsgTxt("Check your inputs!");
        }

    // if we are creating a new instance (and, thus, building the tree), then
    if (!strcmp("new", cmd))
        {
        /* BEGIN: Error Checking */
        if ( (nrhs>2))
            {
            mexPrintf("mexCovertree: 'new'\n");
            mexPrintf("\n");
            mexPrintf("ERROR: need 2 input and 1 output!\n");
            mexPrintf("\n");
            mexPrintf("      INPUTS ORDER \n");
            mexPrintf("      --------------------------------------------------------------------------- \n");
            mexPrintf("      Command String = 'new'                                               0 \n");
            mexPrintf("      Dict (MxN matrix of dictrionary with M atoms)               1 \n");
            mexPrintf("\n");
            mexPrintf("      OUTPUT\n");
            mexPrintf("      --------------------------------------------------------------------------- \n");
            mexPrintf("      'Pointer' (handle) to the C++ quadtree class instance      \n");
            mexPrintf("      Maximum tree levels                                                       \n");
            mexPrintf("\n");
            mexErrMsgTxt("Check the arguments!");
            }
        /* END: Error Checking */

        // declare the object
        covertree*  Cover_Tree_Obj;
        Cover_Tree_Obj = new covertree;

        /* pass the inputs */

        // setup input point data
        Cover_Tree_Obj->Get_Dict_From_Matrix(prhs[INPUT_Points]);

        // build the tree
        Cover_Tree_Obj->Build_Tree();

        // Return a handle to the new C++ instance
        plhs[0] = convertPtr2Mat<covertree>(Cover_Tree_Obj);
        
        return;
        }
    // Delete
    if (!strcmp("delete", cmd))
        {     
        // Destroy the C++ object
        destroyObject<covertree>(prhs[1]);
        // Warn if other commands were ignored
        if (nlhs > 0 || nrhs > 2)
           mexWarnMsgTxt("delete: does not require additional inputs or any outputs.");
        return;
        }  
    // otherwise, the object was created before, so
    //    get the C++ covertree instance pointer from the second input
    covertree* Cover_Tree_Obj = convertMat2Ptr<covertree>(prhs[1]);
    
    if (!strcmp("print_tree_structure", cmd))
    {
        // Check parameters
        if ((nlhs<1) || (nlhs>1) || (nrhs<2) || (nrhs>2))
            {
            mexPrintf("mexCovertree: 'enn_search_batch'\n");
            mexPrintf("\n");
            mexPrintf("ERROR: need 2 inputs and 1 output!\n");
            mexPrintf("\n");
            mexPrintf("      INPUTS ORDER \n");
            mexPrintf("      ---------------------------------------------------------------------------------------- \n");
            mexPrintf("      Command String = 'print_tree_structure'                                          0 \n");
            mexPrintf("      'Pointer' (handle) to the C++ covertree class instance                  1 \n");
            mexPrintf("\n");
            mexPrintf("      OUTPUTS\n");
            mexPrintf("      ------------------------------------------------- \n");
            mexPrintf("      Tree structure (maximum num of levels x num of atoms)     \n");
            mexPrintf("\n");
            mexErrMsgTxt("Check the arguments!");
            }
        int max_tree_level, num_atom;
        max_tree_level = Cover_Tree_Obj->Get_Max_Tree_Level_From_Obj();
        num_atom = Cover_Tree_Obj->Get_Num_Atom_From_Obj();
        plhs[0] = mxCreateDoubleMatrix(max_tree_level,num_atom,mxREAL);        
        Cover_Tree_Obj->print_tree_structure(plhs[0]);      
        return;

    }
 
    // perform multiple e-nearest neighbor batch search in parallel
    if (!strcmp("enn_search_batch", cmd))
        {
        // Check parameters
        if ((nlhs<5) || (nlhs>5) || (nrhs<5) || (nrhs>5))
            {
            mexPrintf("mexCovertree: 'enn_search_batch'\n");
            mexPrintf("\n");
            mexPrintf("ERROR: need 5 inputs and 5 outputs!\n");
            mexPrintf("\n");
            mexPrintf("      INPUTS ORDER \n");
            mexPrintf("      ---------------------------------------------------------------------------------------- \n");
            mexPrintf("      Command String = 'enn_search_batch'                                          0 \n");
            mexPrintf("      'Pointer' (handle) to the C++ covertree class instance                  1 \n");
            mexPrintf("      Enquery Points (d*N matrix, d enqueries)                                     2 \n");
            mexPrintf("      Epsilon                                                                                         3 \n");
            mexPrintf("      Batch size                                                                                    4 \n");          
            mexPrintf("\n");
            mexPrintf("      OUTPUTS\n");
            mexPrintf("      ------------------------------------------------- \n");
            mexPrintf("      Labels of Approximate Nearest Atoms in Covertree (d labels)\n");
            mexPrintf("      Approximate Nearest Atoms                                                 \n");
            mexPrintf("      Mininum distance                                                                 \n");            
            mexPrintf("      Distance calculation counter                                                \n");
            mexPrintf("      Index                                                      ");           
            mexPrintf("\n");
            mexErrMsgTxt("Check the arguments!");
            }
        plhs[0] = mxCreateDoubleMatrix(1,mxGetM(prhs[2]),mxREAL);
        plhs[1] = mxCreateDoubleMatrix(mxGetN(prhs[2]),mxGetM(prhs[2]),mxCOMPLEX);
        plhs[2] = mxCreateDoubleMatrix(1,mxGetM(prhs[2]),mxREAL);
        plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
        plhs[4] = mxCreateDoubleMatrix(1,mxGetM(prhs[2]),mxREAL);
        Cover_Tree_Obj->Get_Epsilon_Aprox_From_Matrix(prhs[3]);
        Cover_Tree_Obj->eNN_bsearch(prhs[4],prhs[2],plhs[0],plhs[1],plhs[2],plhs[3], plhs[4]);
        return;
        }  
    // perform a e-nearest neighbor batch search with current estimate
      if (!strcmp("enn_loop_current", cmd))
        {
        // Check parameters
        if ((nlhs<4) || (nlhs>4) || (nrhs<4) || (nrhs>5))
            {
            mexPrintf("mexCovertree: 'enn_bsearch_current'\n");
            mexPrintf("\n");
            mexPrintf("ERROR: need 5 inputs and 4 outputs!\n");
            mexPrintf("\n");
            mexPrintf("      INPUTS ORDER \n");
            mexPrintf("      --------------------------------------------------------------------------------------- \n");
            mexPrintf("      Command String = 'enn_loop_current'                                          0 \n");
            mexPrintf("      'Pointer' (handle) to the C++ covertree class instance                  1 \n");
            mexPrintf("      Enquery Points (d*N matrix)                                                         2 \n");
            mexPrintf("      Epsilon                                                                                          3 \n");
            mexPrintf("      Upperbound                                                                                  4 \n");
            mexPrintf("\n");
            mexPrintf("      OUTPUTS\n");
            mexPrintf("      ------------------------------------------------- \n");
            mexPrintf("      Labels of Approximate Nearest Atoms in Covertree (d labels)\n");
            mexPrintf("      Approximate Nearest Atoms             \n");
            mexPrintf("      Mininum distance                      \n");            
            mexPrintf("      Distance calculation counter           \n");
            mexPrintf("\n");
            mexErrMsgTxt("Check the arguments!");
            }
        plhs[0] = mxCreateDoubleMatrix(1,mxGetM(prhs[2]),mxREAL);
        plhs[1] = mxCreateDoubleMatrix(mxGetN(prhs[2]),mxGetM(prhs[2]),mxCOMPLEX);
        plhs[2] = mxCreateDoubleMatrix(1,mxGetM(prhs[2]),mxREAL);
        plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
        Cover_Tree_Obj->Get_Epsilon_Aprox_From_Matrix(prhs[3]);
        Cover_Tree_Obj->eNN_loop_current(prhs[4],prhs[2],plhs[0],plhs[1],plhs[2],plhs[3]);
        return;
      }    
    
      if (!strcmp("enn_stop_level", cmd))
        {
        // Check parameters
        if ((nlhs<4) || (nlhs>4) || (nrhs<5) || (nrhs>5))
            {
            mexPrintf("mexCovertree: 'enn_stop_level'\n");
            mexPrintf("\n");
            mexPrintf("ERROR: need 4 inputs and 4 outputs!\n");
            mexPrintf("\n");
            mexPrintf("      INPUTS ORDER \n");
            mexPrintf("      --------------------------------------------------------------------------------------- \n");
            mexPrintf("      Command String = 'enn_stop_level'                                              0 \n");
            mexPrintf("      'Pointer' (handle) to the C++ covertree class instance                  1 \n");
            mexPrintf("      Enquery Points (d*N matrix)                                                         2 \n");
            mexPrintf("      Epsilon                                                                                          3 \n");
            mexPrintf("      Stop level                                                                                      4 \n");
            mexPrintf("\n");
            mexPrintf("      OUTPUTS\n");
            mexPrintf("      ------------------------------------------------- \n");
            mexPrintf("      Labels of Approximate Nearest Atoms in Covertree (d labels)\n");
            mexPrintf("      Approximate Nearest Atoms             \n");
            mexPrintf("      Mininum distance                      \n");            
            mexPrintf("      Distance calculation counter           \n");
            mexPrintf("\n");
            mexErrMsgTxt("Check the arguments!");
            }
        plhs[0] = mxCreateDoubleMatrix(1,mxGetM(prhs[2]),mxREAL);
        plhs[1] = mxCreateDoubleMatrix(mxGetN(prhs[2]),mxGetM(prhs[2]),mxCOMPLEX);
        plhs[2] = mxCreateDoubleMatrix(1,mxGetM(prhs[2]),mxREAL);
        plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
        Cover_Tree_Obj->Get_Epsilon_Aprox_From_Matrix(prhs[3]);               
        Cover_Tree_Obj->Get_Stop_Scale_From_Matrix(prhs[4]);
        Cover_Tree_Obj->eNN_stop_scale(prhs[2],plhs[0],plhs[1],plhs[2],plhs[3]);
        return;
      }      
    
    // if we got here, then the command was not recognized
    mexPrintf("\n");
    mexPrintf("This command was not recognized:  ");
    mexPrintf(cmd);
    mexPrintf("\n\n");
    mexErrMsgTxt("Call this mex file with no inputs to see a list of acceptable command strings.\n");
}

/***/
