/*
============================================================================================
   This class accesses a MATLAB matrix (read only).

   Copyright (c) 01-03-2014,  Shawn W. Walker
============================================================================================
*/

#define  MRO   MATLAB_Matrix_ReadOnly

/* C++ class for accessing data in MATLAB matrices */
template <class myType, unsigned int NUM_COLS>
class MRO
{
public:

    MRO ();  // constructor
    ~MRO (); // DE-structor
    unsigned int Get_Num_Rows ()      { return num_row; }
    unsigned int Get_Num_Cols ()      { return num_col; }
	
    void Setup_Data(const mxArray*);
	void Read(const unsigned int& ri, myType rd[NUM_COLS]);
	const myType Read(const unsigned int&, const unsigned int&);
	const myType* Get_Data_Col_Ptr(const unsigned int&);

private:
    unsigned int   num_row;        // number of rows    (i.e. M)
	unsigned int   num_col;        // number of columns (i.e. C)
    const myType*  Data[NUM_COLS]; // access to MxC data (read-only)
};

/***************************************************************************************/
/* constructor */
template <class myType, unsigned int NUM_COLS>
MRO<myType,NUM_COLS>::MRO ()
{
    num_row = 0; // set later in Setup_Data
    num_col = NUM_COLS;

    // init data information to NULL
    for (unsigned int ii = 0; (ii < num_col); ii++)
        Data[ii] = NULL;
}
/***************************************************************************************/


/***************************************************************************************/
/* DE-structor */
template <class myType, unsigned int NUM_COLS>
MRO<myType,NUM_COLS>::~MRO ()
{
}
/***************************************************************************************/


/***************************************************************************************/
/* put incoming data from MATLAB into a nice struct */
template <class myType, unsigned int NUM_COLS>
void MRO<myType,NUM_COLS>::Setup_Data(const mxArray* mxData) // inputs
{
    // get the number of rows
    num_row = (unsigned int) mxGetM(mxData);
    //num_col = (unsigned int) mxGetN(mxData);
    // get the number of cols
    unsigned int Actual_Num_Col = (unsigned int) mxGetN(mxData);


    /* BEGIN: Simple Error Checking */
    if (typeid(myType) == typeid(unsigned int))
        {
        if (mxGetClassID(mxData)!=mxUINT32_CLASS) mexErrMsgTxt("ERROR: Data must be of type uint32!");
        }
    else if (typeid(myType) == typeid(int))
        {
        if (mxGetClassID(mxData)!=mxINT32_CLASS) mexErrMsgTxt("ERROR: Data must be of type int32!");
        }
    else if (typeid(myType) == typeid(double))
        {
        if (mxGetClassID(mxData)!=mxDOUBLE_CLASS) mexErrMsgTxt("ERROR: Data must be of type double!");
        }
    else if (typeid(myType) == typeid(bool))
        {
        if (mxGetClassID(mxData)!=mxLOGICAL_CLASS) mexErrMsgTxt("ERROR: Data must be of type logical!");
        }
    else
        {
        mexErrMsgTxt("ERROR: Data type not recognized; C++ class must be extended!");
        }
    if (Actual_Num_Col != num_col)
        {
        mexPrintf("ERROR: Actual data has %d columns; was expecting %d columns.\n", Actual_Num_Col, num_col);
        mexErrMsgTxt("ERROR: fix your Data!");
        }
    /* END: Simple Error Checking */


    // split up the columns of the element data
    Data[0] = (const myType*) (mxGetPr(mxData));
	Data[0] = Data[0] - 1; // MATLAB style indexing
    for (unsigned int ii = 1; (ii < num_col); ii++)
        Data[ii] = Data[ii-1] + num_row;
}
/***************************************************************************************/


/***************************************************************************************/
/* read a row of the matrix */
template <class myType, unsigned int NUM_COLS>
void MRO<myType,NUM_COLS>::Read(const unsigned int& row_index,  // input
                                myType row_data[NUM_COLS])     // output
{
	if ((row_index < 1) || (row_index > num_row))
		{
		mexPrintf("The given row index must be between %d and %d.\n",1,num_row);
		mexErrMsgTxt("STOP!\n");
		}

    // read one row of the matrix
    for (unsigned int ii = 0; (ii < num_col); ii++)
        row_data[ii] = Data[ii][row_index];
}
/***************************************************************************************/


/***************************************************************************************/
/* read the i,j element of the matrix */
template <class myType, unsigned int NUM_COLS>
const myType MRO<myType,NUM_COLS>::Read(const unsigned int& row_index,  // input
                                        const unsigned int& col_index)  // input
{
	if ((row_index < 1) || (row_index > num_row))
		{
		mexPrintf("The given row index must be between %d and %d.\n",1,num_row);
		mexErrMsgTxt("STOP!\n");
		}
	if ((col_index < 1) || (col_index > num_col))
		{
		mexPrintf("The given column index must be between %d and %d.\n",1,num_col);
		mexErrMsgTxt("STOP!\n");
		}
		
    // read one element of the matrix
    const myType VAL = Data[col_index-1][row_index];
	return VAL;
}
/***************************************************************************************/


/***************************************************************************************/
/* get pointer to a column of the data */
template <class myType, unsigned int NUM_COLS>
const myType* MRO<myType,NUM_COLS>::Get_Data_Col_Ptr(const unsigned int& col_index)  // input
{
	if ((col_index < 1) || (col_index > num_col))
		{
		mexPrintf("The given column index must be between %d and %d.\n",1,num_col);
		mexErrMsgTxt("STOP!\n");
		}
		
    // get ptr
	return Data[col_index-1];
}
/***************************************************************************************/



#undef MRO

/***/
