/*
============================================================================================
   This class accesses a MATLAB matrix (read and writeable).

   Copyright (c) 01-08-2014,  Shawn W. Walker
============================================================================================
*/

#define  MRW       MATLAB_Matrix_ReadWrite
#define  NUM_COLS  10000

/* C++ class for accessing and writing data to MATLAB matrices */
template <class myType>
class MRW
{
public:

    MRW ();  // constructor
    ~MRW (); // DE-structor
    unsigned int Get_Num_Rows ()      { return num_row; }
    unsigned int Get_Num_Cols ()      { return num_col; }
	
    void Setup_Data(mxArray*);
	void Read(const unsigned int&, myType*);
	const myType Read(const unsigned int&, const unsigned int&);
	myType* Get_Data_Col_Ptr(const unsigned int&);
	void Write(const unsigned int&, const myType*);
	void Write(const unsigned int&, const unsigned int&, const myType&);

private:
    unsigned int   num_row;   // number of rows    (i.e. M)
	unsigned int   num_col;   // number of columns (i.e. C)
    myType*        Data[NUM_COLS]; // access to Mx? data (read-write)
};

/***************************************************************************************/
/* constructor */
template <class myType>
MRW<myType>::MRW ()
{
    num_row = 0; // set later in Setup_Data
    num_col = 0;

    // init data information to NULL
    for (unsigned int ii = 0; (ii < NUM_COLS); ii++)
        Data[ii] = NULL;
}
/***************************************************************************************/


/***************************************************************************************/
/* DE-structor */
template <class myType>
MRW<myType>::~MRW ()
{
}
/***************************************************************************************/


/***************************************************************************************/
/* put incoming data from MATLAB into a nice struct */
template <class myType>
void MRW<myType>::Setup_Data(mxArray* mxData) // inputs
{
    // get the number of rows
    num_row = (unsigned int) mxGetM(mxData);
    // get the number of cols
    num_col = (unsigned int) mxGetN(mxData);


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
    if (num_col > NUM_COLS)
        {
        mexPrintf("ERROR: Actual data has %d columns; cannot have more than %d columns.\n", num_col, NUM_COLS);
        mexErrMsgTxt("ERROR: fix your Data!");
        }
    /* END: Simple Error Checking */


    // split up the columns of the element data
    Data[0] = (myType*) (mxGetPr(mxData));
	Data[0] = Data[0] - 1; // MATLAB style indexing
    for (unsigned int ii = 1; (ii < num_col); ii++)
        Data[ii] = Data[ii-1] + num_row;
}
/***************************************************************************************/


/***************************************************************************************/
/* read a row of the matrix */
template <class myType>
void MRW<myType>::Read(const unsigned int& row_index,  // input
                       myType* row_data)               // output
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
template <class myType>
const myType MRW<myType>::Read(const unsigned int& row_index,  // input
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
template <class myType>
myType* MRW<myType>::Get_Data_Col_Ptr(const unsigned int& col_index)  // input
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


/***************************************************************************************/
/* write to a row of the matrix */
template <class myType>
void MRW<myType>::Write(const unsigned int& row_index,   // input
                        const myType* row_data)          // input
{
	if ((row_index < 1) || (row_index > num_row))
		{
		mexPrintf("The given row index must be between %d and %d.\n",1,num_row);
		mexErrMsgTxt("STOP!\n");
		}

    // write to one row of the matrix
    for (unsigned int ii = 0; (ii < num_col); ii++)
		Data[ii][row_index] = row_data[ii];
}
/***************************************************************************************/


/***************************************************************************************/
/* write to the i,j element of the matrix */
template <class myType>
void MRW<myType>::Write(const unsigned int& row_index,  // input
                        const unsigned int& col_index,  // input
                        const myType&       VAL)        // input
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

    // write to one element of the matrix
	Data[col_index-1][row_index] = VAL;
}
/***************************************************************************************/

#undef MRW
#undef NUM_COLS 

/***/
