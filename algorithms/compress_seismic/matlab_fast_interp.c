#include "mex.h"
#include "matrix.h"
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    const mxArray *xi, *xgrid;
    mxArray *idx;
    size_t nx, m, k, i1, i9, imid;
    double *xiptr, *yptr, *xgridptr, *idxptr;
    double xik;
    mwSize dims[2];
    
    // Get inputs and dimensions
    xgrid = prhs[0];
    nx = mxGetM(xgrid); 
    xi = prhs[2];
    m = mxGetM(xi);
    
    // Create output idx 
    dims[0] = m; dims[1] = 1;
    plhs[0] = idx = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    if (idx==NULL) // Cannot allocate memory
    {
        // Return empty array
        dims[0] = 0; dims[1] = 0;
        plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
        return;
    }
    idxptr = mxGetPr(idx);
    
    // Get pointers
    xiptr = mxGetPr(xi);
    yptr = mxGetPr(prhs[1]);
    xgridptr = mxGetPr(xgrid);
  
    // Loop over the points
    for (k=m; k--;) // Reverse of for (k=0; k<m; k++) {...}
    {
        // Get data value
        xik = xiptr[k];
        
        i1=0;
        i9=nx-1; 
        while (i9>i1+1) // Dichotomy search
        {
            imid = (i1+i9+1)/2;
            if (xgridptr[imid]<xik) i1=imid;
            else i9=imid;
        } // of while loop
        if (i1==i9)
            idxptr[k] = yptr[i1];
        else
            idxptr[k] = yptr[i1] + (yptr[i9]-yptr[i1])*(xik-xgridptr[i1])/(xgridptr[i9]-xgridptr[i1]);
    } // for loop
    
    return;
}
