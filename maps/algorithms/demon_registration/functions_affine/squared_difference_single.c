#include "mex.h"
#include "math.h"

/* squared pixel difference */

/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    float *I1, *I2, *Ierror;
     
    /* Size of input image */
    const mwSize *idims; 
    /*  Dimensions */
    int nsubs;
    /*   number of pixels */
    int npixels=1;
    /* Loop variable  */
    int i;
    /*  Size of output */
    int odims2[2]={1,1};
    /*  Total squared error */
    float ierror=0;
    
    /*  Connect inputs */
    I1=(float *)mxGetData(prhs[0]);
    I2=(float *)mxGetData(prhs[1]);
    
    /*  Initialize output array */
    plhs[0] = mxCreateNumericArray(2, odims2, mxSINGLE_CLASS, mxREAL);
    Ierror=(float *)mxGetData(plhs[0]);
    
    /*  Get the number of dimensions */
    nsubs = mxGetNumberOfDimensions(prhs[0]);
    /* Get the sizes of the grid */
    idims = mxGetDimensions(prhs[0]);   
    for (i=0; i<nsubs; i++) { npixels=npixels*idims[i]; }
    /*  Loop trough all pixels to calculate squared pixel error */
    for (i=0; i<npixels; i++) 
    { 
        ierror+=(I1[i]-I2[i])*(I1[i]-I2[i]);
    }
    /*  Error to output */
    Ierror[0]=ierror/((float)(npixels));;
}
        
