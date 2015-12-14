#include "mex.h"
#include "math.h"
/*   undef needed for LCC compiler  */
#undef EXTERN_C
#ifdef _WIN32
	#include <windows.h>
	#include <process.h>
#else
	#include <pthread.h>
#endif
        
/* This function transforms two pictures with different modalities in to the
 * grey levels of the opposite modality. The function divides the images in
 * small overlapping regions, and calculates mutual histograms from the regions
 * A (gaussian) kernel is used so that pixels far away from the region counts less
 * in the histogram than nearby pixels. The mutual histograms are then used to find the maximum
 * correlation between the intensity of the pixels in image 1 and image 2. This maximum
 * correlation is then used to paint a new picture of the opposite modality.
 *
 * [I1_TF,I2_TF]=mutual_transform_2d_double(I1,I2,I1_moved,I2_moved,range,HistogramBins,Hkernel, SampleDistance)
 *
 * inputs,
 *  I1 : Input image 1
 *  I2 : Input image 2, with different modality then I1
 *  I1_moved : A geometric transformed picture of I2 (from registration), 
 *          used to modality transform I2.
 *  I2_moved : A geometric transformed picture of I2 (from registration), 
 *          used to modality transform I1.
 *  range : The range of the pixel values for instance 0..1
 *  HistogramBins: The number of histogrambins used, for instance 256
 *  Hkernel : A 2D gaussian kernel of a certain size
 *  SampleDistance: The distance between the spacing of the 
 *          kernels  / image regions
 *
 * outputs,
 *  I1_TF : Modality transformed painting of image 1
 *  I2_TF : Modality transformed painting of image 2
 *
 * Function is written by D. Kroon University of Twente(March 2009)
 */

 
int mindex2(int x, int y, int sizx) { return y*sizx+x; }
int mindex3(int x, int y, int z, int sizx, int sizy) { return z*sizy*sizx+y*sizx+x; }
double round(double x) { return floor(x+0.5); }

#ifdef _WIN32
  unsigned __stdcall transformvolume(double **Args) {
#else
  void transformvolume(double **Args) {
#endif
    /* I1 and I2 are both input images with different modality     */
    double *I1, *I2;
    /* I1_moved and I2_moved are position transformed images (registration results) */
    double *I1_moved, *I2_moved;
    /* Range is a vector with the range of the input image */
    double *range;
    /* HistogramBins is the number of bins in the used histogram */
    double *HistogramBins_d;
    double HistogramBins;
    /* Hkernel is the kernel used to multiply with the Histogram update */
    /* for smooth modality transformations */
    double *Hkernel;
    /* Is the distance between the histogram kernels used to determine  */
    /* the modality transformation */
    double *SampleDistance_d;
    double SampleDistance;
    double halfSample;
    double *Histograms_sizes, *kernel_sizes, *image_sizes;
    /* Number of threads */
    double *Nthreadsd;
    int Nthreads;
    double *ThreadID;
    
    /* Pixel index */
    int index;
    
    /* Part done by this thread */
    int start, end;
    
    /* Scale to histogram */
    double scale;
    
    /* Find maximum variables */
    double maxf; 
    /*double maxx;  */
    
    /* For loops variables */
    int i;
    
    /* Linear interpolation variables */
    int xBas[2], yBas[2];
    double perc[4];
    double xCom, yCom;
    
    /* Intensity transformed images */
    double *I1_TF, *I2_TF;
   
    /* Size of input image */
    int image_sizex, image_sizey, image_sizez;
    int image_x, image_y, image_z;
    double image_x_center, image_y_center, image_z_center;
    int image_x_start, image_y_start, image_z_start;
    int image_x_end, image_y_end, image_z_end;
 	
    /* Size of histogram image  */
    int Histograms_size_x, Histograms_size_y, Histograms_size_z;
    int hist_x, hist_y, hist_z;
    
    /* The 2D histograms for I1 on I2 transformed and visa versa */
    double *Histogram_I1, *Histogram_I2;
    double hist2d_x_d, hist2d_y_d;
    int hist2d_x, hist2d_y;
    
	/* Size of kernel */
	int kernel_sizex, kernel_sizey, kernel_sizez;
    int kernel_x, kernel_y, kernel_z;
    int kernel_x_start, kernel_y_start, kernel_z_start;
    int kernel_x_end, kernel_y_end, kernel_z_end;
	int halfKernelx, halfKernely, halfKernelz;

    /* Get the input data */
    I1=Args[0];
    I1_moved=Args[1];
    I2=Args[2];
    I2_moved=Args[3];
    Hkernel=Args[4];
    Histograms_sizes=Args[5]; 
    HistogramBins_d=Args[6];
    SampleDistance_d=Args[7];
    kernel_sizes=Args[8];
    image_sizes=Args[9];
    range=Args[10];
    I1_TF=Args[11];
    I2_TF=Args[12];
    ThreadID=Args[13];
    Nthreadsd=Args[14];
    
    HistogramBins=HistogramBins_d[0];
    SampleDistance=SampleDistance_d[0];
     
    Nthreads=(int)Nthreadsd[0];
    Histograms_size_x=(int)Histograms_sizes[0];
    Histograms_size_y=(int)Histograms_sizes[1];
    Histograms_size_z=(int)Histograms_sizes[2];
   
    kernel_sizex=(int)kernel_sizes[0];
    kernel_sizey=(int)kernel_sizes[1];
    kernel_sizez=(int)kernel_sizes[2];

    image_sizex=(int)image_sizes[0];
    image_sizey=(int)image_sizes[1];
    image_sizez=(int)image_sizes[2];
       
    /* Allocate 2D histogram memory */
    Histogram_I1=(double *)malloc((int)(HistogramBins*HistogramBins*sizeof(double)));  
    Histogram_I2=(double *)malloc((int)(HistogramBins*HistogramBins*sizeof(double)));  

    /* Set intensity to histogram scaling */
    scale=(HistogramBins-1)/(range[1]-range[0]);
    
    /* Half variables */
    halfKernelx=(kernel_sizex-1)/2;
    halfKernely=(kernel_sizey-1)/2;
    halfKernelz=(kernel_sizez-1)/2;
    
    halfSample=(SampleDistance-1)/2;
    
    start=(int)((double)ThreadID[0]*(((double)Histograms_size_z)/((double)Nthreads)));
    end=(int)((double)(ThreadID[0]+1)*(((double)Histograms_size_z)/((double)Nthreads)));
    
    for (hist_z=start; hist_z<end; hist_z++)
    {
        for (hist_y=0; hist_y<Histograms_size_y; hist_y++)
        {
            for (hist_x=0; hist_x<Histograms_size_x; hist_x++)
            {
                /* Calculate the image coordinates of a location histogram */
                image_x_center=hist_x*SampleDistance+halfSample;
                image_y_center=hist_y*SampleDistance+halfSample;
                image_z_center=hist_z*SampleDistance+halfSample;

                /* Check boundary conditions */
                if(image_x_center>(image_sizex-1)) image_x_center = image_sizex-1;
                if(image_y_center>(image_sizey-1)) image_y_center = image_sizey-1;
                if(image_z_center>(image_sizez-1)) image_z_center = image_sizez-1;
                
                /* Calculate the part of the kernel inside the image */
                kernel_x_start=(int)(halfKernelx-image_x_center); 
                kernel_y_start=(int)(halfKernely-image_y_center); 
                kernel_z_start=(int)(halfKernelz-image_z_center); 
                
                kernel_x_end=(int)(image_sizex-1+halfKernelx-image_x_center); 
                kernel_y_end=(int)(image_sizey-1+halfKernely-image_y_center); 
                kernel_z_end=(int)(image_sizez-1+halfKernelz-image_z_center); 

                /* Kernel Boundary conditions */
                if (kernel_x_start<0) kernel_x_start=0;
                if (kernel_y_start<0) kernel_y_start=0;
                if (kernel_z_start<0) kernel_z_start=0;
                
                if (kernel_x_end>(kernel_sizex-1)) kernel_x_end = kernel_sizex-1;
                if (kernel_y_end>(kernel_sizey-1)) kernel_y_end = kernel_sizey-1;
                if (kernel_z_end>(kernel_sizez-1)) kernel_z_end = kernel_sizez-1;

                /* Image Region which uses this kernel */
                image_x_start=(int)(image_x_center-halfSample);
                image_x_end=(int)(image_x_center+halfSample); 
                image_y_start=(int)(image_y_center-halfSample); 
                image_y_end=(int)(image_y_center+halfSample); 
                image_z_start=(int)(image_z_center-halfSample); 
                image_z_end=(int)(image_z_center+halfSample); 

                /* Boundary conditions on the image region */
                if (image_x_start<0) image_x_start=0;
                if (image_x_end>(image_sizex-1)) image_x_end = image_sizex-1;
                if (image_y_start<0) image_y_start=0;
                if (image_y_end>(image_sizey-1)) image_y_end = image_sizey-1;
                if (image_z_start<0) image_z_start=0;
                if (image_z_end>(image_sizez-1)) image_z_end = image_sizez-1;

                /* Empty histogram for I1 on I2 transformed and visa versa */
                for (i=0; i<(HistogramBins*HistogramBins); i++)
                    { Histogram_I1[i]=0; Histogram_I2[i]=0; }    

                /* Loop through the kernel */
                for (kernel_z=kernel_z_start; kernel_z<(kernel_z_end+1); kernel_z++)
                {
                    for (kernel_y=kernel_y_start; kernel_y<(kernel_y_end+1); kernel_y++)
                    {
                        for (kernel_x=kernel_x_start; kernel_x<(kernel_x_end+1); kernel_x++)
                        {
                             /* The current image coordinates */
                             image_x=(int)(image_x_center+kernel_x-halfKernelx);
                             image_y=(int)(image_y_center+kernel_y-halfKernely);
                             image_z=(int)(image_z_center+kernel_z-halfKernelz);

                             /* pixel index  */
                             index = mindex3(image_x,image_y,image_z,image_sizex,image_sizey);

                             /* The current histogram location */
                             hist2d_x_d=scale*(I1[index]-range[0]);
                             hist2d_y_d=scale*(I2_moved[index]-range[0]);

                             /* Pixel Neighbors */
                             xBas[0]=(int) floor(hist2d_x_d); 
                             yBas[0]=(int) floor(hist2d_y_d);
                             xBas[1]=xBas[0]+1; if(xBas[1]>(HistogramBins-1)) xBas[1]=(int)(HistogramBins-1);      
                             yBas[1]=yBas[0]+1; if(yBas[1]>(HistogramBins-1)) yBas[1]=(int)(HistogramBins-1); 

                             /* Linear interpolation constants (percentages) */
                             xCom=hist2d_x_d-floor(hist2d_x_d); yCom=hist2d_y_d-floor(hist2d_y_d);
                             perc[0]=(1-xCom) * (1-yCom);
                             perc[1]=(1-xCom) * yCom;
                             perc[2]=xCom * (1-yCom);
                             perc[3]=xCom * yCom;

                             /* Update the 2D histogram with the kernel value */
                             hist2d_x=xBas[0]; hist2d_y=yBas[0];
                             Histogram_I1[mindex2(hist2d_x,hist2d_y,(int)HistogramBins)]+=perc[0]*Hkernel[mindex3(kernel_x,kernel_y,kernel_z,kernel_sizex,kernel_sizey)];
                             hist2d_x=xBas[0]; hist2d_y=yBas[1];
                             Histogram_I1[mindex2(hist2d_x,hist2d_y,(int)HistogramBins)]+=perc[1]*Hkernel[mindex3(kernel_x,kernel_y,kernel_z,kernel_sizex,kernel_sizey)];
                             hist2d_x=xBas[1]; hist2d_y=yBas[0];
                             Histogram_I1[mindex2(hist2d_x,hist2d_y,(int)HistogramBins)]+=perc[2]*Hkernel[mindex3(kernel_x,kernel_y,kernel_z,kernel_sizex,kernel_sizey)];
                             hist2d_x=xBas[1]; hist2d_y=yBas[1];
                             Histogram_I1[mindex2(hist2d_x,hist2d_y,(int)HistogramBins)]+=perc[3]*Hkernel[mindex3(kernel_x,kernel_y,kernel_z,kernel_sizex,kernel_sizey)];

                             /* The current histogram location */
                             hist2d_x_d=scale*(I1_moved[index]-range[0]);
                             hist2d_y_d=scale*(I2[index]-range[0]);

                             /* Pixel Neighbors */
                             xBas[0]=(int) floor(hist2d_x_d); 
                             yBas[0]=(int) floor(hist2d_y_d);
                             xBas[1]=xBas[0]+1; if(xBas[1]>(HistogramBins-1)) xBas[1]=(int)(HistogramBins-1);      
                             yBas[1]=yBas[0]+1; if(yBas[1]>(HistogramBins-1)) yBas[1]=(int)(HistogramBins-1); 

                             /* Linear interpolation constants (percentages) */
                             xCom=hist2d_x_d-floor(hist2d_x_d); yCom=hist2d_y_d-floor(hist2d_y_d);
                             perc[0]=(1-xCom) * (1-yCom);
                             perc[1]=(1-xCom) * yCom;
                             perc[2]=xCom * (1-yCom);
                             perc[3]=xCom * yCom;

                             /* Update the 2D histogram with the kernel value */
                             hist2d_x=xBas[0]; hist2d_y=yBas[0];
                             Histogram_I2[mindex2(hist2d_x,hist2d_y,(int)HistogramBins)]+=perc[0]*Hkernel[mindex3(kernel_x,kernel_y,kernel_z,kernel_sizex,kernel_sizey)];
                             hist2d_x=xBas[0]; hist2d_y=yBas[1];
                             Histogram_I2[mindex2(hist2d_x,hist2d_y,(int)HistogramBins)]+=perc[1]*Hkernel[mindex3(kernel_x,kernel_y,kernel_z,kernel_sizex,kernel_sizey)];
                             hist2d_x=xBas[1]; hist2d_y=yBas[0];
                             Histogram_I2[mindex2(hist2d_x,hist2d_y,(int)HistogramBins)]+=perc[2]*Hkernel[mindex3(kernel_x,kernel_y,kernel_z,kernel_sizex,kernel_sizey)];
                             hist2d_x=xBas[1]; hist2d_y=yBas[1];
                             Histogram_I2[mindex2(hist2d_x,hist2d_y,(int)HistogramBins)]+=perc[3]*Hkernel[mindex3(kernel_x,kernel_y,kernel_z,kernel_sizex,kernel_sizey)];

                        }
                    }
                }

				/* Loop throught the image region which uses this kernel/histograms */
				for (image_z=image_z_start; image_z<(image_z_end+1);image_z++)
				{
					for (image_y=image_y_start; image_y<(image_y_end+1);image_y++)
					{
						for (image_x=image_x_start; image_x<(image_x_end+1);image_x++)
						{
						   /* pixel index  */
						   index = mindex3(image_x,image_y,image_z,image_sizex,image_sizey);
					   
						   /* The current histogram location */
						   hist2d_x=(int)round(scale*(I1[index]-range[0]));

						   /* Find maximum location */
						   /* maxx=0; maxf=0;for (i=0; i<HistogramBins; i++) { maxf+=i*Histogram_I1[mindex2(hist2d_x,i,(int)HistogramBins)]; maxx=maxx+Histogram_I1[mindex2(hist2d_x,i,(int)HistogramBins)]; } hist2d_y=(int)(maxf/maxx); */
						   hist2d_y=0; maxf=Histogram_I1[mindex2(hist2d_x,0,(int)HistogramBins)]; 
						   for (i=0; i<HistogramBins; i++)
						   {
							   if(Histogram_I1[mindex2(hist2d_x,i,(int)HistogramBins)]>maxf) { maxf=Histogram_I1[mindex2(hist2d_x,i,(int)HistogramBins)]; hist2d_y=i; }
						   }

						   /* Set intensity transformed pixel */
						   I1_TF[index]=(double)hist2d_y/(scale)+range[0];

						   /* The current histogram location */
						   hist2d_y=(int)round(scale*(I2[index]-range[0]));

						   /* Find maximum location */
							/* maxx=0; maxf=0; for (i=0; i<HistogramBins; i++) { maxf+=i*Histogram_I2[mindex2(i,hist2d_y,(int)HistogramBins)]; maxx=maxx+Histogram_I2[mindex2(i,hist2d_y,(int)HistogramBins)]; } hist2d_x=(int)(maxf/maxx); */
						   hist2d_x=0; maxf=Histogram_I2[mindex2(0,hist2d_y,(int)HistogramBins)];
						   for (i=0; i<HistogramBins; i++)
						   {
							   if(Histogram_I2[mindex2(i,hist2d_y,(int)HistogramBins)]>maxf) { maxf=Histogram_I2[mindex2(i,hist2d_y,(int)HistogramBins)]; hist2d_x=i; }
						   }

						   /* Set intensity transformed pixel */
						   I2_TF[index]=(double)hist2d_x/(scale)+range[0];
						}
					}
				}
			}
        }
    }
    free(Histogram_I1);
    free(Histogram_I2);
    

    /*  explicit end thread, helps to ensure proper recovery of resources allocated for the thread */
    #ifdef _WIN32
	_endthreadex( 0 );
    return 0;
	#else
	pthread_exit(NULL);
	#endif
    
}
    
/* The matlab mex function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    /* I1 and I2 are both input images with different modality     */
    double *I1, *I2;
    /* I1_moved and I2_moved are position transformed images (registration results) */
    double *I1_moved, *I2_moved;
    /* Range is a vector with the range of the input image */
    double *range;
    /* HistogramBins is the number of bins in the used histogram */
    double *HistogramBins_d;
    double HistogramBins;
    /* Hkernel is the kernel used to multiply with the Histogram update */
    /* for smooth modality transformations */
    double *Hkernel;
    /* Is the distance between the histogram kernels used to determine  */
    /* the modality transformation */
    double *SampleDistance_d;
    double SampleDistance;

    /* Multi-threading variables */
    mxArray *matlabCallOut[1]={0};
    mxArray *matlabCallIn[1]={0};
    double *Nthreadsd;
    int Nthreads;
	double Histograms_sizes[3]={0,0,0};
    double kernel_sizes[3]={0,0,0};
    double image_sizes[3]={0,0,0};
       
    /* double pointer array to store all needed function variables */
    double ***ThreadArgs;
    double **ThreadArgs1;
    
	/* Handles to the worker threads */
	#ifdef _WIN32
		HANDLE *ThreadList; 
    #else
		pthread_t *ThreadList;
	#endif
	
    
    /* ID of Threads */
    double **ThreadID;              
    double *ThreadID1;
    
    /* For loops variables */
    int i;
    
    /* Intensity transformed images */
    double *I1_TF, *I2_TF;
   
    /* Size of input image */
    const mwSize *image_dims; 
    int image_sizex, image_sizey, image_sizez;
	
    /* Size of histogram image  */
    int Histograms_size_x, Histograms_size_y, Histograms_size_z;
     
	/* Size of kernel */
	const mwSize *kernel_dims; 
    int kernel_sizex, kernel_sizey, kernel_sizez;
   
    /* Check for proper number of arguments. */
    if(nrhs!=8) {
       mexErrMsgTxt("8 inputs are required.");
    } else if(nlhs!=2) {
       mexErrMsgTxt("2 outputs are required");
    }
  
    /* Assign pointers to each input. */
    I1=(double *)mxGetData(prhs[0]);
    I2=(double *)mxGetData(prhs[1]);
    I1_moved=(double *)mxGetData(prhs[2]);
    I2_moved=(double *)mxGetData(prhs[3]);
    range=(double *)mxGetData(prhs[4]);
    HistogramBins_d=(double *)mxGetData(prhs[5]);
    Hkernel=(double *)mxGetData(prhs[6]);
    SampleDistance_d=(double *)mxGetData(prhs[7]); 
    
    /* From array to value */
    HistogramBins = HistogramBins_d[0];
    SampleDistance=SampleDistance_d[0];
    
    /* Get number of threads allowed */
    mexCallMATLAB(1, matlabCallOut, 0, matlabCallIn, "maxNumCompThreads");
	Nthreadsd=mxGetPr(matlabCallOut[0]);
	Nthreads=(int)Nthreadsd[0];
    
	/* Reserve room for handles of threads in ThreadList  */
	#ifdef _WIN32
		ThreadList = (HANDLE*)malloc(Nthreads* sizeof( HANDLE ));
    #else
		ThreadList = (pthread_t*)malloc(Nthreads* sizeof( pthread_t ));
	#endif
	ThreadID = (double **)malloc( Nthreads* sizeof(double *) );
	ThreadArgs = (double ***)malloc( Nthreads* sizeof(double **) );
    
    /* Get the sizes of the input image */
    image_dims = mxGetDimensions(prhs[0]);   
    image_sizex =(int)image_dims[0];
	image_sizey =(int)image_dims[1];
    image_sizez =(int)image_dims[2];
        
    /* Get the sizes of the kernel */
    kernel_dims = mxGetDimensions(prhs[6]);   
    kernel_sizex=(int)kernel_dims[0];
	kernel_sizey=(int)kernel_dims[1];
	kernel_sizez=(int)kernel_dims[2];
	    
    /* Create output mutual information image matrix */
	plhs[0] = mxCreateNumericArray(3, image_dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray(3, image_dims, mxDOUBLE_CLASS, mxREAL);
    
    /* Assign pointers to output. */
    I1_TF=(double *)mxGetData(plhs[0]);
    I2_TF=(double *)mxGetData(plhs[1]);

    /* Calculate number of kernel histograms  */
    Histograms_size_x = (int) ceil(image_sizex/SampleDistance); 
    Histograms_size_y = (int) ceil(image_sizey/SampleDistance); 
    Histograms_size_z = (int) ceil(image_sizez/SampleDistance); 
        
    /* Store sizes in double array format */
    Histograms_sizes[0]=(double)Histograms_size_x;
    Histograms_sizes[1]=(double)Histograms_size_y;
    Histograms_sizes[2]=(double)Histograms_size_z;
    kernel_sizes[0]=(double)kernel_sizex; 
    kernel_sizes[1]=(double)kernel_sizey;
    kernel_sizes[2]=(double)kernel_sizez;
    image_sizes[0]=(double)image_sizex; 
    image_sizes[1]=(double)image_sizey;
    image_sizes[2]=(double)image_sizez;
 
    for (i=0; i<Nthreads; i++)
    {
        /* Make Thread ID */
        ThreadID1= (double *)malloc( 1* sizeof(double) );
        ThreadID1[0]=i;
        ThreadID[i]=ThreadID1;  

        /* Make Thread Structure */
        ThreadArgs1 = (double **)malloc( 15* sizeof( double * ) );  
        ThreadArgs1[0]=I1;
        ThreadArgs1[1]=I1_moved;
        ThreadArgs1[2]=I2;
        ThreadArgs1[3]=I2_moved;
        ThreadArgs1[4]=Hkernel;
        ThreadArgs1[5]=Histograms_sizes; 
        ThreadArgs1[6]=HistogramBins_d;
        ThreadArgs1[7]=SampleDistance_d;
        ThreadArgs1[8]=kernel_sizes;
        ThreadArgs1[9]=image_sizes;
        ThreadArgs1[10]=range;
        ThreadArgs1[11]=I1_TF;
        ThreadArgs1[12]=I2_TF;
        ThreadArgs1[13]=ThreadID[i];
        ThreadArgs1[14]=Nthreadsd;
    
    /* Start a Thread  */
	ThreadArgs[i]=ThreadArgs1;
	#ifdef _WIN32
		ThreadList[i] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume, ThreadArgs[i] , 0, NULL );
	#else
		pthread_create ((pthread_t*)&ThreadList[i], NULL, (void *) &transformvolume, ThreadArgs[i]);
	#endif
  }
   
	#ifdef _WIN32
		for (i=0; i<Nthreads; i++) { WaitForSingleObject(ThreadList[i], INFINITE); }
		for (i=0; i<Nthreads; i++) { CloseHandle( ThreadList[i] ); }
	#else
		for (i=0; i<Nthreads; i++) { pthread_join(ThreadList[i],NULL); }
	#endif


  for (i=0; i<Nthreads; i++) 
  { 
    free(ThreadArgs[i]);
    free(ThreadID[i]);
  }

  free(ThreadArgs);
  free(ThreadID );
  free(ThreadList);
}
        
