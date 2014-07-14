#include "mex.h"

/*#include "sys/queue.h"*/

/*
 * Informally, consider the matrix X2Y as representing a function F that maps 
 * P(N+) to P(N+), where  N+ is the set of positive integers and P() denotes 
 * the power set.
 * X2Y(i,j)==t, if t>0, means that F(i) contains t (t==0 means no mapping).
 * Then the output of this function, Y2X, represents the inverse of F.
 *
 * Input:  matrix X2Y (MxN) with integer values (double format)
 * Output: matrix Y2X (PxQ) where Q is the maximum number of occurences of an 
 *                    element X2Y, and P=max(X(:)). 
 *         vector H   (Px1), where H(i) is the number of occurences of i in X2Y.   
 *
 * If Y2X(i,j)==t, then either:
 *  - t==0, which implies H(t)<j
 *  - t>0,  which implies H(t)>=j and X2Y(t,p)==i for certain p
 *
 * Example: X=round(rand(100,5)*10-5); % generate random numbers
 *          X(X<=0)=0;                 % non-positive values set to zero
 *          X2=SURFING_INVERT_MAPPING(SURFING_INVERT_MAPPING(X))
 *
 *          X and X2 have the same elements in each row
 *
 * NNO May 2010
 */
void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
     int m, n, mn, i, j, maxval, vi, r, nymask;
     double *x2y, *y2x, *y2xrowlength; 
     int *xhist;
     int maxxhist;
     double v;
     double *ymask;
     int hasymask;
     int addtomask;
     mxArray *mx_y2x, *mx_y2xrowlength;
     char txt[100]; /* buffer for warning message */

     
     if (nlhs == 0) {
         /* No output argument, so return immediately */
         return;
     }
     if (nlhs > 2) {
        sprintf(txt, "Expected one or two output arguments, but found %i ",nlhs);
        mexErrMsgTxt(txt);
     }

     if (nrhs == 0 || nrhs > 2) {
         sprintf(txt, "Expected one or two input arguments, but found %i ",nrhs);
         mexErrMsgTxt(txt);
     }
     
     if (!(mxIsDouble(prhs[0]))) {
         mexErrMsgTxt("Only input of type double is supported; use double(X) to convert X to double");
     }


     /* size of input argument*/
     m = mxGetM(prhs[0]);
     n = mxGetN(prhs[0]);
     mn=m*n;
     
     /* pointer to input matrix*/
     x2y = mxGetPr(prhs[0]);
     
     /* NNO Jan 2011: supported for mask in y */
     hasymask=nrhs==2;
     if (hasymask) { 
         if (!(mxIsDouble(prhs[1]))) {
            mexErrMsgTxt("Only input of type double is supported; use double(X) to convert X to double");
         }
         nymask=mxGetM(prhs[1]) * mxGetN(prhs[1]); /* number of elements in mask */
         if (nymask==0) {
             hasymask=0;
         } else {
             ymask = mxGetPr(prhs[1]); /* pointer to mask */
         }
     }
     
     
     
     /* find maximum value */
     maxval=0;
     for (i=0; i<mn; i++) {
         v=x2y[i];
         if ( v > 0 ) {
             vi=(int) v; /* convert to int */
             if (v!=(double)vi) { /* ensure normal values */
                  mexErrMsgTxt("Illegal input value: non-integers not supported");
             }
             if (vi>maxval) {
                 maxval=vi;
             }
         }
     }
         
     /* count how often each element occurs 
      * xhist[k] is the number of occurences of k+1 */
     xhist=calloc(maxval,sizeof(int));
     if (xhist==0) {
         mexErrMsgTxt("Out of memory");
     }
         
     for (i=0; i<mn; i++) {
         if (x2y[i]>0) {
            ++xhist[(int) (x2y[i]-1)]; /* base 0 */
         }
     }
     
     /* maximum number of occurences of an element, i.e. the number of columns in the output */
     maxxhist=0; 
     for (i=0; i<maxval; i++) {
         if (xhist[i]>maxxhist) {
             maxxhist=xhist[i]; 
         }
     }
     
     if (hasymask && nymask != maxval) {
        sprintf(txt,"%i values in mask, but max(x2y)==%i - %s",
                     nymask,maxval,(nymask>maxval?
                                    "ignoring remaining values in mask":
                                    "including remaining values in mask"));
        mexWarnMsgTxt(txt);
     }    
     
     free(xhist); /* we know how big the output is; clear xhist */
     
     /* allocate space for output */ 
     mx_y2x = mxCreateDoubleMatrix(maxval, maxxhist, mxREAL); /* output matrix */ 
     y2x = mxGetPr(mx_y2x);   /* pointer to values in the matrix */
     
     /* free positions of each row in the output */
     mx_y2xrowlength=mxCreateDoubleMatrix(maxval, 1, mxREAL);
     y2xrowlength = mxGetPr(mx_y2xrowlength);
         
     for (i=0; i<m; i++) { /* each column */
         for (j=0; j<n; j++) { /* each row */
             vi=(int) x2y[j*m+i]-1;
             if (vi>=0 && (!hasymask /* accept if: no mask, or ... */
                          || (vi>=nymask       /* ... mask to small, or ... */
                              || (ymask[vi]==ymask[vi] && ymask[vi]!=0)))) /* ... mask value is not zero or NaN */
             { 
                 /* set y2x, and increment y2xrowlength for this value of vi */
                 y2x[vi + ((int)(y2xrowlength[vi]++))*maxval]=(double) i+1;
             }
         }
     }
     
     plhs[0] = mx_y2x; /* set y2x */
     
     if (nlhs==2) {
         plhs[1] = mx_y2xrowlength; /* set h (histogram */
     }
         
}


