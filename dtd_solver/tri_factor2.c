#include <mex.h>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  // setup factoring for solving (spdiags([-b a -b], [-1 0 1]))
  int n = mxGetM(prhs[0]);
  double *a = mxGetPr(prhs[0]);
  double *b = mxGetPr(prhs[1]);
  plhs[0] = mxCreateDoubleMatrix(n,1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(n-1,1, mxREAL);
  double *d = mxGetPr(plhs[0]);
  double *e = mxGetPr(plhs[1]);
  
  d[0] = a[0];
  for (int i = 0; i < n-1; i++) {
    e[i] = -b[i]/d[i];
    d[i+1] = a[i+1] + e[i]*b[i];
  }
}
