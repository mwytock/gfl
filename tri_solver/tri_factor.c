#include <mex.h>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  // setup factoring for solving (D*D' + diag(a))
  int n = mxGetM(prhs[0]);
  double *a = mxGetPr(prhs[0]);
  plhs[0] = mxCreateDoubleMatrix(n,1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(n-1,1, mxREAL);
  double *d = mxGetPr(plhs[0]);
  double *e = mxGetPr(plhs[1]);
  
  d[0] = 2.0+a[0];
  for (int i = 0; i < n-1; i++) {
    e[i] = -1.0/d[i];
    d[i+1] = 2.0 + a[i+1] + e[i];
  }
}
