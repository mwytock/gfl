#include <mex.h>
#include <string.h>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  // solve equations X = B * (ldlt(d,e))^{-1}
  int n = mxGetM(prhs[0]);
  int m = mxGetM(prhs[2]);
  double *d = mxGetPr(prhs[0]);
  double *e = mxGetPr(prhs[1]);
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  double *X = mxGetPr(plhs[0]);
  memcpy(X, mxGetPr(prhs[2]), m*n*sizeof(double));

  for (int i = 1; i < n; i++) {
    for (int j = 0; j < m; j++) {
      X[i*m+j] -= X[(i-1)*m+j]*e[i-1];
    }
  }

  for (int j = 0; j < m; j++) {
    X[(n-1)*m+j] /= d[n-1];
  }
    
  for (int i = n-2; i >= 0; i--) {
    for (int j = 0; j < m; j++) {
      X[i*m+j] /= d[i];
      X[i*m+j] -= X[(i+1)*m+j] * e[i];
    }
  }
}
