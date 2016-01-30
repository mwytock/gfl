#include <mex.h>
#include <string.h>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  // setup factoring for solving (D*D' + diag(a))
  int n = mxGetM(prhs[0]);
  double *d = mxGetPr(prhs[0]);
  double *e = mxGetPr(prhs[1]);
  double *X = mxGetPr(prhs[2]);
  int m = mxGetM(prhs[2]);

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
