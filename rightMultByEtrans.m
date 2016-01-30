function V=rightMultByEtrans(U)
n=size(U,2);
K=U(:,2:n)-U(:,1:n-1);
V=[U(:,1),K,-U(:,n)];
