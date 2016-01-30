T = 1000000;
D = spdiags([-ones(T,1) ones(T,1)], [0 1], T-1, T);

u = randn(T-1,100);
z = rand(T-1,1);

tic;
R = chol(D*D' + spdiags(z, 0, T-1, T-1));
y = R \ (R' \ (u));
toc;

ut = u';
tic;
[d,e] = tri_factor(z);
tri_solve(d,e,ut);
toc;
