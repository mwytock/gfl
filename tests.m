T = 10000;
n = 100;

addpath mtv;
addpath dtd_solver;
addpath(genpath('./GFLseg-1.0'))
Y = [repmat(randn(n,1), 1, T/2) + 0.1*randn(n,T/2) ...
     repmat(randn(n,1), 1, T/2) + 0.1*randn(n,T/2)];

lam = 10*ones(T-1,1);
opt = options('pn');
opt.k_max = 400;
opt.max_iters = 1000;
opt.ls_quadratic=1;
profile on; out = pnTV(Y,lam,opt); profile viewer;
profile on; out = aspnTV(Y,lam,opt); profile viewer;

gflopt.weights = ones(T-1,1);
tic; res = gflasso(Y', lam(1), gflopt); toc;



figure(1); plot(1:T,Y');
figure(2); plot(1:T,out.X');


%Y = load('GFLseg-1.0/data/bladderData.txt')';
Y = load('GFLseg-1.0/data/lungData.txt')';
Y = Y(4:end,:);
Y(isnan(Y)) = 0;
lam = 300*ones(size(Y,2)-1,1);
opt = options('aspn');
opt.k_max = 200;
opt.max_iters = 1000;
opt.ls_quadratic=1;
profile on; out = aspnTV(Y,lam,opt); profile viewer;

gflopt.weights = lam/lam(1);
gflopt.max_time = 10;
tic; res = gflasso(Y', lam(1), gflopt); toc;


%%% weighted tests
T = 10000;
n = 100;
Y = [repmat(randn(n,1), 1, T/2) + 0.1*randn(n,T/2) ...
     repmat(randn(n,1), 1, T/2) + 0.1*randn(n,T/2)];
lam = 20*ones(T-1,1);
tic; out = pnwTV(Y,ones(T,1),lam,opt); toc
X = out.X;

[d,e] = dtd_factor(zeros(T-1,1));
U = dtd_solve(d,e,diff(Y,1,2));


% construct reduced version
ch = find(abs(out.z) > 0)';
t = [0 ch T];
T0 = length(ch)+1;
Y0 = zeros(n,T0); w = zeros(T0,1);
Ycs = [zeros(n,1) cumsum(Y,2)];
Y2cs = [zeros(n,1) cumsum(Y.^2,2)];
for i=1:T0,
  w(i) = t(i+1) - t(i);
  Y0(:,i) = (Ycs(:,t(i+1)+1) - Ycs(:,t(i)+1))/w(i);
  f(i) = sum(Y2cs(:,t(i+1)+1) - Y2cs(:,t(i)+1)) - sum(Y0(:,i).^2*w(i));
end

tic; out2 = pnwTV(Y0,w,lam(ch),opt); toc;


% primal
D0 = spdiags([-ones(T0,1) ones(T0,1)], [-1 0], T0, T0-1);
cvx_begin
variable X0(n,T0);
minimize(0.5*square_pos(norm((X0-Y0)*diag(sqrt(w)),'fro')) + norms(X0*D0,2)*lam(ch))
cvx_end



% first dual
DtD = D0'*diag(1./w)*D0;
cvx_begin
variable U0(n,T0-1);
maximize(-0.5*square_pos(norm(U0*sqrtm(DtD),'fro')) + trace(U0*D0'*Y0'))
norms(U0,2) <= lam(ch)'
cvx_end


%%% test new tri solvers
T = 10;
w = rand(T,1);
a = w(1:end-1)+w(2:end);
b = w(2:end-1);
K = full(spdiags([-[b;0] a -[0;b]], [-1 0 1], T-1,T-1));
B = randn(2,T-1);

cslogb = cumsum(log(b));


Ki = tri_invert(a,b,cslogb)
[rd,re] = tri_factor(a,b);
(tri_solve(rd0,re0,B) - B/K)



D0'*diag(w)*D0