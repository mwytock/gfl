function D = lattice_operator(m,n)
% First order differences on m x n lattice

A = zeros(m,n); A(:,1:end-1) = 1;
B = zeros(m,n); B(:,2:end) = -1;
ai = find(A); bi = find(B);
D1 = [sparse(1:length(ai), ai, 1, length(ai), m*n) + ...
      sparse(1:length(bi), bi, -1, length(bi), m*n)];

A = zeros(m,n); A(1:end-1,:) = 1;
B = zeros(m,n); B(2:end,:) = -1;
ai = find(A); bi = find(B);
D2 = [sparse(1:length(ai), ai, 1, length(ai), m*n) + ...
      sparse(1:length(bi), bi, -1, length(bi), m*n)];

D = [D1; D2];
