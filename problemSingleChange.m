function Y = problemSingleChange(T,n)
% Generate problem instance with single change point
Y = [repmat(randn(n,1), 1, T/2) + 0.1*randn(n,T/2) ...
     repmat(randn(n,1), 1, T/2) + 0.1*randn(n,T/2)];
