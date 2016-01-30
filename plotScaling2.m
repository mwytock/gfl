
load scaling_lambda;

% Fit k^3 line to data
n = length(nnzs);
X = [(nnzs/1e2).^3];
y = time;
yh = X*((X'*X)\(X'*y));

figure;
loglog(nnzs,time,'LineWidth',1);
hold on;
loglog(nnzs,yh,'Color',[0.5 0.5 0.5],'LineWidth',1);
set(gca, 'FontSize', 8);
axis([1 1e4 1e-2 1e3])
legend('Actual', 'O(k^3)', 'Location', 'NorthWest');
prepare_figure('scaling_lambda3.pdf', [3.5 1.7], 'Number of change points (k)', 'Time (seconds)');

load scaling;

n = length(Ts);
X = [Ts'];
y = time2;
yh = X*((X'*X)\(X'*y));

figure;
loglog(Ts,time2,'LineWidth',1);
hold on;
loglog(Ts,yh,'Color',[0.5 0.5 0.5],'LineWidth',1);
set(gca, 'FontSize', 6);
legend('Actual', 'O(T)', 'Location', 'NorthWest');
prepare_figure('scaling_length2.pdf', [3.5 1.7], 'Number of time points (T)', ...
               'TIme (seconds)');
