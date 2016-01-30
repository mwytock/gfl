
randn('seed',1);
rand('seed',1);
problems = {problemTimeSeries(1000,20,[0.1 0.15],0.05)
            problemTimeSeries(10000,20,[0.1 0.1],0.05)};

load('avecx_tests');

axlim = [0 5 1e-7 1e2
         0 60 1e-7 1e2];

% Convergence timing
for i=1:length(problems)
  min_f = outs{i}{4}.obj;
  figure;
  semilogy(outs{i}{1}.time, outs{i}{1}.obj - min_f, 'b', 'LineWidth', 1);
  hold on;
  semilogy(outs{i}{2}.time, outs{i}{2}.objval - min_f, 'r', 'LineWidth', 1);
  axis(axlim(i,:));
  legend('ADMM w/ TV', 'ADMM');
  prepare_figure(['arx' num2str(i) '_time.pdf'], [3 1.5], 'Time (seconds)', ...
                 'f - f^*');
end

set(0, 'DefaultAxesFontSize', 6);

% Results
figure;
plot(problems{2}.y,'LineWidth', 1);
box off;
ylim([-0.5 0.5]);
prepare_figure('arx1e4_y.pdf', [3 1.5], 'T', 'y');

figure;
X = [repmat(problems{2}.A1(1,:), 5000, 1);
     repmat(problems{2}.A2(1,:), 5000, 1)];
plot(X);
box off;
ylim([-0.25 0.25]);
prepare_figure('arx1e4_X.pdf', [3 1.5], 'T', 'X');

figure;
plot(outs{2}{1}.X');
box off;
ylim([-0.25 0.25]);
prepare_figure('arx1e4_admm_aspn.pdf', [3 1.5], 'T', 'X');

figure;
plot(reshape(outs{2}{2}.x, [], 1e4)');
box off;
ylim([-0.25 0.25]);
prepare_figure('arx1e4_admm.pdf', [3 1.5], 'T', 'X');
