
load scaling_lambda;
nnzs = zeros(length(outs),1);
time = zeros(length(outs),1);
for i=1:length(outs)
  time(i) = outs{i}.time(end);
  nnzs(i) = length(find(outs{i}.z>0));
end

figure;
loglog(nnzs,time,'LineWidth',1);
set(gca, 'FontSize', 6);
prepare_figure('scaling_lambda.pdf', [3.5 1.5], 'Number of change points', 'Time (seconds)');

load scaling_length;
time2 = zeros(length(outs),1);
for i=1:length(outs)
  time2(i) = outs{i}.time(end);
end

figure;
loglog(Ts,time,'LineWidth',1);
set(gca, 'FontSize', 6);
prepare_figure('scaling_length.pdf', [3.5 1.5], 'Number of time points (T)', ...
               'TIme (seconds)');
