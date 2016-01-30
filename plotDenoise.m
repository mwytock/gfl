
img = 'colors.png';
load(['denoise_' img '.mat']);
f0 = min(out3.f);
f0 = f0 - f0*1e-4;
figure;
semilogy((out3.f-f0)/f0, 'LineWidth', 1);
hold on;
semilogx((out2.objval-f0)/f0, 'r', 'LineWidth', 1);
set(gca, 'FontSize', 6);
legend('ADMM', 'Dykstra+ASPN');
ylim([1e-3 1]);
prepare_figure([img '_iter.pdf'], [3.5 1.7], 'Iteration', '(f-f^*)/f^*');

fprintf('D+ASPN: %f\n', (out2.time(end)-out2.time0) / length(out2.time));
fprintf('ADMM: %f\n', (out3.time(end)-out3.time0) / length(out3.time));
