load('tests');
colors = get(0, 'DefaultAxesColorOrder');
colors = [colors;colors];
names = {'syn1e3', 'syn1e4', 'bladder', 'lung', 'bladder2', 'lung2'};

for i=1:length(outs)
  % Convergence in primal objective
  figure;
  solvers = {};
  for j=1:length(outs{i})
    if ~isfield(outs{i}{j}, 'obj')
      continue
    end
    h = semilogy(outs{i}{j}.time, outs{i}{j}.obj - min_f(i) + min_gap(i), 'Color', colors(j,:), ...
                 'LineWidth', 1);
    if j >= 8
      set(h, 'Marker', 'o');
    end
    hold on;
    solvers{end+1} = opts{i}{j}.solver;
  end
  legend(solvers);
  prepare_figure(['obj_' names{i} '.pdf'], [8 6], 'Time (seconds)', 'f - f^*');

  figure;
  solvers = {};
  for j=1:length(outs{i})
    if ~isfield(outs{i}{j}, 'gap')
      continue
    end
    semilogy(outs{i}{j}.time, outs{i}{j}.gap, 'Color', colors(j,:), ...
             'LineWidth', 1);
    hold on;
    solvers{end+1} = opts{i}{j}.solver;
  end
  legend(solvers);
  prepare_figure(['gap_' names{i} '.pdf'], [8 6], 'Time (seconds)', ...
                 'Duality gap');
end