function plotChromosome(chr,Y,style,ylim)
% Pretty-plot profiles on chromosomes with chromosome separation
%
% plotUpDown(chr,Y)
% plotUpDown(chr,Y,style)
% plotUpDown(chr,Y,style,ylim)
%
% INPUT
% chr : a description of the chromosome number and boundaries, as created
% by the function chrInfo. It should contain:
%   chr.chr_nums
%   chr.chrs
%   chr.nchrs
%   chr.chr_data_len
% Y : a matrix of profiles to plot
% style : a cell array of strings describing the styles to plot each
% profile (default: all '-k' to have black lines)
% ylim : the limits on the Y axis (default: automaticly adjusted)
%
%
%    This file is part of GFLseg
%    Copyright (C) 2009-2011 Kevin Bleakley and Jean-Philippe Vert


[n,p] = size(Y);

% Define the Y limits
if nargin<4
    ylim = [min(min(Y))-0.1 ; max(max(Y))+0.1];
end

% Define default styles
if nargin<3
    style = cell(p,1);
    for i=1:p
        style{i} = 'k';
    end
end

% Draw a vertical bar at the end of a chromosome to indicate the border
x_vbar = repmat(chr.chr_nums, 2, 1);
y_vbar = repmat(ylim, 1, chr.nchrs);

% Label the autosomes with their chromosome numbers, and the sex chromosome
% with X.
x_label = chr.chr_nums - ceil(chr.chr_data_len/2);
y_label = zeros(1, length(x_label)) + ylim(1) - (ylim(2)-ylim(1))*0.12;
chr_labels=num2str(chr.chrs);
chr_labels = cellstr(chr_labels);
% Replace chromosome 24 by 'X'
chr24 = find(strcmp('24',chr_labels));
if length(chr24)>0
    chr_labels{chr24} = 'X';
end

hold on
% Plot the profile
for i=1:p
    h_ratio = plot(Y(:,i) , style{i});
end

% Separation between chromosomes
h_vbar = line(x_vbar, y_vbar, 'color', [0.8 0.8 0.8]);
% Chromosome names
h_text = text(x_label, y_label, chr_labels,...
    'fontsize', 8, 'HorizontalAlignment', 'center');

h_axis = get(h_ratio, 'parent');
set(h_axis, 'xtick', [], 'ygrid', 'on', 'box', 'on',...
    'xlim', [0 chr.chr_nums(end)], 'ylim', ylim)

xlabel({'', 'Chromosome'})
ylabel('Log-ratio')
hold off