function [xmin, xmax, ymin, ymax, ytext] = plotParam_beta(speed_matrix, n_treat)

%PLOTPARAM This function computes the ideal plot parameters for cell speed
%          vs. substrate concentration 

% Unpack speed matrix
speeds = speed_matrix(:,2);
errors = speed_matrix(:,3);
concentrations = speed_matrix(:,5);
% offset the 0 ug/ml concentration for logarithmic plot (also because there
% is some substrate in media)
concentrations(concentrations == 0) = 1; 

% compute peaks and troughs
peaks = speeds + errors/2;
troughs = speeds - errors/2;

% compute extra space above and below plot
extra_space = 0.25*(max(peaks) - min(troughs));

% return plot limits and y-coordinate of text
xmin = min(concentrations)-min(concentrations)/2;
xmax = max(concentrations)+max(concentrations)/2;
ymax = max(peaks) + extra_space;
ymin = min(troughs) - extra_space;
ytext = ymax - linspace(0.2, 0.8, n_treat)*extra_space;
% ytext1 = ymax - 0.75*extra_space;
% ytext2 = ymax - 0.25*extra_space;