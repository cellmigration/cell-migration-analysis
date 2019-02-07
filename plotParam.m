function [ymin, ymax, ytext1, ytext2] = plotParam(speed_matrix1, speed_matrix2)

%PLOTPARAM This function computes the ideal plot parameters for cell speed
%          vs. substrate concentration 

% combine matrices
speed_matrix = [speed_matrix1; speed_matrix2];

% Unpack speed matrix
speeds = speed_matrix(:,1);
errors = speed_matrix(:,2);

% compute peaks and troughs
peaks = speeds + errors/2;
troughs = speeds - errors/2;

% compute extra space above and below plot
extra_space = 0.2*(max(peaks) - min(troughs));

% return plot limits and y-coordinate of text
ymax = max(peaks) + extra_space;
ymin = min(troughs) - extra_space;
ytext1 = ymax - 0.75*extra_space;
ytext2 = ymax - 0.25*extra_space;