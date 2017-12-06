function [GN, meanGN, sigmaGN] = Gn(kernel, timeline, x, T, R, varargin)
%GN simulates G_N = S_{1}(t)/N
% IN : 
% kernel                scalar              1-4 kernel reference number
% timeline              vector 1xS          time linspace for sampling
% x                     vector Nx1          initial particles sizes
% T                     scalar              sample on [0, T]
% R                     scalar              number of samples
% varargin              (optinal)           parameter a for kernel 2
% OUT : 
% GN                    vector
% meanGN
% sigmaGN
% ---------------------------------------------------------------------
GN = zeros(R, length(timeline));
N = length(x); 
for i = 1: R 
    %tic
    [JC, J] = CTDSMC(kernel, x, T, false, varargin{:});                     % simulate process
    out = TableProcess(timeline, J, JC);                                    % interpolate the process
    GN(i, :) = max(out, [], 1)./N;                                          % define GN
    %t = toc; 
    %clc
    %fprintf('Computing for N = %.0f\nEst. time remaining : %.03f sec\n', N, t*(R-i))
end

meanGN = mean(GN,1);                                                        % empiric mean of samples
sigmaGN = 1/(N-1).*sum((GN - meanGN).^2, 1);                                % empiric standard deviation of samples
end