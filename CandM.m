function [meanC, sigmaC, meanM, sigmaM] = CandM(kernel, timeline, R, x, T, k, p, varargin)
 %CANDM computes approximates for average particle concentration C and
 %moments M
 % IN : 
 % kernel           scalar          1-4 kernel reference number
 % timeline         vector 1xS      time linspace on which we want to
 %                                  sample M and C
 % R                scalar          number of simulation repeats
 % x                vector Nx1      sizes of the particles
 % T                scalar          sample on [0, T]
 % k                vector 1xK      values for which to compute c(t,k)
 % p                vector 1xP      values for which to compute m_p
 % varargin         (optional)      contains values for a (kernel2)
 % OUT : 
 % meanC            matrix 1xSxK    mean of average particle concentration
 % sigmaC           matrix 1xSxK    sd of average particle concentration
 % meanM            matrix 1xSxP    mean of moments over R simulations
 % sigmaM           matrix 1xSxP    sd of moments over R simulations
 % ---------------------------------------------------------------
 %N = length(x); 
 C = zeros(R, length(timeline), length(k)); 
 M = zeros(R, length(timeline), length(p)); 
 
 for i = 1:R
    %tic
    [JC, J] = CTDSMC(kernel, x, T, false, varargin{:}); % simulate process
    out = TableProcess(timeline, J, JC); % interpolate on linspace
    
    for j = 1:length(k) 
        C(i , :, j) = mean(out==k(j),1); % define average particle concentration
    end
    for j = 1:length(p) 
        M(i , :, j) = mean(out.^p(j),1); % define moments
    end
    %t = toc;
    %clc
    %fprintf('Running for N = %.0f\nest. time remaining : %0.3f sec\n', N, (R-i)*t)
 end

meanC = mean(C, 1);  % compute estimate average on R simulations
meanM = mean(M, 1); 

sigmaC = sum((C - meanC).^2, 1); % compute sd estimate on R simulations
sigmaM = sum((M - meanM).^2, 1); 
end