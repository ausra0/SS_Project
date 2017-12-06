function [IndxI, IndxJ] = pidx2(kernel, x, n, test, varargin)
% PIDX samples according to the index distribution 
% IN : 
% kernel        scalar          along which kernel we want to sample
% x             vector Nx1      the sizes of the particles
% n             scalar          number of samples wanted
% test          boolean         if test script has to be executed
% varargin      (optional)      parameter a for kernel 2
% OUT : 
% Indx          matrix Nx2      coodinates of colliding particules
% ------------------------------------------------------------------

N = length(x); 
                                                                            % --- generate kernel
switch kernel
    case 1 
        K = single(x>0)*single(x>0)';                                       % we compute K1
    case 2
        a = varargin{:};
        K = (x.^a)*(x.^a)';                                                 % we compute K2
    case 3
        I = single(x>0)*single(x>0)';                                       % we compute 1{xy>0}
        K = (repmat(x, 1, N) + x').*I;                                      % we compute K3
    case 4
        I = single(x>0)*single(x>0)';                                       % we compute 1{xy>0}
        K = 0.25*((repmat(x.^(1/3), 1, N) + (x.^(1/3))').^3).*I;            % we compute K4
    otherwise 
        error('input kernel in range 1-4');
end
K(logical(speye(N)))=0;                                                     % set the diagonal to zero

% sample along j
margj = sum(K, 1); 
NonZeroj = find(margj);
cumweigths = cumsum(margj(NonZeroj));
cumweigths = cumweigths./(sum(margj));

Uj = rand(n,1); 

compte = histc(Uj, [0,cumweigths]);                                             
RawValues = repelem(NonZeroj, compte(1:length(compte)-1)); 
IndxJ = [mod(RawValues-1,N)+1]; 

% sample along i


% randomly generate U, map to j




                                                                            % --- prepare cdf for inverison
weights = K(K~=0);                                                          % weights of non-zero values
cumweigths = cumsum(weights);                                               % build cumsum
S = cumweigths(length(cumweigths));                                         % extract maximal value of cumweights
cumweigths = cumweigths./S;                                                 % compute the discrete CDF
NonZeroIndx = find(K);                                                      % extract non-zero indexes                                                                          % --- generate U~unif([1:S])
while(isempty(Indx))
U = rand(1, n);                                                             % vector 1xn, we generate all samples wanted                                                                       
                                                                            % --- map U back to the matrix 
compte = histc(U, cumweigths);                                              % convert U to indexes
RawValues = repelem(NonZeroIndx, compte); 
Indx = [ceil(RawValues./N),mod(RawValues-1,N)+1];                           % convert indexes to coordinates
end
if(test)                                                                    % --- run test script
   % cumweights = real CDF
   [compte, IndxKhat] = hist(RawValues, unique(RawValues));                 % count collisions
   Khat = zeros(N);                                                         % initialize matrix Khat, our approximate of the kernel
   Khat(IndxKhat) = compte;                                                 % assign count of collisions to the matrix
   Khat= Khat./sum(Khat(:));                                                % normalize Khat to obtain frquencies of collisions
   
   weigthsKhat = Khat(NonZeroIndx);                                         % estimate weights of non-zero values
   cumweigthsKhat = cumsum(weigthsKhat);                                    % build estimate cumsum
   Shat = cumweigthsKhat(length(cumweigthsKhat));                           % extract estimate maximal value of cumweights
   cumweigthsKhat = cumweigthsKhat./Shat;                                   % compute estimate CDF
   
   % plot QQ-plot
   figure
   plot(cumweigths, cumweigthsKhat); 
   xlabel('CDF(Khat)'); 
   ylabel('CDF(K)'); 
   title(sprintf('Q-Q plot for kernel %.0f and n = %.0f', kernel, n))
   
   D = max(abs(cumweigths - cumweigthsKhat));                               % Perform Kolmogorov-Smirnov test
   k = 1.36/sqrt(n);                                                        % K-S critical value
   fprintf('Testing for n=%.0f, kernel %.0f \n', n, kernel);
   if (D>k)
       fprintf('By Kolmogorov-Smirnov test, H0 accepted (alpha = 0.05)\n');
   else
       fprintf('By Kolmogorov-Smirnov test, H0 not accepted (alpha = 0.05)\n');
   end
end
   
end