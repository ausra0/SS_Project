function [Indx,RawValues, NonZeroIndx, cumweigths] = pidx3(kernel, x, n, K, varargin)
% PIDX samples according to the index distribution
% IN : 
% kernel        scalar          along which kernel we want to sample
% x             vector Nx1      the sizes of the particles
% n             scalar          number of samples wanted
% K             
% varargin      (optional)      parameter a for kernel 2
% OUT : 
% Indx          matrix Nx2      coodinates of colliding particules
% ------------------------------------------------------------------

y = x(x~=0); 
Ny = length(y); 

% --- generate kernel
switch kernel
    case 0
        K = K; 
    case 1 
        K = single(y>0)*single(y>0)'-eye(Ny).*y.^2;                                       % we compute K1
    case 2
        a = varargin{:};
        K = (x.^a)*(x.^a)' - speye(N).x.^2;                                                 % we compute K2
    case 3
        I = single(x>0)*single(x>0)';                                       % we compute 1{xy>0}
        K = x.*I + x'.*I;                                                   % we compute K3
    case 4
        I = single(x>0)*single(x>0)';                                       % we compute 1{xy>0}
        K = 0.25*((x.^(1/3)).*I +((x.^(1/3))').*I).^3;                      % we compute K4
        
    otherwise 
        error('input kernel in range 1-4');
end
K(logical(speye(N)))=0;                                                     % set the diagonal to zero

                                                                            % --- prepare cdf for inverison
weights = K(K~=0);                                                          % weights of non-zero values
cumweigths = cumsum(weights);                                               % build cumsum
S = cumweigths(length(cumweigths));                                         % extract maximal value of cumweights
cumweigths = cumweigths./S;                                                 % compute the discrete CDF
NonZeroIndx = find(K);                                                      % extract non-zero indexes                                                                          
                                                                            % --- generate U~unif([1:S])
U = rand(1, n);                                                             % vector 1xn, we generate all samples wanted                                                                       
                                                                            % --- map U back to the matrix 
compte = histc(U, [0; cumweigths]);                                         % convert U to indexes
RawValues = repelem(NonZeroIndx, compte(1:length(compte)-1)); 
Indx = [ceil(RawValues./N),mod(RawValues-1,N)+1];                           % convert indexes to coordinates
end