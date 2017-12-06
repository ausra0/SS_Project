function [JumpChain, Jumps] = CTDSMC(kernel, x0, T, process, varargin)
% CTDSMC computes continuous time/discrete state Markov chain for our
% assignments transition probabilities 

% IN : 
% kernel        scalar 1-4          kernel to be used for simulation
% x0            vector Nx1          inintial sizes of the molecules
% T             scalar              generate process for t in [0,T]
% process       boolean             processing to get nice ploting data
% varargin      (optional)          parameter a for kernel 2
% OUT : 
% JumpChain     vector mx1          history of values Xjn n in {1, m}
% Jumps         vector mx1          history of jumps (time values when jump occurs)
% -----------------------------------------------
X = x0; % initialize 
J = 0; 
Y = X; 
N = length(x0); 

JumpChain = X; % initialize output variables 
Jumps = J; 
while(J<=T && length(find(X))>1)
                                                                            % we compute K, the kernel
    switch kernel                                                           % choose which kernel to compute
        case 1 
            K = single(X>0)*single(X>0)';                                   % we compute K1
            K(logical(speye(N)))=0;                                         % make sure that the diagonal is zero
        case 2
            a = varargin{:};
            K = (X.^a)*(X.^a)';                                             % we compute K2
            K(logical(speye(N)))=0;                                         % make sure that the diagonal is zero
        case 3
            I = single(X>0)*single(X>0)';                                   % we compute 1{xy>0}
            K = X.*I + X'.*I;                                               % we compute K3
            K(logical(speye(N)))=0;                                         % make sure that the diagonal is zero
        case 4
            I = single(X>0)*single(X>0)';                                   % we compute 1{xy>0}
            K = 0.25*((X.^(1/3)).*I + (X.^(1/3))'.*I).^3;        % we compute K4
            K(logical(speye(N)))=0;                                         % make sure that the diagonal is zero
        otherwise 
            error('input kernel in range 1-4');
    end
    
    q = double(sum(sum(K))/(2*N));                                          % compute q
    S = InverseExp(q, 1);                                                   % generate S from exp(q)
    J = J + S;                                                              % define time at which the jump occured
    Jumps = [Jumps, J];                                                     % store jump history
    
    couple = pidx(0, X, 1, K, varargin{:});                        % generate Y ~Pidx
    Y(couple(1)) = X(couple(1)) + X(couple(2)); 
    Y(couple(2)) = 0; 
    
    X = Y;                                                                  % update Jump Chain
    JumpChain = [JumpChain, Y];                                             % Store Jump Chain
end

if(length(find(X))<=1)
    disp('Process Atteined Stable State');
end

if (process)
   OldJumpChain = JumpChain; 
   OldJumps = Jumps; 
   
   %JumpChain = zeros(size(OldJumpChain, 1), 2*size(OldJumpChain)
   JumpChain(:,2:2:2*size(OldJumpChain, 2)-2) = OldJumpChain(:,1:size(OldJumpChain, 2)-1); 
   JumpChain(:,1:2:2*size(OldJumpChain, 2)-1) = OldJumpChain; 
   
   Jumps(1:2:2*length(OldJumps)-1) = OldJumps; 
   Jumps(2:2:2*length(OldJumps)-2) = OldJumps(2:size(OldJumps, 2)) - eps;
end
end