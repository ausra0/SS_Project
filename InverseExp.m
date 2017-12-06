function x = InverseExp(lambda, n)
%INVERSEEXP generates values from exponential distribution using the
%inverse CDF method 
% IN : 
% lambda            scalar              distribution parameter
% n                 scalar              the number of samples wanted
% OUT : 
% x                 1xn vector          values generated
% -----------------------------------------------------------------
u = rand(1, n); % generate u ~ unif[0,1]
x = -1/lambda .* log(1-u); 
end