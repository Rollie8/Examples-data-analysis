% We write the LogLikelihood function, which we can append 
% to the main file or put in the same folder
function [LL] = LL_Poisson(theta,Y,X); 
% LL = -sum(Y.*(X*theta)-exp(X*theta)-log(factorial(Y))) % correct but unstable
LL = -sum(Y.*(X*theta) - exp(X*theta) - gammaln(Y+1))
end