function [logL] = LL_negbin(params, X, y)
    beta = params(1:(end-1));
    alpha = params(end);  
    mu = exp(X * beta);
    
    logL = -sum(gammaln(y + inv(alpha)) - gammaln(y + 1) - gammaln(inv(alpha)) + ...
           inv(alpha) .* log(inv(alpha) ./ (inv(alpha) + mu)) + ...
           y .* log(mu ./ (inv(alpha) + mu)));
end