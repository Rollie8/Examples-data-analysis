% (1) Load the dataset PoissonDATA.mat.
load('PoissonDATA.mat')

% Take the regressors from the slide to understand the context:
% 1) Health insurance measures:
% LC (log(coinsrate+1) where coinsurance rate is 0 to 100);
% IDP (1 if individual has deductible plan, 0 otherwise); 
% LPI (log(annual participation incentive payment) or 0 if no payment); 
% FMDE (log(medical deductible expenditure) if IDP=1, 0 otherwise) ->
% should we remove one of them?
% 2) Health status measures:
% DISEASE (1 if chronic disease);
% NDISEASE (number of chronic diseases);
% PHYSLIM (1 if physical limitation); 
% HLTHG (1 if good health);
% HLTHF (1 if fair health)
% 3) Socioeconomic characteristics:
% LINC (log of annual family income);
% LFAM (log of family size); 
% EDUCDEC (years of schooling of decision maker); 
% AGE (age); 
% BLACK (1 if black); 
% FEMALE (1 if female); 
% CHILD (1 if child); 
% FEMCHILD (1 if female child)

X=[ones(length(Y),1), x]; 
[N,k]=size(X);

% (2) Write a MATLAB function for the log-likelihood as a function of parameters
% and data. 

function [logL] = LL_negbin(params, X, y)
    beta = params(1:(end-1));
    alpha = params(end);  
    mu = exp(X * beta);
    
    logL = -sum(gammaln(y + inv(alpha)) - gammaln(y + 1) - gammaln(inv(alpha)) + ...
           inv(alpha) .* log(inv(alpha) ./ (inv(alpha) + mu)) + ...
           y .* log(mu ./ (inv(alpha) + mu)));
end

% (3) Estimate the parameters of the NBRM using this data.
% starting point, add 0.05 for the alpha overdispersion parameter
rng(22);
parameters_init=[randn(k,1)/50; 0.05];
% lower and upper bound, let's put a number for LB != 0 but very close
% otherwise it is not possible to do inv(alpha)
LB=[-50*ones(k,1); 0.0000000000000000000000001]; 
UB=[50*ones(k,1); 5];

options = optimoptions('fmincon','Display','iter',...
'MaxIterations', 1000,...
'OptimalityTolerance',1e-6, ...
'StepTolerance', 1e-8, ...
'MaxFunctionEvaluations',1e4);

% let's now maximize the function using fmincon and estimate the parameters
% such that the function is maximized

[parameters_nb,LL_nb] = fmincon(@(parameters) LL_negbin(parameters,X,Y), parameters_init,[],[],[],[],LB,UB,[],options); 
% Extract constant and betas parameters from the parameters_nb
betas = parameters_nb(1:end-1);
% Extract the overdispersion parameter alpha
alpha = parameters_nb(end);
disp("Alpha =")
disp(alpha)
W = diag(exp(X*betas)./(1+exp(X*betas)./alpha));
% This is the formula for negative binomial. We search it on the web
% because in the slides there is only the one for the Poisson.
std_err = sqrt(diag(inv(X'* W *X)));

t_stat = betas./std_err;

regressor_names = ['Constant'; 'LC'; "IDP"; "LPI"; "FMDE"; "DISEASE1"; "NDISEASE"; "PHYSLIM"; "HLTHG"; "HLTHF"; "LINC"; "LFAM"; "EDUCDEC"; "AGE"; "BLACK"; "FEMALE"; "CHILD"; "FEMCHILD"];

%(4) Report and briefly interpret the results (coefficients, standard errors, t-statistics).
disp("-------------------- Negative binomial --------------------");
disp(table(regressor_names, betas, std_err, t_stat, betas-1.96*std_err, betas+1.96*std_err,...
    'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));


% INTERPRETATION OF RESULTS
% All variables, except for PHYSLIM are significant, so let's comment them:
 
% • LC, so logharitm of coinsurance rate is significant in decreasing the
% number of doctor visits, this result is quite explainable by the fact
% that the higher the mount of health cost that person has to pay, the
% lower the care they are willing to receive.
% • IDP, A deductible plan is a health insurance plan that requires you to pay a
% set amount of money before your insurance pays for covered service,
% so passing from 0 to 1 in deductible plan will decrease the number of
% doctor visits, this results can be explained with the same concept of the
% previous regressor, if a person has to pay in advance for their service
% before the insurance, probably she will not require care very often.
% • LPI, defined as: An annual participation incentive payment is a reward
%  given to employees for meeting specific performance goals within a year,
%  increased the number of doctor visits, explaining that if people have
%  incentives in money terms, they will increase their expected number of doctor
%  visits.
% • FMDE, defined the amount of money you pay before your health insurance pays for covered health care
% services, decrease the number of doctor visits, all these variables are
% referring to the same concept of money-quantity of care relationship.
% • DISEASE, as it is obvious to think, having a chronic disease will
% increase the amount of visits to the doctor in a year.
% • NDISEASE Of course, in case of multiple number of diseases, the number
% of visits increase even more.
% • HLTHG, HLTHF having good health surprisingly increase the number of
% doctor visits in a significant manner, this probably is due to prevention
% behavior.
% • LINC, the higher the family income, the higher the expected count as
% expected.
% • LFAM, on the other hand, the higher the family size, the lower the
% number of visits, this is probably due to the fact that as family size
% increases, large number of priority overcome personal care in the
% priority list.
% • EDUCDEC, years of schooling increase the expected count, the higher the
% level of education the higher the attention to personal care.
% • AGE, the higher the age the more the number of doctor visits as
% expected.
% • BLACK and FEMALE increase the number of doctor visits.
% • On the other hand having a child, especiallly if she is female decrease
% the expected count.

%(5) Compute the sample average marginal effect of age on the expected number
% of doctor visits.
disp('---AME Negative Binomial---')
% As it sis possible to notice, the AME in this case tells us that
% after a one unit increase in age increases the expectation of y
% (number of doctor visits) of almost 1 unit, this make sense since it is
% trivial that the greater it is a person, the most likely the person need
% to be visited by a doctor.
AME = mean(exp(X(:, 14) *betas(14) .* betas(14)));

% (6) Test the overdispersion hypothesis (i.e., test whether α = 0).
% To test the overdispersion hypotesis, a likelihood ratio test can be
% conducted.

% So let's do first a poisson model and take the value of the function from
% there.

function [LL] = LL_Poisson(theta,Y,X); 
% LL = -sum(Y.*(X*theta)-exp(X*theta)-log(factorial(Y))) % correct but unstable
LL = -sum(Y.*(X*theta) - exp(X*theta) - gammaln(Y+1));
end

% remove the alpha values from the initializers and the lower and upper bound 
parameters_init = randn(k,1)/50;
LB=-50*ones(k,1); 
UB=50*ones(k,1);
[parameters_poisson,LL_poiss] = fmincon(@(parameters) LL_Poisson(parameters,Y,X), parameters_init,[],[],[],[],LB,UB,[],options); 

% Let's now perform the LR test, usually this test should be
% loglik_unrestricted - loglik_restricted(poisson put a restriction on alpha = 0)but in this case we do the
% opposite because we changed sign in fmincon, in fact either LL_poiss and
% LL_nb have very large positive values, while usual loglik is negative.
LR = 2*(LL_poiss-LL_nb);

% Let's calculate the critical value to determine the significance level
% for the LR test
critical_value = chi2inv(0.95, 1)
display(LR)

% As it is posible to see, the LR value is very large with respect to the
% critical value, and so the hypotesis that the two models are equal, and
% so alpha = 0, it's rejected, negative binomial will better fit the data
% in case of overdispersion presence.

%(7) Plot the empirical distribution of counts in the sample,
% along with the distribution of predicted counts from the Poisson and NBRM models.

% Let's calculate the predicted values of the counts from both poisson and
% negbin, each expected value of P(y|x) = exp(X * beta), so:
y_hat_poisson = exp(X * parameters_poisson);
y_hat_nb = exp(X * betas);

% 3. Predicted counts – NB (using theta_hat from Q2)

%Q7
% 1. Observed counts
figure;
hold on;
histogram(Y, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 2);

% 2. Predicted counts – Poisson (using glmfit)
mu_poisson = exp(X * parameters_poisson);
histogram(round(mu_poisson), 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 2);

mu_nb = exp(X * betas);
histogram(round(mu_nb), 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'LineWidth', 2);

% Add legend and title
legend('Observed Y', 'Predicted Poisson', 'Predicted NB');
title('Q7 – Comparison of Count Distributions');
xlabel('Number of doctor visits');
ylabel('Probability');
grid on;
hold off;

% let's now plot the final plot with the comparisons
figure;
hold on;
plot(sort(Y), 'ko', 'LineWidth', 0.5);
plot(sort(y_hat_poisson), 'go', 'LineWidth', 0.5);
plot(sort(y_hat_nb), 'bo', 'LineWidth', 0.5);
legend({'Real Y', 'Y_hat of Poisson', 'Y_hat of NegBinomial'}, 'Location', 'best');
hold off;

% As it is possible to notice from the plot, poisson and negbin follows the real curve
% larger part of the values for count < 10, however for large values of real
% Y, the predicted value of both poisson and negbin are very far from real one, but negbin has larger max values.

% From the plot we can see that the two models are similar but the negbin is a bit better. 
% Both the models predict less zeros compared to the observed counts and
% overestimate the counts for small values of doctor visits; 
% underestimating, instead, the dependent value of larger values of Y.

sum(Y == 0)/length(Y) % = 0.3125
% There is a big quantity of zeros, almost 1/3 of the data are zeros.
% We can consider using a zero-inflated negative binomial model, 
% that uses a logit model to detect the presence of structural zeros in data.