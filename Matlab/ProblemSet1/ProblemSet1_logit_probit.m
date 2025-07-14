% BINARY REGRESSION

% (1) Import the dataset Data Logit.xlsx.
data = xlsread('Dati Logit.xlsx', "Sheet1");

% (2) Estimate a logit regression of the FOSD violation dummy on all other vari-
% ables and a constant. Maximize the log-likelihood using both fminunc and
% fmincon in MATLAB. (This is primarily for practice with both optimizers;
% fminunc is recommended.)

Y = data(:, 1);
X = [ones(length(data), 1), data(:,2:end)];
[N,k]=size(X);

% Since we will use fmincon (or fminunc) to optimize our function and find
% the parameters, we have to define the starting points and the boundaries
% in order to not make the process too heavy in computational terms.
 
% Define starting point, let's set a seed in order to have the possibility
% to have always the same result even though there is a random start
rng(22);
parameters_init=randn(k,1)/50; 
% Define lower and upper bound, let's put +- 50 given that these values are
% already very huge in a logit model, so it is useless to allow the alghorithm to have larger
% values.
LB=-50*ones(k,1); UB=50*ones(k,1); 

% Let's define the options for alghorithm convergence, in particular we
% fixed an optimality tolerance of 1e-6 and a steptolerance of 1e-8, that means that if the alghoritm, 
% from an iteration to the next one (specified in max iterations = 1000),
% found that parameters change is within a value of 1e-8, it can stop since 
% the "convergence" is specified, at least for the options that the analyst
% put.
options = optimoptions('fmincon','Display','iter',...
'MaxIterations', 1000,...
'OptimalityTolerance',1e-6, ...
'StepTolerance', 1e-8, ...
'MaxFunctionEvaluations',1e4);

% Let's now define the objective function by expliciting the loglik of the
% logistic function.
function [LL] = LL_logit(parameters,Y,X) 
p = 1 ./ (1 + exp(-X * parameters));
% put -1 because we want to maximize log-lik while both fmincon and fminunc
% minimize the function.
LL = -1* (sum(Y .* log(p) + (1 - Y) .* log(1 - p)));
end
objective_fun = @(parameters) LL_logit(parameters, Y, X);

% in fmincon the first two [] are the inequality constraints,
% third and fourth [] are equality constraints,
% last one are the non linear constraints, in these case there are no
% constraints in fact we can also use fminunc.
tic
[parameters_fmincon, LL, ~,~,~,GRAD,HESSIAN] = fmincon(objective_fun, parameters_init, [], [], ...
[], [], LB, UB, [], options);
toc

P_logit = 1 ./ (1 + exp(-X * parameters_fmincon));
% We should use the inverse of the hessian matrix. But, as stated in slide
% 24, in the logit model the estimated variance matrix of the MLE simplifies
% as this formula
std_err_fmincon = sqrt(diag(inv(X'*diag(P_logit .* (1 - P_logit))*X)));

options = optimoptions('fminunc','Display','iter',...
'MaxIterations', 1000,...
'OptimalityTolerance',1e-6, ...
'StepTolerance', 1e-8, ...
'MaxFunctionEvaluations',1e4);

% Let's try fminunc
tic
[parameters_fminunc, LL] = fminunc(objective_fun, parameters_init, options);
toc

P_logit = 1 ./ (1 + exp(- X * parameters_fminunc));
std_err_fminunc = sqrt(diag(inv(X'*diag(P_logit .* (1 - P_logit))*X)));

regressor_names = ['Constant'; 'Cognitive ability score'; "Time_ability_test";
    "Understanding_experiment"; "Understanding_experiment2"; "Gender"; "Education";
    "Response_time_experiment"; "Exp_value_risk_neutrality"; "High risk aversion";
    "Randomization_behavior"];

% fmincon is slower (elapsed time is 0.792732 seconds in our laptop) but more precise than fminunc for standard errors.
% fminunc is faster (0.475166 seconds) but less precise than fmincon for coefficients and standard errors.
% fsolve is both faster (0.109685 seconds) and precise. So, using the gradient is recommended
% if there is the possibility to compute it (gradient computation can be very time consuming in real life).
disp("-------------------- FMINCON --------------------");
disp(table(regressor_names, parameters_fmincon, std_err_fmincon, parameters_fmincon./std_err_fmincon, parameters_fmincon-1.96*std_err_fmincon, parameters_fmincon+1.96*std_err_fmincon,...
    'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));
disp("-------------------- FMINUNC --------------------");
disp(table(regressor_names, parameters_fminunc, std_err_fminunc, parameters_fminunc./std_err_fminunc, parameters_fminunc-1.96*std_err_fminunc, parameters_fminunc+1.96*std_err_fminunc,...
    'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));
 

% (3) Re-estimate the parameters by solving the first-order conditions using fsolve.
% The FOC is:
% Sum(yi − Λ(x′θ)xik = 0, k = 1, . . . , K .

% So, let's write the gradient function and give it in input to fsolve:
% Definition of the function that returns the gradient of the
% log-likelihood, in which we substitute Λ(x′θ) with the prob(y = 1|x)
gradLogLik = @(beta) X' * (Y - 1 ./ (1 + exp(-X * beta))); 

% From the help, this is the syntax of fsolve: x = fsolve(fun,x0,options),
% so:
options = optimoptions('fsolve','Display','iter',...
'MaxIterations', 1000,...
'OptimalityTolerance',1e-6, ...
'StepTolerance', 1e-8, ...
'MaxFunctionEvaluations',1e4);

tic
[parameters_fsolve, fval] = fsolve(gradLogLik, parameters_init, options);
toc

P_logit = 1 ./ (1 + exp(- X * parameters_fsolve));
std_err_fsolve = sqrt(diag(inv(X'*diag(P_logit .* (1 - P_logit))*X)));

% As it is possible to see, parameters of fsolve are different from the ones
% of the fmincon or fminunc function, this is because actually giving to
% matlab te gradient improves both the speed(only 0.059 seconds!) and the accuracy of the
% estimates.
disp("-------------------- FMINCON --------------------");
disp(table(regressor_names, parameters_fmincon, std_err_fmincon, parameters_fmincon./std_err_fmincon, parameters_fmincon-1.96*std_err_fmincon, parameters_fmincon+1.96*std_err_fmincon,...
    'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));
disp("-------------------- FMINUNC --------------------");
disp(table(regressor_names, parameters_fminunc, std_err_fminunc, parameters_fminunc./std_err_fminunc, parameters_fminunc-1.96*std_err_fminunc, parameters_fminunc+1.96*std_err_fminunc,...
    'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));
disp("-------------------- FSOLVE --------------------");
disp(table(regressor_names, parameters_fsolve, std_err_fsolve, parameters_fsolve./std_err_fsolve, parameters_fsolve-1.96*std_err_fsolve, parameters_fsolve+1.96*std_err_fsolve,...
    'VariableNames',{'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}))

%(4) Re-estimate the parameters and standard errors using MATLAB’s glmfit
%function with a logit link.

% remove the constant from X since this function already add them by themselves
X = [data(:,2:11)];

[b_logit, ~, stats_logit] = glmfit(X, Y, 'binomial', 'link', 'logit');

% (5) Estimate a probit model using glmfit with a probit link.
[b_probit, ~, stats_probit] = glmfit(X, Y, 'binomial', 'link', 'probit');

% Standard Errors estimation
std_err_logit = stats_logit.se;  
std_err_probit = stats_probit.se;  

% T-ratios
t_ratio_logit = b_logit ./ std_err_logit;
t_ratio_probit = b_probit ./ std_err_probit;

% Confidence intervals
Cmin_logit = b_logit - 1.96 * std_err_logit;
Cmax_logit = b_logit + 1.96 * std_err_logit;

Cmin_probit = b_probit - 1.96 * std_err_probit;
Cmax_probit = b_probit + 1.96 * std_err_probit;

% (6) Compare numerically the estimated logit and probit coefficients. Do you
% notice any patterns?

% As it is possible to notice, logit coefficient are equal to the previous
% estimated one with fsolve (gradient is very precise!), std_errors are very similar but not identical
% this may be due to the option tolerance chosen before.

disp('-------------------- LOGIT --------------------');
disp(table(regressor_names  ,b_logit, std_err_logit, t_ratio_logit, Cmin_logit, Cmax_logit, ...
     'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));

disp('-------------------- PROBIT --------------------');
disp(table(regressor_names, b_probit, std_err_probit, t_ratio_probit, Cmin_probit, Cmax_probit, ...
    'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));

% As it is possible to see, logit coefficients are almost 1.6 times larger
% than probit ones. This is because probit model uses a gaussian
% distribution and the logit a logistic distribution which has fatter tail.
% Error variances of the logit model are usually normaized to π^2/6, while
% for probit they are normalized to 1.

% The logit coefficients are usually a scaled up version of the probit one
% by a 1.6 factor.
disp('--Logit/Probit Coefficient Ratios--')
disp(b_logit ./ b_probit);

% (7)Estimate the average marginal effect of having high risk aversion on the prob-
% ability of violating FOSD in both the logit and probit models. What do you
% observe?

% to estimate the average marginal effect we only have to take the mean
% of the marginal effect, so take the mean of:
% ∂E[y|xi]/∂xj = p(y = 1|xi)(1−p(y = 1|xi))θj for the logit model.


% Average marginal effect is computed, so in this case we can say that, 
% passing from 0 to 1 in the high risk aversion variable will decrease the
% P(Y = 1|x) of, on average, 25 percentage points. This means that for high risk
% aversion people, violating the FOSD is less likely, as expected.
P_logit = 1 ./ (1 + exp(-[ones(size(X,1),1), X]* b_logit));
AME_logit = mean(P_logit.* (1-P_logit) .* b_logit');
disp('---AME Logit---')
disp(AME_logit)

% Now let's do it with probit, we obtain a very similar result both in the
% sign and in the magnitude. 
% Thus, even though coefficients of logit and probit are scaled up, the
% average marginal effect provides very similar values, and so the same
% information provided, regardless the choice related to the distibution of the errors.
AME_probit = mean(normpdf([ones(size(X,1),1), X] * b_probit) .* b_probit');
disp('---AME Probit---')
disp(AME_probit)

 % (8) Comment on your results.
 % • Response_time_experiment is found to have a significant negative effect
 %   on the probability of FOSD violations. This may
 %   be due to the fact that the greater the time spent, the more the accuracy
 %   in the decision they made, and so less likely will be the violation.

 % • High risk aversion is found to have a significant negative effect too.
 % The AME explains that, according to the data we
 % have, a one unit increase in the variable will decrease the prob of
 % violation of, on average, almost 25 percentage points.
 % Also in here, probably the behavior of the respondent influence the analysis that he makes, and so
 % the violation of FOSD are less likely, leading to a decision
 % in which, say, a lottery "A" always yields higher outcomes than a lottery "B". 
 % Indeed an high risk aversion person will want to minimize its own risk,
 % and so will probably violate the FOSD less likely since he wants the
 % lottery with best outcomes without taking the risk to gamble.
