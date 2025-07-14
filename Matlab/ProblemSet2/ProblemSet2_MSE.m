%(1) Generate n = 1000 independent samples of size m = 10, where each obser-
% vation Yij ∼N(θ,σ2), with known true mean θ = 1 and standard deviation
% σ = 1. For each sample (i.e., each row of your simulated matrix), compute
% the sample mean estimator θˆi
rng(22)

n =1000;
m = 10;
theta_true = 1;
sd = 1;

Y = theta_true + sd^2*randn(1000, 10);

% First estimator is the sample mean
unbiased_estimator = mean(Y, 2);

% (2) Implement two alternative estimators:
% •Shrinkage estimator:
% θi = λˆθi with λ= 0.8
% •Constant estimator: θ^const_i = 0 for all i

lambda = 0.8;
% Let's introduce a little of bias to reduce variance 
shrinkage_estimator = lambda .* unbiased_estimator;

% The constant estimator will have always the same value for each
% observation, as we know this estimator may be seen as a stupid one.
% However, this estimator is actually admissible, according to the
% definition in fact the estimator is not outperformed by the others, 
% i.e. there will exist at least one region in which the constant estimator
% provides better MSE (thus lower) than the others (In reasonable terms,
% probably the region is the one closest to the constant itself).
constant_estimator = zeros(1000, 1);

%(3) For all three estimators (sample mean, shrinkage, and constant): Compute
% the bias, variance, and mean squared error (MSE)

% Let's start from the first estimator: 

% UNBIASED ESTIMATOR

% Bias of the unbiased estimator, should be very close to 0
% AS it is possible to notice, bias is computed as 0.0037 so the results
% are aligned with our expectations.
bias_unbiased_estimator = mean(unbiased_estimator - theta_true);
% Variance is instead computed as follows
var_unbiased_estimator = var(unbiased_estimator) % sigma/m = 0.10
% Variance decreases when m increases.

% SHRINKAGE ESTIMATOR

% Bias of the shrinkage estimator should be close to 0.8 -1 = -0.2
% As it is possible to notice, bias is computed as -0.1970 so the results
% are aligned with our expectations.
bias_shrinkage_estimator = mean(shrinkage_estimator - theta_true);

% As expected, the variance if the shrinkage estimator is lower than the 
% unbiased estimator, having introduced bias has decreased the variance.
var_shrinkage_estimator = var(shrinkage_estimator); % 0.64 * 0.10 = 0.064

% According to the shrinkage estimator properties, the variance and bias
% can also be computed as follows:
% As it is possible to notice, results are very aligned with the ones
% computed previously.
var_shrinkage_estimator_prop = lambda^2*var_unbiased_estimator;
bias_shrinkage_estimator_prop = (1-lambda)*theta_true;

% CONSTANT ESTIMATOR

% In this case the estimator may be seen as a shrinkage with lambda = 0 

% Results are aligned with expectation, bias_constant = -1
bias_constant_estimator = mean(constant_estimator - theta_true);
% Of course, being the estimator constant, the variance is 0
var_constant_estimator = var(constant_estimator);

MSE_unbiased = bias_unbiased_estimator^2 + var_unbiased_estimator;
MSE_shrinkage = bias_shrinkage_estimator^2 + var_shrinkage_estimator; 
MSE_constant = bias_constant_estimator^2 + var_constant_estimator; 


% We can check the correct estimation of the MSE by using the following
% equality, and computing the first term:
% E[( θˆ − θ)^2] = Bias2(θˆ) + Var(θˆ)

MSE_unbiased_prop = mean((unbiased_estimator - theta_true).^2)
MSE_shrinkage_prop = mean((shrinkage_estimator - theta_true).^2)
MSE_constant_prop = mean((constant_estimator - theta_true).^2)

% Results are very close to the ones computed previously, as expected. 

% (4) Present the results in a table with three rows and three columns: Bias2
% , Variance, MSE

estimator_names = ["Unbiased";
                   "Shrinkage";
                   "Constant"];

bias = [bias_unbiased_estimator;
        bias_shrinkage_estimator;
        bias_constant_estimator];

variance = [var_unbiased_estimator;
            var_shrinkage_estimator;
            var_constant_estimator];

MSE = [MSE_unbiased;
       MSE_shrinkage;
       MSE_constant];

% As it is possible to notice from the following table, the best estimator
% is the unbiased one, however this is not necessarily true in general.
% Sometimes it may be useful to introduce a little of bias to reduce
% variance and thus minimize MSE.

disp("-------------------- Summary Bias-Variance-MSE --------------------");
disp(table(estimator_names, bias, variance, MSE,...
    'VariableNames', {'Estimator', 'Bias', 'Variance', 'MSE'}));

%(5) Plot the empirical distribution (histogram or density) of the three estimators,
% labeling the estimators in the plot

% Let's compute density using Kernel smoothing function estimate

x = -0.5:0.001:1.5;
unbiased_pdf = normpdf(x, mean(unbiased_estimator), var_unbiased_estimator^.5);
shrinkage_pdf = normpdf(x, mean(shrinkage_estimator), var_shrinkage_estimator^.5);

figure;
plot(x, unbiased_pdf, 'b-', 'LineWidth', 2); 
hold on;
plot(x, shrinkage_pdf, 'r--', 'LineWidth', 2);
xline(0, 'g--', 'LineWidth', 2); 

xline(theta_true, 'k', 'LineWidth', 2, 'DisplayName', 'True \theta');

title('Empirical Distribution of Three Estimators');
xlabel('Estimated \theta');
ylabel('Density');
legend({'Unbiased estimator', ...
        'Shrinkage estimator', ...
        'Constant estimator', ...
        'True \theta = 1'}, 'Location', 'northwest');
grid on;

%(6) Which estimator has the smallest MSE? Repeat the analysis with m= 5 and
% m = 50. What do you observe about the bias-variance tradeoff as sample
% size increases?
%(7) Be sure to comment your code and your results

% As already said, the smallest MSE belong to the unbiased estimator for m = 10.

% M = 5

% Now, let's do all the things we did before but with m = 5, given that the
% constant estimator doe not change over sample size, we compute only the
% new unbiased and shrinkage estimators.

rng(22);

n =1000;

Y_5 = theta_true + sd^2*randn(1000, 5);

unbiased_estimator_5 = mean(Y_5, 2);
shrinkage_estimator_5 = lambda .* unbiased_estimator_5;


bias_unbiased_estimator_5 = mean(unbiased_estimator_5 - theta_true);
var_unbiased_estimator_5 = var(unbiased_estimator_5);
bias_shrinkage_estimator_5 = mean(shrinkage_estimator_5 - theta_true);
var_shrinkage_estimator_5 = var(shrinkage_estimator_5);

MSE_unbiased_5 = bias_unbiased_estimator_5^2 + var_unbiased_estimator_5;
MSE_shrinkage_5 = bias_shrinkage_estimator_5^2 + var_shrinkage_estimator_5;

bias_5 = [bias_unbiased_estimator_5;
        bias_shrinkage_estimator_5;
        bias_constant_estimator];

variance_5 = [var_unbiased_estimator_5;
            var_shrinkage_estimator_5;
            var_constant_estimator];

MSE_5 = [MSE_unbiased_5;
       MSE_shrinkage_5;
       MSE_constant];

disp("-------------------- Summary Bias-Variance-MSE with m = 5 --------------------");
disp(table(estimator_names, bias_5, variance_5, MSE_5,...
    'VariableNames', {'Estimator', 'Bias', 'Variance', 'MSE'}));

unbiased_pdf_5 = normpdf(x, mean(unbiased_estimator_5), var_unbiased_estimator_5^.5);
shrinkage_pdf_5 = normpdf(x, mean(shrinkage_estimator_5), var_shrinkage_estimator_5^.5);

figure;
plot(x, unbiased_pdf_5, 'b-', 'LineWidth', 2); 
hold on;
plot(x, shrinkage_pdf_5, 'r--', 'LineWidth', 2);
xline(0, 'g--', 'LineWidth', 2); 


xline(theta_true, 'k', 'LineWidth', 2, 'DisplayName', 'True \theta');

title('Empirical Distribution of Three Estimators, m = 5');
xlabel('Estimated \theta');
ylabel('Density');
legend({'Unbiased estimator', ...
        'Shrinkage estimator', ...
        'Constant estimator', ...
        'True \theta = 1'}, 'Location', 'northwest');
grid on;

% M = 50

rng(22);

n =1000;

Y_50 = theta_true + sd^2*randn(1000, 50);

unbiased_estimator_50 = mean(Y_50, 2);
shrinkage_estimator_50 = lambda .* unbiased_estimator_50;


bias_unbiased_estimator_50 = mean(unbiased_estimator_50 - theta_true);
var_unbiased_estimator_50 = var(unbiased_estimator_50);
bias_shrinkage_estimator_50 = mean(shrinkage_estimator_50 - theta_true);
var_shrinkage_estimator_50 = var(shrinkage_estimator_50);

MSE_unbiased_50 = bias_unbiased_estimator_50^2 + var_unbiased_estimator_50;
MSE_shrinkage_50 = bias_shrinkage_estimator_50^2 + var_shrinkage_estimator_50;

bias_50 = [bias_unbiased_estimator_50;
        bias_shrinkage_estimator_50;
        bias_constant_estimator];

variance_50 = [var_unbiased_estimator_50;
            var_shrinkage_estimator_50;
            var_constant_estimator];

MSE_50 = [MSE_unbiased_50;
       MSE_shrinkage_50;
       MSE_constant];

disp("-------------------- Summary Bias-Variance-MSE with m = 50 --------------------");
disp(table(estimator_names, bias_50, variance_50, MSE_50,...
    'VariableNames', {'Estimator', 'Bias', 'Variance', 'MSE'}));
unbiased_pdf_50 = normpdf(x, mean(unbiased_estimator_50), var_unbiased_estimator_50^.5);
shrinkage_pdf_50 = normpdf(x, mean(shrinkage_estimator_50), var_shrinkage_estimator_50^.5);

figure;
plot(x, unbiased_pdf_50, 'b-', 'LineWidth', 2); 
hold on;
plot(x, shrinkage_pdf_50, 'r--', 'LineWidth', 2);
xline(0, 'g--', 'LineWidth', 2); 

xline(theta_true, 'k', 'LineWidth', 2, 'DisplayName', 'True \theta');

title('Empirical Distribution of Three Estimators, m = 50');
xlabel('Estimated \theta');
ylabel('Density');
legend({'Unbiased estimator', ...
        'Shrinkage estimator', ...
        'Constant estimator', ...
        'True \theta = 1'}, 'Location', 'northwest');
grid on;

% Analyzing the obtained results, it is possible to notice that in the
% different analysis the best estimator (thus, the one that minimizes MSE) 
% is the Shrinkage one for m = 5 and again the unbiased one for m = 50.
% Regardless the sample size, the Shrinkage estimator provides always more
% bias but less variance than the unbiased one, and in particular for m = 5 
% the gain in variance reduction is greater 
% than the loss in bias increment, resulting in lower MSE than Unbiased
% estimator.
% Finally, as it is possible to see, the higher the sample size the lower
% the bias, resulting in values closer to 0 as m increases. Also the
% variance of the sample mean decreases of course with an increment in
% sample size, as expected by its formula, resulting in lower MSE.

% It is possible to get the same conclusions also looking at the plots.
% For which regards instead the plots, it is possible to notice that as m
% increases, the first two estimators presents tight tails,
% decreasing their variance as already argued before, the constant estimator
% instead provides always the same result since it is independent to the
% sample mean.