%(1) Consider the following data from the Pfizer COVID-19 vaccine trial:
% •Vaccinated group: 8 out of 18,198 participants were infected
% •Placebo group: 162 out of 18,325 participants were infected
% You may treat these outcomes as binomial. The vaccine efficacy (VE) is
% defined as: VE = 1−θv/θp

y_v = 8;
n_v = 18198;

y_p = 162;
n_p = 18325;

%(2) Consider an uninformative prior: Beta(1,1) for the binomial success proba-
%bility θ:
%•Compute the posterior distributions of θv and θp analytically and plot
% them
a_prior = 1;
b_prior = 1;

% posterior distribution is computed as a Beta with parameters:
% • y + alpha and 
% • n -y + b

x_v = [0:0.001:0.3];
a_v = y_v + a_prior;
b_v = n_v - y_v + b_prior;


posterior_v = betapdf(x_v, a_v, b_v);

% Almost 40% of mass lies below the MLE 
betacdf(y_v/n_v, a_v, b_v)

x_p = [0:0.001:0.3];
a_p = y_p + a_prior;
b_p = n_p - y_p + b_prior;


posterior_p = betapdf(x_p, a_p, b_p);

% Almost 48% of mass lies below the MLE 
betacdf(y_p/n_p, a_p, b_p)

figure;
plot(x_v, posterior_v, 'b-', 'LineWidth', 2);
hold on;
plot(x_p, posterior_p, 'r-', 'LineWidth', 2);
xlabel('\theta (infection probability)');
ylabel('Posterior Density');
title('Posterior Distributions: Vaccinated vs. Placebo');
legend('Vaccinated', 'Placebo');
grid on;

% •Compute the posterior means, variances, and 95% credible intervals for
% these two posterior distributions

% The means are very close to the MLE estimators y/n = 8/18198 and 162/18325, respectively. This
% result is expected since the posterior estimate is the weighted average
% of the MLE and the prior one, and the prior estimate has always less
% importance as sample size n increases, like in this case.
[M_v,V_v] = betastat(a_v, b_v)

[M_p,V_p] = betastat(a_p, b_p)

CI_v = [x_v(sum(betacdf(x_v, a_v, b_v)<.025)), x_v(sum(betacdf(x_v, a_v, b_v)>.975))]

CI_p = [x_p(sum(betacdf(x_p, a_p, b_p)<.025)), x_p(sum(betacdf(x_p, a_p, b_p)>.975))]


% •Simulate from the posterior distributions of θv and θp, compute the de-
% rived distribution of VE, and plot its histogram.

rng(22);
theta_v_samples = betarnd(a_v, b_v, [1000, 1]);
theta_p_samples = betarnd(a_p, b_p, [1000, 1]);

VE_samples = 1 - (theta_v_samples ./ theta_p_samples);

figure;
histogram(VE_samples);

VE_mean = mean(VE_samples)
VE_CI = quantile(VE_samples, [0.025 0.975])

% (3) Redo the analysis with an informative prior: Beta(0.5,99.5)
% The informative prior should change the results, howver the sample size
% is too large. posterior is more influenced by data (and so more weight to
% the MLE than the prior estimate).

a_prior_inf = 0.5;
b_prior_inf = 99.5;

x_v = [0:0.001:0.3];
a_v = y_v + a_prior_inf;
b_v = n_v - y_v + b_prior_inf;


posterior_v = betapdf(x_v, a_v, b_v);

% Almost 48% of mass lies below the MLE 
betacdf(y_v/n_v, a_v, b_v)

x_p = [0:0.001:0.3];
a_p = y_p + a_prior_inf;
b_p = n_p - y_p + b_prior_inf;

posterior_p = betapdf(x_p, a_p, b_p);

% Almost 52% of mass lies below the MLE 
betacdf(y_p/n_p, a_p, b_p)

figure;
plot(x_v, posterior_v, 'b-', 'LineWidth', 2);
hold on;
plot(x_p, posterior_p, 'r-', 'LineWidth', 2);
xlabel('\theta (infection probability)');
ylabel('Posterior Density');
title('Posterior Distributions: Vaccinated vs. Placebo');
legend('Vaccinated', 'Placebo');
grid on;

% •Compute the posterior means, variances, and 95% credible intervals for
% these two posterior distributions

 % Posterior means and variances 
[M_v_inf,V_v_inf] = betastat(a_v, b_v)
[M_p_inf,V_p_inf] = betastat(a_p, b_p)

CI_v = [x_v(sum(betacdf(x_v, a_v, b_v)<.025)), x_v(sum(betacdf(x_v, a_v, b_v)>.975))]

CI_p = [x_p(sum(betacdf(x_p, a_p, b_p)<.025)), x_p(sum(betacdf(x_p, a_p, b_p)>.975))]

rng(22);
theta_v_samples = betarnd(a_v, b_v, [1000, 1]);
theta_p_samples = betarnd(a_p, b_p, [1000, 1]);

VE_samples = 1 - (theta_v_samples ./ theta_p_samples);

figure;
histogram(VE_samples);

% By using iformative prior, the mean of the estimated vaccine efficacy
% decreased a little, however the intepretation is quite the same
VE_mean_inf = mean(VE_samples)
VE_CI_inf = quantile(VE_samples, [0.025 0.975])

% (4) Interpret the results:
% •What is the estimated vaccine efficacy, and what is its 95% credible
% interval?

% Vaccine efficacy is on average around the 94% for both the analysis,
% with a confidence interval that goes from almost 90% to almost 97%.

% This means that the vaccine has a very high efficacy and very little
% person who vaccinates will be infected.

% •How does the choice of prior affect the inference?
% Actually the posterior mean in our case should be the weighted average
% among the MLE and the prior mean, computed respectively as y/n and
% a/(a+b). However the weights of the average depends upon the sample size
% n, given lambda = n/(n+a +b). So even though prior change, in this case
% it is almost ininfluent since sample size n is very large compared to a
% and b values. 
% The team indeed expect that results may change in the next answer, where the sample size is strongly 
% decreased. 

% (5) Repeat the analysis assuming a smaller sample (10% of the original size)
% (6) Be sure to comment your code and your results

a_prior_inf = 0.5;
b_prior_inf = 99.5;
y_v_smallerSample = 0.1*  y_v;
y_p_smallerSample = 0.1*  y_p;
n_v_smallerSample = 0.1*  n_v;
n_p_smallerSample = 0.1 * n_p;

x_v = [0:0.001:0.3];
a_v = y_v_smallerSample + a_prior_inf;
b_v = n_v_smallerSample - y_v_smallerSample + b_prior_inf;


posterior_v = betapdf(x_v, a_v, b_v);
% Almost 53% of mass lies below MLE
betacdf(y_v_smallerSample/n_v_smallerSample, a_v, b_v)

x_p = [0:0.001:0.3];
a_p = y_p_smallerSample + a_prior_inf;
b_p = n_p_smallerSample - y_p_smallerSample + b_prior_inf;

posterior_p = betapdf(x_p, a_p, b_p);


% Almost 75% of mass lies below MLE, we can start to see that results are changing
betacdf(y_p_smallerSample/n_p_smallerSample, a_p, b_p)

figure;
plot(x_v, posterior_v, 'b-', 'LineWidth', 2);
hold on;
plot(x_p, posterior_p, 'r-', 'LineWidth', 2);
xlabel('\theta (infection probability)');
ylabel('Posterior Density');
title('Posterior Distributions: Vaccinated vs. Placebo');
legend('Vaccinated', 'Placebo');
grid on;

% •Compute the posterior means, variances, and 95% credible intervals for
% these two posterior distributions

 % Posterior means and variances are close to the previous results,
 % this is due to the fact that on one hand we decreasded the sample size
 % and so the weights are different but they are already very unbalanced).

[M_v_smallerSample,V_v_smallerSample] = betastat(a_v, b_v)
[M_p_smallerSample,V_p_smallerSample] = betastat(a_p, b_p)

% Let's check the weights using the two samples 
weight_data = n_v/(n_v + 0.5 + 99.5);
weight_prior = 1-weight_data;

weight_data_smallerSample = n_v_smallerSample/(n_v_smallerSample + 0.5 + 99.5);
weight_prior_smallerSample = 1-weight_data_smallerSample;

% As it is possible to notice, even though the sample size decreases the
% weights are still unbalanced: before the weights were almost 99.5% and
% 0.5%, now with smaller sample size the weights are 95% and 5%. Therefore
% the results are closer because the weights are stil unbalanced. Hoewever,
% it is identified that weights tends to be more balanced for decreasing
% values of the sample size.

% Let's check the mean computation 
Posterior_mean_v = weight_data * y_v_smallerSample/n_v_smallerSample + weight_prior * 0.5/100
Posterior_mean_p = weight_data * y_p_smallerSample/n_p_smallerSample + weight_prior * 0.5/100
% the results are equal to "M_v_smallerSample" and "M_p_smallerSample" as expected, check satisfied.

CI_v = [x_v(sum(betacdf(x_v, a_v, b_v)<.025)), x_v(sum(betacdf(x_v, a_v, b_v)>.975))]

CI_p = [x_p(sum(betacdf(x_p, a_p, b_p)<.025)), x_p(sum(betacdf(x_p, a_p, b_p)>.975))]

rng(22);
theta_v_samples = betarnd(a_v, b_v, [0.1*1000, 1]);
theta_p_samples = betarnd(a_p, b_p, [0.1*1000, 1]);

VE_samples = 1 - (theta_v_samples ./ theta_p_samples);

figure;
histogram(VE_samples);

% By using iformative prior and smaller sample size, the mean of the estmated vaccine effcacy
% decreased a little, however the intepretaion is quite the same
VE_mean_smallerSample = mean(VE_samples)
VE_CI_smallerSample = quantile(VE_samples, [0.025 0.975])


% The results are close to the previous, this is mainly
% due to the fact that posterior distribution is influenced by both the
% data and the prior parameters, but the weights are again unbalanced.

% Looking at the plots, it is possible to highlight that, as already
% argued in the analysis, the VE distribution does not change too much
% among analysis, therefore we are in presence of a Vaccine with almost 
% on average 95% of efficacy. 
% Instead, regarding the analysis of the beta distributions coming from
% the bayesian procedure, it is possible to analyze that mean of the posterior are very close to
% their respective MLE, however the Vaccinated group for the entire sample size presents
% higher density than the one analyzed in the smallest sample size. 
