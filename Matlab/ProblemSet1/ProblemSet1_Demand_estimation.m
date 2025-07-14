% DEMAND ESTIMATION

%The data set Data Broiler.mat contains aggregate data on quantity, price, cost,
%and demographic variables related to edible meat of young chickens (”broilers”) in
%the United States from 1962 to 1999.
%There are 12 variables in the dataset:
%(1) Year
%(2) Quantity of broiler chicken
%(3) Income
%(4) Price of broiler chicken
%(5) Price of beef
% (6) Price of corn
% (7) Price of chicken feed
% (8) Consumer Price Index
% (9) Aggregate production of chicken in pounds
% (10) US population
% (11) Exports of beef, veal, and pork in pounds
% (12) Time trend

%load the dataset
load('Data_Broiler.mat')


%(1) Regress the quantity of broiler chicken on the price of broiler chicken (and
%a constant), excluding all other variables. Interpret and report your results
%(coefficients, standard errors, and t-statistics). What exactly are you esti-
% mating by running OLS on this model? What is the interpretation of the
% price coefficient, and what do you make of its sign?

% We are estimating the coefficient b of the price such that the squared
% residuals are minimized, the b coefficient can be interpreted as how much
% the demand changes after an increase in price of 1€ (or 1$, depends on the
% unit of the dataset).


% Assign to y the Quantity of broiler chicken and
% to x the Price of broiler chicken (and a constant)
y = Data_broiler(:, 2);
N = length(Data_broiler);
const = ones(N,1);
X = [const, Data_broiler(:, 4)];
k = size(X, 2);
% Let's estimates beta_hat from the least squares:
% in particular let's put d(MSE)/dB = 0
% -> X'(-2y + 2XB) = 0 -> B = (X'X)^-1 * X'y
B =(X'*X)\X'*y; 

% calculate the residuals ad the difference beyween the true value and the
% estimated one
e = y - X*B;

% Let's now calculate std. errors, we will use both homoskedastic
% std.errors and heteroskedastic ones thanks to thw white formula for
% variance estimation

% Ressidual variance calculated as the mean of the squared residuals
% over N-k degrees of freedom.
s_2=e'*e/(N-k); 

% Homoskedastc std.error calculated as the square root of the diagonal
% of the variance-covariance matrix, we take the diagonal since we are
% interested in the variance and then take the square root
Var_hom=s_2*eye(k)/(X'*X); 
std_err_hom=sqrt(diag(Var_hom));

% Variance estimatation using robust std.errors, acording to white formula
Var_het=(X'*X)\(X'*diag(e.^2)*X)/(X'*X);
std_err_white=sqrt(diag(Var_het));

% Let's now visualize the values of the parameters estimation, their
% confidence interval and the t-statistic

% confidence intervals calculation
Cmin=B - 1.96*std_err_white;
Cmax=B + 1.96*std_err_white;

% Let's create a summary table in which we can observe what are we
% estimating -> the parameters of the regression such that least squares
% are obtained from data, its std.err and the significance level coming
% from t-statistic.

% std.error estimates are of course very important for two reasons:
% • first of all they allow the reader to look at the confidence interval
% of the estimate B, and so consider the degree of variability.
% • furthermore, std_error are fundamental to calculate t-statistic,
% allowing the analyst to understand if the B estimation is significantly
% different from zero (reject null hypotesis H0) or not (not reject H0).


% In particular we expect that the coefficient of price will be negative, since it is quite obvious that
% price and quantity are charachterized by a negative relationship 
disp('-------------------OLS-Response: Quantity of broiler chicken ------------------');
regressor_names = ["Intercept"; "Price of broiler chicken"];
disp(table(regressor_names, B, Cmin, Cmax, std_err_hom, B./std_err_hom, std_err_white, B./std_err_white, ...
    'VariableNames', {'Regressor', 'B','Cmin','Cmax','std_err_hom','t_stat_hom',  'std_err_white', 't-statisic_white'}));

SS_tot = sum((y - mean(y)).^2);      
SS_res = sum((e).^2);      
R2 = 1 - (SS_res / SS_tot);
R2_adj = 1 - ((SS_res / (N - k)) / (SS_tot / (N - 1)));
disp("R2 = ")
disp(R2)
disp("Adjusted R2 = ")
disp(R2_adj)

% As it is possible to notice from the table, the results found are not
% as expected, in particular for one unit increase in the price, demand of broiler chicken will increase of 
% 0.25, this is a particular weird result since usually in a demand curve the equation can be like: q = a - b*q, 
% while instead in here we have a positive sign.


%(2) Re-estimate the OLS model, now including all available demand covariates.
% That is, include x, which contains: income, price of beef, consumer price
% index, aggregate production of chicken in pounds, US population, exports of
% beef/veal/pork in pounds, and time trend. Report the estimated coefficients
% and standard errors, and comment on what happened after including these
% additional variables. What happened to the fit of the model overall? What
% about the precision of the individual estimates? Do the coefficients have the
%expected signs?

% Let's start by assigning each variable with the related column of the dataset
income = Data_broiler(:, 3);
price_broiler_chicken = Data_broiler(:, 4);
price_beef = Data_broiler(:, 5);
consumer_price_index = Data_broiler(:, 8);
aggr_prod_chicken = Data_broiler(:, 9);
us_population = Data_broiler(:, 10);
exports_bvp = Data_broiler(:, 11);
time_trend = Data_broiler(:, 12);

X = [const, income, price_broiler_chicken, price_beef, consumer_price_index, aggr_prod_chicken, us_population, exports_bvp, time_trend];
k = size(X, 2);

% Basically, let's re-determine all the things that we did before but for
% more than one covariates.
B =(X'*X)\X'*y; 

e = y - X*B;

s_2=e'*e/(N-k); 

Var_hom=s_2*eye(k)/(X'*X); 
std_err_hom=sqrt(diag(Var_hom));

Var_het=(X'*X)\(X'*diag(e.^2)*X)/(X'*X);
std_err_white=sqrt(diag(Var_het));

Cmin=B - 1.96*std_err_white;
Cmax=B + 1.96*std_err_white;

disp('-------------------OLS-Response:Quantity of broiler chicken------------------');
regressor_names = ["Intercept"; "income"; "price_broiler_chicken"; "price_beef"; "consumer_price_index"; "aggr_prod_chicken";  "us_population"; "exports_bvp"; "time_trend"];
disp(table(regressor_names, B, Cmin, Cmax, std_err_hom, B./std_err_hom, std_err_white, B./std_err_white, ...
    'VariableNames', {'Regressor', 'B','Cmin','Cmax','std_err_hom','t_stat_hom',  'std_err_white', 't-statisic_white'}));

% Let's now determine R^2 to calculate how much variance is explained by
% the model and so the data and the parameters of the model fit the
% response variable
SS_tot = sum((y - mean(y)).^2);      
SS_res = sum((e).^2);      
R2 = 1 - (SS_res / SS_tot)
R2_adj = 1 - ((SS_res / (N - k)) / (SS_tot / (N - 1)))

% As is is possible to notice, now both the R^2 and the R^2 adj increased,
% so adding parameters to the regression has increased the goodness of fit.

% We can say that, adding variables, the precision of the estimate for price decreased(and also the t-ratio, thus also significance).
% However this result is more reliable since we know that, in the context
% of the analysis, the price should be negative related with demand.

% INTERPRETATION OF THE RESULTS
% Let's now interpret the results by taking into account the significant
% variables (to understand the significance we can either see whether 
% the confidence interval of the parameters include the 0 or not or look at the t-ratio,
% (usually values close to |2| explains significance levels of the parameter)
% in particular:

% • As income increases, demand increases as it is quite obvious. Consumers
% will have more money to spend for their goods.

% • price_broiler_chicken now provides better interpretation than before,
% suggesting a negative relationship with the response, this is the result that he team expected to see.
% A marginal increase in price will negatively influence the demand, decreasing of almost 0.058 points. 

% • price_beef presents a positive significant sign, this result is
% expected since beef and broiler_chicken are substitute products, and so
% increasing the price of one of the two products will increase the quantity of the
% other one -> positive cross-elasticity.

% • consumer_price_index presents instead a negative significant sign, this
% is explainable by the fact that increasing the cost of the 'basket of goods' at the
% current period, and so increase the CPI, will result in a lower demand as
% the results tell us.

% • aggregate_production_chicken in pound has a significant positive sign of the estimate 
% and very large parameter beta. This can be explained by the fact that the
% higher the supply of chicken available in the market, the lower will be
% the price and so of course the demand for chicken will increase.

% • finally, exports_bvp has a negative significant sign. This may be
% explained by the fact that the higher the level of exports, 
% the lower the availability of substitutes product in the market, this may
% cause an higher price of the chicken, leading to a decrease of the demand.



%(3) To address the potential endogeneity of price, use as instrumental variables
% the price of corn, the price of chicken feed, and their interaction. Explain
% the rationale behind this choice. Estimate the model using two-stage least
% squares (2SLS) with these three instruments. Comment on your results.

% When we have to handle endogeneity, usually we have to find an 
% instrumental variable(IV) such that the cov(IV, X) != 0 
% and cov (IV, e) = 0, this means that we should find a variable
% that actually affects X (relevance condition) but not directly affect the
% dependent variable, for this purpose often a suitable strategy is to use a linear combination of the variables
% as IV. However this will result in creating an *overidentified* model, i.e. we
% have more instruments than endogenous regressors.

% Endogeneity occurs because of simultaneity and omitted variables, since
% it is very likely to have variables not present in the dataset that are
% correlated with price, violating exogeneity conditions. Simultaneity
% occurs because supply and demand are determined jointly adn simultaneously.

% Since we are in a context of overdentified model, we should use a two
% stage least squares model.

price_corn = Data_broiler(:, 6);
price_chicken_feed = Data_broiler(:, 7);

% So, the first step is to regress X on X_IV and check
% for the relevance condition, instruments should be significant on X

% Let's start initializing the matrix containing the instruments, called
% X_IV

% UNRESTRICTED
X_IV = [const,income, price_corn, price_chicken_feed, price_corn.*price_chicken_feed, price_beef, consumer_price_index, aggr_prod_chicken, us_population, exports_bvp, time_trend];
beta_iv = (X_IV'*X_IV)^-1 * X_IV'*price_broiler_chicken;
e_unres= price_broiler_chicken - X_IV*beta_iv;

s_2=e'*e/(N-k); 

Var_het=(X_IV'*X_IV)\(X_IV'*diag(e.^2)*X_IV)/(X_IV'*X_IV);
std_err_white=sqrt(diag(Var_het));

Cmin=beta_iv - 1.96*std_err_white;
Cmax=beta_iv + 1.96*std_err_white;

% unrestricted residual sum of squares
URSS = sum(e_unres.^2);
% 
k_unres = size(X_IV, 2);

% RESTRICTED
X_res = [const,income, price_beef, consumer_price_index, aggr_prod_chicken, us_population, exports_bvp, time_trend];
beta_res = (X_res'*X_res)^-1 * X_res'*price_broiler_chicken;
e_res= price_broiler_chicken - X_res*beta_res;
RRSS = sum(e_res.^2);

% numer of restriction q = 3 = k_unres - k_res
q = k_unres - size(X_res, 2);

% As a rule of thumb, usually F_stat should be at least 10 to have strong
% instruments, in this case they may be weak.
F_stat=((RRSS-URSS)/q)/(URSS/(N-k_unres))
P = fcdf(F_stat,q,N-k_unres)

disp('-------------------2SLS-1 Stage: Relevance condition of the instruments on Price of broiler chicken ------------------');
regressor_names = ["Intercept"; "income"; "price_corn"; "price_chicken_feed"; "price_corn*price_chicken_feed"; "price_beef"; "consumer_price_index"; "aggr_prod_chicken"; "us_population"; "exports_bvp"; "time_trend"];
disp(table(regressor_names, beta_iv, Cmin, Cmax, std_err_white, beta_iv./std_err_white, ...
    'VariableNames', {'Regressor', 'beta_iv','Cmin','Cmax', 'std_err_white', 't-statisic_white'}));

% We don't pass the rule of thumb of an F-stat higher than 10, but we get
% 5.8. This means that the instruments are significant but they may be
% not so strong.
% Since the instruments parameters(both the interaction and the single instruments) are significantly different from zero in the
% X^hat estimation, so influence the endogenous regressor, in particular they positive influence taking into account 
% both the direct effect and the interaction effect.
% The relevance condition is now tested and so it is possible to do 2SLS.

% Obviously, this model works only with the exogeneity assumption,
% instruments are in fact taken in order to have an E(X_IV, e) = 0.

% Now let's move forward into the 2SLS:
X = [const, income, price_broiler_chicken,price_beef, consumer_price_index, aggr_prod_chicken,  us_population, exports_bvp, time_trend];
beta_1_stage = (X_IV'*X_IV)^(-1)*X_IV'*X;
Xhat= X_IV * beta_1_stage;

beta_2_stage =(Xhat'*Xhat)^(-1)*(Xhat'*y);
e=y-X*beta_2_stage;
s_2 = mean(e.^2)';

% std.error determination procedure
QX_IVX=X_IV'*X/N; QXX_IV=X'*X_IV/N;QX_IVX_IV=X_IV'*X_IV/N;

% homoskedasticity
V=inv(QXX_IV*inv(QX_IVX_IV)*QX_IVX)*s_2;
se2SLS=sqrt(diag(V/N));

% White covariance matrix estimate
O=(X_IV'*(diag(e'*e)/N)*X_IV)/N;
V_white=inv(QXX_IV*inv(QX_IVX_IV)*QX_IVX)*(QXX_IV*inv(QX_IVX_IV)*O*inv(QX_IVX_IV)*QX_IVX);

se2SLS_white=sqrt(diag(V_white/N));
tstat2SLS_white=beta_2_stage./se2SLS_white;

Cmin2SLS=beta_2_stage - 1.96*se2SLS;
Cmax2SLS=beta_2_stage + 1.96*se2SLS;


disp('-------------------2SLS-2Stage ------------------');
regressor_names = ["Intercept"; "income"; "IV_price_broiler_chicken"; "price_beef"; "consumer_price_index"; "aggr_prod_chicken"; "us_population"; "exports_bvp"; "time_trend"];
disp(table(regressor_names, beta_2_stage, se2SLS, se2SLS_white, tstat2SLS_white, Cmin2SLS, Cmax2SLS, ...
    'VariableNames', {'Regressor', 'beta2_stage', 'se2SLS', 'se2SLS_white', 'tstat2SLS_white', ' Cmin2SLS', 'Cmax2SLS'}));

% INTERPRETATION OF THE RESULTS
% Using 2SLS, it is possible to see that results changed a bit, in
% particular:
% • IV_price_broiler_chicken now significantly and negative influence the demand
% of chicken as it is normal to think. 
% • price_beef, aggregate_production_chicken and
% exports_bpv maintained the same results of question (2)
% • CPI lost its significance
% • US_population has a significant negative sign in decreasing demand,
% this may be due to the fact that as population increases, there are food
% alternatives becoming more popular.

% (4) Plot the price elasticity of demand for broiler chicken over the sample years.

price_elast_IV=beta_2_stage(3)*(Xhat(1:N,3))./(y(1:N));
plot([1960:1999],price_elast_IV);

% As we know, elasticity is not constant for this models but must be
% calculated for each observation, constant elasticty can in fact be seen
% only in log-log models by construction.

% As it is possible to notice, price elasticity is quite rigid until 1972.
% After that year it starts to be more elastic, swinging in the following
% years.


