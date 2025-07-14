%(1) Import data from file RevealedPreferenceData.mat, which contains panel
% scanner data on purchases of margarine. This data set was originally ana-
% lyzed in Allenby and Rossi (Quality Perceptions and Asymmetric Switching
% Between Brands, Marketing Science, 1991), and contains purchases on ten
% brands of margarine by 561 panelists, making a total of 4470 purchases. We
% also have information on various demographic characteristics of each house-
% hold, including household income, family size, education, and if the head of
% the household is retired. The data has 44700 rows and 8 columns:
% •Column 1: Person = Individual identifier.
% •Column 2: Number of occasions.
% •Column 3: Choice = Choice indicator, 0/1; for each individual and each
% choice occasion, it is is a multinomial indicator of one of the 10 brands,
% namely PPkStk,PBBSt, PFlStk, PHseStk, PGenStk, PImpStk, PSSTub,
% PPkTub, PFlTub, PHseTub, where Pk is Parkay; BB is BlueBonnett, Fl
% is Fleischmanns, Hse is house, Gen is generic, Imp is Imperial, SS is
% Shed Spread. Stk indicates stick, and Tub indicates Tub form.
% •Column 4: Income
% •Column 5: Family size
% •Column 6: College educated
% •Column 7: Retired
% •Column 8: Price in dollar

load("RevealedPreferenceData.mat");

% Individual identifier
PERSON=RevealedPreferenceData(:,1); 
% Number of occasions
NOCC=RevealedPreferenceData(:,2); 
% Choice, binary variable
CHOICE=RevealedPreferenceData(:,3); 
% Income
INCOME = RevealedPreferenceData(:,4); 
% Family size
FAMILYSIZE = RevealedPreferenceData(:,5); 
% College
COLLEGE = RevealedPreferenceData(:,6); 
% Retired
RETIRED = RevealedPreferenceData(:,7); 
% Price (dollars $)
PRICE = RevealedPreferenceData(:,8); 

% (2) Write the Log-Likelihood function and the Gradient. Hint: you have a sample
% of N individuals who make repeated choices. Assuming independence, the
% individual likelihood is a product of the densities for each choice; assuming
% i.i.d. samples, the overall likelihood to be maximized is the product of the
% individual ones. The formula for the individual likelihood and the gradient
% (score equations) are given in the lecture notes. To write the Log-Likelihood
% you need to define the utility of option j for individual i in occasion t as
% a function of the covariates and the parameters. The structure of the LL
% function should be as follows: as inputs you should have the parameters to be
% estimated, the set of regressors (covariates) and the set of individual choices;
% as outputs you should have the value of the LL and the value of the gradient.

% The Log-likelihood function should look like something as [LL,G] =...
% loglik($x 0$,covariates,choices), so you have to set not only LL as
% output but also the gradient G as described in the slides. Refer to the Lecture
% Notes 6, example with Stated Preference data for inspiration...

% Since we are in presence of alternative-varying regressors (for example price), the group
% finds a conditional logit model to be more suitable to the given data.

function [LL,G] = loglik(par,COV,CHOICE)
    U=COV * par; % [44700,6] * [6,1] = [44700,1]
    W=reshape(U,10,size(COV,1) / 10); % size 10,4470
    P=exp(W)./sum(exp(W)); % size 10,4470
    p=reshape(P,size(COV,1) ,1); % size 44700,1
    L=log(p(CHOICE==1)); % size 4470,1
    LL=-sum(L); % float
    G=COV'*(p - CHOICE); % Gradient [6,44700] * ([44700,1]-[44700,1]) = [6,1]
end

%(3) Consider two cases:
%(a) CASE 1: The utility of subject i for option j = 1,...,10 in choice task
% t is a function of price and a set of Alternative Specific Constants βj
% Uijt= αPricejt + βj + ϵijt
% where ϵijt are i.i.t. Type 1 Extreme Value (Gumbel) distributed. (Hint:
%To create the ASC’s, you might want to appropriately use the Matlab
%function eye).
%(b) CASE 2: The utility of subject i for option j = 1,...,10 in choice task
%t is a function of price which depends on te set of observed individ-
%ual demographic characteristics (Income, FamilySize, CollegeEducated,
%Retired), and again on the set of Alternative Specific Constants βj
%Uijt= α0Pricejt + αIPricejtIncomei + αFPricejtFamilySizei +....
%αCPricejtCollegeEducatedi + αRPricejtRetiredi + βj + ϵijt
%where ϵijt are again i.i.t. Type 1 Extreme Value (Gumbel) distributed.

%(4) Maximize the LL function using the MATLAB command fmincon, and find
%the point estimates of the parameters, their standard errors and 95% con-
%fidence intervals in both cases (Hint: to find the standard errors, use the
%estimated Hessian, an output of the fmincon function). Report and comment
%your estimates.

% Let's start with CASE 1, for the ASC's we use identity matrix as follows:
ASC = [eye(9); zeros(1,9)];

ASC = repmat(ASC, 4470,1);

COV=[ASC, PRICE]; 
k = size(COV,2); 

% starting point
rng(22);
parameters_init=randn(k,1)/50; 

% lower and upper bound
LB=-50*ones(k,1); UB=50*ones(k,1); 

% options
options=optimoptions('fmincon','FiniteDifferenceType', 'central'); 

% FMINCON  WITHOUT GRADIENT
tic
[betas_withoutGradient,FVAL,~,~,~,~,HESSIAN]=fmincon(@(x)loglik(x,COV,CHOICE),parameters_init,[],[],[],[],LB,UB,[],options); 
toc
std_err_fmincon_withoutGradient = sqrt(diag(inv(HESSIAN)));

% FMINCON WITH GRADIENT
options=optimoptions('fmincon','FiniteDifferenceType', 'central','GradObj','on'); 

tic
[betas_withGradient,FVAL,~,~,~,~,HESSIAN]=fmincon(@(x) loglik(x,COV,CHOICE),parameters_init,[],[],[],[],LB,UB,[],options); 
toc
std_err_fmincon_withGradient = sqrt(diag(inv(HESSIAN)));

regressor_names = ["Const_Parkay_stick", 
                    "Const_BlueBonnett_stick";
                    "Const_Fleischmanns_stick";
                    "Const_house_stick";
                    "Const_generic_stick";
                    "Const_Imperial_stick";
                    "Const_ShedSpread_TubForm";
                    "Const_Parkay_TubForm";
                    "Const_Fleischmanns_TubForm";
                    "Price"];

% As it is possible to notice, both of the methods provides identical
% results for coefficients but slightly different for std_error and thus
% also T-ratios. Furthermore, as expected, the computation with the
% gradient is faster (0.183 seconds vs 1.074 on our laptop)
disp("-------------------- CL-FMINCON WITHOUT GRADIENT --------------------");
disp(table(regressor_names, betas_withoutGradient, std_err_fmincon_withoutGradient, betas_withoutGradient./std_err_fmincon_withoutGradient, betas_withoutGradient-1.96*std_err_fmincon_withoutGradient, betas_withoutGradient+1.96*std_err_fmincon_withoutGradient,...
    'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));
disp("-------------------- CL-FMINCON WITH GRADIENT --------------------");
disp(table(regressor_names, betas_withGradient, std_err_fmincon_withGradient, betas_withGradient./std_err_fmincon_withGradient, betas_withGradient-1.96*std_err_fmincon_withGradient, betas_withGradient+1.96*std_err_fmincon_withGradient,...
    'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));

% Interpretation of results.
% Products with higher constants have higher brand reputation, better
% packaging or other latent variables that aren't possible to measure
% numerically. From this regression, it is possible to notice some insights.
% First, we can notice that named products are more popular, and thus have higher
% coefficients, than the generic "generic_stick" margarine. This may be due to the
% fact that the branded products are more likely to invest in advertising
% or other typologies of differentiation investments. Therefore, insights
% from this regression suggest that consumers are more likely to choose
% products with a known brand, compared to the generic ones.
% Secondly, it's possible to recognize that ceteris paribus the brand,
% consumers are more likely to choose TubForm packaging compared to the
% stick. 
% Price, as expected has a negative coefficient, since an increment of it
% reduces the likelihood of choosing that type of margarine.
% All the coefficients are statistically significant, as it is possible to
% see from the t-ratios and the confidence intervals.


% CASE 2:

% Since the new variables are fixed over the alternatives, and only the
% differences in utility matter, we proceeded by multiplying each new
% regressor with price.
COV_2=[ASC, PRICE, PRICE.*INCOME, PRICE.*FAMILYSIZE, PRICE.*COLLEGE, PRICE.*RETIRED]; 
k_2 = size(COV_2,2); 

% starting point
rng(22);
parameters_init_2=randn(k_2,1)/50; 
 
% lower and upper bound
LB=-50*ones(k_2,1); UB=50*ones(k_2,1); 

% options
options=optimoptions('fmincon','FiniteDifferenceType', 'central'); 

% FMINCON  WITHOUT GRADIENT
tic
[betas_withoutGradient_2,FVAL,~,~,~,~,HESSIAN]=fmincon(@(x)loglik(x,COV_2,CHOICE),parameters_init_2,[],[],[],[],LB,UB,[],options); 
toc
std_err_fmincon_withoutGradient_2 = sqrt(diag(inv(HESSIAN)));

% FMINCON WITH GRADIENT
options=optimoptions('fmincon','FiniteDifferenceType', 'central','GradObj','on'); 

tic
[betas_withGradient_2,FVAL,~,~,~,~,HESSIAN]=fmincon(@(x) loglik(x,COV_2,CHOICE),parameters_init_2,[],[],[],[],LB,UB,[],options); 
toc
std_err_fmincon_withGradient_2 = sqrt(diag(inv(HESSIAN)));

regressor_names = ["Const_Parkay_stick", 
                    "Const_BlueBonnett_stick";
                    "Const_Fleischmanns_stick";
                    "Const_house_stick";
                    "Const_generic_stick";
                    "Const_Imperial_stick";
                    "Const_ShedSpread_TubForm";
                    "Const_Parkay_TubForm";
                    "Const_Fleischmanns_TubForm";
                    "Price";
                    "Price*Income";
                    "Price*FamilySize";
                     "Price*College";
                     "Price*Retired"];

% Also in this case, the function with gradient is faster (0.335 seconds vs 2.209 on our laptop)
disp("-------------------- CL-FMINCON WITHOUT GRADIENT CASE 2 --------------------");
disp(table(regressor_names, betas_withoutGradient_2, std_err_fmincon_withoutGradient_2, betas_withoutGradient_2./std_err_fmincon_withoutGradient_2, betas_withoutGradient_2-1.96*std_err_fmincon_withoutGradient_2, betas_withoutGradient_2+1.96*std_err_fmincon_withoutGradient_2,...
    'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));
disp("-------------------- CL-FMINCON WITH GRADIENT CASE 2--------------------");
disp(table(regressor_names, betas_withGradient_2, std_err_fmincon_withGradient_2, betas_withGradient_2./std_err_fmincon_withGradient_2, betas_withGradient_2-1.96*std_err_fmincon_withGradient_2, betas_withGradient_2+1.96*std_err_fmincon_withGradient_2,...
    'VariableNames', {'Regressor', 'Coefficient', 'Std_Error', 'T_Ratio', 'Cmin', 'Cmax'}));

% Interpretation of results
% Constants are similar to the previous analysis.
% It is possibile to notice from the positive coefficient of price*income
% that this term will mitigate the sensibility of th consumer to the price,
% suggesting that for increasing values of income the consumer will care
% less about the price. On the other hand, family size emphasizes the
% sensibility of price, stating that for larger families, consumers will be
% more influenced by the price in their choice.
% From the price*retired coefficient we can see that price positively
% influence the choice for retired people. 
% We suppose that retired people usually have high savings, therefore they
% may be less sensitive to the price as for the income term. 


% (5) Using the estimated parameter calculate the averaged own-price and cross-
% price elasticities for the 10 goods (a 10 ×10 matrix) in the two cases, and
% comment on your results.

% CASE 1:

% Data preparation
I = eye(10);
N = size(COV,1);
E = zeros(10, 10, N);

U=COV*betas_withGradient; 
W=reshape(U,10,size(COV,1)/10);  
P=exp(W)./sum(exp(W));

for j = 1:10
    for k = 1:10
        for i = 1:size(P,2)
            index = find(NOCC == i);
            E(j,k,i) = (betas_withGradient(end) .* (I(j,k) - P(k,i)) .* PRICE(index(k)));
        end
    end
end

% 

% All the own price elasticities are negative as expected, since the
% probability tho choose a product is decreasing when the price of that
% product is increasing, and viceversa.
% Cross-price elasticities between j and k are equal among them
% (moving j and fixing k), as we saw in slide 24 of lecture 7. 
% The reason behind this results is that Logit cross elasticities are not
% product specific, it is only based on the k characteristics and the
% possible relationship among j and k is not included.
% Cross price elasticities are positive since the products are substitutes.
% We can state this because we are analysing a dataset of margarines and we
% know from microeconomic that this kind of products are substitutes
% If the price of a competitor k increaseas it's more likely that the
% product j will be chosen ( j != k).
E_mean = mean(E, 3)


% CASE 2:
% Basically we did the same things of before but add the other covariates
% Data preparation
N = size(COV_2,1);
E_2 = zeros(10, 10, N);

U=COV_2*betas_withGradient_2; 
W=reshape(U,10,size(COV_2,1)/10);  
P=exp(W)./sum(exp(W));

for j = 1:10
    for k = 1:10
        for i = 1:size(P,2)
            index = find(NOCC == i);
            E_2(j,k,i) = (betas_withGradient_2(10) .* (I(j,k) - P(k,i)) .* PRICE(index(k)));
        end
    end
end

E_mean_2 = mean(E_2, 3)