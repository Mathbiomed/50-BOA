function Lambda = CV_Inhibition(file_name, L, isSSRE)
%% Default values of input variable
if nargin < 3
    isSSRE = true;
end

%% Load data
data = readmatrix(file_name);
%Data check
if (height(data) < 2 || width(data) < 4)
    error('Invalid Input')
end

%Load the formatted data
Vmax = data(1,1); Ks = data(1,2); IC50 = data(1,3); St_IC50 = data(1,4);
St_setup = data(2:end,1); It_setup = data(2:end,2); V0 = data(2:end,3); a = data(2,4);

X_setup = [St_setup It_setup]; C = [Vmax Ks a]; IC50s = [St_IC50 IC50];

%Cross-validation to select regularization constant
cv_value = zeros(1, numel(L));

for i = 1:numel(L)
    r = L(i);
    cv_value(i) = CV_error(X_setup, V0, C, IC50s, r, isSSRE);
end

best_r_idx = cv_value == min(cv_value);
best_r = L(best_r_idx);
if numel(best_r) > 1
    best_r = best_r(1);
end

Lambda = best_r;
end

%% Inhibition model
function v = Inhibition(K, X, C)
v = C(1)*(C(2)*X(:,1).*(1+X(:,2)/K(2))+(1+X(:,2)/(K(2)^2/K(1))).*X(:,1).^2/C(3))...
    ./(C(2)^2*(1+X(:,2)/K(1))+2*C(2)*X(:,1).*(1+X(:,2)/K(2))+(1+X(:,2)/(K(2)^2/K(1))).*X(:,1).^2/C(3));
end

%% Cheng-Prusoff equation
function v = Cheng_Prusoff(K, X, C)
h1 = @(S, a, x) a*S*x + S^2;
h2 = @(S, a, x) a*x^2 + 2*a*S*x + S^2;

v = K(1)*(K(2)/K(1))^2*(h1(X,C(2),C(1))*h2(X,C(2),C(1)))...
    /(h1(X,C(2),C(1))*h2(X,C(2),C(1)*K(2)/K(1))-2*h1(X,C(2),C(1)*K(2)/K(1))*h2(X,C(2),C(1)));
end

%% Error structure
function s = Error_Structure(X, Y, isSSRE)
if isSSRE
    s = Y;
else
    % User can put his or her own model of stndard deviation here.

    % Model 1
    % K1 = 0.011; K2 = 0.04;
    % s = K1 + K2*Y;

    % Model 2
    % K1 = 0.036; K2 = 0.008;
    % s = K1 + K2*Y.^2;

    % Model 3
    % K1 = 0.05; K2 = 0.889;
    % s = K1*Y.^K2;

    % Model 4
    K1 = 0.259; K2 = 0.529; K3 = 0.699; K4 = 0.561; K5 = 0.505;
    s = K1*X(:,1).^K2./(X(:,1).^K3 + K4*X(:,2).^K5);

    % Model 5
    % K1 = 0.308; K2 = 0.063; K3 = 0.273; K4 = 0.095;
    % s = K1*X(:,1)./(K2 + X(:,1) + K3*X(:,1).^2 + K4*X(:,2));

    % Model 6
    % K1 = 0.308; K2 = 0.063; K3 = 0.273; K4 = 0.095; K5 = 0.007;
    % s = (K1*X(:,1) + K5*X(:,2))./(K2 + X(:,1) + K3*X(:,1).^2 + K4*X(:,2));
end
end
%% Loss function with lambda
function loss = CV_loss(K, X, Y, C, IC50s, lambda, isSSRE)
Y_predict = Inhibition(K, X, C);
loss = mean(((Y-Y_predict)./Error_Structure(X, Y, isSSRE)).^2) +...
                 lambda*mean(((IC50s(:,2)-Cheng_Prusoff(K, IC50s(:,1), [C(2) C(3)]))./IC50s(:,2)).^2);
end

%% Fitting
function params = Fit_inhibition(X, Y, C, IC50s, lambda, isSSRE)
K0 = [max(IC50s(:,2)) max(IC50s(:,2))];
objFun = @(K)CV_loss(K, X, Y, C, IC50s, lambda, isSSRE);

options = optimset('Display', 'off');
params = fminsearch(objFun, K0, options);
end

%% Test error
function loss = Test_error(Xtrain, Ytrain, Xtest, Ytest, C, IC50s, lambda, isSSRE)
params = Fit_inhibition(Xtrain, Ytrain, C, IC50s, lambda, isSSRE);
Ypredict = Inhibition(params, Xtest, C);
loss = mean((Ytest - Ypredict).^2);
end

%% Cross-validation error
function loss = CV_error(X, Y, C, IC50s, lambda, isSSRE)
%Leave-one-out
cv = cvpartition(height(X), 'LeaveOut');

loss_set = zeros(1,cv.NumTestSets);

for i = 1:cv.NumTestSets
    trainIdx = training(cv, i);
    testIdx = test(cv, i);

    %Split data
    Xtrain = X(trainIdx, :);
    Ytrain = Y(trainIdx);
    Xtest = X(testIdx, :);
    Ytest = Y(testIdx);

    %Train & Test
    loss = Test_error(Xtrain, Ytrain, Xtest, Ytest, C, IC50s, lambda, isSSRE);
    loss_set(i) = loss;
end

meanLoss = mean(loss_set);
loss = meanLoss;
end