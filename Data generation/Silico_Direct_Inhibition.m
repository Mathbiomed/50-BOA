clear
clc

%% Model
% Models with direct inhibition
S = @(X, K) 0.5*(sqrt((X(:,2)-X(:,1)+K).^2+4*K*X(:,1))-(X(:,2)-X(:,1)+K));
v = @(C,X,K) C(1)*S(X,K)./(C(2)+S(X,K));

% Equation for estimating IC50
f = @(alpha, IC50, d, I) d - d./(1+(I/IC50).^alpha);

% Object function
loss = @(X, Y, E) sum((Y - f(E(2),E(1),E(3),X)).^2);
%% Parameters and variables
% Observed parameters
Vmax = 0.1; Ks = 1;
C = [Vmax Ks];

% True value of inhibition parameter
K = 1;
%% Generate %Inhibition data
% Substrate concentration for estimating IC50
St_IC50 = Ks;
It_ranged = logspace(-3,3,10); It_ranged = [0 It_ranged];
control = [St_IC50*ones(height(It_ranged'),1) zeros(height(It_ranged'),1)];
experimental = [St_IC50*ones(height(It_ranged'),1) It_ranged'];
noise = normrnd(0, 1, [height(It_ranged) 1]);
p_inhibition = (1-v(C,experimental,K)./(v(C,control,K)))*100;

% Estimate IC50
E0 = [1 8 50];
objFun = @(E)loss(It_ranged', p_inhibition, E);

options = optimset('Display', 'off');
IC50s = fminsearch(objFun, E0, options);
ABS_IC50 = IC50s(1)*(50/(IC50s(3)-50))^(1/IC50s(2));

IC50 = ABS_IC50;

St_set = [0.2 1 5]*Ks;
It_set = [0.1 1 10]*ABS_IC50;
%% Make in silico data with 5 replicates
raw_data = zeros(numel(St_set)*numel(It_set), 5);
raw_X = zeros(numel(St_set)*numel(It_set),2);
for i = 1:5
    n = 1;
    for j = 1:numel(St_set)
        for k = 1:numel(It_set)
            e = normrnd(0, 1, [1 1]);
            X = [St_set(j) It_set(k)];
            raw_data(n, i) = v(C, X, K)*(1 + 0.1*e);
            raw_X(n, :) = X;
            n = n + 1;
        end
    end
end

% Raw data: 1 column - Substrate, 2 column - Inhibitor, 3 column - V0
data = [raw_X raw_data];