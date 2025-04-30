clear
clc

%% Random bireactant system with inhibitor competing with one substrate
% Model with multiple inhibition
v1 = @(C,X,K) C(3)*(1+X(:,3)/K(1)) + X(:,2).*(1+X(:,3)/K(2));
v2 = @(C,X,K) X(:,2) + C(4)*C(3);
v = @(C,X,K) C(1)*X(:,1).*X(:,2)./(C(4)*C(2)*v1(C,X,K) + X(:,1).*v2(C,X,K));

% Equation for estimating IC50
f = @(alpha, IC50, d, I) d - d./(1+(I/IC50).^alpha);

% IC50 equation
r1 = @(K, X, C) C(4)*C(2)*C(3)/(C(4)*C(2)*(C(3)+X(2)) + X(1)*(X(2)+C(4)*C(3)));
r2 = @(K, X, C) C(4)*C(2)*X(2)/(C(4)*C(2)*(C(3)+X(2)) + X(1)*(X(2)+C(4)*C(3)));

H = @(C,X,K) K(1)*K(2)/(r1(K,X,C)*K(2)+r2(K,X,C)*K(1));

% Object function
loss = @(X, Y, E) sum((Y - f(E(2),E(1),E(3),X)).^2);
%% Parameters and variables
% Observed parameters
Vmax = 0.1; Ka = 1; Kb = 1; a = 0.5;
C = [Vmax Ka Kb a];

% True values for inhibition constants
Kic = 1; Kiu = 10;
K = [Kic Kiu];

%% Generate %Inhibition data
% Substrate concentrations for estimating IC50
At_IC50 = Ka; Bt_IC50 = Kb;

It_ranged = logspace(-3,3,10); It_ranged = [0 It_ranged];
control = [At_IC50*ones(height(It_ranged'),1) Bt_IC50*ones(height(It_ranged'),1) zeros(height(It_ranged'),1)];
experimental = [At_IC50*ones(height(It_ranged'),1) Bt_IC50*ones(height(It_ranged'),1) It_ranged'];
noise = normrnd(0, 1, [height(It_ranged) 1]);
p_inhibition = (1-v(C,experimental,K)./(v(C,control,K)))*100;
p_inhibition = p_inhibition + noise;
for i = 1:height(p_inhibition)
    if p_inhibition(i) < 0
        p_inhibition(i) = 0;
    end
end

% Estimate IC50
E0 = [1 8 50];
objFun = @(E)loss(It_ranged', p_inhibition, E);

options = optimset('Display', 'off');
IC50s = fminsearch(objFun, E0, options);
ABS_IC50 = IC50s(1)*(50/(IC50s(3)-50))^(1/IC50s(2));

IC50 = ABS_IC50;

% Experimental setup
At_set = [0.2 1 5]*Ka; Bt_set = [0.2 1 5]*Kb; It_set = [0.1 1 10]*ABS_IC50;
%% Make in silico data with 5 replicates
% Replicate
r = 5;
raw_data = zeros(numel(At_set)*numel(Bt_set)*numel(It_set)*r,1);
raw_X = zeros(numel(At_set)*numel(Bt_set)*numel(It_set)*r,3);

n = 1;
for i = 1:r
    e = normrnd(0, 1, [1 1]);
    for j = 1:numel(At_set)
        for k = 1:numel(Bt_set)
            for l = 1:numel(It_set)
                X = [At_set(j) Bt_set(k) It_set(l)];
                raw_data(n, 1) = v(C, X, K)*(1 + 0.1*e);
                raw_X(n, :) = X;
                n = n + 1;
            end
        end
    end
end

% Raw data: 1 column - Substrate, 2 column - Inhibitor, 3 column - V0
data = [raw_X raw_data];