clear
clc

%% Model
% Model with partial inhibition
v = @(C,X,K) C(1)*(C(2)*X(:,1).*(1+X(:,2)/K(2))+(1+X(:,2)/(K(2)^2/K(1))).*X(:,1).^2/C(3))...
    ./(C(2)^2*(1+X(:,2)/K(1))+2*C(2)*X(:,1).*(1+X(:,2)/K(2))+(1+X(:,2)/(K(2)^2/K(1))).*X(:,1).^2/C(3));

% Equation for estimating IC50
f = @(alpha, IC50, d, I) d - d./(1+(I/IC50).^alpha);

% IC50 equation
h1 = @(S, a, x) a*S*x + S^2;
h2 = @(S, a, x) a*x^2 + 2*a*S*x + S^2;
H = @(Kic, Kiu, S, Ks, a) Kic*(Kiu/Kic)^2*(h1(S,a,Ks)*h2(S,a,Ks))...
    /(h1(S,a,Ks)*h2(S,a,Ks*Kiu/Kic)-2*h1(S,a,Ks*Kiu/Kic)*h2(S,a,Ks));

% Object function
loss = @(X, Y, E) sum((Y - f(E(2),E(1),E(3),X)).^2);
%% Parameters and variables
% Observed parameters
Vmax = 0.1; Ks = 1; a = 0.04;
C = [Vmax Ks a];

% True values for inhibition constants
Kic = 1; Kiu = 10;
K = [Kic Kiu];
%% Generate %Inhibition data
% Substrate concentration for estimating IC50
St_IC50 = Ks; 

It_ranged = logspace(-3,3,10); It_ranged = [0 It_ranged];
control = [St_IC50*ones(height(It_ranged'),1) zeros(height(It_ranged'),1)];
experimental = [St_IC50*ones(height(It_ranged'),1) It_ranged'];
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

St_set = [0.1 0.2 1 2.5 5]*Ks;
It_set = [0.01 0.1 1 10]*ABS_IC50;
%% Make in silico data with 5 replicates
replicate = 5;
raw_data = zeros(numel(St_set)*numel(It_set)*replicate,1);
raw_X = zeros(numel(St_set)*numel(It_set)*replicate,2);

n = 1;
for i = 1:replicate
    for j = 1:numel(St_set)
        for k = 1:numel(It_set)
            e = normrnd(0, 1, [1 1]);
            X = [St_set(j) It_set(k)];
            raw_data(n, 1) = v(C, X, K)*(1 + 0.1*e);
            raw_X(n, :) = X;
            n = n + 1;
        end
    end
end

% Raw data: 1 column - Substrate, 2 column - Inhibitor, 3 column - V0
data = [raw_X raw_data];