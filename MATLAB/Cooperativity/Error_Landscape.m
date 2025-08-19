function Error_Landscape(file_name, L, isReg, isSSRE, isMatched)
%% Load data
data = readmatrix(file_name);

%% 1. BOA_Condition
BOA_Condition(file_name);

%% 2. CV_Inhibition
%Range of regularization constant
if isReg
    lambda = CV_Inhibition(file_name, L, isSSRE);
else
    lambda = 0;
end
fprintf('The regularization constant is %.2f.\n', lambda);

%% 3. Estimate Kic and Kiu
%Load the formatted data
Vmax = data(1,1); Ks = data(1,2); IC50 = data(1,3); St_IC50 = data(1,4);
St_setup = data(2:end,1); It_setup = data(2:end,2); V0 = data(2:end,3); a = data(2,4);

X_setup = [St_setup It_setup]; C = [Vmax Ks a]; IC50s = [St_IC50 IC50];

%Estimation
K = Fit_inhibition(X_setup, V0, C, IC50s, lambda, isSSRE);

%Compute 95% confidence interval via bootstrapping
bootSample = 1000;
betaSample = zeros(bootSample, 2);

for i = 1:bootSample
    %Resampling
    indices = randsample(1:length(V0), length(V0), true);
    X_boot = X_setup(indices, :);
    V0_boot = V0(indices);

    %Re-estimate parameter using bootstrap sample
    betaSample(i,:) = Fit_inhibition(X_boot, V0_boot, C, IC50s, lambda, isSSRE);
end

%Compute confidence interval
CI_lower = prctile(betaSample, 2.5);
CI_upper = prctile(betaSample, 97.5);

%Estimates of Kic and Kiu
fprintf('Kic: %.4f, (%.4f, %.4f)\n', K(1), CI_lower(1), CI_upper(1));
fprintf('Kiu: %.4f, (%.4f, %.4f)\n', K(2), CI_lower(2), CI_upper(2));

%% Generate error landscape
%Range of parameter
K_round = round(log10(IC50));
Kic_min = K_round - 2; Kic_max = K_round + 2;
Kiu_min = K_round - 2; Kiu_max = K_round + 2;

%Set up the ranges
Kicrange = logspace(Kic_min,Kic_max,100);
Kiurange = logspace(Kiu_min,Kiu_max,100);

S1 = zeros(numel(Kicrange)*numel(Kiurange),1);
S2 = zeros(numel(Kicrange)*numel(Kiurange),1);
total_error = zeros(numel(Kicrange)*numel(Kiurange),1);

%Calculate total error
idx = 1;
for i = 1:numel(Kicrange)
    for j = 1:numel(Kiurange)
        Kicr = Kicrange(i);
        Kiur = Kiurange(j);
        K = [Kicr Kiur];
        total_error(idx) = CV_loss(K, X_setup, V0, C, IC50s, lambda, isSSRE);
        S1(idx) = Kicr; S2(idx) = Kiur;
        idx = idx + 1;
    end
end

error_table = table;
error_table.Kic = S1; error_table.Kiu = S2; error_table.Error = total_error;
min_error = min(total_error);

%Generate heatmap
figure

%Axis label
CustomXLabels = string(Kicrange); CustomYLabels = string(Kiurange);
CustomXLabels(mod(Kicrange,10^Kic_min) ~= 0) = " ";
CustomYLabels(mod(Kiurange,10^Kiu_min) ~= 0) = " ";

%X-axis label
for i = linspace(Kic_min, Kic_max, Kic_max - Kic_min + 1)
    diff = abs(Kicrange - 10^i*ones(1,numel(Kicrange)));
    min_diff = min(diff); idx_X = diff==min_diff;
    CustomXLabels(idx_X) = string(i);
end

%Y-axis label
for i = linspace(Kiu_min, Kiu_max, Kiu_max - Kiu_min + 1)
    diff = abs(Kiurange - 10^i*ones(1,numel(Kiurange)));
    min_diff = min(diff); idx_Y = diff==min_diff;
    CustomYLabels(idx_Y) = string(i);
end

%Draw heatmap
h= heatmap(error_table,'Kic','Kiu','ColorVariable','Error');
grid off
if isMatched
    h.ColorLimits = [min_error 2*min_error];
else
    h.ColorLimits = [0 0.05];
end
h.Colormap = hot;
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomYLabels;
h.XLabel = 'log[K_{ic}]'; h.YLabel = 'log[K_{iu}]';
h.FontSize = 20;
h.Title = '';
h.NodeChildren(3).YDir = 'normal';
h.NodeChildren(3).DataAspectRatio = [1 1 1];
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
    % User can assign his or her own model of standard deviation to "s" here.

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
