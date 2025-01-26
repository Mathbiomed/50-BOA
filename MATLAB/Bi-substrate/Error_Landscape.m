function Error_Landscape(file_name)
%% Load data
data = readmatrix(file_name);

% Change the heatmap color range
isMatched = true;

%% 1. BOA_Condition
BOA_Condition(file_name);

%% 2. CV_Inhibition
%Range of regularization constant
L = logspace(-3,3,100);
lambda = CV_Inhibition(file_name, L);
fprintf('The regularization constant is %.2f.\n', lambda);

%% 3. Estimate Kic and Kiu
%Load the formatted data
Vmax = data(1,1); Ka = data(1,2); Kb = data(1,3); a = data(1,4); 
IC50 = data(1,5); At_IC50 = data(1,6); Bt_IC50 = data(1,7);
At_setup = data(2:end,1); Bt_setup = data(2:end,2); It_setup = data(2:end,3);
V0 = data(2:end,4);

X_setup = [At_setup Bt_setup It_setup]; C = [Vmax Ka Kb a]; IC50s = [At_IC50 Bt_IC50 IC50];

%Estimation
K = Fit_inhibition(X_setup, V0, C, IC50s, lambda);

%Compute 95% confidence interval via bootstrapping
bootSample = 1000;
betaSample = zeros(bootSample, 2);

for i = 1:bootSample
    %Resampling
    indices = randsample(1:length(V0), length(V0), true);
    X_boot = X_setup(indices, :);
    V0_boot = V0(indices);

    %Re-estimate parameter using bootstrap sample
    betaSample(i,:) = Fit_inhibition(X_boot, V0_boot, C, IC50s, lambda);
end

%Compute confidence interval
CI_lower = prctile(betaSample, 2.5);
CI_upper = prctile(betaSample, 97.5);
CI = [CI_lower' CI_upper'];

%Estimates of Kic and Kiu
Estimate = [K' CI];
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
        total_error(idx) = CV_loss(K, X_setup, V0, C, IC50s, lambda);
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
v1 = @(C,X,K) C(3)*(1+X(:,3)/K(1)) + X(:,2).*(1+X(:,3)/K(2));
v2 = @(C,X,K) X(:,2) + C(4)*C(3);
v = C(1)*X(:,1).*X(:,2)./(C(4)*C(2)*v1(C,X,K) + X(:,1).*v2(C,X,K));
end

%% Cheng-Prusoff equation
function v = Cheng_Prusoff(K, X, C)
r1 = @(K, X, C) C(4)*C(2)*C(3)/(C(4)*C(2)*(C(3)+X(2)) + X(1)*(X(2)+C(4)*C(3)));
r2 = @(K, X, C) C(4)*C(2)*X(2)/(C(4)*C(2)*(C(3)+X(2)) + X(1)*(X(2)+C(4)*C(3)));

v = K(1)*K(2)/(r1(K,X,C)*K(2)+r2(K,X,C)*K(1));
end

%% Loss function with lambda
function loss = CV_loss(K, X, Y, C, IC50s, lambda)
Y_predict = Inhibition(K, X, C);
loss = mean(((Y-Y_predict)./Y).^2) +...
                 lambda*mean(((IC50s(:,3)-Cheng_Prusoff(K, [IC50s(:,1) IC50s(:,2)], C))./IC50s(:,3)).^2);
end

%% Fitting
function params = Fit_inhibition(X, Y, C, IC50s, lambda)
K0 = [max(IC50s(:,2)) max(IC50s(:,2))];
objFun = @(K)CV_loss(K, X, Y, C, IC50s, lambda);

options = optimset('Display', 'off');
params = fminsearch(objFun, K0, options);
end