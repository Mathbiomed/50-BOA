clear
clc

% Rate constant in inverse direction
ki1 = 100; ki3 = 100; ki4 = 100; ki5 = 100;
k = [ki1 ki3 ki4 ki5];
%% Hyperparameters
%Time points
tpts = 5;

% Number of replicates
replicate = 5;
replicate_IC50 = 5;

% Error structure model
model = 7;
%% Parameters

% Observed parameter
Vmax = 0.1; Km = 1;

% Substrate concentration for estimating IC50
St_IC50 = Km;

% True values for inhibition constants
Kic = 1; Kiu = 1;

%Enzyme, substrate concentration
Et = 0.01; St_set = [0.2 1 5]*Km;
%% Compute IC50
% Equation for estimating IC50
f = @(alpha, IC50, d, I) d - d./(1+(I/IC50).^alpha);

% Object function
loss = @(X, Y, E) sum((Y - f(E(2),E(1),E(3),X)).^2);

% Initial velocity without IT
IC50_cand = logspace(-3, 3, 10);
IC50_cand = [0 IC50_cand];

% Get the initial velocity time
T = logspace(-3,8,5000);
St_set_IC50 = Km;
T_row = zeros(numel(IC50_cand),1); T_set = zeros(numel(St_set_IC50),numel(IC50_cand));
T_row_SS = zeros(numel(IC50_cand),1); T_set_SS = zeros(numel(St_set_IC50),numel(IC50_cand));
n = 1; m = 1;
for x = St_set_IC50
    for y = IC50_cand
        Tend = 0; Tend_SS = 0; diff = 0;
        St = x; It = y;
        [~,full_SS] = ode15s(@(t,p)fullDynamics(t,p,[Vmax;Km;Kic;Kiu], k, Et),...
            T,[St;Et;0;0;0;It;0]);
        P = full_SS(:,7);
        initial = false;
        for z = 1:numel(P)
            p = P(z);
            % Initial velocity
            if initial == false
                if p/St > 0.01
                    Tend = T(z);
                    break
                end
            end
        end
        T_row(n) = Tend;
        T_row_SS(n) = Tend_SS;
        n = n + 1;
    end
    T_set(m,:) = T_row;
    T_set_SS(m,:) = T_row_SS;
    n = 1;
    m = m + 1;
end

% Compute initial velocity
IC50_vel_raw = zeros(numel(St_set_IC50)*numel(IC50_cand), replicate_IC50);
n = 1;
for r = 1:replicate_IC50
    for t1 = 1:numel(St_set_IC50)
        St = St_set_IC50(t1);
        for t2 = 1:numel(IC50_cand)
            It = IC50_cand(t2);
            
            % Select the appropriate T
            T_end = T_set(t1, t2);
            T = linspace(0, T_end, tpts);

            % In silico data with noise
            [~, full] = ode15s(@(t,p)fullDynamics(t,p,[Vmax;Km;Kic;Kiu], k, Et),...
                T,[St;Et;0;0;0;It;0]);
            Pdata_no_noise = full(end,end)/T_end;
            Pdata = introduceNoise(Pdata_no_noise, [St It], model);
            if Pdata < 0
                Pdata = 10^-3;
            end
            IC50_vel_raw(n, r) = Pdata;
            n = n + 1;
        end
    end
    n = 1;
end
IC50_vel = zeros(height(IC50_vel_raw), 1);
for i = 1:numel(IC50_vel)
    IC50_vel(i) = mean(IC50_vel_raw(i, :));
end
p_Inh = (1 - IC50_vel/IC50_vel(1))*100;

% Estimation and plot the prediction
E0 = [10 1 100];
objFun = @(E)loss(IC50_cand', p_Inh, E);

options = optimset('Display', 'off');
IC50s = fminsearch(objFun, E0, options);

% Absolute IC50
ABS_IC50 = IC50s(1)*(50/(IC50s(3)-50))^(1/IC50s(2));

% Inhibitor concentration
It_set = ABS_IC50*[0.1 1 10];

%% Get the initial velocity time
T = logspace(-3,8,5000);
T_row = zeros(numel(It_set),1); T_set = zeros(numel(St_set),numel(It_set));
T_row_SS = zeros(numel(It_set),1); T_set_SS = zeros(numel(St_set),numel(It_set));
n = 1; m = 1;
for x = St_set
    for y = It_set
        Tend = 0; Tend_SS = 0; diff = 0;
        St = x; It = y;
        [~,full_SS] = ode15s(@(t,p)fullDynamics(t,p,[Vmax;Km;Kic;Kiu], k, Et),...
            T,[St;Et;0;0;0;It;0]);
        P = full_SS(:,7);
        initial = false;
        for z = 1:numel(P)
            p = P(z);
            % Initial velocity: until product concentration reaches to 1 %
            % of St
            if initial == false
                if p/St > 0.01
                    Tend = T(z);
                    break
                end
            end
        end
        T_row(n) = Tend;
        T_row_SS(n) = Tend_SS;
        n = n + 1;
    end
    T_set(m,:) = T_row;
    T_set_SS(m,:) = T_row_SS;
    n = 1;
    m = m + 1;
end
%% Set the range of parameters and concentrations
% Select the indices of concentration to use for data generation
Et_type_set = 1; St_type_set = linspace(1, numel(St_set), numel(St_set)); It_type_set = [1 2 3];

inputSt = zeros(numel(St_set)*numel(It_type_set)*replicate,1);
inputIt = zeros(numel(St_set)*numel(It_type_set)*replicate,1);
whole_P = zeros(numel(St_set)*numel(It_type_set),1);
whole_P_no_noise = zeros(numel(St_set)*numel(It_type_set),1);
raw_data = zeros(numel(St_set)*numel(It_type_set)*replicate, 1);
n = 1;

for r = 1:replicate
    for t1 = 1:numel(St_type_set)
        St_type = St_type_set(t1);
        for t2 = 1:numel(It_type_set)
            It_type = It_type_set(t2);
            St = St_set(St_type); It = It_set(It_type);
            % Select the appropriate T
            T_end = T_set(St_type, It_type);
            T = linspace(0, T_end, tpts);

            % In silico data with noise
            [~, full] = ode15s(@(t,p)fullDynamics(t,p,[Vmax;Km;Kic;Kiu], k, Et),T,[St;Et;0;0;0;It;0]);
            Pdata_no_noise = full(end,end)/T_end;
            Pdata = introduceNoise(Pdata_no_noise, [St It], model);
            % If the product concentration < 0 : assign zero value.
            if Pdata < 0
                Pdata = 0;
                disp("Alert")
            end
            inputSt(n + numel(St_set)*numel(It_type_set)*(r-1)) = St;
            inputIt(n + numel(St_set)*numel(It_type_set)*(r-1)) = It;
            whole_P(n) = Pdata;
            whole_P_no_noise(n) = Pdata_no_noise;
            n = n + 1;
        end
    end
    raw_data(1+numel(St_set)*numel(It_type_set)*(r-1):numel(St_set)*numel(It_type_set)*r,1) = whole_P;
    n = 1;
end

% Raw data: 1 column - Substrate, 2 column - Inhibitor, 3 column - V0
data = [inputSt inputIt raw_data];
%% Full model
function dpdt = fullDynamics(t, p, b, k, Et)
Vmax = b(1); Km = b(2); Kic = b(3); Kiu = b(4);
ki1 = k(1); ki3 = k(2); ki4 = k(3); ki5 = k(4);

k2 = Vmax/Et; A = Kiu/Kic;
k1 = (ki1+k2)/Km; k3 = ki3/Kic; k4 = ki4/Kiu; k5 = ki5/(A*Km);

S = p(1); E = p(2); C = p(3); B = p(4); Y = p(5); I = p(6);

odeS = -k1*E*S-k5*Y*S+ki1*C+ki5*B;
odeE = -k1*E*S-k3*E*I+(ki1+k2)*C+ki3*Y;
odeC = -(ki1+k2)*C-k4*C*I+k1*E*S+ki4*B;
odeB = -ki4*B-ki5*B+k4*C*I+k5*Y*S;
odeY = -ki3*Y-k5*Y*S+k3*E*I+ki5*B;
odeI = ki3*Y+ki4*B-k3*E*I-k4*C*I;
odeP = k2*C;

dpdt = [odeS;odeE;odeC;odeB;odeY;odeI;odeP];
end
%% Noise function
function Pdata = introduceNoise(P, X, model)
if model == 1
    K1 = 0.011; K2 = 0.04;
    s = K1 + K2*P;
elseif model == 2
    K1 = 0.036; K2 = 0.008;
    s = K1 + K2*P.^2;
elseif model == 3
    K1 = 0.05; K2 = 0.889;
    s = K1*P.^K2;
elseif model == 4
    K1 = 0.15; K2 = 0.529; K3 = 0.699; K4 = 0.561; K5 = 0.505;
    s = K1*X(:,1).^K2./(X(:,1).^K3 + K4*X(:,2).^K5);
elseif model == 5
    K1 = 0.15; K2 = 0.063; K3 = 0.273; K4 = 0.095;
    s = K1*X(:,1)./(K2 + X(:,1) + K3*X(:,1).^2 + K4*X(:,2));
elseif model == 6
    K1 = 0.05; K2 = 0.063; K3 = 0.273; K4 = 0.095; K5 = 0.007;
    s = (K1*X(:,1) + K5*X(:,2))./(K2 + X(:,1) + K3*X(:,1).^2 + K4*X(:,2));
else
    s = 0.1*P;
end
noise = zeros(numel(P), 1);
for i = 1:numel(P)
    noise(i) = normrnd(0, s(i), [1 1]);
end
Pdata = P + noise;
end