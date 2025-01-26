function BOA_Condition(file_name)
%% Load the data
data = readmatrix(file_name);
%Data check
if (height(data) < 2 || width(data) < 4)
    error('Invalid Input')
end

%Load the formatted data
Km = data(1,2); IC50 = data(1,3);
St_setup = data(2:end,1); It_setup = data(2:end,2);

%% 1. Check whether It >= IC50
reject_1 = "";
check_1 = false;
check_It = any(It_setup < IC50 * ones(numel(It_setup),1));
if check_It
    reject_1 = "inhibitor concentration < IC50";
    check_1 = true;
end

%% 2. Check whether St varies sufficiently
reject_2 = "";
check_2 = false;

%Check whether St does not vary
matrix_vary = St_setup - St_setup(1)*ones(numel(St_setup),1);
check_vary = isequal(matrix_vary,zeros(height(matrix_vary), width(matrix_vary)));

%Check whether St ranges sufficiently over 0.2Km ~ 5Km
check_sufficient = min(St_setup) > 0.2*Km | max(St_setup) < 5*Km;

if check_vary || check_sufficient
    reject_2 = "Substrate concentration should vary";
    check_2 = true;
end

%% Output
if check_1 || check_2
    disp("Estimation may be insufficient for precise results:");
    disp(reject_1);
    disp(reject_2);
end
end
