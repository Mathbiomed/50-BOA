function BOA_Condition(file_name)
%% Load the data
data = readmatrix(file_name);
%Data check
if (height(data) < 2 || width(data) < 4)
    error('Invalid Input')
end

%Load the formatted data
Ka = data(1,2); Kb = data(1,3); IC50 = data(1,4);
At_setup = data(2:end,1); Bt_setup = data(2:end,2); It_setup = data(2:end,3);

%% 1. Check whether It >= IC50
reject_1 = "";
check_1 = false;
check_It = any(It_setup < IC50 * ones(numel(It_setup),1));
if check_It
    reject_1 = "inhibitor concentration < IC50";
    check_1 = true;
end

%% 2. Check whether At varies sufficiently
reject_2 = "";
check_2 = false;

%Check whether At does not vary
matrix_vary_A = At_setup - At_setup(1)*ones(numel(At_setup),1);
check_vary_A = isequal(matrix_vary_A,zeros(height(matrix_vary_A), width(matrix_vary_A)));

%Check whether At ranges sufficiently over 0.2Ka ~ 5Ka
check_sufficient_A = min(At_setup) > 0.2*Ka | max(At_setup) < 5*Ka;

if check_vary_A || check_sufficient_A
    reject_2 = "First type substrate concentration should vary";
    check_2 = true;
end

%% 3. Check whether Bt varies sufficiently
reject_3 = "";
check_3 = false;

%Check whether St does not vary
matrix_vary_B = Bt_setup - Bt_setup(1)*ones(numel(Bt_setup),1);
check_vary_B = isequal(matrix_vary_B,zeros(height(matrix_vary_B), width(matrix_vary_B)));

%Check whether Bt ranges sufficiently over 0.2Kb ~ 5Kb
check_sufficient_B = min(Bt_setup) > 0.2*Kb | max(Bt_setup) < 5*Kb;

if check_vary_B || check_sufficient_B
    reject_3 = "Second type substrate concentration should vary";
    check_3 = true;
end

%% Output
if check_1 || check_2 || check_3
    disp("Estimation may be insufficient for precise results:");
    disp(reject_1);
    disp(reject_2);
    disp(reject_3);
end
end
