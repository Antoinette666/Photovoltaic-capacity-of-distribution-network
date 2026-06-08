function result = NoQControlcapacity(cap_kW)
%NOQCONTROLCAPACITY Exact DistFlow simulation without reactive-power control.

format short

NT = 24;
SB = 10;                       % MVA
VBase = 12.66;                 % kV
ZB = VBase^2 / SB;
Sbase_kVA = SB * 1000;

% Use the same exact power-flow implementation as the droop model.
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(thisDir), 'droop'));

[PDN_Data, Pd_total, Qd_total, loc2, solar_data] = get_PDN_Data();

NL = size(PDN_Data, 1);
NB = NL + 1;
L_i = PDN_Data(:, 2);
L_j = PDN_Data(:, 3);
P_load = PDN_Data(:, 4);
Q_load = PDN_Data(:, 5);
R = PDN_Data(:, 6) / ZB;
X = PDN_Data(:, 7) / ZB;

% Column 8 is the branch apparent-power rating in MVA.
Smax = PDN_Data(:, 8) / SB;
if any(Smax <= 0)
    error('Branch apparent-power ratings must be positive.');
end

Bc = zeros(NL, 1);
U_min = 0.90;
U_max = 1.10;
I2_max = (Smax / U_min).^2;

Pd = P_load * Pd_total / 1e6;
Qd = Q_load * Qd_total / 1e6;
Pd = [zeros(1, NT); Pd];
Qd = [zeros(1, NT); Qd];

Psolar = zeros(NB, NT);
Qsolar = zeros(NB, NT);
nPV = length(loc2);

if nPV == 0 || cap_kW <= 1e-6
    pv_cap_each = zeros(nPV, 1);
else
    pv_cap_each = cap_kW / nPV * ones(nPV, 1);
end

for k = 1:nPV
    profile = solar_data(k, :);
    if max(profile) > 1e-6
        Psolar(loc2(k), :) = pv_cap_each(k) * profile / max(profile);
    end
end

U_all = zeros(NB, NT);
U2_all = zeros(NB, NT);
I2_all = zeros(NL, NT);
P_all = zeros(NL, NT);
Q_all = zeros(NL, NT);
pf_converged_record = false(NT, 1);

for t = 1:NT
    pf = radial_DistFlow_powerflow( ...
        Pd(:, t), Qd(:, t), ...
        Psolar(:, t) / Sbase_kVA, Qsolar(:, t) / Sbase_kVA, ...
        L_i, L_j, R, X, Bc, NB, NL);

    U2_all(:, t) = pf.U2;
    U_all(:, t) = sqrt(max(pf.U2, 0));
    I2_all(:, t) = pf.I2;
    P_all(:, t) = pf.P;
    Q_all(:, t) = pf.Q;
    pf_converged_record(t) = pf.converged;
end

line_loading_ratio = hypot(P_all, Q_all) ./ repmat(Smax, 1, NT);
I2_ratio = I2_all ./ repmat(I2_max, 1, NT);

result = struct();
result.cap_total_kW = cap_kW;
result.is_feasible = false;
result.violation_type = "";
result.sol_info = "Exact radial DistFlow simulation without Q control.";
result.converged = all(pf_converged_record);
result.pf_converged_record = pf_converged_record;

result.max_U = max(U_all(:));
result.min_U = min(U_all(:));
result.max_U_all = result.max_U;
result.min_U_all = result.min_U;
result.max_I2 = max(I2_all(:));
result.max_I2_ratio = max(I2_ratio(:));
result.max_S_MVA = max(hypot(P_all, Q_all), [], 'all') * SB;
result.max_line_loading_ratio = max(line_loading_ratio(:));

result.Obj_value = NaN;
result.soc_gap_max = 0;
result.U_opt_pu = U_all;
result.U2_opt = U2_all;
result.I2_opt = I2_all;
result.P_opt_pu = P_all;
result.Q_opt_pu = Q_all;
result.Psolar_kW = Psolar;
result.Qsolar_kvar = Qsolar;

tol = 1e-5;
if ~result.converged
    result.violation_type = "潮流计算未收敛";
elseif result.max_U > U_max + tol
    result.violation_type = "节点电压越上限";
elseif result.min_U < U_min - tol
    result.violation_type = "节点电压越下限";
elseif result.max_line_loading_ratio > 1 + tol
    result.violation_type = "线路视在功率越限";
elseif result.max_I2_ratio > 1 + tol
    result.violation_type = "线路电流越限";
else
    result.is_feasible = true;
    result.violation_type = "无越限";
end
end
