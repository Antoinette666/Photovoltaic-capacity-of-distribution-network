function result = Centralizedcapacity(cap_kW)
%CENTRALIZEDCAPACITY Hybrid centralized Q control with exact DistFlow.
% SOCP and droop control generate initial points. The final dispatch and
% feasibility decision always use exact radial DistFlow power flow.

format short

NT = 24;
SB = 10;
VBase = 12.66;
ZB = VBase^2 / SB;
Sbase_kVA = SB * 1000;

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
Smax = PDN_Data(:, 8) / SB;

Bc = zeros(NL, 1);
U_min = 0.90;
U_max = 1.10;
I2_max = (Smax / U_min).^2;

Pd = P_load * Pd_total / 1e6;
Qd = Q_load * Qd_total / 1e6;
Pd = [zeros(1, NT); Pd];
Qd = [zeros(1, NT); Qd];

nPV = length(loc2);
Psolar = zeros(NB, NT);
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

inv_oversize = 1.1;
Qmax_kvar = zeros(NB, NT);
for k = 1:nPV
    node = loc2(k);
    S_node = inv_oversize * pv_cap_each(k);
    Qmax_kvar(node, :) = sqrt(max(0, S_node^2 - Psolar(node, :).^2));
end

% Two independent initial-point generators.
droop_init = droop_control(cap_kW);
socp_init = Centralized_SOCP_initial(cap_kW);
Qdroop = droop_init.Qsolar_opt_kvar;
Qsocp = zeros(NB, NT);
if ~isempty(socp_init.Qsolar_opt_kvar)
    Qsocp = socp_init.Qsolar_opt_kvar;
end

U_all = zeros(NB, NT);
U2_all = zeros(NB, NT);
I2_all = zeros(NL, NT);
P_all = zeros(NL, NT);
Q_all = zeros(NL, NT);
Qsolar = zeros(NB, NT);
exitflag_record = zeros(NT, 1);
pf_converged_record = false(NT, 1);
objective_record = zeros(NT, 1);

options = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display', 'off', ...
    'MaxIterations', 200, ...
    'MaxFunctionEvaluations', 5000, ...
    'ConstraintTolerance', 1e-7, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-8);

for t = 1:NT
    qmax = Qmax_kvar(loc2, t);
    candidates = [Qdroop(loc2, t), Qsocp(loc2, t), zeros(nPV, 1)];

    best = exact_candidate(candidates(:, 1), t, Psolar, Pd, Qd, loc2, ...
        Sbase_kVA, L_i, L_j, R, X, Bc, NB, NL, Smax, I2_max, U_min, U_max);

    for j = 2:size(candidates, 2)
        trial = exact_candidate(candidates(:, j), t, Psolar, Pd, Qd, loc2, ...
            Sbase_kVA, L_i, L_j, R, X, Bc, NB, NL, Smax, I2_max, U_min, U_max);
        best = choose_candidate(best, trial);
    end

    x0 = min(max(best.q, -qmax), qmax);
    objective = @(q) exact_objective(q, t, Psolar, Pd, Qd, loc2, ...
        Sbase_kVA, L_i, L_j, R, X, Bc, NB, NL);
    nonlcon = @(q) exact_constraints(q, t, Psolar, Pd, Qd, loc2, ...
        Sbase_kVA, L_i, L_j, R, X, Bc, NB, NL, Smax, I2_max, U_min, U_max);

    try
        [qopt, ~, exitflag] = fmincon(objective, x0, [], [], [], [], ...
            -qmax, qmax, nonlcon, options);
        optimized = exact_candidate(qopt, t, Psolar, Pd, Qd, loc2, ...
            Sbase_kVA, L_i, L_j, R, X, Bc, NB, NL, Smax, I2_max, U_min, U_max);
        best = choose_candidate(best, optimized);
    catch
        exitflag = -999;
    end

    Qsolar(loc2, t) = best.q;
    U_all(:, t) = best.U;
    U2_all(:, t) = best.pf.U2;
    I2_all(:, t) = best.pf.I2;
    P_all(:, t) = best.pf.P;
    Q_all(:, t) = best.pf.Q;
    pf_converged_record(t) = best.pf.converged;
    objective_record(t) = best.objective;
    exitflag_record(t) = exitflag;
end

line_loading_ratio = hypot(P_all, Q_all) ./ repmat(Smax, 1, NT);
I2_ratio = I2_all ./ repmat(I2_max, 1, NT);

result = struct();
result.cap_total_kW = cap_kW;
result.is_feasible = false;
result.violation_type = "";
result.sol_info = "SOCP/droop initialization followed by exact DistFlow fmincon.";
result.socp_initial_info = socp_init.sol_info;
result.socp_initial_gap = socp_init.soc_gap_max;
result.exact_optimization_exitflag = exitflag_record;
result.pf_converged_record = pf_converged_record;
result.converged = all(pf_converged_record);

result.U_opt_pu = U_all;
result.U2_opt = U2_all;
result.I2_opt = I2_all;
result.P_opt_pu = P_all;
result.Q_opt_pu = Q_all;
result.Psolar_kW = Psolar;
result.Qsolar_opt_kvar = Qsolar;

result.Obj_value = sum(objective_record);
result.soc_gap_max = 0;
result.max_U = max(U_all(:));
result.min_U = min(U_all(:));
result.max_I2 = max(I2_all(:));
result.max_I2_ratio = max(I2_ratio(:));
result.max_S_MVA = max(hypot(P_all, Q_all), [], 'all') * SB;
result.max_line_loading_ratio = max(line_loading_ratio(:));

tol = 1e-5;
if ~result.converged
    result.violation_type = "精确潮流计算未收敛";
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

function value = exact_objective(q, t, Psolar, Pd, Qd, loc2, Sbase_kVA, ...
    L_i, L_j, R, X, Bc, NB, NL)
pf = run_exact_pf(q, t, Psolar, Pd, Qd, loc2, Sbase_kVA, ...
    L_i, L_j, R, X, Bc, NB, NL);
U = sqrt(max(pf.U2, 0));
loss = sum(R .* pf.I2);
voltage_deviation = sum((U(2:end) - 1).^2);
value = loss + 10 * voltage_deviation + 1e-8 * sum(q.^2);
if ~pf.converged
    value = value + 1e6;
end
end

function [c, ceq] = exact_constraints(q, t, Psolar, Pd, Qd, loc2, ...
    Sbase_kVA, L_i, L_j, R, X, Bc, NB, NL, Smax, I2_max, U_min, U_max)
pf = run_exact_pf(q, t, Psolar, Pd, Qd, loc2, Sbase_kVA, ...
    L_i, L_j, R, X, Bc, NB, NL);
U = sqrt(max(pf.U2, 0));
c = [U_min - U; U - U_max; hypot(pf.P, pf.Q) - Smax; pf.I2 - I2_max];
if ~pf.converged
    c = c + 1;
end
ceq = [];
end

function candidate = exact_candidate(q, t, Psolar, Pd, Qd, loc2, ...
    Sbase_kVA, L_i, L_j, R, X, Bc, NB, NL, Smax, I2_max, U_min, U_max)
pf = run_exact_pf(q, t, Psolar, Pd, Qd, loc2, Sbase_kVA, ...
    L_i, L_j, R, X, Bc, NB, NL);
U = sqrt(max(pf.U2, 0));
violations = [U_min - U; U - U_max; hypot(pf.P, pf.Q) - Smax; pf.I2 - I2_max];
max_violation = max([0; violations]);
candidate.q = q;
candidate.pf = pf;
candidate.U = U;
candidate.objective = sum(R .* pf.I2) + 10 * sum((U(2:end) - 1).^2) ...
    + 1e-8 * sum(q.^2);
candidate.feasible = pf.converged && max_violation <= 1e-5;
candidate.merit = candidate.objective + 1e6 * max_violation;
end

function best = choose_candidate(best, trial)
if trial.feasible && ~best.feasible
    best = trial;
elseif trial.feasible == best.feasible && trial.merit < best.merit
    best = trial;
end
end

function pf = run_exact_pf(q, t, Psolar, Pd, Qd, loc2, Sbase_kVA, ...
    L_i, L_j, R, X, Bc, NB, NL)
Qpv = zeros(NB, 1);
Qpv(loc2) = q / Sbase_kVA;
pf = radial_DistFlow_powerflow(Pd(:, t), Qd(:, t), ...
    Psolar(:, t) / Sbase_kVA, Qpv, L_i, L_j, R, X, Bc, NB, NL);
end
