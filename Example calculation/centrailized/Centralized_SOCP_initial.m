function result = Centralized_SOCP_initial(cap_kW)
%CENTRALIZEDCAPACITY 集中式光伏无功优化的 DistFlow-SOCP 模型
% cap_kW 为全网光伏总装机容量，单位 kW。

NT = 24;
SB = 10;                       % MVA
VBase = 12.66;                 % kV
ZB = VBase^2 / SB;             % ohm
Sbase_kVA = SB * 1000;         % kVA

%% 系统数据
[PDN_Data, Pd_total, Qd_total, loc2, solar_data] = get_PDN_Data();

NL = size(PDN_Data, 1);
NB = NL + 1;
L_i = PDN_Data(:, 2);
L_j = PDN_Data(:, 3);
P_load = PDN_Data(:, 4);
Q_load = PDN_Data(:, 5);
R = PDN_Data(:, 6) / ZB;
X = PDN_Data(:, 7) / ZB;

% 数据第 8 列作为线路视在功率额定值，单位 MVA。
Smax = PDN_Data(:, 8) / SB;
if any(Smax <= 0)
    error('线路视在功率额定值必须为正数。');
end

% 数据中没有线路对地电纳，不能使用未定义或人为假设的 Bc。
Bc = zeros(NL, 1);

U_min = 0.90;
U_max = 1.10;
U2_min = U_min^2;
U2_max = U_max^2;

% 由线路容量和最低允许电压得到保守的电流平方上限。
I2_max = (Smax / U_min).^2;

%% 径向网络拓扑
children = cell(NL, 1);
for ell = 1:NL
    children{ell} = find(L_i == L_j(ell));
end

%% 负荷标幺化
Pd = P_load * Pd_total / 1e6;
Qd = Q_load * Qd_total / 1e6;
Pd = [zeros(1, NT); Pd];
Qd = [zeros(1, NT); Qd];

%% 光伏有功及逆变器无功能力
nPV = length(loc2);
Psolar = zeros(NB, NT);        % kW

if nPV == 0 || cap_kW <= 1e-6
    pv_cap_each = zeros(nPV, 1);
else
    pv_cap_each = cap_kW / nPV * ones(nPV, 1);
end

pv_shape = zeros(nPV, NT);
for k = 1:nPV
    profile = solar_data(k, :);
    if max(profile) > 1e-6
        pv_shape(k, :) = profile / max(profile);
    end
    Psolar(loc2(k), :) = pv_cap_each(k) * pv_shape(k, :);
end

% 逆变器适度超配，使光伏满有功输出时仍保留无功调节能力。
% 若取 1.0，峰值时段 Qmax 恰好为 0，会形成大量退化的零宽约束。
inv_oversize = 1.1;
Qmax_kvar = zeros(NB, NT);
for k = 1:nPV
    node = loc2(k);
    S_node = inv_oversize * pv_cap_each(k);
    Qmax_kvar(node, :) = sqrt(max(0, S_node^2 - Psolar(node, :).^2));
end
Qmax_pu = Qmax_kvar / Sbase_kVA;

%% 决策变量
P = sdpvar(NL, NT, 'full');       % 线路发送端有功，p.u.
Q = sdpvar(NL, NT, 'full');       % 线路发送端无功，p.u.
U2 = sdpvar(NB, NT, 'full');      % 节点电压平方，p.u.
I2 = sdpvar(NL, NT, 'full');      % 线路电流平方，p.u.
Qsolar = sdpvar(NB, NT, 'full');  % 光伏无功，p.u.
Vdev = sdpvar(NB - 1, NT, 'full');

Cons_P = [];
Cons_Q = [];
Cons_U = [];
Cons_I = [];
Cons_SOC = [];
Cons_Line = [];
Cons_PV = [];
Cons_Vdev = [0 <= Vdev, ...
    Vdev >= U2(2:NB, :) - 1, ...
    Vdev >= 1 - U2(2:NB, :)];

Ind_NoPV = setdiff(1:NB, loc2);
Cons_PV = [Cons_PV, Qsolar(Ind_NoPV, :) == 0];
if ~isempty(loc2)
    Cons_PV = [Cons_PV, ...
        -Qmax_pu(loc2, :) <= Qsolar(loc2, :), ...
         Qsolar(loc2, :) <= Qmax_pu(loc2, :)];
end

%% DistFlow 功率平衡
for ell = 1:NL
    child = children{ell};
    if isempty(child)
        P_down = zeros(1, NT);
        Q_down = zeros(1, NT);
    else
        P_down = sum(P(child, :), 1);
        Q_down = sum(Q(child, :), 1);
    end

    from = L_i(ell);
    to = L_j(ell);
    Qc_line = Bc(ell) * (U2(from, :) + U2(to, :));

    Cons_P = [Cons_P, ...
        P(ell, :) - R(ell) * I2(ell, :) + Psolar(to, :) / Sbase_kVA ...
        == P_down + Pd(to, :)];
    Cons_Q = [Cons_Q, ...
        Q(ell, :) - X(ell) * I2(ell, :) + Qsolar(to, :) + Qc_line ...
        == Q_down + Qd(to, :)];
end

%% 电压降约束
Cons_U = [Cons_U, ...
    U2(L_j, :) == U2(L_i, :) ...
    - 2 * (repmat(R, 1, NT) .* P + repmat(X, 1, NT) .* Q) ...
    + repmat(R.^2 + X.^2, 1, NT) .* I2];
Cons_U = [Cons_U, U2(1, :) == 1, U2_min <= U2, U2 <= U2_max];

%% 线路电流、潮流容量约束及二阶锥松弛
Cons_I = [Cons_I, 0 <= I2, I2 <= repmat(I2_max, 1, NT)];

for t = 1:NT
    for ell = 1:NL
        % P^2 + Q^2 <= U_from^2 * I^2 的旋转二阶锥形式。
        Cons_SOC = [Cons_SOC, ...
            cone([2 * P(ell, t); 2 * Q(ell, t); ...
                  I2(ell, t) - U2(L_i(ell), t)], ...
                 I2(ell, t) + U2(L_i(ell), t))];

    end
end

% 线路发送端视在功率约束。使用凸二次形式，避免旧版 YALMIP 在同一
% 支路上叠加多个锥约束时引发 Gurobi 数值错误。
% 当前线路额定容量为 9.9 MVA，而算例最大潮流远低于该值。
% 暂不将冗余热稳定约束放入优化模型，求解后仍计算并检查负载率。

Cons = [Cons_P, Cons_Q, Cons_U, Cons_I, Cons_SOC, Cons_Line, Cons_PV, Cons_Vdev];

%% 目标函数
% 网损项用于收紧 SOCP；线性电压项抑制高光伏时段的过电压。
% 二次电压偏差目标会使当前旧版 YALMIP/Gurobi 组合出现数值问题。
voltage_weight = 0.05;
Obj = sum(sum(I2)) + voltage_weight * sum(sum(Vdev));

%% 求解
ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
sol = optimize(Cons, Obj, ops);

%% 输出初始化
result = struct();
result.cap_total_kW = cap_kW;
result.is_feasible = false;
result.violation_type = "";
result.sol_info = sol.info;
result.max_U = NaN;
result.min_U = NaN;
result.max_I2 = NaN;
result.max_I2_ratio = NaN;
result.max_S_MVA = NaN;
result.max_line_loading_ratio = NaN;
result.Obj_value = NaN;
result.soc_gap_max = NaN;
result.U_opt_pu = [];
result.U2_opt = [];
result.I2_opt = [];
result.P_opt_pu = [];
result.Q_opt_pu = [];
result.Psolar_kW = Psolar;
result.Qsolar_opt_kvar = [];

if sol.problem ~= 0
    result.violation_type = "优化不可行或求解失败：" + string(sol.info);
    return;
end

%% 提取并检查结果
U2_opt = value(U2);
I2_opt = value(I2);
P_opt_pu = value(P);
Q_opt_pu = value(Q);
Qsolar_opt = value(Qsolar);
U_opt_pu = sqrt(max(U2_opt, 0));
S_opt_pu = hypot(P_opt_pu, Q_opt_pu);
soc_gap = I2_opt .* U2_opt(L_i, :) - P_opt_pu.^2 - Q_opt_pu.^2;
line_loading_ratio = S_opt_pu ./ repmat(Smax, 1, NT);
I2_ratio = I2_opt ./ repmat(I2_max, 1, NT);

result.U_opt_pu = U_opt_pu;
result.U2_opt = U2_opt;
result.I2_opt = I2_opt;
result.P_opt_pu = P_opt_pu;
result.Q_opt_pu = Q_opt_pu;
result.Qsolar_opt_kvar = Qsolar_opt * Sbase_kVA;
result.Obj_value = value(Obj);
result.soc_gap_max = max(abs(soc_gap(:)));
result.max_U = max(U_opt_pu(:));
result.min_U = min(U_opt_pu(:));
result.max_I2 = max(I2_opt(:));
result.max_I2_ratio = max(I2_ratio(:));
result.max_S_MVA = max(S_opt_pu(:)) * SB;
result.max_line_loading_ratio = max(line_loading_ratio(:));

tol = 1e-5;
soc_gap_tol = 1e-4;

if result.soc_gap_max > soc_gap_tol
    result.violation_type = "SOCP松弛间隙过大";
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
