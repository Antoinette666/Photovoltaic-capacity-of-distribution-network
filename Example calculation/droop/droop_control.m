function result = droop_control(cap_kW)
% 给定光伏总容量下的分散式无功优化模型
%内层控制
% 输入：
% cap_total_kW  当前测试的光伏总接入容量，kW
%
% 输出：
% result        当前容量下的求解结果和可行性判断
format short

NT=24;%优化时间
SB=10;%基准功率MVA标幺化基准值
VBase=12.66;%基准电压KV
ZB=VBase^2/SB;%基准电抗kv

%% ===== 从数据文件加载系统参数 =====
[PDN_Data, Pd_total, Qd_total, loc2, solar_data] = get_PDN_Data();

NL=size(PDN_Data,1);
NB=NL+1;%NL为线路数
L_i=PDN_Data(:,2);
L_j=PDN_Data(:,3);
P_load=PDN_Data(:,4);
Q_load=PDN_Data(:,5);
R=PDN_Data(:,6)/ZB;
X=PDN_Data(:,7)/ZB;
Smax=PDN_Data(:,8)/SB;

Bc = zeros(NL, 1); % 与无无功电压控制模型保持一致

U_min = 0.90;
U_max = 1.10;
I2_max = (Smax / U_min).^2;

%% ===== 负荷转换为 p.u.，并映射到节点 =====
Pd=P_load*Pd_total/1000000;%转换成PU
Qd=Q_load*Qd_total/1000000;%转换成PU
Pd=[zeros(1,24);Pd];%加了一行0
Qd=[zeros(1,24);Qd];

%% ========== PV 光伏变量 ==========
Psolar = zeros(NB, NT);              % kW，固定预测值   Psolar=sdpvar(33,24);
nPV = length(loc2);

%% ========== 按承载力容量缩放光伏有功 ==========
% cap_kW 是当前外层循环测试的总光伏容量。
% solar_data 只作为日内出力形状曲线，先归一化，再按容量缩放。

if nPV == 0 || cap_kW <= 1e-6
    pv_cap_each = zeros(nPV, 1);
else
    pv_cap_each = cap_kW / nPV * ones(nPV, 1);%平均分配到每个光伏节点
end

pv_shape = zeros(nPV, NT);

for k = 1:nPV

    temp = solar_data(k, :);

    if max(temp) <= 1e-6
        pv_shape(k, :) = zeros(1, NT);
    else
        pv_shape(k, :) = temp / max(temp);%归一化为 0~1 的日内出力曲线
    end

end

for k = 1:nPV
    node = loc2(k);
    Psolar(node, :) = pv_cap_each(k) * pv_shape(k, :);
end

%% ========== 可用无功容量 ==========
Qmax = zeros(NB, NT);
inv_oversize = 1.1;   % 逆变器容量超配至光伏额定容量的 1.1 倍

for k = 1:nPV
    node = loc2(k);
    S_node = inv_oversize * pv_cap_each(k);   % kVA
    Qmax(node, :) = sqrt(max(0, S_node^2 - Psolar(node, :).^2));
end

%% ===== 分散式下垂控制迭代参数 =====
maxIter = 500;
tol_U = 1e-5;
tol_Q = 1e-3;     % kvar
alpha_init = 0.20;
alpha_min = 0.01;
alpha_max = 0.30;

Ind_NoPV = setdiff(1:NB, loc2);

%% ===== 结果变量 =====
U_all = zeros(NB, NT);
U2_all = zeros(NB, NT);
I2_all = zeros(NL, NT);
P_all = zeros(NL, NT);
Q_all = zeros(NL, NT);

Qsolar = zeros(NB, NT);       % kvar

iter_record = zeros(NT, 1);
droop_converged_record = false(NT, 1);
pf_converged_record = false(NT, 1);

%% ===== 逐时段进行分散式下垂控制仿真 =====
for t = 1:NT

    if t == 1
        Qpv_old = zeros(NB, 1);
    else
        % 以上一时段的收敛结果作为初值
        Qpv_old = Qsolar(:, t - 1);
    end

    U_old = ones(NB, 1);
    alpha = alpha_init;
    prev_err_Q = inf;

    droop_converged = false;

    for iter = 1:maxIter

        %% 当前 Qpv 下进行潮流计算
        Ppv_pu = Psolar(:, t) / (SB * 1000);   % kW -> p.u.
        Qpv_pu = Qpv_old / (SB * 1000);        % kvar -> p.u.
        pf = radial_DistFlow_powerflow( ...
            Pd(:, t), Qd(:, t), ...
            Ppv_pu, Qpv_pu, ...
            L_i, L_j, R, X, Bc, NB, NL);

        U_now = sqrt(max(pf.U2, 0));

        %% 各 PV 节点根据本地电压进行 Q-U 下垂更新
        Qcmd = QU_Droop(U_now, Qmax(:, t), loc2, NB);

        %% 阻尼更新，避免振荡
        Qpv_new = Qpv_old;

        if ~isempty(loc2)
            Qpv_new(loc2) = (1 - alpha) * Qpv_old(loc2) + alpha * Qcmd(loc2);
        end

        %% 非 PV 节点无功出力为 0
        Qpv_new(Ind_NoPV) = 0;

        %% 收敛判断
        err_U = max(abs(U_now - U_old));
        err_Q = max(abs(Qpv_new - Qpv_old));

        if err_U <= tol_U && err_Q <= tol_Q && pf.converged
            droop_converged = true;
            Qpv_old = Qpv_new;
            break;
        end

        if err_Q > 1.02 * prev_err_Q
            alpha = max(alpha_min, 0.5 * alpha);
        elseif err_Q < 0.70 * prev_err_Q
            alpha = min(alpha_max, 1.05 * alpha);
        end

        prev_err_Q = err_Q;
        U_old = U_now;
        Qpv_old = Qpv_new;

    end

    %% 用最终 Qpv 再计算一次潮流，保证存储结果和最终控制量一致
    Ppv_pu = Psolar(:, t) / (SB * 1000);
    Qpv_pu = Qpv_old / (SB * 1000);

    pf = radial_DistFlow_powerflow( ...
        Pd(:, t), Qd(:, t), ...
        Ppv_pu, Qpv_pu, ...
        L_i, L_j, R, X, Bc, NB, NL);

    U_now = sqrt(max(pf.U2, 0));

    %% 存储结果
    Qsolar(:, t) = Qpv_old;

    U_all(:, t) = U_now;
    U2_all(:, t) = pf.U2;
    I2_all(:, t) = pf.I2;
    P_all(:, t) = pf.P;
    Q_all(:, t) = pf.Q;

    iter_record(t) = iter;
    droop_converged_record(t) = droop_converged;
    pf_converged_record(t) = pf.converged;

end

%% ===== 初始化输出 =====
result = struct();

result.cap_total_kW = cap_kW;
result.is_feasible = false;
result.violation_type = "";

result.max_U = max(U_all(:));
result.min_U = min(U_all(:));
result.max_I2 = max(I2_all(:));
line_loading_ratio = hypot(P_all, Q_all) ./ repmat(Smax, 1, NT);
I2_ratio = I2_all ./ repmat(I2_max, 1, NT);
result.max_line_loading_ratio = max(line_loading_ratio(:));
result.max_I2_ratio = max(I2_ratio(:));

result.Obj_value = NaN;
result.sol_info = "No optimization. Distributed Q-U droop simulation.";

result.U_opt_pu = U_all;
result.U2_opt = U2_all;
result.I2_opt = I2_all;

result.P_branch_pu = P_all;
result.Q_branch_pu = Q_all;

result.Psolar_kW = Psolar;
result.Qsolar_opt_kvar = Qsolar;

result.iter_record = iter_record;
result.droop_converged_record = droop_converged_record;
result.pf_converged_record = pf_converged_record;

result.converged = all(droop_converged_record) && all(pf_converged_record);

%% ===== 越限判断 =====
tol = 1e-5;

if ~result.converged

    result.is_feasible = false;
    result.violation_type = "分散式下垂控制或潮流计算未收敛";

elseif result.max_U > U_max + tol

    result.is_feasible = false;
    result.violation_type = "节点电压越上限";

elseif result.min_U < U_min - tol

    result.is_feasible = false;
    result.violation_type = "节点电压越下限";

elseif result.max_line_loading_ratio > 1 + tol

    result.is_feasible = false;
    result.violation_type = "线路视在功率越限";

elseif result.max_I2_ratio > 1 + tol

    result.is_feasible = false;
    result.violation_type = "线路电流越限";

else
    result.is_feasible = true;
    result.violation_type = "无越限";

end

end
