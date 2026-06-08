function pf = radial_DistFlow_powerflow(Pd, Qd, Ppv, Qpv, L_i, L_j, R, X, Bc, NB, NL)
%RADIAL_DISTFLOW_POWERFLOW
% 径向配电网 DistFlow 前推回代潮流计算
%
% 输入：
% Pd, Qd   NB × 1，负荷，p.u.
% Ppv,Qpv  NB × 1，光伏注入，p.u.
%
% 输出：
% pf.U2    NB × 1，节点电压平方
% pf.I2    NL × 1，支路电流平方
% pf.P     NL × 1，支路首端有功
% pf.Q     NL × 1，支路首端无功

maxPFIter = 100;
tolPF = 1e-8;

%% ===== 构造各线路的子线路编号 =====
% Ind_Subline(i,:) 存放第 i 条线路末端节点的下游线路编号
% 本程序假设每条线路最多只有 2 条子线路
Ind_Subline = zeros(NL, 2);

for i = 1:NL

    temp = find(L_i == L_j(i));

    if length(temp) > 2
        error('第 %d 条线路存在超过两个子支路，当前 Ind_Subline 维度不适用。', i);
    end

    if ~isempty(temp)
        Ind_Subline(i, 1:length(temp)) = temp;
    end

end

%% ===== 在函数内部生成线路拓扑顺序 =====
% order 是从上游到下游的线路顺序，用于前推计算电压
% reverse_order 是从下游到上游的线路顺序，用于后推计算功率

knownBus = false(NB, 1);
knownBus(1) = true;          % 默认 1 号节点为平衡节点

usedLine = false(NL, 1);
order = zeros(NL, 1);

cnt = 0;

while cnt < NL

    candidate = find(~usedLine & knownBus(L_i));

    if isempty(candidate)
        error('网络拓扑错误：无法从 1 号节点生成径向拓扑顺序。');
    end

    for m = 1:length(candidate)

        ell = candidate(m);

        cnt = cnt + 1;
        order(cnt) = ell;

        usedLine(ell) = true;
        knownBus(L_j(ell)) = true;

    end

end

reverse_order = flipud(order);

%% ===== 初始化 =====
U2 = ones(NB, 1);
I2 = zeros(NL, 1);
P = zeros(NL, 1);
Q = zeros(NL, 1);

converged = false;

%% ===== 前推回代迭代 =====
for iter = 1:maxPFIter

    U2_old = U2;

    %% ===== 后推：由下游向上游计算支路功率 =====
    for idx = 1:length(reverse_order)

        ell = reverse_order(idx);

        from = L_i(ell);
        to = L_j(ell);

        child = Ind_Subline(ell, :);
        child = child(child > 0);

        if isempty(child)
            P_down = 0;
            Q_down = 0;
        else
            P_down = sum(P(child));
            Q_down = sum(Q(child));
        end

        % 线路对地电纳无功注入
        QcLine = Bc(ell) * U2(from) + Bc(ell) * U2(to);

        % 有功功率平衡：
        % P - R*I2 + Ppv = P_down + Pd
        P(ell) = P_down + Pd(to) - Ppv(to) + R(ell) * I2(ell);

        % 无功功率平衡：
        % Q - X*I2 + Qpv + QcLine = Q_down + Qd
        Q(ell) = Q_down + Qd(to) - Qpv(to) - QcLine + X(ell) * I2(ell);

    end

    %% ===== 前推：由上游向下游计算节点电压 =====
    U2 = zeros(NB, 1);

    % 对应原优化模型中的 U2(1,:) == 1.0^2
    U2(1) = 1.0;

    invalid_voltage = false;

    for idx = 1:length(order)

        ell = order(idx);

        from = L_i(ell);
        to = L_j(ell);

        U2(to) = U2(from) ...
            - 2 * (R(ell) * P(ell) + X(ell) * Q(ell)) ...
            + (R(ell)^2 + X(ell)^2) * I2(ell);

        if U2(to) <= 0
            invalid_voltage = true;
            break;
        end

    end

    if invalid_voltage
        converged = false;
        break;
    end

    %% ===== 更新支路电流平方 =====
    for ell = 1:NL

        from = L_i(ell);

        if U2(from) <= 0
            invalid_voltage = true;
            break;
        end

        % 对应 DistFlow 电流关系：
        % I2 * U2 = P^2 + Q^2
        I2(ell) = (P(ell)^2 + Q(ell)^2) / U2(from);

    end

    if invalid_voltage
        converged = false;
        break;
    end

    %% ===== 潮流收敛判断 =====
    if max(abs(U2 - U2_old)) <= tolPF
        converged = true;
        break;
    end

end

%% ===== 输出 =====
pf = struct();

pf.U2 = U2;
pf.I2 = I2;
pf.P = P;
pf.Q = Q;
pf.iter = iter;
pf.converged = converged;

end