%% 清空环境
clear; clc;

%% ===== 从数据文件加载系统参数 =====
[PDN_Data, Pd_total, Qd_total, loc2, solar_data] = get_PDN_Data();

%% 2. 将您的自定义数据转换为 MATPOWER 兼容的标准 mpc 结构体
SB=10;%基准功率MVA标幺化基准值
VBase=12.66;%基准电压KV
ZB=VBase^2/SB;%基准电抗kv

num_buses = 33;              % 节点数
num_branches = size(PDN_Data, 1); % 线路数

mpc.version = '2';
mpc.baseMVA = SB;

% 2.1 构建 bus 矩阵
% MATPOWER格式: [bus_i, type, Pd, Qd, Gs, Bs, area, Vm, Va, baseKV, zone, maxV, minV]
mpc.bus = zeros(num_buses, 13);
mpc.bus(:, 1) = 1:num_buses;    % 节点编号
mpc.bus(:, 2) = 1;              % 默认设为 PQ 节点 (type=1)
mpc.bus(1, 2) = 3;              % 1号节点设为平衡节点 (type=3)
mpc.bus(:, 8) = 1.0;            % 初始电压幅值 1.0 p.u.
mpc.bus(:, 9) = 0.0;            % 初始电压相角 0.0
mpc.bus(:, 10) = VBase;        % 基准电压
mpc.bus(:, 12) = 1.1;           % 电压上限
mpc.bus(:, 13) = 0.9;           % 电压下限

% 将您的 kW 和 kVar 负荷转换为 MW 和 MVAr，并映射到线路的末节点(to_bus)
for i = 1:num_branches
    to_bus = PDN_Data(i, 3);
    mpc.bus(to_bus, 3) = PDN_Data(i, 4) / 1000; % 转换为 MW
    mpc.bus(to_bus, 4) = PDN_Data(i, 5) / 1000; % 转换为 MVAr
end

% 2.2 构建 branch 矩阵
% MATPOWER格式: [fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax]
mpc.branch = zeros(num_branches, 13);
mpc.branch(:, 1) = PDN_Data(:, 2);       % 首节点 (From)
mpc.branch(:, 2) = PDN_Data(:, 3);       % 末节点 (To)
mpc.branch(:, 3) = PDN_Data(:, 6) / ZB; % 阻抗 R 转换为标幺值(p.u.)
mpc.branch(:, 4) = PDN_Data(:, 7) / ZB; % 阻抗 X 转换为标幺值(p.u.)
mpc.branch(:, 11) = 1;                   % 投入运行状态
mpc.branch(:, 12) = -360;                % 相角差下限
mpc.branch(:, 13) = 360;                 % 相角差上限

% 2.3 构建 gen 矩阵 (平衡节点发电机)
% MATPOWER格式: [bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, ...]
mpc.gen = [1, 0, 0, 999, -999, 1.0, SB, 1, 999, 0];

%% 3. 运行潮流计算
opt = mpoption('verbose', 0, 'out.all', 0); % 关闭多余的窗口打印
results = runpf(mpc, opt);

%% 4. 提取全网雅可比矩阵，并只截取无功-电压一阶导数分块 (dQ/dV)
J_full = makeJac(results, 1); % 获取全网完整雅可比矩阵 (2nb x 2nb)
nb = size(results.bus, 1);    % 节点总数 (33)

% 提取出 33x33 的全网纯无功-电压雅可比矩阵
J22_full = J_full(nb+1:2*nb, nb+1:2*nb); 

%% 5. 剔除平衡节点（1号节点）
noslack = 2:nb; % 2 到 33 号节点
J22 = J22_full(noslack, noslack); % 得到 32x32 的一阶雅可比矩阵

%% 6. 直接求逆得到无功-电压灵敏度矩阵
S_VQ_sparse = inv(J22); % 此时得到的仍是稀疏矩阵

% 【核心修正】将稀疏矩阵转换为常规的全矩阵（稠密矩阵），避免 fprintf 报错
S_VQ = full(S_VQ_sparse); 

%% 5. 结果输出与绘图
self_sensitivity = diag(S_VQ); % 提取自灵敏度
[max_val, max_node_idx] = max(self_sensitivity);
max_node = noslack(max_node_idx);

fprintf('--- 33节点无功-电压灵敏度分析 ---\n');
fprintf('电压对无功最敏感的薄弱节点是: %d 号节点\n', max_node);
fprintf('该节点的自灵敏度值为: %f (p.u./MVar)\n', max_val);

% 绘制曲线图
figure;
plot(noslack, self_sensitivity, '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on;
set(gca, 'XTick', noslack); % 显示所有节点坐标
xlabel('节点编号 (Bus Number)');
ylabel('无功-电压自灵敏度 \partialV / \partialQ');
title('自定义33节点配电网无功-电压自灵敏度分析');