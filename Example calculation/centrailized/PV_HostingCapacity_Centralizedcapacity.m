function result_HC = PV_HostingCapacity_Centralizedcapacity(Cap_init_kW, Cap_step_kW, Cap_max_kW)
%PV_HOSTINGCAPACITY_CENTRALIZEDCAPACITY
% 集中式无功电压控制下的光伏承载力评估
%
% 外层循环：逐步增加光伏容量
% 内层模型：给定光伏容量，求解集中式无功优化下的配电网运行状态

%% ========== 初始化 ==========
cap_kW = Cap_init_kW;

last_feasible_cap = 0;
last_feasible_result = [];

record = struct([]);
iter = 0;

%% ========== 外层循环 ==========
while cap_kW <= Cap_max_kW

    iter = iter + 1;

    fprintf('\n========== 第 %d 次集中式无功控制承载力评估 ==========\n', iter);
    fprintf('当前测试光伏总容量 = %.2f kW\n', cap_kW);

    %% 内层：集中式无功电压控制
    result = Centralizedcapacity(cap_kW);

    %% 记录结果
    record(iter).iter = iter;
    record(iter).cap_kW = cap_kW;
    record(iter).is_feasible = result.is_feasible;
    record(iter).max_U = result.max_U;
    record(iter).min_U = result.min_U;
    record(iter).max_I2 = result.max_I2;
    record(iter).max_line_loading_ratio = result.max_line_loading_ratio;
    record(iter).max_I2_ratio = result.max_I2_ratio;
    record(iter).violation_type = result.violation_type;

    %% 判断是否继续增加容量
    if result.is_feasible

        last_feasible_cap = cap_kW;
        last_feasible_result = result;

        fprintf('当前容量可行，继续增加光伏容量。\n');

    else

        fprintf('当前容量不可行，继续扫描后续容量点。\n');
        fprintf('越限类型：%s\n', result.violation_type);

    end

    cap_kW = cap_kW + Cap_step_kW;
end

%% ========== 输出 ==========
result_HC = struct();

result_HC.control_mode = "Centralized_Q_Control";
result_HC.HC_kW = last_feasible_cap;
result_HC.HC_MW = last_feasible_cap / 1000;
result_HC.best_result = last_feasible_result;

if isempty(record)
    result_HC.record = table();
else
    result_HC.record = struct2table(record);
end

end
