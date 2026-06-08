clc;
clear;
close all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(fileparts(scriptDir));

Cap_init_kW = 6000;
Cap_step_kW = 100;
Cap_max_kW  = 12000;

result_HC_NoQ = PV_HostingCapacity_NoQControl(Cap_init_kW, Cap_step_kW, Cap_max_kW);

disp('========== 无无功电压控制下光伏承载力 ==========');
disp(['光伏承载力 = ', num2str(result_HC_NoQ.HC_kW), ' kW']);
disp(['光伏承载力 = ', num2str(result_HC_NoQ.HC_MW), ' MW']);

disp(result_HC_NoQ.record);
