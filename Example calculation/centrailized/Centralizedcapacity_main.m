clc;
clear;
close all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fileparts(scriptDir));

Cap_init_kW = 12000;
Cap_step_kW = 100;
Cap_max_kW  = 14000;

result_HC_C = PV_HostingCapacity_Centralizedcapacity( ...
    Cap_init_kW, Cap_step_kW, Cap_max_kW);

disp('========== 集中式无功控制下光伏承载力 ==========');
disp(['光伏承载力 = ', num2str(result_HC_C.HC_kW), ' kW']);
disp(['光伏承载力 = ', num2str(result_HC_C.HC_MW), ' MW']);

disp(result_HC_C.record);
