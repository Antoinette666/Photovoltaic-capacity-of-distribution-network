clc;
clear;
close all;

addpath(fileparts(fileparts(mfilename('fullpath'))));

Cap_init_kW = 8000;
Cap_step_kW = 100;
Cap_max_kW  = 14000;

result_HC_C = PV_HostingCapacity_QUD(Cap_init_kW, Cap_step_kW, Cap_max_kW);

disp('========== 分散式 Q-U 控制下光伏承载力 ==========');
disp(['光伏承载力 = ', num2str(result_HC_C.HC_kW), ' kW']);
disp(['光伏承载力 = ', num2str(result_HC_C.HC_MW), ' MW']);

disp(result_HC_C.record);
