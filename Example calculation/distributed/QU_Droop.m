function Qcmd = QU_Droop(U, Qmax, loc2, NB)
%QU_Droop
% 分散式 Q-U 下垂控制数值更新函数
%
% 输入：
% U     NB × 1，节点电压幅值，p.u.
% Qmax  NB × 1，当前时刻 PV 可用无功容量，kvar
% loc2  PV 接入节点
% NB    节点数
%
% 输出：
% Qcmd  NB × 1，PV 无功指令，kvar
%
% 符号约定：
% Qcmd > 0：逆变器发出无功，支撑低电压
% Qcmd < 0：逆变器吸收无功，抑制高电压

Qcmd = zeros(NB, 1);

%% Q-U 下垂曲线参数
U1  = 0.90;
U2d = 0.98;
U3d = 1.02;
U4  = 1.10;

for k = 1:length(loc2)

    node = loc2(k);

    u = U(node);
    Qm = Qmax(node);

    if Qm <= 1e-6

        Qcmd(node) = 0;

    elseif u <= U1

        Qcmd(node) = Qm;

    elseif u > U1 && u < U2d

        Qcmd(node) = Qm * (U2d - u) / (U2d - U1);

    elseif u >= U2d && u <= U3d

        Qcmd(node) = 0;

    elseif u > U3d && u < U4

        Qcmd(node) = -Qm * (u - U3d) / (U4 - U3d);

    else

        Qcmd(node) = -Qm;

    end

end

end
