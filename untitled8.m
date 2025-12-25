function experiment1_attitude_compare()
% 实验一：无扰动下的快速收敛性验证（基准性能）
% 对比：PID / 常规 SMC / FFTSMC（基于四元数）
%
% 输出图：
% 1) 欧拉角响应对比（roll/pitch/yaw）
% 2) 误差四元数矢量部 e_v 收敛曲线
% 3) 相平面图：phi vs phidot（以 roll 为例）

clc; clear; close all;

%% ========== 仿真参数 ==========
Tend = 8;                 % 仿真时长(s)
dt_plot = 0.001;          % 插值采样用于画图
t_grid = (0:dt_plot:Tend)';

% 转动惯量（可按你的机体改）
J = diag([0.02, 0.02, 0.04]);

% 期望姿态：单位四元数（对齐）
qd = [1; 0; 0; 0];

% 初始姿态（给定欧拉角，按 ZYX: yaw-pitch-roll 转四元数）
roll0  = deg2rad(30);
pitch0 = deg2rad(-20);
yaw0   = deg2rad(45);
q0 = eulZYX_to_quat(roll0, pitch0, yaw0);  % [q0;q1;q2;q3]
q0 = quat_normalize(q0);

% 初始角速度（机体系）
w0 = [0.1; -0.1; 0.05];

% 外部扰动（本实验：无扰动）
d = [0; 0; 0]; %#ok<NASGU>

%% ========== 控制器参数 ==========
% --- PID ---
pid.Kp = diag([6, 6, 4]);
pid.Kd = diag([2.0, 2.0, 1.5]);
pid.Ki = diag([1.0, 1.0, 0.8]);
pid.int_limit = [0.5; 0.5; 0.5];     % 积分限幅（避免发散）

% --- 常规 SMC ---
smc.c   = 4.0;                        % 滑模面系数
smc.K   = diag([3.0, 3.0, 2.5]);      % 线性到达项
smc.eta = diag([0.6, 0.6, 0.5]);      % 开关项
smc.phi = 0.05;                       % 饱和边界层宽度（抑制抖振）

% --- FFTSMC（快速有限时间终端滑模） ---
% 采用 s = e_dot + λ1*sig(e)^(p/q) + λ2*sig(e)^(q/p)（结构见文档(3.3)）
% 到达律采用双幂次项（形式见文档(3.13)）
ffts.p = 1; ffts.q = 2;               % p/q=0.5, q/p=2
ffts.lambda1 = 3.0;
ffts.lambda2 = 1.2;

ffts.alpha = 0.6;                     % 0<alpha<1
ffts.beta  = 1.2;                     % beta>1
ffts.K1 = diag([6.0, 6.0, 5.0]);
ffts.K2 = diag([1.5, 1.5, 1.2]);
ffts.K3 = diag([1.2, 1.2, 1.0]);      % 线性阻尼项（增强数值稳定性）

%% ========== 初值 ==========
x0_pid  = [q0; w0; zeros(3,1)];       % PID: [q(4); w(3); int_e(3)]
x0_smc  = [q0; w0];                   % SMC: [q(4); w(3)]
x0_ffts = [q0; w0];                   % FFTSMC: [q(4); w(3)]

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

%% ========== 数值积分 ==========
[t_pid,  x_pid]  = ode45(@(t,x) dyn_pid(t,x,J,qd,pid), [0 Tend], x0_pid,  opts);
[t_smc,  x_smc]  = ode45(@(t,x) dyn_smc(t,x,J,qd,smc), [0 Tend], x0_smc,  opts);
[t_ffts, x_ffts] = ode45(@(t,x) dyn_fftsmc(t,x,J,qd,ffts), [0 Tend], x0_ffts, opts);

%% ========== 插值到统一时间轴（方便画对比图） ==========
Xp = interp1(t_pid,  x_pid,  t_grid, 'pchip');
Xs = interp1(t_smc,  x_smc,  t_grid, 'pchip');
Xf = interp1(t_ffts, x_ffts, t_grid, 'pchip');

qp = normalize_rows_quat(Xp(:,1:4));
qs = normalize_rows_quat(Xs(:,1:4));
qf = normalize_rows_quat(Xf(:,1:4));

wp = Xp(:,5:7);
ws = Xs(:,5:7);
wf = Xf(:,5:7);

% 误差四元数 qe = qd^{-1} ⊗ q
[qe_p, ev_p] = quat_error_series(qp, qd);
[qe_s, ev_s] = quat_error_series(qs, qd);
[qe_f, ev_f] = quat_error_series(qf, qd);

% 欧拉角（用于展示：roll/pitch/yaw）
eul_p = quat_to_eulZYX_series(qp);
eul_s = quat_to_eulZYX_series(qs);
eul_f = quat_to_eulZYX_series(qf);

%% ========== 图1：欧拉角响应对比 ==========
figure('Name','Euler Angle Response','Color','w');
subplot(3,1,1);
plot(t_grid, rad2deg(eul_p(:,1)),'LineWidth',1.2); hold on;
plot(t_grid, rad2deg(eul_s(:,1)),'LineWidth',1.2);
plot(t_grid, rad2deg(eul_f(:,1)),'LineWidth',1.2);
grid on; ylabel('\phi (deg)'); title('Roll Response');
legend('PID','SMC','FFTSMC','Location','best');

subplot(3,1,2);
plot(t_grid, rad2deg(eul_p(:,2)),'LineWidth',1.2); hold on;
plot(t_grid, rad2deg(eul_s(:,2)),'LineWidth',1.2);
plot(t_grid, rad2deg(eul_f(:,2)),'LineWidth',1.2);
grid on; ylabel('\theta (deg)'); title('Pitch Response');

subplot(3,1,3);
plot(t_grid, rad2deg(eul_p(:,3)),'LineWidth',1.2); hold on;
plot(t_grid, rad2deg(eul_s(:,3)),'LineWidth',1.2);
plot(t_grid, rad2deg(eul_f(:,3)),'LineWidth',1.2);
grid on; ylabel('\psi (deg)'); xlabel('Time (s)'); title('Yaw Response');

%% ========== 图2：姿态误差四元数矢量部 e_v 收敛曲线 ==========
figure('Name','Quaternion Vector Error','Color','w');
subplot(3,1,1);
plot(t_grid, ev_p(:,1),'LineWidth',1.2); hold on;
plot(t_grid, ev_s(:,1),'LineWidth',1.2);
plot(t_grid, ev_f(:,1),'LineWidth',1.2);
grid on; ylabel('e_1'); title('Error Quaternion Vector Part');
legend('PID','SMC','FFTSMC','Location','best');

subplot(3,1,2);
plot(t_grid, ev_p(:,2),'LineWidth',1.2); hold on;
plot(t_grid, ev_s(:,2),'LineWidth',1.2);
plot(t_grid, ev_f(:,2),'LineWidth',1.2);
grid on; ylabel('e_2');

subplot(3,1,3);
plot(t_grid, ev_p(:,3),'LineWidth',1.2); hold on;
plot(t_grid, ev_s(:,3),'LineWidth',1.2);
plot(t_grid, ev_f(:,3),'LineWidth',1.2);
grid on; ylabel('e_3'); xlabel('Time (s)');

%% ========== 图3：相平面图（以 roll 误差 phi vs phidot 为例） ==========
% 这里 phi 使用欧拉角 roll，phidot 用数值微分（rad/s）
phi_p = eul_p(:,1); phi_s = eul_s(:,1); phi_f = eul_f(:,1);
phidot_p = gradient(phi_p, dt_plot);
phidot_s = gradient(phi_s, dt_plot);
phidot_f = gradient(phi_f, dt_plot);

figure('Name','Phase Plot (Roll)','Color','w');
plot(phi_p, phidot_p, 'LineWidth',1.2); hold on;
plot(phi_s, phidot_s, 'LineWidth',1.2);
plot(phi_f, phidot_f, 'LineWidth',1.2);
grid on; xlabel('\phi (rad)'); ylabel('\dot{\phi} (rad/s)');
title('Phase Plot: \phi vs \dot{\phi}');
legend('PID','SMC','FFTSMC','Location','best');

%% ========== 可选：打印收敛信息 ==========
fprintf('Done.\n');
fprintf('Final |e_v| (PID ) = %.3e\n', norm(ev_p(end,:)));
fprintf('Final |e_v| (SMC ) = %.3e\n', norm(ev_s(end,:)));
fprintf('Final |e_v| (FFTS) = %.3e\n', norm(ev_f(end,:)));

end

%% ======================= 动力学：PID =======================
function dx = dyn_pid(~, x, J, qd, pid)
q  = quat_normalize(x(1:4));
w  = x(5:7);
ie = x(8:10);

% 误差四元数 qe = qd^{-1} ⊗ q
qe = quat_mul(quat_conj(qd), q);
qe = shortest_quat(qe);
e  = qe(2:4);                 % 误差矢量部

% 积分（限幅）
ie = ie + e * 0; %#ok<NASGU> % 留给 ODE 积分（此行无效，仅占位说明）

% 控制：期望角加速度命令
w_dot_cmd = -pid.Kp*e - pid.Kd*w - pid.Ki*ie;

% 刚体动力学：J*w_dot + w×(Jw) = tau
tau = cross(w, J*w) + J*w_dot_cmd;
w_dot = J \ (tau - cross(w, J*w));  %#ok<NASGU> (等于 w_dot_cmd)

% 四元数运动学
q_dot = quat_kinematics(q, w);

% 积分误差动态（抗积分饱和）
ie_dot = clamp_vec(e, pid.int_limit);

dx = [q_dot; w_dot_cmd; ie_dot];
end

%% ======================= 动力学：SMC =======================
function dx = dyn_smc(~, x, J, qd, smc)
q = quat_normalize(x(1:4));
w = x(5:7);

qe = quat_mul(quat_conj(qd), q);
qe = shortest_quat(qe);
e  = qe(2:4);

% e_dot = 0.5*(qe0*I + skew(e))*w   （四元数误差矢量部导数）
qe0 = qe(1);
A = 0.5*(qe0*eye(3) + skew3(e));
e_dot = A*w;

% 滑模面：s = w + c*e
s = w + smc.c*e;

% 饱和函数（边界层）
sat_s = sat_vec(s, smc.phi);

% 令 s_dot = w_dot + c*e_dot = -(K*s + eta*sat)
w_dot_cmd = -smc.c*e_dot - smc.K*s - smc.eta*sat_s;

tau = cross(w, J*w) + J*w_dot_cmd; %#ok<NASGU>
q_dot = quat_kinematics(q, w);

dx = [q_dot; w_dot_cmd];
end

%% ======================= 动力学：FFTSMC =======================
function dx = dyn_fftsmc(~, x, J, qd, ffts)
q = quat_normalize(x(1:4));
w = x(5:7);

qe = quat_mul(quat_conj(qd), q);
qe = shortest_quat(qe);
e  = qe(2:4);
qe0 = qe(1);

% e_dot = A*w
A = 0.5*(qe0*eye(3) + skew3(e));
e_dot = A*w;

% s = e_dot + λ1*sig(e)^(p/q) + λ2*sig(e)^(q/p)   （结构见文档(3.3)）
r1 = ffts.p/ffts.q;   % <1
r2 = ffts.q/ffts.p;   % >1

f1 = ffts.lambda1 * sig_pow(e, r1);
f2 = ffts.lambda2 * sig_pow(e, r2);
s  = e_dot + f1 + f2;

% 为实现 s_dot ≈ -K1*sig(s)^alpha - K2*sig(s)^beta - K3*s
% 近似处理：忽略 A_dot*w，并令
% s_dot = e_ddot + df/de * e_dot
% 选 e_ddot_cmd = - dfde_e_dot - K1*sig(s)^alpha - K2*sig(s)^beta - K3*s

eps0 = 1e-6;
dfde1 = ffts.lambda1 * r1 * (abs(e)+eps0).^(r1-1);
dfde2 = ffts.lambda2 * r2 * (abs(e)+eps0).^(r2-1);
dfde_e_dot = dfde1 .* e_dot + dfde2 .* e_dot;

e_ddot_cmd = -dfde_e_dot ...
            - ffts.K1*sig_pow(s, ffts.alpha) ...
            - ffts.K2*sig_pow(s, ffts.beta)  ...
            - ffts.K3*s;

% 近似：e_ddot = A*w_dot_cmd  => w_dot_cmd = pinv(A)*e_ddot_cmd
w_dot_cmd = pinv(A) * e_ddot_cmd;

tau = cross(w, J*w) + J*w_dot_cmd; %#ok<NASGU>
q_dot = quat_kinematics(q, w);

dx = [q_dot; w_dot_cmd];
end

%% ======================= 工具函数 =======================
function qdot = quat_kinematics(q, w)
% q = [q0;qv], w in body frame
q0 = q(1); qv = q(2:4);
qdot0 = -0.5 * (qv.' * w);
qdotv =  0.5 * (q0*eye(3) + skew3(qv)) * w;
qdot = [qdot0; qdotv];
end

function S = skew3(v)
S = [  0   -v(3)  v(2);
      v(3)  0    -v(1);
     -v(2) v(1)   0  ];
end

function qn = quat_normalize(q)
qn = q ./ max(norm(q), 1e-12);
end

function qcon = quat_conj(q)
qcon = [q(1); -q(2:4)];
end

function qc = quat_mul(a, b)
% Hamilton product, both as [q0;q1;q2;q3]
a0=a(1); av=a(2:4);
b0=b(1); bv=b(2:4);
c0 = a0*b0 - dot(av,bv);
cv = a0*bv + b0*av + cross(av,bv);
qc = [c0; cv];
end

function qe = shortest_quat(qe)
% 取最短旋转（避免 q 与 -q 的等价导致控制抖动）
if qe(1) < 0
    qe = -qe;
end
end

function y = sig_pow(x, r)
% sig(x)^r := |x|^r * sign(x)（逐元素）
y = (abs(x)).^r .* sign(x);
end

function y = sat_vec(x, phi)
% 饱和：sat(x/phi) * phi? 这里输出为 sat(x/phi)（保持量纲与 sign 一致）
y = x ./ max(phi, 1e-12);
y = max(-1, min(1, y));
end

function y = clamp_vec(x, lim)
% 逐元素限幅到 [-lim, lim]
y = min(max(x, -lim), lim);
end

function q = eulZYX_to_quat(roll, pitch, yaw)
% ZYX: yaw(Z) - pitch(Y) - roll(X)
cy = cos(yaw*0.5);  sy = sin(yaw*0.5);
cp = cos(pitch*0.5);sp = sin(pitch*0.5);
cr = cos(roll*0.5); sr = sin(roll*0.5);

q0 = cr*cp*cy + sr*sp*sy;
q1 = sr*cp*cy - cr*sp*sy;
q2 = cr*sp*cy + sr*cp*sy;
q3 = cr*cp*sy - sr*sp*cy;

q = [q0; q1; q2; q3];
end

function eul = quat_to_eulZYX_series(Q)
% Q: Nx4, each row [q0 q1 q2 q3]
N = size(Q,1);
eul = zeros(N,3);
for i=1:N
    q = Q(i,:)';
    q = quat_normalize(q);
    q0=q(1); q1=q(2); q2=q(3); q3=q(4);

    % roll
    sinr = 2*(q0*q1 + q2*q3);
    cosr = 1 - 2*(q1*q1 + q2*q2);
    roll = atan2(sinr, cosr);

    % pitch
    sinp = 2*(q0*q2 - q3*q1);
    sinp = max(-1, min(1, sinp));
    pitch = asin(sinp);

    % yaw
    siny = 2*(q0*q3 + q1*q2);
    cosy = 1 - 2*(q2*q2 + q3*q3);
    yaw = atan2(siny, cosy);

    eul(i,:) = [roll, pitch, yaw];
end
end

function Qn = normalize_rows_quat(Q)
Qn = Q;
for i=1:size(Q,1)
    Qn(i,:) = (Q(i,:)./max(norm(Q(i,:)),1e-12));
end
end

function [QE, EV] = quat_error_series(Q, qd)
% Q: Nx4 (current), qd: 4x1
N = size(Q,1);
QE = zeros(N,4);
EV = zeros(N,3);
qd_conj = quat_conj(qd);
for i=1:N
    q = Q(i,:)';
    qe = quat_mul(qd_conj, q);
    qe = shortest_quat(qe);
    QE(i,:) = qe.';
    EV(i,:) = qe(2:4).';
end
end
