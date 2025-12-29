function main_sim1()
clc; clear; close all;

%% ========== 仿真工况：无扰动 + 大初始误差自稳 ==========
sim.Tend = 8;
sim.dt   = 0.001;

% inertia
sim.J = diag([0.054, 0.053, 0.089]);

% desired attitude & rate
sim.qd = [1;0;0;0];
% sim.wd = [0;0;0]; % 本仿真不需要用到

% initial attitude: axis-angle -> quaternion (大初始误差)
axis0  = [1; -0.3; 0.5];
theta0 = deg2rad(100);                 % 大误差（可改 60~140deg）
sim.q0  = axis_angle_to_quat_local(axis0, theta0);

% initial angular rate (非零中等幅值)
sim.w0 = [0.2; -0.15; 0.10];

% torque limit（可选，防止数值爆掉；不想限幅可设 inf）
sim.tau_limit = [8; 8; 8];             % N·m（按你机体调整）

%% ========== 控制器参数 ==========
% PID
sim.pid.Kp = diag([6, 6, 4]);
sim.pid.Kd = diag([2.0, 2.0, 1.5]);
sim.pid.Ki = diag([1.0, 1.0, 0.8]);
sim.pid.int_limit = [0.5; 0.5; 0.5];

% SMC
sim.smc.c   = 4.0;
sim.smc.K   = diag([3.0, 3.0, 2.5]);
sim.smc.eta = diag([0.6, 0.6, 0.5]);
sim.smc.phi = 0.05;

% FFTSMC
sim.ffts.p = 1; sim.ffts.q = 2;        % p/q=0.5
sim.ffts.lambda1 = 3.0;
sim.ffts.lambda2 = 1.2;
sim.ffts.alpha = 0.6;                  % 0<alpha<1
sim.ffts.beta  = 1.2;                  % beta>1
sim.ffts.K1 = diag([6.0, 6.0, 5.0]);
sim.ffts.K2 = diag([1.5, 1.5, 1.2]);
sim.ffts.K3 = diag([1.2, 1.2, 1.0]);

% FTDO（给 FFTSMC+FTDO 用）
sim.ftdo.L     = diag([15,15,15]);
sim.ftdo.Ko1   = diag([6,6,6]);
sim.ftdo.Ko2   = diag([3,3,3]);
sim.ftdo.alpha = 0.6;

%% ========== 一键仿真（对比：FFTSMC+FTDO / PID / SMC） ==========
res_pid  = attitude_sim('PID',          sim);
res_smc  = attitude_sim('SMC',          sim);
res_ffts = attitude_sim('FFTSMC_FTDO',  sim);

t = res_pid.t;

%% ========== 必画图1：姿态误差角 theta_e(t) ==========
figure('Name','Attitude Error Angle','Color','w');
plot(t, rad2deg(res_ffts.theta_e),'LineWidth',1.3); hold on;
plot(t, rad2deg(res_pid.theta_e),'LineWidth',1.3);
plot(t, rad2deg(res_smc.theta_e),'LineWidth',1.3);
grid on; xlabel('Time (s)'); ylabel('\theta_e (deg)');
title('Attitude Error Angle \theta_e(t)');
legend('FFTSMC+FTDO','PID','SMC','Location','best');

%% ========== 必画图2：误差四元数矢量部 e=[e1,e2,e3] ==========
figure('Name','Error Quaternion Vector Part','Color','w');
for i=1:3
    subplot(3,1,i);
    plot(t, res_ffts.ev(:,i),'LineWidth',1.2); hold on;
    plot(t, res_pid.ev(:,i),'LineWidth',1.2);
    plot(t, res_smc.ev(:,i),'LineWidth',1.2);
    grid on; ylabel(sprintf('e_%d',i));
    if i==1, title('Error Quaternion Vector Part e(t)'); end
    if i==3, xlabel('Time (s)'); end
end
legend('FFTSMC+FTDO','PID','SMC','Location','best');

%% ========== 必画图3：角速度 wx,wy,wz ==========
figure('Name','Angular Rate','Color','w');
labels = {'\omega_x','\omega_y','\omega_z'};
for i=1:3
    subplot(3,1,i);
    plot(t, res_ffts.w(:,i),'LineWidth',1.2); hold on;
    plot(t, res_pid.w(:,i),'LineWidth',1.2);
    plot(t, res_smc.w(:,i),'LineWidth',1.2);
    grid on; ylabel([labels{i} ' (rad/s)']);
    if i==1, title('Angular Rate \omega(t)'); end
    if i==3, xlabel('Time (s)'); end
end
legend('FFTSMC+FTDO','PID','SMC','Location','best');

%% ========== 必画图4：控制力矩 tau_x,tau_y,tau_z ==========
figure('Name','Control Torque','Color','w');
tau_labels = {'\tau_x','\tau_y','\tau_z'};
for i=1:3
    subplot(3,1,i);
    plot(t, res_ffts.tau(:,i),'LineWidth',1.2); hold on;
    plot(t, res_pid.tau(:,i),'LineWidth',1.2);
    plot(t, res_smc.tau(:,i),'LineWidth',1.2);
    grid on; ylabel([tau_labels{i} ' (N·m)']);
    if i==1, title('Control Torque \tau(t)'); end
    if i==3, xlabel('Time (s)'); end
end
legend('FFTSMC+FTDO','PID','SMC','Location','best');

%% ========== 打印末端误差（可选） ==========
fprintf('Final theta_e (deg): FFTS=%.4f, PID=%.4f, SMC=%.4f\n', ...
    rad2deg(res_ffts.theta_e(end)), rad2deg(res_pid.theta_e(end)), rad2deg(res_smc.theta_e(end)));

end

%% -------- local axis-angle to quat (main内) --------
function q = axis_angle_to_quat_local(axis, angle)
axis = axis(:);
axis = axis ./ max(norm(axis),1e-12);
q = [cos(angle/2); axis*sin(angle/2)];
q = q ./ max(norm(q),1e-12);
end
