function out = attitude_sim(ctrlName, sim)
% out = attitude_sim(ctrlName, sim)
% ctrlName: 'PID' | 'SMC' | 'FFTSMC_FTDO'
% sim: struct with fields
%   .Tend, .dt, .J(3x3), .qd(4x1), .q0(4x1), .w0(3x1)
%   .pid, .smc, .ffts, .ftdo (params)
%   .tau_limit (optional, scalar or 3x1), default inf
%
% out fields:
%   .t, .q, .w, .tau, .ev, .theta_e (rad)
%   .dhat (Nx3, only for FFTSMC_FTDO)

arguments
    ctrlName (1,:) char
    sim struct
end

% defaults
if ~isfield(sim,'tau_limit'); sim.tau_limit = inf; end

Tend = sim.Tend; dt = sim.dt;
t_grid = (0:dt:Tend)';

% build initial state
q0 = quat_normalize(sim.q0(:));
w0 = sim.w0(:);

switch upper(ctrlName)
    case 'PID'
        x0 = [q0; w0; zeros(3,1)];          % [q; w; int_e]
    case 'SMC'
        x0 = [q0; w0];                      % [q; w]
    case 'FFTSMC_FTDO'
        what0 = w0;                         % observer init
        dhat0 = zeros(3,1);
        x0 = [q0; w0; what0; dhat0];        % [q; w; what; dhat]
    otherwise
        error('Unknown ctrlName: %s', ctrlName);
end

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

[t, X] = ode45(@(t,x) dyn_all(t,x,sim,ctrlName), [0 Tend], x0, opts);

% interpolate to uniform grid for plotting
Xg = interp1(t, X, t_grid, 'pchip');

% parse
q = normalize_rows_quat(Xg(:,1:4));
w = Xg(:,5:7);

% compute error quaternion vector part and angle error
[qe, ev, theta_e] = quat_error_series(q, sim.qd(:));

% compute tau on grid (post-process via same control law)
N = numel(t_grid);
tau = zeros(N,3);
dhat_hist = [];

for k = 1:N
    xk = Xg(k,:)';
    [tau_k, dbg] = control_law(ctrlName, xk, sim);
    tau(k,:) = tau_k(:).';
    if strcmpi(ctrlName,'FFTSMC_FTDO')
        if isempty(dhat_hist); dhat_hist = zeros(N,3); end
        dhat_hist(k,:) = dbg.dhat(:).';
    end
end

out.t = t_grid;
out.q = q;
out.w = w;
out.tau = tau;
out.ev = ev;
out.theta_e = theta_e;

if strcmpi(ctrlName,'FFTSMC_FTDO')
    out.dhat = dhat_hist;
end

end

%% ===================== Dynamics Wrapper =====================
function dx = dyn_all(~, x, sim, ctrlName)
J = sim.J;
q = quat_normalize(x(1:4));
w = x(5:7);

% disturbance for this sim: d = 0
d = [0;0;0];

% compute control torque
[tau, dbg] = control_law(ctrlName, x, sim);

% rigid body rotational dynamics
w_dot = J \ (tau - cross(w, J*w) + d);

% quaternion kinematics
q_dot = quat_kinematics(q, w);

% extra states
switch upper(ctrlName)
    case 'PID'
        int_e_dot = dbg.int_e_dot;
        dx = [q_dot; w_dot; int_e_dot];

    case 'SMC'
        dx = [q_dot; w_dot];

    case 'FFTSMC_FTDO'
        % observer states exist: what(3), dhat(3)
        what = x(8:10);
        dhat = x(11:13);

        % observer: w_dot = f(w) + J^{-1}tau + J^{-1}d
        % f(w) = -J^{-1} (w x (Jw))
        f = -(J \ cross(w, J*w));
        Jinv = inv(J);

        L  = sim.ftdo.L;
        Ko1 = sim.ftdo.Ko1;
        Ko2 = sim.ftdo.Ko2;
        alpha = sim.ftdo.alpha;

        e_obs = (w - what);

        what_dot = f + Jinv*tau + Jinv*dhat + L*e_obs;
        dhat_dot = Ko1*sig_pow(e_obs, alpha) + Ko2*e_obs;

        dx = [q_dot; w_dot; what_dot; dhat_dot];

    otherwise
        error('Unknown ctrlName in dyn_all.');
end

end

%% ===================== Control Laws =====================
function [tau, dbg] = control_law(ctrlName, x, sim)
% returns torque tau(3x1) and dbg fields for extra states
J = sim.J;
qd = sim.qd(:);
tau_limit = sim.tau_limit;

q = quat_normalize(x(1:4));
w = x(5:7);

% error quaternion qe = qd^{-1} ⊗ q
qe = quat_mul(quat_conj(qd), q);
qe = shortest_quat(qe);
e  = qe(2:4);
qe0 = qe(1);

% error vector derivative e_dot = 0.5*(qe0 I + skew(e))*w
A = 0.5*(qe0*eye(3) + skew3(e));
e_dot = A*w;

dbg = struct();

switch upper(ctrlName)
    case 'PID'
        int_e = x(8:10);

        Kp = sim.pid.Kp; Kd = sim.pid.Kd; Ki = sim.pid.Ki;
        int_lim = sim.pid.int_limit;

        tau = -Kp*e - Kd*w - Ki*int_e;

        % integral state derivative with anti-windup clamp
        dbg.int_e_dot = clamp_vec(e, int_lim);

    case 'SMC'
        c   = sim.smc.c;
        K   = sim.smc.K;
        eta = sim.smc.eta;
        phi = sim.smc.phi;

        s = w + c*e;
        sat_s = sat_vec(s, phi); % sat(s/phi) in [-1,1]

        % w_dot_cmd = -c*e_dot - K*s - eta*sat(s/phi)
        w_dot_cmd = -c*e_dot - K*s - eta*sat_s;

        tau = cross(w, J*w) + J*w_dot_cmd;

    case 'FFTSMC_FTDO'
        % get observer estimate dhat from state
        dhat = x(11:13);
        dbg.dhat = dhat;

        p = sim.ffts.p; qn = sim.ffts.q;
        r1 = p/qn;                 % <1
        r2 = qn/p;                 % >1
        lambda1 = sim.ffts.lambda1;
        lambda2 = sim.ffts.lambda2;

        alpha = sim.ffts.alpha;    % 0<alpha<1
        beta  = sim.ffts.beta;     % beta>1
        K1 = sim.ffts.K1;
        K2 = sim.ffts.K2;
        K3 = sim.ffts.K3;

        % s = e_dot + λ1*sig(e)^(r1) + λ2*sig(e)^(r2)
        f1 = lambda1 * sig_pow(e, r1);
        f2 = lambda2 * sig_pow(e, r2);
        s  = e_dot + f1 + f2;

        % df/de * e_dot  (elementwise approximation)
        eps0 = 1e-6;
        dfde1 = lambda1 * r1 * (abs(e)+eps0).^(r1-1);
        dfde2 = lambda2 * r2 * (abs(e)+eps0).^(r2-1);
        dfde_e_dot = (dfde1 + dfde2) .* e_dot;

        % reaching law (dual power) + linear damping
        s_reach = K1*sig_pow(s, alpha) + K2*sig_pow(s, beta) + K3*s;

        % command e_ddot
        e_ddot_cmd = -dfde_e_dot - s_reach;

        % robust mapping: e_ddot = A*w_dot  => w_dot_cmd = DLS(A)*e_ddot_cmd
        mu = 1e-4;
        w_dot_cmd = (A.'*A + mu*eye(3)) \ (A.'*e_ddot_cmd);

        % torque with disturbance compensation
        tau = cross(w, J*w) + J*w_dot_cmd - dhat;

    otherwise
        error('Unknown ctrlName in control_law.');
end

% torque limit (optional)
tau = clamp_vec(tau, tau_limit);

end

%% ===================== Quaternion / Math Utils =====================
function qdot = quat_kinematics(q, w)
q0 = q(1); qv = q(2:4);
qdot0 = -0.5 * (qv.' * w);
qdotv =  0.5 * (q0*eye(3) + skew3(qv)) * w;
qdot  = [qdot0; qdotv];
end

function S = skew3(v)
S = [  0   -v(3)  v(2);
      v(3)   0   -v(1);
     -v(2)  v(1)   0  ];
end

function qn = quat_normalize(q)
qn = q ./ max(norm(q), 1e-12);
end

function Qn = normalize_rows_quat(Q)
Qn = Q;
for i=1:size(Q,1)
    Qn(i,:) = Q(i,:)./max(norm(Q(i,:)),1e-12);
end
end

function qcon = quat_conj(q)
qcon = [q(1); -q(2:4)];
end

function qc = quat_mul(a, b)
a0=a(1); av=a(2:4);
b0=b(1); bv=b(2:4);
c0 = a0*b0 - dot(av,bv);
cv = a0*bv + b0*av + cross(av,bv);
qc = [c0; cv];
end

function qe = shortest_quat(qe)
if qe(1) < 0
    qe = -qe;
end
end

function y = sig_pow(x, r)
y = (abs(x)).^r .* sign(x);
end

function y = sat_vec(x, phi)
% returns sat(x/phi) elementwise
z = x ./ max(phi, 1e-12);
y = max(-1, min(1, z));
end

function y = clamp_vec(x, lim)
% lim can be scalar or 3x1
if isscalar(lim)
    y = min(max(x, -lim), lim);
else
    lim = lim(:);
    y = min(max(x, -lim), lim);
end
end

function q = axis_angle_to_quat(axis, angle)
axis = axis(:);
axis = axis ./ max(norm(axis),1e-12);
q = [cos(angle/2); axis*sin(angle/2)];
q = quat_normalize(q);
end

function [QE, EV, TH] = quat_error_series(Q, qd)
% Q: Nx4 rows [q0 q1 q2 q3], qd: 4x1
N = size(Q,1);
QE = zeros(N,4);
EV = zeros(N,3);
TH = zeros(N,1);

qd_conj = quat_conj(qd);
for i=1:N
    q = Q(i,:)';
    qe = quat_mul(qd_conj, q);
    qe = shortest_quat(qe);
    QE(i,:) = qe.';
    EV(i,:) = qe(2:4).';

    q0 = max(-1,min(1,qe(1)));
    qv = qe(2:4);
    TH(i) = 2*atan2(norm(qv), q0); % robust angle error [0,pi]
end
end
