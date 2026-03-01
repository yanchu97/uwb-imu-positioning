clc
clear
close all

%% deepseek1 (improved v4) - Minimal factor-graph style improvements (KEEP EKF STRUCTURE)
% Goals:
%   - Reduce tight/loose EKF error WITHOUT heavily changing baseline algorithm
%   - Apply SIMPLE but EFFECTIVE factor-graph ideas:
%       (A) Per-anchor outlier rejection (switch-like) using MAD (robust statistics)
%       (B) Mild robust weighting (Huber) only for relative outliers, NOT for common-mode state offset
%       (C) Optional 2-iteration IEKF (relinearize once) ~ factor-graph relinearization
%       (D) Loose-coupling: gate bad LS position updates (avoid LS spikes), no smoothing lag
%       (E) Joseph covariance update (numerical stability)
%
% Notes:
%   - This file is self-contained (helper functions at the end)
%   - Saved as UTF-8 to avoid garbled text

%% Load data
load datas2
dataset = datas;
N = length(dataset.imu.time);

MC_GS   = [];
MC_EKF  = [];
MC_JEKF = [];

times = 10;
count = 0;

%% ===== Switches (minimal-change defaults) =====
ENABLE_IEKF_TIGHT        = true;   % 2-iter IEKF for tight UWB update (simple relinearization)
IEKF_ITERS               = 2;      % 2 is usually enough

ENABLE_ROBUST_TIGHT      = true;   % tight: per-anchor outlier rejection + mild Huber
ENABLE_JOSEPH_UPDATE     = true;   % tight/loose: Joseph form covariance update

ENABLE_SMOOTH_LS         = false;  % IMPORTANT: disable LS smoothing to avoid lag error
ENABLE_POS_GATE_LOOSE    = true;   % loose: gate LS spikes


ENABLE_GEOM_R_LOOSE        = false;  % default OFF: can over-trust LS and create spikes   % loose: use geometry-based position covariance from anchors+ranges (factor-graph linearization idea)
SIGMA_POS_MIN              = 0.05;   % loose: minimum position sigma (m) for Rk (prevents over-trust)
SIGMA_POS_MAX              = 0.80;   % loose: maximum position sigma (m) for Rk (prevents under-trust)

ENABLE_COMMON_BIAS_CANCEL  = false;  % default OFF: can increase error if no common bias
ENABLE_COMMON_BIAS_AUTO    = true;   % if ON, cancel common bias ONLY when detected (safe)   % tight: cancel common range bias by removing median residual before robust weighting
ENABLE_RESET_LOOSE         = false;  % default OFF (keep baseline behavior)   % loose: if innovation explodes, softly reset position towards LS (prevents divergence lock-in)
RESET_GATE_M               = 3.0;    % loose reset gate in meters (innovation norm)
RESET_BLEND                = 0.7;    % loose reset blend (0..1), 0.7 means mostly trust LS when reset triggers
% Robust parameters (tune lightly if needed)
MAD_TAU                  = 3.0;    % outlier threshold in MAD scale (3~4 typical)
HUBER_K_REL              = 2.0;    % Huber threshold (relative to MAD scale)
W_MIN                    = 0.00;   % do NOT drop; just downweight (safer)
INFL_MAX                 = 25;     % max variance inflation factor (keep update from becoming too weak)
Z_ABS_MAX                = 12;     % absolute standardized residual hard cap

L_GATE                   = 200;    % keep original global NIS gate
LOOSE_GATE_K             = 3.0;    % tighter gate to suppress LS outlier spikes    % gate LS position update if innovation too large

USE_FAKE_ERROR           = false;  % keep OFF

ENABLE_PLANAR_Z_CONSTRAINT = false;  % set true if motion is planar (z is constant)
PLANAR_Z_VALUE             = 0.0;    % meters

for M_1 = 4:4
    %#ok<NASGU>
    count = count + 1;

    for t = 1:times %#ok<NASGU>

        dt = mean(diff(dataset.imu.time));

        m_div_cntr = 0;
        m_div = 1;

        UWB_LS_MODE = 3;           % baseline setting
        UWB_EKF_UPDATE_MODE = 3;   % baseline setting

        %% Output containers
        out_data = struct();
        out_data.uwb = struct();
        out_data.uwb.time = dataset.uwb.time;
        out_data.imu = struct();
        out_data.imu.time = dataset.imu.time;
        out_data.uwb.anchor = dataset.uwb.anchor;

        %% Filter init
        settings = uwb_imu_example_settings();
        R = diag(ones(dataset.uwb.cnt, 1) * settings.sigma_uwb^2);

        noimal_state       = init_navigation_state(settings);
        noimal_state_loose = init_navigation_state(settings);

        % Initial state from reference (baseline behavior)
        noimal_state(1:3)       = dataset.pos(:, 1);
        noimal_state_loose(1:3) = dataset.pos(:, 1);

        % Bias states (tight / loose)
        du       = zeros(6, 1);
        du_loose = zeros(6, 1);

        [P, Q1, Q2] = init_filter(settings);
        P_loose = P;

        fprintf("%d frames, IMU rate: %.0f Hz, duration: %.2f s", N, 1/dt, N*dt);
        fprintf("UWB anchors: %d", dataset.uwb.cnt);
        fprintf("UWB update rate: ~%.0f Hz", (1/dt)/m_div);
        fprintf("UWB EKF mode: %dD", UWB_EKF_UPDATE_MODE);
        fprintf("UWB LS mode: %dD", UWB_LS_MODE);
        fprintf("Start filtering...");

        out_data.x(1,:)       = noimal_state;
        out_data.delta_u(1,:) = du';
        out_data.diag_P(1,:)  = trace(P);

        % Loose coupling measurement model: z = p + n
        I3 = eye(3);
        H_Loose = [I3 zeros(3,12)];

        % Keep baseline-scale loose measurement noise (do NOT over-trust LS)
        sigma_pos_loose = max(0.3, settings.sigma_uwb);
        R_Loose_base = (sigma_pos_loose^2) * eye(3);

        %% Pre-compute UWB LS positions (baseline)
        fprintf("Computing UWB LS position...");
        j = 1;
        uwb_pos = [1 1 1]';
        Nuwb = length(dataset.uwb.time);

        for i = 1:Nuwb
            pr = dataset.uwb.tof(:, i);
            uwb_pos = ch_multilateration(dataset.uwb.anchor, uwb_pos, pr', UWB_LS_MODE);

            if ENABLE_SMOOTH_LS
                if i == 1
                    uwb_pos_f = uwb_pos;
                else
                    ALPHA_LS = 0.2; % if enabled, keep it small to avoid lag
                    uwb_pos_f = ALPHA_LS * uwb_pos_f + (1-ALPHA_LS) * uwb_pos;
                end
                uwb_pos = uwb_pos_f;
            end

            out_data.uwb.pos(:, j) = uwb_pos;
            j = j + 1;
        end

        %% Main loop
        x_loose = zeros(N, 10);

        % pre-allocate error arrays (avoid carry-over when times>1)
        E_tightlycoupled = nan(1,N);
        E_loose          = nan(1,N);
        E_ls             = nan(1,N);

        for k = 2:N

            %% ------------ IMU propagation (tight/loose) ------------
            acc_t = dataset.imu.acc(:,k) - du(1:3);
            gyr_t = dataset.imu.gyr(:,k) - du(4:6);

            acc_l = dataset.imu.acc(:,k) - du_loose(1:3);
            gyr_l = dataset.imu.gyr(:,k) - du_loose(4:6);

            % State unpack
            p = noimal_state(1:3);
            v = noimal_state(4:6);
            q = noimal_state(7:10);

            p_loose = noimal_state_loose(1:3);
            v_loose = noimal_state_loose(4:6);
            q_loose = noimal_state_loose(7:10);

            % Propagate
            [p, v, q] = ch_nav_equ_local_tan(p, v, q, acc_t, gyr_t, dt, [0, 0, -9.8]');
            [p_loose, v_loose, q_loose] = ch_nav_equ_local_tan(p_loose, v_loose, q_loose, acc_l, gyr_l, dt, [0, 0, -9.8]');

            % Constrain vertical velocity (baseline)
            v(3) = 0;
            v_loose(3) = 0;
            if ENABLE_PLANAR_Z_CONSTRAINT
                p(3) = PLANAR_Z_VALUE;
                p_loose(3) = PLANAR_Z_VALUE;
            end
            % Pack back
            noimal_state(1:3) = p;
            noimal_state(4:6) = v;
            noimal_state(7:10) = q;

            noimal_state_loose(1:3) = p_loose;
            noimal_state_loose(4:6) = v_loose;
            noimal_state_loose(7:10) = q_loose;

            out_data.eul(k,:) = ch_q2eul(q);

            %% ------------ Covariance propagation ------------
            [F, G] = state_space_model(noimal_state, acc_t, dt);
            [F_loose, G_loose] = state_space_model(noimal_state_loose, acc_l, dt);

            P       = F * P * F' + G * blkdiag(Q1, Q2) * G';
            P_loose = F_loose * P_loose * F_loose' + G_loose * blkdiag(Q1, Q2) * G_loose';

            %% ------------ UWB updates ------------
            m_div_cntr = m_div_cntr + 1;
            if m_div_cntr == m_div
                m_div_cntr = 0;

                pr   = dataset.uwb.tof(1:dataset.uwb.cnt, k);
                anch = dataset.uwb.anchor;
                R1   = R;

                % ===== Tight coupling: UWB range EKF update =====
                if ENABLE_IEKF_TIGHT
                    [noimal_state, du, P, Lval] = tight_iekf_update( ...
                        noimal_state, du, P, pr, anch, R1, UWB_EKF_UPDATE_MODE, ...
                        ENABLE_ROBUST_TIGHT, MAD_TAU, HUBER_K_REL, W_MIN, INFL_MAX, Z_ABS_MAX, ...
                        IEKF_ITERS, ENABLE_COMMON_BIAS_CANCEL, ENABLE_COMMON_BIAS_AUTO, ENABLE_JOSEPH_UPDATE, L_GATE);
                    out_data.L(k,:) = Lval;
                else
                    [noimal_state, du, P, Lval] = tight_ekf_update_once( ...
                        noimal_state, du, P, pr, anch, R1, UWB_EKF_UPDATE_MODE, ...
                        ENABLE_ROBUST_TIGHT, MAD_TAU, HUBER_K_REL, W_MIN, INFL_MAX, Z_ABS_MAX, ...
                        ENABLE_COMMON_BIAS_CANCEL, ENABLE_COMMON_BIAS_AUTO, ENABLE_JOSEPH_UPDATE, L_GATE);
                    out_data.L(k,:) = Lval;
                end

                % ===== Loose coupling: LS position EKF update =====

                z = out_data.uwb.pos(:, k);

                V_Loose = z - noimal_state_loose(1:3);


                Rk = R_Loose_base;

                % Geometry-based covariance (factor-graph linearization idea):
                %   For ranges r_i with sigma_r, the position covariance can be approximated by:
                %     Cov_p approx (J' * (1/sigma_r^2) * J)^(-1),  J_i = (p-a_i)'/||p-a_i||
                if ENABLE_GEOM_R_LOOSE
                    Rk = uwb_geo_pos_cov(anch, z, settings.sigma_uwb, SIGMA_POS_MIN, SIGMA_POS_MAX, R_Loose_base);
                end
                % Robustify LS position factor (factor-graph style): do NOT skip, just inflate R when innovation is huge

                if ENABLE_POS_GATE_LOOSE
                    % NOTE: In EKF, if we inflate R too much when innovation is large,
                    % the filter may "lock-in" and fail to recover. So we:
                    %   (1) prefer a soft reset towards LS when innovation is huge (optional),
                    %   (2) otherwise only mildly inflate R.

                    dpos = norm(V_Loose);
                    if isfinite(dpos)

                        % soft reset (prevents divergence lock-in)
                        if ENABLE_RESET_LOOSE && dpos > RESET_GATE_M
                            noimal_state_loose(1:3) = (1-RESET_BLEND) * noimal_state_loose(1:3) + RESET_BLEND * z;
                            % increase uncertainty a bit so measurement keeps pulling
                            P_loose(1:3,1:3) = P_loose(1:3,1:3) + (RESET_BLEND^2) * Rk;
                        end

                        % mild inflation only (keep update effective)
                        sigpos = sqrt(max(1e-12, trace(Rk)/3));
                        gate = LOOSE_GATE_K * sigpos;
                        if dpos > gate
                            infl = min(100, (dpos/gate)^2); % allow stronger inflation to suppress LS outliers
                            Rk = Rk * infl;
                        end
                    end
                end

                S_Loose = H_Loose * P_loose * H_Loose' + Rk;

                K_Loose = (P_loose * H_Loose') / S_Loose;

                deltaL  = K_Loose * V_Loose;


                noimal_state_loose(1:6) = noimal_state_loose(1:6) + deltaL(1:6);


                ql = noimal_state_loose(7:10);

                ql = ch_qmul(ch_rv2q(deltaL(7:9)), ql);

                ql = ql ./ max(1e-12, norm(ql));

                noimal_state_loose(7:10) = ql;


                du_loose = du_loose + deltaL(10:15);


                if ENABLE_JOSEPH_UPDATE

                    I15 = eye(15);

                    P_loose = (I15 - K_Loose*H_Loose) * P_loose * (I15 - K_Loose*H_Loose)' + K_Loose * Rk * K_Loose';

                else

                    P_loose = (eye(15) - K_Loose*H_Loose) * P_loose;

                end
            end

            % log
            out_data.x(k,:)       = noimal_state;
            x_loose(k,:)          = noimal_state_loose;
            out_data.delta_u(k,:) = du';
            out_data.diag_P(k,:)  = trace(P);

            % 2D errors (XY)
            E_tightlycoupled(1,k) = hypot(dataset.pos(1,k) - out_data.x(k,1), dataset.pos(2,k) - out_data.x(k,2)); %#ok<AGROW>
            E_loose(1,k)          = hypot(dataset.pos(1,k) - x_loose(k,1),    dataset.pos(2,k) - x_loose(k,2));    %#ok<AGROW>
            E_ls(1,k)             = hypot(dataset.pos(1,k) - out_data.uwb.pos(1,k), dataset.pos(2,k) - out_data.uwb.pos(2,k)); %#ok<AGROW>

            if USE_FAKE_ERROR
                % keep OFF
            end
        end

        Mean_E_tightlycoupled(1,t) = mean(E_tightlycoupled, 'omitnan'); %#ok<AGROW>
        Mean_E_loose(1,t)          = mean(E_loose, 'omitnan');          %#ok<AGROW>
    end

    MC_GS   = [MC_GS,   E_ls];   %#ok<AGROW>
    MC_EKF  = [MC_EKF,  E_loose]; %#ok<AGROW>
    MC_JEKF = [MC_JEKF, E_tightlycoupled]; %#ok<AGROW>

    RMSE_Mean_E_tightlycoupled(count) = mean(Mean_E_tightlycoupled, 'omitnan'); %#ok<AGROW>
end

%% CDF plot
figure;
cdfplot(MC_GS); hold on
cdfplot(MC_EKF); hold on
cdfplot(MC_JEKF); hold on
legend("UWB (LS)", "Loose EKF", "Tight EKF", "Location", "southeast");
xlabel('2D error / m'); ylabel('CDF');
axis([0,5,0,1]);

%% Error time series
figure;
plot(E_ls,'k'); hold on;
plot(E_loose,'b'); hold on;
plot(E_tightlycoupled,'r'); hold on;
xlabel('Frame'); ylabel('2D error / m');
legend("UWB (LS)", "Loose EKF", "Tight EKF");

%% Trajectory plot
out_data.uwb.fusion_pos = out_data.x(:,1:3)';

figure;
plot(dataset.pos(1,:), dataset.pos(2,:), 'm--'); hold on;
plot(out_data.uwb.pos(1,:), out_data.uwb.pos(2,:), 'k.');
plot(x_loose(:,1), x_loose(:,2), 'b');
plot(out_data.uwb.fusion_pos(1,:), out_data.uwb.fusion_pos(2,:), 'r');

anch = out_data.uwb.anchor;
scatter(anch(1, :),anch(2, :),'filled');
for i=1:size(anch,2)
    text(anch(1, i),anch(2, i), "A" + (i-1));
end
legend("Reference", "UWB (LS)", "Loose EKF", "Tight EKF", "Anchors", "Location", "best");
xlabel('X / m'); ylabel('Y / m');

%% ===================== Helper functions =====================

function Rk = uwb_geo_pos_cov(anch, p, sigma_r, sigma_pos_min, sigma_pos_max, R_fallback)
% Geometry-based position covariance from range measurement linearization.
% anch: 3xM, p:3x1 (LS position), sigma_r: scalar range std.
% Returns diagonal covariance (3x3) with floors/caps to remain well-conditioned.

Rk = R_fallback;

try
    if any(~isfinite(p)) || isempty(anch)
        return;
    end

    M = size(anch,2);
    if M < 3
        return;
    end

    J = zeros(M,3);
    for i = 1:M
        d = (p - anch(:,i));
        dn = norm(d);
        if dn < 1e-6
            return;
        end
        J(i,:) = (d./dn)';
    end

    W = (1/max(1e-6, sigma_r)^2) * eye(M);
    H = J' * W * J;

    % regularize
    H = H + 1e-6 * eye(3);

    Cp = inv(H);

    % Use only diagonal (simple, stable)
    s2 = diag(Cp);
    s2 = max((sigma_pos_min^2)*ones(3,1), min((sigma_pos_max^2)*ones(3,1), s2));
    Rk = diag(s2);

catch
    Rk = R_fallback;
end
end


function x = init_navigation_state(~)
q = ch_eul2q(deg2rad([0 0 0]));
x = [zeros(6,1); q];
end

function [P, Q1, Q2] = init_filter(settings)
P = zeros(15);
P(1:3,1:3)     = settings.factp(1)^2 * eye(3);
P(4:6,4:6)     = settings.factp(2)^2 * eye(3);
P(7:9,7:9)     = diag(settings.factp(3:5)).^2;
P(10:12,10:12) = settings.factp(6)^2 * eye(3);
P(13:15,13:15) = settings.factp(7)^2 * eye(3);

Q1 = zeros(6);
Q1(1:3,1:3) = (settings.sigma_acc^2) * eye(3);
Q1(4:6,4:6) = (settings.sigma_gyro^2) * eye(3);

Q2 = zeros(6);
Q2(1:3,1:3) = settings.sigma_acc_bias^2  * eye(3);
Q2(4:6,4:6) = settings.sigma_gyro_bias^2 * eye(3);
end

function [F,G] = state_space_model(x, acc, dt)
Cb2n = ch_q2m(x(7:10));
sk = ch_askew(Cb2n * acc);

O = zeros(3);
I = eye(3);
F = [
    O I   O     O      O;
    O O  -sk  -Cb2n    O;
    O O   O     O   -Cb2n;
    O O   O     O      O;
    O O   O     O      O];

F = eye(15) + dt * F;

G = dt * [
    O     O     O   O;
    Cb2n  O     O   O;
    O   -Cb2n   O   O;
    O     O     I   O;
    O     O     O   I];
end

function [Y, H] = uwb_hx(x, anchor, dim)
N = size(anchor,2);
position = x(1:3);
if dim == 2
    position = position(1:2);
    anchor = anchor(1:2, 1:N);
end

Y = zeros(N,1);
H = zeros(N,15);

perd_pr = repmat(position,1,N) - anchor(:,1:N);
for i = 1:N
    d = norm(perd_pr(:,i));
    if d < 1e-9
        d = 1e-9;
    end
    if dim == 2
        H(i,:) = [ (perd_pr(:,i)'/d) zeros(1,13) ];
    else
        H(i,:) = [ (perd_pr(:,i)'/d) zeros(1,12) ];
    end
    Y(i) = d;
end
end

function [x, du, P, Lval] = tight_iekf_update(x, du, P, pr, anch, R1, dim, ...
    enable_robust, mad_tau, huber_k_rel, w_min, infl_max, z_abs_max, iters, enable_common_bias_cancel, enable_common_bias_auto, enable_joseph, L_gate)
% Simple IEKF: relinearize UWB range update a few times
Lval = inf;
for it = 1:iters
    [x, du, P, Lval] = tight_ekf_update_once(x, du, P, pr, anch, R1, dim, ...
        enable_robust, mad_tau, huber_k_rel, w_min, infl_max, z_abs_max, enable_common_bias_cancel, enable_common_bias_auto, enable_joseph, L_gate);
    if ~isfinite(Lval) || Lval > L_gate
        break;
    end
end
end

function [x, du, P, Lval] = tight_ekf_update_once(x, du, P, pr, anch, R1, dim, ...
    enable_robust, mad_tau, huber_k_rel, w_min, infl_max, z_abs_max, enable_common_bias_cancel, enable_common_bias_auto, enable_joseph, L_gate)

% Linearize
[Y, H] = uwb_hx(x, anch, dim);
residual = pr - Y;

% Cancel common (global) range bias (simple nuisance elimination, factor-graph style):
% Many datasets have a shared positive bias (e.g., clock/NLOS). Removing the median
% keeps relative geometry constraints while avoiding "all residuals large" lock-in.
if enable_common_bias_cancel
    rr = residual(isfinite(residual));
    if numel(rr) >= 3
        b0 = median(rr);
        if enable_common_bias_auto
            % Apply only when residuals show a CONSISTENT common-mode offset:
            %   (1) most residuals share the same sign
            %   (2) dispersion is small compared to the offset magnitude
            sgn = sign(rr);
            sameSign = (mean(sgn > 0) > 0.75) || (mean(sgn < 0) > 0.75);
            dispStd  = 1.4826 * median(abs(rr - b0)) + 1e-12; % MAD -> std
            if sameSign && abs(b0) > 2.0*dispStd
                residual = residual - b0;
            end
        else
            residual = residual - b0;
        end
    end
end

H_eff = H;
res_eff = residual;
R_eff = R1;

if enable_robust
    [H_eff, res_eff, R_eff] = uwb_build_mad_robust_system(H, residual, R1, mad_tau, huber_k_rel, w_min, infl_max, z_abs_max);
end

if isempty(res_eff)
    Lval = inf;
    return;
end

S = H_eff * P * H_eff' + R_eff;
Lval = (res_eff' / S * res_eff);

if Lval > L_gate
    return;
end

K = (P * H_eff') / S;
delta = K * res_eff;

% position/velocity
x(1:6) = x(1:6) + delta(1:6);

% attitude (multiplicative)
q = x(7:10);
q = ch_qmul(ch_rv2q(delta(7:9)), q);
q = q ./ max(1e-12, norm(q));
x(7:10) = q;

% biases
if numel(delta) >= 15
    du = du + delta(10:15);
end

% covariance update
if enable_joseph
    I15 = eye(15);
    P = (I15 - K*H_eff) * P * (I15 - K*H_eff)' + K * R_eff * K';
else
    P = (eye(15) - K*H_eff) * P;
end
end

function [H_eff, res_eff, R_eff] = uwb_build_mad_robust_system(H, residual, R1, mad_tau, huber_k_rel, w_min, infl_max, z_abs_max)
% Robust per-anchor weighting using MAD on standardized residuals.
% This is a SIMPLE factor-graph idea: downweight only anchors that deviate
% from the majority (relative outliers), instead of suppressing all measurements
% when the state has a common-mode offset.
%
% Output system keeps (almost) all finite measurements:
%   - Inliers: weight ~ 1
%   - Relative outliers: weight floored to 1/infl_max (large R)

sigma = sqrt(max(1e-12, diag(R1)));
sigma = sigma(:);
r = residual(:);

finite = isfinite(r) & isfinite(sigma) & (sigma > 0);
if ~any(finite)
    H_eff = []; res_eff = []; R_eff = [];
    return;
end

H = H(finite,:);
r = r(finite);
sigma = sigma(finite);

z = r ./ sigma;                  % standardized residuals
medz = median(z);
dz = abs(z - medz);

mad = median(dz) + 1e-12;
scale = 1.4826 * mad;            % robust std estimate in z-space
scale = max(scale, 1e-3);

% Relative inlier mask (only used to decide if anchor is an outlier)
inlier = (dz <= mad_tau * scale) & (abs(z) <= z_abs_max);

% Mild Huber weights on relative deviation (only meaningful for inliers)
a = huber_k_rel * scale;
w = ones(size(z));
idx = (dz > a) & inlier;
w(idx) = a ./ max(dz(idx), 1e-12);

% Outliers: do NOT drop; just downweight strongly
w_floor = 1 / max(1, infl_max);
w(~inlier) = w_floor;

% Optional switch-like drop (default w_min=0: keep all)
if w_min > 0
    keep = w >= w_min;
    H = H(keep,:);
    r = r(keep);
    sigma = sigma(keep);
    w = w(keep);
end

if isempty(r)
    H_eff = []; res_eff = []; R_eff = [];
    return;
end

% Convert weights to variance inflation, cap by infl_max
w = max(w, w_floor);
R_eff = diag((sigma.^2) ./ w);

H_eff = H;
res_eff = r;
end
