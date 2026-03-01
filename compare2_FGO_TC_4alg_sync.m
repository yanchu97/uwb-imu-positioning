%% IMM-PDAF based INS/UWB Integration
% simulation tool: navego
% by lizhengjian
% 3.2.2020

%%
clc
close all
clear
matlabrc

addpath FUNCTIONS/
addpath FUNCTIONS/simulation
addpath FUNCTIONS/conversions
addpath OTHER/
% freq_imu = 100; freq_UWB = 5;

%% CONVERSION CONSTANTS

G =  9.80665;       % Gravity constant, m/s^2
G2MSS = G;          % g to m/s^2
MSS2G = (1/G);      % m/s^2 to g

D2R = (pi/180);     % degrees to radians
R2D = (180/pi);     % radians to degrees

KT2MS = 0.514444;   % knot to m/s
MS2KMH = 3.6;       % m/s to km/h

%% LOAD REFERENCE DATA
fprintf('NaveGo: loading reference dataset from a trajectory generator... \n')
load ref2.mat
% ref=load('shiyan_1.mat'); %新加实验
%% ADIS16405 IMU error profile
ADIS16488.arw      = 0.125  .* ones(1,3);     % Angle random walks [X Y Z] (deg/root-hour)
ADIS16488.arrw     = zeros(1,3);            % Angle rate random walks [X Y Z] (deg/root-hour/s)
ADIS16488.vrw      = 0.02.* ones(1,3);     % Velocity random walks [X Y Z] (m/s/root-hour)
ADIS16488.vrrw     = zeros(1,3);            % Velocity rate random walks [X Y Z] (deg/root-hour/s)
ADIS16488.gb_sta   = 0.1  .* ones(1,3);     % Gyro static biases [X Y Z] (deg/s)
ADIS16488.ab_sta   = 6   .* ones(1,3);     % Acc static biases [X Y Z] (mg)
ADIS16488.gb_dyn   = 4.5/3600  .* ones(1,3);% Gyro dynamic biases [X Y Z] (deg/s)
ADIS16488.ab_dyn   = 0.07  .* ones(1,3);     % Acc dynamic biases [X Y Z] (mg)
ADIS16488.gb_corr  = 100*ones(1,3);     % Gyro correlation times [X Y Z] (seconds)
ADIS16488.ab_corr  = 100*ones(1,3);     % Acc correlation times [X Y Z] (seconds)
ADIS16488.freq     = ref.freq;              % IMU operation frequency [X Y Z] (Hz)
% ref dataset will be used to simulate IMU sensors.
ADIS16488.t = ref.t;                        % IMU time vector
dt = mean(diff(ADIS16488.t));               % IMU sampling interval

imu2 = imu_si_errors(ADIS16488, dt);        % Transform IMU manufacturer error units to SI units.

imu2.ini_align_err = [0 0 0] .* D2R;                     % Initial attitude align errors for matrix P in Kalman filter, [roll pitch yaw] (radians)
imu2.ini_align = [ref.roll(1) ref.pitch(1) ref.yaw(1)];  % Initial attitude align at t(1) (radians).

%% MONTECARLO SIMULATION
rmse_EKF = zeros(1,8);
rmse_IMMPDAF = zeros(1,8);
rmse_PRPF = zeros(1,8);
rmse_SHFAEKF = zeros(1,8);
rmse_FGO = zeros(1,8);
% num_anchor = 5
% p_nlos = 0.4 % 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8
% u_nlos = 3 % 1,2,3,4,5,6
% for u_nlos = 1: 1 :8
%     u_nlos
tic
%采样次数
times = 1;
T = times;
T1 = times;
T2= times;
T_FGO = times;
EKF = 'ON';
SHFAEKF = 'ON';
PRPF = 'ON';
IMMPDAF = 'OFF';
FGO = 'ON';%%%%较IMMPDAF新增?
GDOP = 'ON';
if strcmp( IMMPDAF , 'ON' )||strcmp( FGO , 'ON' )%%%%较IMMPDAF新增?
    fq.freq1 = 50;
end
if strcmp( EKF , 'ON' ) || strcmp( SHFAEKF , 'ON' ) || strcmp( PRPF , 'ON' )
    fq.freq2 = 5;
end
uwbmode = 'random';
num_particle = 100; % number of particles (PF)
num_anchor = 5;
p_nlos = 0.4; % proportion of nlos
std_los = 0.5; % standard deviation of los
u_nlos = 1; % shape parameter of nlos gamma distribution
std_nlos = 1; % scale parameter of nlos gamma distribution
std_mode1 = 1; % R1 of imm model1 (IMMPDAF)
std_mode2 = 2; % R2 of imm model2 (IMMPDAF)
gamma1 = 16; % gate value of model1 (IMMPDAF)
gamma2 = 25; % gate value of model2 (IMMPDAF)
LW = 85; % 滑动窗口长度 (SHFAKF)
kxi = 3; % 异常值检测门�??? (SHFAKF)
alpha = 4; % scale factor of (SHFAKF)
zzz_flag=zeros(times,1); % algorithm divergence flag (IMMPDAF)
alpha_gam = (u_nlos^2)/(std_nlos^2);
betta_gam = (std_nlos^2)/u_nlos;
for mc = 1:times
    fprintf('%dth montecarlo run. \n', mc)

    %% 10Hz UWB error profile
    uwb.std = [std_los, alpha_gam, betta_gam];                  % UWB positions standard deviations [lat lon h] (meters)
    uwb.stdm = [std_mode1 * ones(1,3);
        std_mode2 * ones(1,3)];
    uwb.stdv = [0.3 0.3 0.3];   % UWV velocities standard deviations [Vn Ve Vd] (meters/s)
    % Parameters for ZUPT detection algorithm
    uwb.zupt_th = 0.05;   % ZUPT threshold (m/s).
    uwb.zupt_win = 4;    % ZUPT time window (seconds).
    uwb.eps = 1E-3;
    %% UWB SYNTHETIC DATA
    rng('shuffle')                  % Reset pseudo-random seed
    if isfield(fq,'freq1')
        uwb.freq = fq.freq1;                          % UWB operation frequency (Hz)
        fprintf('NaveGo: generating %dHz UWB synthetic data... \n', fq.freq1)
        uwb = uwb_gen (ref, uwb,5, p_nlos, uwbmode,'em');
        save uwb50.mat uwb
    end
    if isfield(fq,'freq2')
        uwb.freq = fq.freq2;
        fprintf('NaveGo: generating %dHz UWB synthetic data... \n', fq.freq2)
        uwb2 = uwb_gen (ref, uwb,5, p_nlos, uwbmode,'em');
        save uwb5.mat uwb2
    end
    if isfield(fq,'freq3')
        uwb.freq = fq.freq3;
        fprintf('NaveGo: generating %dHz UWB synthetic data... \n', fq.freq3)
        uwb3 = uwb_gen (ref, uwb,5, p_nlos, uwbmode,'em');
        save uwb50.mat uwb3
    end
% ---- Common evaluation time grid (sync all methods on the same timestamps) ----
% Prefer the 5Hz UWB grid (used by EKF/SHFAEKF/PRPF). If not available, fall back to UWB grid.
if exist('uwb2','var') && isfield(uwb2,'t')
    t_eval = uwb2.t(:);
else
    t_eval = uwb.t(:);
end
L = numel(t_eval);
% Reference on evaluation grid
x_ref_eval = interp1(ref.t, ref.x, t_eval, 'linear', 'extrap');
y_ref_eval = interp1(ref.t, ref.y, t_eval, 'linear', 'extrap');
z_ref_eval = interp1(ref.t, ref.z, t_eval, 'linear', 'extrap');

rmse11=zeros(L,1);
rmse12=zeros(L,1);
rmse2=zeros(L,1);
rmse3=zeros(L,1);
rmse4=zeros(L,1);
rmse5=zeros(L,1);
rmset11=zeros(L,1);
rmset12=zeros(L,1);
rmset2=zeros(L,1);
rmset3=zeros(L,1);
rmset4=zeros(L,1);
rmset5=zeros(L,1);
    %% IMU2 SYNTHETIC DATA

    rng('shuffle')					% Reset pseudo-random seed

    fprintf('NaveGo: generating IMU2 synthetic data... \n')

    fb = acc_gen (ref, imu2);   % Generate acc in the body frame
    imu2.fb = fb;

    wb = gyro_gen (ref, imu2);  % Generate gyro in the body frame
    imu2.wb = wb;

    save imu2.mat imu2

    %% Print navigation time

    to = (ref.t(end) - ref.t(1));

    fprintf('NaveGo: navigation time is %.2f minutes or %.2f seconds. \n', (to/60), to)
    %% INS/UWB integration using IMU2

    if strcmp( EKF , 'ON' )
        fprintf('NaveGo: INS/UWB navigation estimates for IMU2, EKF employed... \n')
        tic
        [nav2_e11] = uwb_ins_ekf(imu2, uwb2, ref);
        toc
        %         err_uiekf = [ ref.x - nav2_e11.x  ref.y - nav2_e11.y  ref.z - nav2_e11.z ];
        %         for i = 1:length(ref.t)
        %             rmse11(i) = norm( err_uiekf(i,:) );
        %         end        % --- rmse11 on common grid t_eval (robust time/shape handling) ---
        [t_nav, x_series, y_series, z_series] = nav_pos_series(nav2_e11, ref, uwb2, imu2);
        x_nav = interp1(t_nav, x_series, t_eval, 'linear', 'extrap');
        y_nav = interp1(t_nav, y_series, t_eval, 'linear', 'extrap');
        z_nav = interp1(t_nav, z_series, t_eval, 'linear', 'extrap');
        err = [x_ref_eval - x_nav, y_ref_eval - y_nav, z_ref_eval - z_nav];
        rmse11 = sqrt(sum(err.^2, 2));
        rmset11 = rmset11 + rmse11;
    end

    if strcmp( SHFAEKF , 'ON' )
       
        fprintf('NaveGo: INS/UWB navigation estimates for IMU2, SHFAKF employed... \n')
        tic
        [nav2_e4] = uwb_ins_shfakf2D(imu2, uwb2, ref, LW, kxi, alpha);
        toc
        %         err_uiakf_shf = [ ref.x - nav2_e4.x  ref.y - nav2_e4.y  ref.z - nav2_e4.z ];
        %         for i = 1:length(ref.t)
        %             rmse4(i) = norm( err_uiakf_shf(i,:) );
        %         end        % --- rmse4 on common grid t_eval (robust time/shape handling) ---
        [t_nav, x_series, y_series, z_series] = nav_pos_series(nav2_e4, ref, uwb2, imu2);
        x_nav = interp1(t_nav, x_series, t_eval, 'linear', 'extrap');
        y_nav = interp1(t_nav, y_series, t_eval, 'linear', 'extrap');
        z_nav = interp1(t_nav, z_series, t_eval, 'linear', 'extrap');
        err = [x_ref_eval - x_nav, y_ref_eval - y_nav, z_ref_eval - z_nav];
        rmse4 = sqrt(sum(err.^2, 2));
        if mean(rmse4) < 5
            rmset4 = rmset4 + rmse4;
        else
            T1 = T1 - 1;
        end

    end
    if strcmp( PRPF , 'ON' )
        fprintf('NaveGo: INS/UWB navigation estimates for IMU2, PRPF employed... \n')
        tic
        [nav2_e3] = uwb_ins_nspf2(imu2, uwb2, 1, num_particle, alpha_gam, betta_gam, ref);
        toc
        %         err_uipf = [ ref.x - nav2_e3.x  ref.y - nav2_e3.y  ref.z - nav2_e3.z ];
        %         for i = 1:length(ref.t)
        %             rmse3(i) = norm( err_uipf(i,:) );
        %         end
        %         rmset3 = rmset3 + rmse3;        % --- rmse3 on common grid t_eval (robust time/shape handling) ---
        [t_nav, x_series, y_series, z_series] = nav_pos_series(nav2_e3, ref, uwb2, imu2);
        x_nav = interp1(t_nav, x_series, t_eval, 'linear', 'extrap');
        y_nav = interp1(t_nav, y_series, t_eval, 'linear', 'extrap');
        z_nav = interp1(t_nav, z_series, t_eval, 'linear', 'extrap');
        err = [x_ref_eval - x_nav, y_ref_eval - y_nav, z_ref_eval - z_nav];
        rmse3 = sqrt(sum(err.^2, 2));
        if mean(rmse3) < 5
            rmset3 = rmset3 + rmse3;
        else
            T2 = T2 - 1;
        end

    end
    if strcmp( IMMPDAF , 'ON' )
        fprintf('NaveGo: INS/UWB navigation estimates for IMU2, IMMPDAF employed... \n')
        [nav2_e2] = uwb_ins_immpdaf2D(imu2, uwb, num_anchor, ref, gamma1, gamma2);
        err_uiekf_immpda = [ ref.x - nav2_e2.x  ref.y - nav2_e2.y  ref.z - nav2_e2.z ];
        for i = 1:length(ref.t)
            rmse2(i) = norm( err_uiekf_immpda(i,:) );
        end
        if mean(rmse2)<5
            rmset2 = rmset2 + rmse2;
        else
            T = T-1;
        end
    end
    %%%%%%%%%FGO%%%%%%%%%%%%
    if strcmp( FGO , 'ON' )
        fprintf('NaveGo: INS/UWB navigation estimates for IMU2, FGO employed... \n')
        tic
        [nav2_e5] = uwb_ins_FGO2_TC(imu2, uwb, num_anchor, ref);
        toc        % --- FGO error on common grid t_eval (robust time/shape handling) ---
        [t_fgo, x_series, y_series, z_series] = nav_pos_series(nav2_e5, ref, uwb, imu2);
        x_fgo = interp1(t_fgo, x_series, t_eval, 'linear', 'extrap');
        y_fgo = interp1(t_fgo, y_series, t_eval, 'linear', 'extrap');
        z_fgo = interp1(t_fgo, z_series, t_eval, 'linear', 'extrap');
        err_uiekf_FGO = [ x_ref_eval - x_fgo, y_ref_eval - y_fgo, z_ref_eval - z_fgo ];
        rmse5 = sqrt(sum(err_uiekf_FGO.^2, 2));
        if mean(rmse5) < 5
            rmset5 = rmset5 + rmse5;
        else
            T_FGO = T_FGO - 1;
        end

    end

end

%%
% rmset12 = rmset12./times;
if strcmp( EKF , 'ON' )
    rmset11 = rmset11./times;
    armse.rmset_ekf = rmset11;
    nav_e.nav_ekf = nav2_e11;
end
if strcmp( IMMPDAF , 'ON' )
    rmset2 = rmset2./T;
    armse.rmset_immpdaf = rmset2;
    nav_e.nav_immpdaf = nav2_e2;
end
%%%%%%%%FGO%%%%%%%%%%
if strcmp( FGO , 'ON' )
    rmset5 = rmset5./T_FGO;
    armse.rmset_FGO = rmset5;
    nav_e.nav_FGO = nav2_e5;
end

if strcmp( PRPF , 'ON' )
    rmset3 = rmset3./T2;
    armse.rmset_pf = rmset3;
    nav_e.nav_pf = nav2_e3;
end
if strcmp( SHFAEKF , 'ON' )
    rmset4 = rmset4./T1;
    armse.rmset_shfaekf = rmset4;
    nav_e.nav_shfaekf = nav2_e4;
end

toc
save armse.mat armse
save nav_e.mat nav_e

rmse_EKF(u_nlos) = mean(rmset11);
rmse_PRPF(u_nlos) = mean(rmset3);
rmse_SHFAEKF(u_nlos) = mean(rmset4);
rmse_FGO(u_nlos) = mean(rmset5);

% end
% rmse_p_nlos.EKF=rmse_EKF;
% rmse_p_nlos.PRPF=rmse_PRPF;
% rmse_p_nlos.SHFAEKF=rmse_SHFAEKF;
% rmse_p_nlos.FGO=rmse_FGO;
% save rmse_p_nlos1.mat rmse_p_nlos
% rmse_u_nlos.EKF=rmse_EKF;
% rmse_u_nlos.PRPF=rmse_PRPF;
% rmse_u_nlos.SHFAEKF=rmse_SHFAEKF;
% rmse_u_nlos.FGO=rmse_FGO;
% save rmse_u_nlos1.mat rmse_u_nlos
%% PLOT
%非视距阈?
% figure;
% plot(rmse_EKF,'r-+',0.1:0.1:0.8,...
%     0.1:0.1:0.8,rmse_PRPF,'g-*',0.1:0.1:0.8, rmse_SHFAEKF,'c-p',0.1:0.1:0.8,rmse_FGO,'b-^');
% xlabel('Different probability of NLOS errors');
% grid;
% % axis([0.1,1,0,10]);
% ylabel('RMSE/m');
% legend('EKF', 'SHFAF', 'PR-PF','FGO')
% 非视距均?
% figure;
% plot(rmse_EKF,'r-+')%
% hold on
% plot(rmse_PRPF,'g-*')
% hold on
% plot(rmse_SHFAEKF,'c-p')
% hold on
% plot(rmse_FGO,'b-^');
% xlabel('The mean value of NLOS errors');
% grid;
% % axis([0.1,1,0,10]);
% ylabel('RMSE/m');
% legend('EKF', 'SHFAF', 'PR-PF','FGO')
% TRAJECTORY
figure;
hold on
h = gobjects(0); lbl = {};
h(end+1) = plot(ref.x, ref.y, '--k','LineWidth', 1); lbl{end+1} = 'Planned';
if strcmp( EKF , 'ON' )
    [~, x_tr, y_tr, ~] = nav_pos_series(nav2_e11, ref, uwb2, imu2);
    h(end+1) = plot(x_tr, y_tr, 'r','LineWidth', 1); lbl{end+1} = 'EKF';
end
if strcmp( SHFAEKF , 'ON' )
    [~, x_tr, y_tr, ~] = nav_pos_series(nav2_e4, ref, uwb2, imu2);
    h(end+1) = plot(x_tr, y_tr, 'Color', [0.4660 0.6740 0.1880],'LineWidth', 1); lbl{end+1} = 'SHFAEKF';
end
if strcmp( PRPF , 'ON' )
    [~, x_tr, y_tr, ~] = nav_pos_series(nav2_e3, ref, uwb2, imu2);
    h(end+1) = plot(x_tr, y_tr, 'Color',[0.9290 0.6940 0.1250],'LineWidth', 1); lbl{end+1} = 'PR-PF';
end
if strcmp( IMMPDAF , 'ON' )
    [~, x_tr, y_tr, ~] = nav_pos_series(nav2_e2, ref, uwb, imu2);
    h(end+1) = plot(x_tr, y_tr, 'b','LineWidth', 1); lbl{end+1} = 'IMMPDAF';
end
if strcmp( FGO , 'ON' )
    [~, x_tr, y_tr, ~] = nav_pos_series(nav2_e5, ref, uwb, imu2);
    h(end+1) = plot(x_tr, y_tr, 'b-','LineWidth', 1); lbl{end+1} = 'FGO';
end
h(end+1) = plot(ref.x(1), ref.y(1), 'ok', 'MarkerSize', 8, 'LineWidth', 3); lbl{end+1} = 'Start Point';
h(end+1) = plot(ref.x(end), ref.y(end), 'sk', 'MarkerSize', 8, 'LineWidth', 3); lbl{end+1} = 'End Point';
% Anchors (if available)
if exist('uwb2','var') && isfield(uwb2,'anchor')
    ha = scatter(uwb2.anchor(:,1), uwb2.anchor(:,2), 36, 'k^', 'filled');
elseif isfield(uwb,'anchor')
    ha = scatter(uwb.anchor(:,1), uwb.anchor(:,2), 36, 'k^', 'filled');
else
    ha = [];
end
if ~isempty(ha)
    h(end+1) = ha; lbl{end+1} = 'Anchors';
end
axis tight
axis equal
title('TRAJECTORIES')
xlabel('east [m]')
ylabel('north [m]')
legend(h, lbl, 'Location', 'best')
grid on
% TRACKING ERROR
% uwb.t(L)=[];
figure;
hold on
h = gobjects(0); lbl = {};
if strcmp( EKF , 'ON' )
    h(end+1) = plot(t_eval, rmse11, 'r'); lbl{end+1} = 'EKF';
end
if strcmp( SHFAEKF , 'ON' )
    h(end+1) = plot(t_eval, rmse4, 'Color', [0.3010 0.7450 0.9330]); lbl{end+1} = 'SHFAEKF';
end
if strcmp( PRPF , 'ON' )
    h(end+1) = plot(t_eval, rmse3, 'Color', [0.9290 0.6940 0.1250]); lbl{end+1} = 'PR-PF';
end
if strcmp( IMMPDAF , 'ON' )
    h(end+1) = plot(t_eval, rmse2, 'b'); lbl{end+1} = 'IMMPDAF';
end
if strcmp( FGO , 'ON' )
    h(end+1) = plot(t_eval, rmse5, 'b'); lbl{end+1} = 'FGO';
end
axis tight
title('TRACKING ERROR')
xlabel('time [s]')
ylabel('armse [m]')
legend(h, lbl, 'Location', 'best')
grid on


%% ---- Local helpers (robust nav time/shape extraction) ----
function [t_nav, x, y, z] = nav_pos_series(nav, ref, uwb, imu)
%NAV_POS_SERIES Robustly extract (t,x,y,z) column vectors from NaveGo-style nav structs.
% Handles cases where nav.t is missing or where nav.x/nav.y are stored as Nx3 or 3xN.

    xraw = getfield_safe(nav, 'x', []);
    yraw = getfield_safe(nav, 'y', []);
    zraw = getfield_safe(nav, 'z', []);

    n = series_nominal_len(xraw);
    if n == 0, n = series_nominal_len(yraw); end
    if n == 0, n = series_nominal_len(zraw); end
    if n == 0, n = numel(ref.t); end

    % Choose time vector
    if isfield(nav,'t') && numel(nav.t) == n
        t_nav = nav.t(:);
    elseif nargin >= 3 && ~isempty(uwb) && isfield(uwb,'t') && numel(uwb.t) == n
        t_nav = uwb.t(:);
    elseif numel(ref.t) == n
        t_nav = ref.t(:);
    elseif nargin >= 4 && ~isempty(imu) && isfield(imu,'t') && numel(imu.t) == n
        t_nav = imu.t(:);
    else
        if isfield(nav,'t') && ~isempty(nav.t)
            m = min(numel(nav.t), n);
            t_nav = nav.t(1:m);
            n = m;
        else
            t_nav = linspace(ref.t(1), ref.t(end), n).';
        end
    end

    x = series_to_col(xraw, n, 1);
    y = series_to_col(yraw, n, 1);
    if isempty(zraw)
        z = zeros(n,1);
    else
        z = series_to_col(zraw, n, 1);
    end

    % Remove duplicate timestamps for interp1 stability
    [t_nav, ia] = unique(t_nav, 'stable');
    x = x(ia); y = y(ia); z = z(ia);
end

function v = getfield_safe(s, f, defaultVal)
    if isstruct(s) && isfield(s,f)
        v = s.(f);
    else
        v = defaultVal;
    end
end

function n = series_nominal_len(v)
    if isempty(v)
        n = 0;
    elseif isvector(v)
        n = numel(v);
    else
        n = max(size(v));
    end
end

function out = series_to_col(v, n, colIdx)
    if isempty(v)
        out = zeros(n,1); return;
    end
    if isvector(v)
        out = v(:);
    else
        if size(v,1) == n
            out = v(:, min(colIdx, size(v,2)));
        elseif size(v,2) == n
            out = v(min(colIdx, size(v,1)), :).';
        else
            out = v(:);
        end
    end

    if numel(out) < n
        out(end+1:n,1) = out(end);
    elseif numel(out) > n
        out = out(1:n);
    end
end
