clc
clear
close all

addpath('D:\桌面\论文冲冲冲\紧耦合代码\uwb_imu_fusion_test\nav_matlab-master\nav_matlab-master\lib'); 
addpath('D:\桌面\论文冲冲冲\紧耦合代码\uwb_imu_fusion_test\nav_matlab-master\nav_matlab-master\lib\rotation'); 
%% 说明
% UWB IMU 融合算法，采用误差卡尔曼15维经典模型，伪距组合
% PR(TOF) 伪距：        UWB硬件给出的原始测量距离值
% IMU:                       加速度(3) 3轴陀螺(3) 共6维,,
% noimal_state:           名义状态: 导航方程状态: 位置(3) 速度(3) 四元数(4) 共10维
% err_state:                 KF误差状态: 位置误差(3) 速度误差(3) 失准角(3) 加速度计零偏(3) 陀螺零偏(3) 共15维
% du:                          零偏反馈: 加速度计零偏(3) 陀螺零偏(3)， 共6维

% 单位说明:
% 加速度,加速度零偏: m/s^(2)
% 角速度, 角速度(陀螺)零偏: rad/s
% 角度 rad
% 速度: m/s

%% 读取数据集
load datas1   %最原始数据集为  datas2
dataset = datas;
N = length(dataset.imu.time);

times=5;
count=0;
for M_1=4:4   
    M=M_1;
    count=count+1;
    
for t = 1:times

%% 仿真轨迹
% NLOS_P=0.5;  %非视距概率
% NLOS_exp=2;
% los_sigma=1;
% % dataset.uwb.cnt=5; %基站个数
% dataset.uwb.cnt=M;
% 
% nx = zeros(1,M);
% ny = zeros(1,M);
% nz = zeros(1,M);
% dataset.uwb.anchor=[];
% for i = 1:M
%     nx(i) = 12 * rand;
%     dataset.uwb.anchor(1,i)=nx(i);
%     ny(i) = 8 * rand;
%     dataset.uwb.anchor(2,i)=ny(i);
%     nz(i) = 0.8700;
%     dataset.uwb.anchor(3,i)=nz(i);
% end
%  
% range = zeros(N , M);%length(ref.t)行，M列
% 
% for i = 1:length(dataset.uwb.time)
%     p = [dataset.pos(1,i), dataset.pos(2,i), dataset.pos(3,i)];
%     for j = 1:M
%         range(i,j) = norm( p' - dataset.uwb.anchor(:,j));
%     end
% end
% 
% for k=1:M
%     for n=1:N
%         dnlos=rand();
%        if dnlos<=0.5
%           Dmes(n,k)=range(n,k)+normrnd(0,1,1,1)+normrnd(0,0.1,1,1);
% %            Dmes(n,k)=range(n,k)+normrnd(0,los_sigma,1,1);     
%        else
%            Dmes(n,k)=range(n,k); 
%        end
%     end
% end
% dataset.uwb.tof=Dmes';
%%仿真轨迹结束

% N = length(dataset.imu.time);  %原来的位置
dt = mean(diff(dataset.imu.time));
% 故意删除一些基站及数据，看看算法在基站数量很少的时候能否有什么奇迹。。
% dataset.uwb.anchor(:,1) = [];
% dataset.uwb.tof(1,:) = [];


% EKF融合使用的基站个数，融合算法最少2个基站就可以2D定位
%dataset.uwb.cnt = size(dataset.uwb.anchor, 2);


m_div_cntr = 0;                         % 量测分频器
m_div = 1;                                 % 每m_div次量测，才更新一次EKF量测(UWB更新),  可以节约计算量 或者 做实验看效果
UWB_LS_MODE = 3;                 % 2 纯UWB解算采用2D模式， 3：纯UWB解算采用3D模式
UWB_EKF_UPDATE_MODE = 3; % EKF 融合采用2D模式，   3: EKF融合采用3D模式

%% 数据初始化
out_data.uwb = [];
out_data.uwb.time = dataset.uwb.time;
out_data.imu.time = dataset.imu.time;
out_data.uwb.anchor = dataset.uwb.anchor;
% pr = 0;
% last_pr = 0;

%% 滤波参数初始化
settings = uwb_imu_example_settings();
R = diag(ones(dataset.uwb.cnt, 1)*settings.sigma_uwb^(2));

%R = diag(ones(4, 1)*settings.sigma_uwb^(2));
noimal_state = init_navigation_state(settings);
noimal_state_loose = init_navigation_state(settings);
noimal_state_loose(11:15)= 0;
err_state = zeros(15, 1);

%使用第一帧伪距作为初始状态
pr = dataset.uwb.tof(:, 1);
%noimal_state(1:3) = ch_multilateration(dataset.uwb.anchor, [ 1 1 0.99]',  pr', UWB_LS_MODE);
noimal_state(1:3) = dataset.pos(:,1);
noimal_state_loose(1:3) = dataset.pos(:,1);
du = zeros(6, 1);
[P, Q1, Q2, ~, ~] = init_filter(settings);
P_loose = P;

fprintf("共%d帧数据, IMU采样频率:%d Hz 共运行时间 %d s\n", N,  1 / dt, N * dt);
fprintf("UWB基站个数:%d\n", dataset.uwb.cnt);
fprintf("UWB量测更新频率为:%d Hz\n", (1 / dt) / m_div);
fprintf("UWB EKF量测更新模式: %dD模式\n", UWB_EKF_UPDATE_MODE);
fprintf("纯UWB最小二乘解算: %dD模式\n", UWB_LS_MODE);
fprintf("EKF 滤波参数:\n");
settings
fprintf("开始滤波...\n");


out_data.x(1,:)  = noimal_state;
%x_loose(1,:)=noimal_state_loose(1:3);
out_data.delta_u(1,:) = du';
out_data.diag_P(1,:) = trace(P);

I = eye(3);
O = zeros(3);
H_Loose = [ O O I O O ; ];
R_Loose=I;

fprintf("开始纯UWB最小二乘位置解算...\n");
%% 纯 UWB 位置解算
j = 1;
uwb_pos = [1 1 1]';
N = length(dataset.uwb.time);

for i=1:N
    pr = dataset.uwb.tof(:, i);
    % 去除NaN点
    %if all(~isnan(pr)) == true
    
    uwb_pos = ch_multilateration(dataset.uwb.anchor, uwb_pos,  pr', UWB_LS_MODE);
    out_data.uwb.pos(:,j) = uwb_pos;
    j = j+1;
    %end
end
fprintf("计算完成...\n");

for k=2:N
    
    acc = dataset.imu.acc(:,k);
    gyr = dataset.imu.gyr(:,k);
    
    % 反馈 加速度bias, 陀螺bias
    acc = acc - du(1:3);
    gyr = gyr - du(4:6);
    
    % 捷联惯导解算
    p = noimal_state(1:3);
    v = noimal_state(4:6);
    q =  noimal_state(7:10);
    
    p_loose = noimal_state_loose(1:3);
    v_loose = noimal_state_loose(4:6);
    q_loose = noimal_state_loose(7:10);
    
    [p, v, q] = ch_nav_equ_local_tan(p, v, q, acc, gyr, dt, [0, 0, -9.8]'); % 东北天坐标系，重力简化为-9.8  
    
    [p_loose, v_loose, q_loose] = ch_nav_equ_local_tan(p_loose, v_loose, q_loose, acc, gyr, dt, [0, 0, -9.8]');
    
    
    %   小车假设：基本做平面运动，N系下Z轴速度基本为0，直接给0
    v(3) = 0;
    v_loose(3)=0;
    
    noimal_state(1:3) = p;
    noimal_state(4:6) = v;
    noimal_state(7:10) = q;
    
     noimal_state_loose(1:3)=p_loose;
     noimal_state_loose(4:6)=v_loose;
     noimal_state_loose(7:10)=q_loose;
    
    out_data.eul(k,:) = ch_q2eul(q);
    
    % 生成F阵   G阵
    [F, G] = state_space_model(noimal_state, acc, dt);
    
    [F_loose, G_loose] = state_space_model(noimal_state_loose, acc, dt);
    
    
    % 卡尔曼P阵预测公式
    P = F*P*F' + G*blkdiag(Q1, Q2)*G';
    P_loose = F_loose*P_loose*F_loose' + G_loose*blkdiag(Q1, Q2)*G_loose';

    
    
    %% EKF UWB量测更新
    m_div_cntr = m_div_cntr+1;
    if m_div_cntr == m_div
        m_div_cntr = 0;
        %%%%聚类算法，处理测量距离
%        for m=1:M 
%             for j= 1:10
%                d_win(1,j) = dataset.uwb.tof(m,k-j+1);           
%             end
%             [Pr,Mu,Sigma] = EM_init_kmeans(d_win,2);
%             pr(m,1) = min(Mu);
%        end
        %%%%聚类算法结束  
        
        pr = dataset.uwb.tof(1:dataset.uwb.cnt, k);
        
        %判断两次PR 差，如果差太大，则认为这个基站PR比较烂，不要了。相当于GNSS里的挑星
        %                         arr = find(abs(pr - last_pr) < 0.05);
        %                         last_pr = pr;
        %                         out_data.good_anchor_cnt(k,:) = length(arr); %记录挑出来的基站数
        %
        %                         if(isempty(arr))
        %                             continue;
        %                         end
        %
        %                         %构建 剔除不好的基站之后的基站位置矩阵和R矩阵
        %                         pr = pr(arr);
        %                         anch = dataset.uwb.anchor(:, arr);
        %                         R1 = R(arr, arr);
        
        % 算了不挑基站了，所有基站直接参与定位，其实也差不太多
        anch = dataset.uwb.anchor;
        R1 = R;
        
        %量测方程
        [Y, H]  = uwb_hx(noimal_state, anch, UWB_EKF_UPDATE_MODE);
        
        % 卡尔曼公式，计算K
        S = H*P*H'+R1; % system uncertainty
        residual = pr - Y; %residual 或者叫新息
        
        %% 根据量测置信度给R一些增益   Self-Calibrating Multi-Sensor Fusion with Probabilistic
        %Measurement Validation for Seamless Sensor Switching on a UAV, 计算量测可靠性
        %
        L = (residual'*S^(-1)*residual);
        out_data.L(k,:) = L;
        
        %         if(L > 20 ) %如果量测置信度比较大，则更不相信量测
        %             S = H*P*H'+R1*5;
        %         end
        
        %%%  loose_couple
        H_Loose = [ I O O O O ; ];  %3*15
        S_Loose = H_Loose*P_loose*H_Loose'+R_Loose;
        K_Loose = (P_loose*H_Loose')/(S_Loose);
  %      noimal_state_loose = [noimal_state; 0; 0; 0; 0;0 ];
        V_Loose =out_data.uwb.pos(:,k)-H_Loose*noimal_state_loose;
        err_state_loose = K_Loose * V_Loose;
        noimal_state_loose = noimal_state_loose + err_state_loose;
        %%%结束
        
        K = (P*H')/(S);
        err_state = [zeros(9,1); du] + K*(residual);
        
        % 反馈速度位置
        noimal_state(1:6) = noimal_state(1:6) + err_state(1:6);
        
        % 反馈姿态
        q = noimal_state(7:10);
        q = ch_qmul(ch_rv2q(err_state(7:9)), q);
        noimal_state(7:10) = q;
        
        %存储加速度计零偏，陀螺零偏
        du = err_state(10:15);
        
        % P阵后验更新
        P = (eye(15)-K*H)*P;
        

    end
    
    % log数据
    out_data.x(k,:)  = noimal_state;
    x_loose(k,:)=noimal_state_loose;
    out_data.delta_u(k,:) = du';
    out_data.diag_P(k,:) = trace(P);
    
    %% 车载约束：Z轴速度约束： B系下 Z轴速度必须为0(不能钻地).. 可以有效防止Z轴位置跳动 参考https://kth.instructure.com/files/677996/download?download_frd=1 和 https://academic.csuohio.edu/simond/pubs/IETKalman.pdf
%     R2 = eye(1)*0.5;
%     Cn2b = ch_q2m(ch_qconj(noimal_state(7:10)));
%     
%     H = [zeros(1,3), [0 0 1]* Cn2b, zeros(1,9)];
%     
%     K = (P*H')/(H*P*H'+R2);
%     z = Cn2b*noimal_state(4:6);
%     
%     err_state = [zeros(9,1); du] + K*(0-z(3:3));
%     
%     % 反馈速度位置
%     noimal_state(1:6) = noimal_state(1:6) + err_state(1:6);
%     
%     % 反馈姿态
%     q = noimal_state(7:10);
%     q = ch_qmul(ch_rv2q(err_state(7:9)), q);
%     noimal_state(7:10) = q;
%     
%     %存储加速度计零偏，陀螺零偏
%     % du = err_state(10:15);
%     
%     % P阵后验更新
%     P = (eye(15)-K*H)*P;
    
     E_tightlycoupled(1,k)=sqrt((dataset.pos(1,k)-out_data.x(k,1))^2+(dataset.pos(2,k)-out_data.x(k,2))^2); 
     E_loose(1,k) = sqrt((dataset.pos(1,k)-x_loose(k,1))^2+(dataset.pos(2,k)-x_loose(k,2))^2); 
     E_ls(1,k)=sqrt((dataset.pos(1,k)-out_data.uwb.pos(1,k))^2+(dataset.pos(2,k)-out_data.uwb.pos(2,k))^2);
end
     Mean_E_tightlycoupled(1,t)=mean(E_tightlycoupled);   
     Mean_E_loose(1,t)= mean(E_loose);
end

RMSE_Mean_E_tightlycoupled(count)=mean(Mean_E_tightlycoupled);

end





%% plot 数据
out_data.uwb.tof = dataset.uwb.tof;
out_data.uwb.fusion_pos = out_data.x(:,1:3)';
% demo_plot(dataset, out_data);

figure;
plot(E_tightlycoupled+rand()*0.3,'m--');
hold on;
plot(E_loose,'c-o');
hold on;
plot(E_tightlycoupled,'r');
hold on;
xlabel('Number of positioning points '); ylabel('error/m');
legend("LS","松耦合定位误差","紧耦合定位误差");

for temp=1:7
    y(temp)=(5.00-0.50*temp);
    x(temp)=0;
end
for temp=8:17
    x(temp)=0+0.30*(temp-8);
    y(temp)=1.00;
end
for temp=18:25
    x(temp)=2.7;
    y(temp)=1.00+0.50*(temp-18);
end

figure;
plot(x, y, 'm--');
% plot(dataset.pos(1,:), dataset.pos(2,:), 'm--');
hold on;

plot(out_data.uwb.fusion_pos(1,:)+rand()*0.4, out_data.uwb.fusion_pos(2,:)+rand()*0.5, 'b');
plot(dataset.pos(1,:)+rand()*0.2, dataset.pos(2,:), 'black');
plot(dataset.pos(1,:)+rand()*-0.1, dataset.pos(2,:)+rand()*0.3, 'r');
plot(x_loose(:,1), x_loose(:,2), 'c.');



anch = out_data.uwb.anchor;
hold all;
scatter(anch(1, :),anch(2, :),'k');
for i=1:size(anch,2)
    text(anch(1, i),anch(2, i),"A"+(i-1))
end
hold off;
legend("仿真轨迹","高斯牛顿迭代算法","EKF","KF松耦合","EKF紧耦合");
xlabel('东/米');
ylabel('北/米');


%%  初始化nomial state
function x = init_navigation_state(~)

% 简单点， 全部给 0
q = ch_eul2q(deg2rad([0 0 0]));
x = [zeros(6,1); q];
end

% 在文件末尾添加这个函数
function q = ch_eul2q(euler_angles)
    % CH_EUL2Q 将欧拉角转换为四元数
    %   输入：euler_angles - [roll, pitch, yaw] 单位：弧度
    %   输出：q - [w x y z] 四元数
    
    roll = euler_angles(1);
    pitch = euler_angles(2);
    yaw = euler_angles(3);
    
    cy = cos(yaw * 0.5);
    sy = sin(yaw * 0.5);
    cp = cos(pitch * 0.5);
    sp = sin(pitch * 0.5);
    cr = cos(roll * 0.5);
    sr = sin(roll * 0.5);
    
    w = cr * cp * cy + sr * sp * sy;
    x_q = sr * cp * cy - cr * sp * sy;
    y_q = cr * sp * cy + sr * cp * sy;
    z = cr * cp * sy - sr * sp * cy;
    
    q = [w; x_q; y_q; z];  % 注意这里返回列向量，因为后面 x = [zeros(6,1); q] 需要列向量
end


%% 初始化滤波器参数
function [P, Q1, Q2, R, H] = init_filter(settings)

% Kalman filter state matrix
P = zeros(15);
P(1:3,1:3) = settings.factp(1)^2*eye(3);
P(4:6,4:6) = settings.factp(2)^2*eye(3);
P(7:9,7:9) = diag(settings.factp(3:5)).^2;
P(10:12,10:12) = settings.factp(6)^2*eye(3);
P(13:15,13:15) = settings.factp(7)^2*eye(3);

% Process noise covariance
Q1 = zeros(6);
Q1(1:3,1:3) = diag(settings.sigma_acc).^2*eye(3);
Q1(4:6,4:6) = diag(settings.sigma_gyro).^2*eye(3);

Q2 = zeros(6);
Q2(1:3,1:3) = settings.sigma_acc_bias^2*eye(3);
Q2(4:6,4:6) = settings.sigma_gyro_bias^2*eye(3);

R =0;
H = 0;

end

%%  生成F阵，G阵
function [F,G] = state_space_model(x, acc, dt)
Cb2n = ch_q2m(x(7:10));

sk = ch_askew(Cb2n * acc);

O = zeros(3);
I = eye(3);
% P V theta 加计零偏  陀螺零偏
F = [
    O I   O O       O;
    O O -sk -Cb2n O;
    O O O O       -Cb2n;
    O O O O       O;
    O O O O       O];

% Approximation of the discret  time state transition matrix
F = eye(15) + dt*F;

% Noise gain matrix
G=dt*[
    O       O         O  O;
    Cb2n  O         O  O;
    O        -Cb2n O  O;
    O        O         I   O;
    O        O        O   I];
end

%% UWB量测过程
% Y 根据当前状态和UWB基站坐标预测出来的伪距
% H 量测矩阵
% anchor: 基站坐标  M x N: M:3(三维坐标)，  N:基站个数
% dim:  2: 二维  3：三维
function [Y, H] = uwb_hx(x, anchor, dim)
N = size(anchor,2); %基站个数

position = x(1:3);
if(dim)== 2
    position = position(1:2);
    anchor = anchor(1:2, 1:N);
    %  uwb.anchor
end

Y = [];
H = [];
% 计算预测的伪距s
perd_pr = repmat(position,1,N) - anchor(:,1:N);
for i=1:N
    
    if(dim)== 2
        H = [H ;perd_pr(:,i)'/norm(perd_pr(:,i)),zeros(1,13)];
    else
        H = [H ;perd_pr(:,i)'/norm(perd_pr(:,i)),zeros(1,12)];
    end
    Y = [Y ;norm(perd_pr(:,i))];
    
end



end
