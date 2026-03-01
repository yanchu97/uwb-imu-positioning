function [nav2_e5] = uwb_ins_FGO2(imu2, uwb,~, ref)

import gtsam.*;

%% Read metadata and compute relative sensor Position transforms 读取元数据并计算相对传感器位姿变换
% IMU metadata
L=length(uwb.t);
IMU_metadata = importdata(findExampleDataFile('KittiEquivBiasedImu_metadata.txt'));
IMU_metadata = cell2struct(num2cell(IMU_metadata.data), IMU_metadata.colheaders, 2);
IMUinBody = Pose3.Expmap([IMU_metadata.BodyPtx; IMU_metadata.BodyPty; IMU_metadata.BodyPtz;
    IMU_metadata.BodyPrx; IMU_metadata.BodyPry; IMU_metadata.BodyPrz; ]);
if ~IMUinBody.equals(Pose3, 1e-5)
    error 'Currently only support IMUinBody is identity, i.e. IMU and body frame are the same';
end
%% IMU2 SYNTHETIC DATA
rng('shuffle')					% Reset pseudo-random seed
% fprintf('NaveGo: generating IMU2 synthetic data... \n') % 打印文字
imu2.dt=0.01*ones(length(ref.t),1);
% 整合IMU数据
IMU_data.data=cat(2,imu2.t,imu2.dt,imu2.fb,imu2.wb);
IMU_data.colheaders={'Time','dt','accelX','accelY','accelZ','omegaX','omegaY','omegaZ'};
IMU_data = cell2struct(num2cell(IMU_data.data), IMU_data.colheaders, 2);
imum = cellfun(@(x) x', num2cell([ [IMU_data.accelX]' [IMU_data.accelY]' [IMU_data.accelZ]' [IMU_data.omegaX]' [IMU_data.omegaY]' [IMU_data.omegaZ]' ], 2), 'UniformOutput', false);
[IMU_data.acc_omega] = deal(imum{:});
clear imum
%% uwb PROCESS
BS=uwb.anchor(1:5,1:2);
[uwb.X,uwb.Y] =fun_PDA(ref.x(1,1),ref.y(1,1),uwb.range_30,BS,0.1);

%% uwb data
uwb.Z(1:size(uwb.t,1),1)=2;
uwb_data.data=cat(2,uwb.t,uwb.X',uwb.Y',uwb.Z);
uwb_data.colheaders={'Time','X','Y','Z'};
uwb_data = cell2struct(num2cell(uwb_data.data), uwb_data.colheaders, 2);
for i = 1:numel(uwb_data)
    uwb_data(i).Position = gtsam.Point3(uwb_data(i).X, uwb_data(i).Y, uwb_data(i).Z);
end
noiseModeluwb = noiseModel.Diagonal.Precisions([ [0.5;0.5;0.5]; [0.5;0.5;0.5] ]);
firstuwbPosition = 2;
%% Get initial conditions for the estimated trajectory 获取估计轨迹的初始条件
currentPositionGlobal = Pose3(Rot3, uwb_data(firstuwbPosition).Position); % initial Position is the reference frame (navigation frame) 初始姿势是参考框架（导航框架）
currentVelocityGlobal = LieVector([0;0;0]); % the vehicle is stationary at the beginning 车辆一开始是静止的
currentBias = imuBias.ConstantBias([0.1;0.1;0.1], [0.1;0.1;0.1]);
sigma_init_x = noiseModel.Isotropic.Precisions([ 0.0; 0.0; 0.0; 1; 1; 1 ]);
sigma_init_v = noiseModel.Isotropic.Sigma(3, 1000.0);
sigma_init_b = noiseModel.Isotropic.Sigmas([ 0.100; 0.100; 0.100; 5.00e-05; 5.00e-05; 5.00e-05 ]);
sigma_between_b = [ IMU_metadata.AccelerometerBiasSigma * ones(3,1); IMU_metadata.GyroscopeBiasSigma * ones(3,1) ];
g = [0;0;-9.8];
w_coriolis = [0;0;0]; %coriolis 对旋转体系中进行直线运动的质点由于惯性相对于旋转体系产生的直线运动的偏移的一种描述

%% Solver object 设置求解器相关 isam
isamParams = ISAM2Params;
isamParams.setFactorization('CHOLESKY'); %Cholesky 分解是把一个对称正定的矩阵表示成一个下三角矩阵L和其转置的乘积的分解
isamParams.setRelinearizeSkip(5); %设置重新线性化的间隔
isam = gtsam.ISAM2(isamParams);
newFactors = NonlinearFactorGraph;
newValues = Values;

%% Main loop: 主循环 FGO
IMUtimes = [IMU_data.Time];
% 启动主循环：在每个时间步执行推理
for measurementIndex = firstuwbPosition:L

    % At each non=IMU measurement we initialize a new node in the graph 在每个非 IMU 测量中，我们在图中初始化一个新节点
    currentPositionKey = symbol('x',measurementIndex);
    currentVelKey =  symbol('v',measurementIndex);
    currentBiasKey = symbol('b',measurementIndex);
    t = uwb_data(measurementIndex, 1).Time;

    if measurementIndex == firstuwbPosition
        %% Create initial estimate and prior on initial Position, velocity, and biases
        % 在初始姿势、速度和偏差上创建初始估计和先验
        newValues.insert(currentPositionKey, currentPositionGlobal);
        newValues.insert(currentVelKey, currentVelocityGlobal);
        newValues.insert(currentBiasKey, currentBias);
        newFactors.add(PriorFactorPose3(currentPositionKey, currentPositionGlobal, sigma_init_x));
        newFactors.add(PriorFactorLieVector(currentVelKey, currentVelocityGlobal, sigma_init_v));
        newFactors.add(PriorFactorConstantBias(currentBiasKey, currentBias, sigma_init_b));
    else
        t_previous = uwb_data(measurementIndex-1, 1).Time;
        %% Summarize IMU data between the previous uwb measurement and now
        IMUindices = find(IMUtimes >= t_previous & IMUtimes < t);

        currentSummarizedMeasurement = gtsam.ImuFactorPreintegratedMeasurements( ...
            currentBias, IMU_metadata.AccelerometerSigma.^2 * eye(3), ...
            IMU_metadata.GyroscopeSigma.^2 * eye(3), IMU_metadata.IntegrationSigma.^2 * eye(3));

        for imuIndex = IMUindices
            accMeas = [ IMU_data(imuIndex).accelX; IMU_data(imuIndex).accelY; IMU_data(imuIndex).accelZ ];
            omegaMeas = [ IMU_data(imuIndex).omegaX; IMU_data(imuIndex).omegaY; IMU_data(imuIndex).omegaZ ];
            deltaT = IMU_data(imuIndex).dt;
            currentSummarizedMeasurement.integrateMeasurement(accMeas, omegaMeas, deltaT);
        end

        % Create IMU factor
        newFactors.add(ImuFactor( ...
            currentPositionKey-1, currentVelKey-1, ...
            currentPositionKey, currentVelKey, ...
            currentBiasKey, currentSummarizedMeasurement, g, w_coriolis));

        % Bias evolution as given in the IMU metadata IMU 元数据中给出的偏差演化 偏置因子
        newFactors.add(BetweenFactorConstantBias(currentBiasKey-1, currentBiasKey, imuBias.ConstantBias(zeros(3,1), zeros(3,1)), ...
            noiseModel.Diagonal.Sigmas(sqrt(numel(IMUindices)) * sigma_between_b)));

        % Create uwb factor
        uwbPosition = Pose3(currentPositionGlobal.rotation, uwb_data(measurementIndex).Position);
        newFactors.add(PriorFactorPose3(currentPositionKey, uwbPosition, noiseModeluwb));

        % Add initial value
        newValues.insert(currentPositionKey, uwbPosition);
        newValues.insert(currentVelKey, currentVelocityGlobal);
        newValues.insert(currentBiasKey, currentBias);
        isam.update(newFactors, newValues);
%         measurementIndex
        result = isam.calculateEstimate();
        result_Position = utilities.extractPose3(result);
        newFactors = NonlinearFactorGraph;
        newValues = Values;
    end
end
nav2_e5.x(1)=ref.x(1,1);
nav2_e5.y(1)=ref.y(1,1);
nav2_e5.z(1)=ref.z(1,1);
for i=1:L-1
    nav2_e5.x(i+1)=result_Position(i,10);
    nav2_e5.y(i+1)=result_Position(i,11);
    nav2_e5.z(i+1)=2;
end
end