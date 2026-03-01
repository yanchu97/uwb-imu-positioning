# Paper Notes

## Goal
这篇论文主要解决 UWB + IMU 融合定位问题。

## Key ideas
- 用 IMU 做短时运动预测
- 用 UWB 做位置校正
- 通过滤波/优化减少漂移

## Variables / outputs
- 输入：UWB ranges, IMU accel, gyro
- 输出：position, velocity, orientation

## Map to this repo
- `src/...`：传感器预处理
- `...`：融合主逻辑
- `...`：评估/可视化

## What Codex should focus on
请重点检查：
1. 代码里的状态量定义是否和论文一致
2. 观测模型是否对应论文公式
3. 初始化和坐标系处理是否合理
