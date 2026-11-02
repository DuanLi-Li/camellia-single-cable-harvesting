# Camellia oleifera single-cable vibratory harvesting — data and MATLAB code

油茶单索振枝振动采收：试验数据、MATLAB 处理脚本与论文插图源文件。

Data and code supporting the paper *"Piecewise nonlinear dynamics and inverse design of
excitation parameters for single-cable vibratory harvesting of Camellia oleifera fruit"*.

## 1. 目录说明 / Contents

全部文件平铺在仓库根目录。数据文件命名规则：`sdd_<位置><频率><幅值>a<偏置>b.mat`，
其中位置 `d` = 底端（近根部）、`z` = 中间、`g` = 高端（近自由端），`a` 前为交变激励幅值，
`b` 前为预紧激励偏置，单位均为伺服电机额定力矩（2.39 N·m）的倍数。

### 1.1 自由衰减标定（4.1 节 / 表 3-2 / 图 15）

| 文件 | 夹持位置 | 等效刚度 k (N/m) |
|---|---|---|
| `sdd_928d_1~3.mat` | 近根部 | 8667.1 |
| `sdd_928z_1~3.mat` | 中间位置 | 3917.0 |
| `sdd_928g_1~3.mat` | 近自由端 | 2290.35 |

### 1.2 典型工况响应（4.3 节 / 表 3-3 / 图 17）

| 位置 | 7 Hz | 10 Hz | 13 Hz |
|---|---|---|---|
| 近根部 | `sdd_d7hz0.25a0.4b.mat` | `sdd_d10hz0.25a0.9b.mat` | `sdd_d13hz0.3a1b.mat` |
| 中间位置 | `sdd_z7hz0.15a0.4b.mat` | `sdd_z10hz0.35a1.1b.mat` | `sdd_z13hz0.35a0.95b.mat` |
| 近自由端 | `sdd_g7hz0.2a1b.mat` | `sdd_g10hz0.7a1b.mat` | `sdd_g13hz0.2a1b.mat` |

### 1.3 逆向设计验证（4.4 节 / 表 3-4 / 图 18）

| 位置 | 10 Hz | 13 Hz |
|---|---|---|
| 近根部 | `sdd_nitd10hz0.235a1b.mat` | `sdd_nitd13hz0.34a1b.mat` |
| 中间位置 | `sdd_nitz10hz0.29a1.1b.mat` | `sdd_nitz13hz0.34a1b.mat` |
| 近自由端 | `sdd_g10hz0.3584a1b.mat` | `sdd_13hz0.4090a1b.mat` |

### 1.4 田间原位试验（4.5 节 / 表 3-5 / 图 21）

| 文件 | 说明 |
|---|---|
| `sdd_10_8_1~3.mat` | 田间枝条三次自由衰减标定 |
| `sdd_2cha7hz0.2a0.7b.mat` | 田间 7 Hz（Amp = 0.20，Bia = 0.70） |
| `sdd_2cha10hz0.45a0.85b.mat` | 田间 10 Hz（Amp = 0.45，Bia = 0.85） |

### 1.5 庞加莱映射插图源文件（4.2 节 / 图 16）

| 文件 | 说明 |
|---|---|
| `experiment.fig` | 实测映射点 |
| `modeAl.fig` | 仿真映射点 |
| `kkmodel2zhouqi.fig` | 仿真 2-周期映射序列 |
| `kkmodelandexp1zhouqi.fig` | 仿真与实测 1-周期映射对比 |

### 1.6 代码

| 文件 | 说明 | 对应论文 |
|---|---|---|
| `bianshi.m` | 单自由度树枝等效参数辨识：对数衰减率法 + 非线性最小二乘，双法交叉校验并输出拟合优度 R² | 4.1 节；表 3-2；图 15 |

## 2. 运行 / How to run

把本目录设为 MATLAB 当前工作目录，直接运行：

```matlab
bianshi
```

- 环境：MATLAB R2020b 或更高版本；需要 Signal Processing Toolbox（`findpeaks`）。
- `fminsearch` 属 MATLAB 基础函数，无需 Optimization Toolbox。
- 脚本默认 `k = 2290.35`（近自由端）、载入 `sdd_928g_3.mat`；辨识其他位置时改这两处即可。

## 3. 数据格式 / Data format

每个 `.mat` 为 Simulink 记录的 `Simulink.SimulationData.Dataset`：

- `X.Data` —— 时间序列 (s)
- `Y(1).Data` / `Y(2).Data` —— 树枝末端加速度 (m/s²)，部分文件含仿真与实测两路

采样率由时间序列差分得到；原始数据未做滤波，仅去直流零点。

## 4. 引用 / Citation

若使用本仓库的代码或数据，请引用：

> （论文信息待补：作者、题名、期刊、年、DOI）

## 5. 许可 / License

代码采用 MIT 许可，见 `LICENSE`。
