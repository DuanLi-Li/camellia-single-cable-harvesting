%% ===================== 二阶线弹簧系统参数辨识（修复版）=====================
% 模型：m*ddx + b*dx + k*x = 0
% 说明：全程无滤波，仅去除直流零点；修复峰值识别逻辑，保证两方法结果一致
clear; clc; close all;

%% ===================== 1. 数据导入与系统参数 =====================
load('sdd_928g_3.mat'); 
t_full = sdd_928g_3.X.Data';       % 完整时间数据
a_full = sdd_928g_3.Y.Data';       % 完整加速度数据(m/s²)

% 已知线刚度
k = 2290.35;  

fprintf('==================== 系统基础参数 ====================\n');
fprintf('线刚度 k = %.2f N/m\n', k);
fprintf('======================================================\n\n');

%% ===================== 2. 截取自由振动段与预处理（无滤波）=====================
start_time = 0.097;  
end_time = 1.8;     % 截取前1.5秒高信噪比数据，舍弃后期噪声段

start_idx = find(t_full >= start_time, 1, 'first');
end_idx = find(t_full >= end_time, 1, 'first');

t = t_full(start_idx:end_idx) - t_full(start_idx);
a_raw = a_full(start_idx:end_idx);
n = length(t);
fs = 1/(t(2)-t(1));

% 仅去直流零点，不做任何滤波
offset = mean(a_raw(round(n*0.8):end));
a = a_raw - offset;

fprintf('已截取数据：%.2f ~ %.2f s，共 %d 个采样点\n', start_time, end_time, n);
fprintf('采样率 fs = %.1f Hz\n\n', fs);

%% ===================== 方法1：对数衰减率法（修复峰值识别）=====================
fprintf('==================== 方法1：对数衰减率法 ====================\n');

% ---------- 核心修复：限制最小峰间距，过滤噪声毛刺 ----------
% 预估频率约10Hz，周期0.1s，最小峰间距设为周期的70%，避免一个周期识别多个峰
min_peak_dist = round(0.07 * fs);  
% 最小峰值高度设为最大值的20%，过滤后期小噪声峰
[pks, locs] = findpeaks(a, 'MinPeakHeight', max(a)*0.2, 'MinPeakDistance', min_peak_dist);
peak_num = length(pks);
t_peaks = t(locs);

if peak_num < 5
    error('有效峰值过少，请调整截取时间或峰值阈值');
end

% 1. 计算振动周期与频率
T_list = diff(t_peaks);
Td_avg = mean(T_list);          % 平均有阻尼周期
fd = 1 / Td_avg;                 % 有阻尼固有频率
omega_d = 2 * pi * fd;           % 有阻尼固有角频率

% 2. 计算对数衰减率与阻尼比（取相邻峰值的平均值，更稳健）
delta_list = log(pks(1:end-1) ./ pks(2:end));
delta_avg = mean(delta_list);
zeta1 = delta_avg / sqrt(4*pi^2 + delta_avg^2);

% 3. 无阻尼固有频率
omega_n1 = omega_d / sqrt(1 - zeta1^2);
fn1 = omega_n1 / (2*pi);

% 4. 计算等效质量与线阻尼系数
m1 = k / (omega_n1^2);
b1 = 2 * zeta1 * sqrt(m1 * k);

fprintf('识别有效峰值数量：%d 个\n', peak_num);
fprintf('平均振动周期 Td = %.4f s\n', Td_avg);
fprintf('有阻尼固有频率 fd = %.2f Hz\n', fd);
fprintf('阻尼比 ζ = %.4f (%.2f%%)\n', zeta1, zeta1*100);
fprintf('无阻尼固有频率 fn = %.2f Hz\n', fn1);
fprintf('等效质量 m = %.4f kg\n', m1);
fprintf('线阻尼系数 b = %.4f N·s/m\n', b1);
fprintf('==============================================================\n\n');

%% ===================== 方法2：非线性最小二乘法（全局拟合）=====================
fprintf('==================== 方法2：非线性最小二乘法 ====================\n');

% 欠阻尼加速度解析模型
acc_model = @(p, t) p(1) .* exp(-p(2)*p(3)*t) .* cos( p(3)*sqrt(1-p(2)^2)*t + p(4) );

% 自动提取初始参数（用修复后的峰值识别逻辑，保证初始值准确）
A0 = pks(1);
zeta0 = zeta1;
omega_n0 = omega_n1;
phi0 = -omega_d * t_peaks(1);  % 根据第一个峰值位置反推初相位
params0 = [A0, zeta0, omega_n0, phi0];

% 拟合设置
options = optimset(...
    'MaxIter', 30000, ...
    'MaxFunEvals', 100000, ...
    'TolX', 1e-12, ...
    'TolFun', 1e-12, ...
    'Display', 'off');

% 执行全局最优拟合
params_fit = fminsearch(@(p) sum((a - acc_model(p, t)).^2), params0, options);

% 提取拟合结果
A_fit = params_fit(1);
zeta2 = abs(params_fit(2));
omega_n2 = params_fit(3);
phi_fit = params_fit(4);
omega_d2 = omega_n2 * sqrt(1 - zeta2^2);
fn2 = omega_n2 / (2*pi);

% 计算物理参数
m2 = k / (omega_n2^2);
b2 = 2 * zeta2 * sqrt(m2 * k);

fprintf('固有频率 fn = %.2f Hz\n', fn2);
fprintf('阻尼比 ζ = %.4f (%.2f%%)\n', zeta2, zeta2*100);
fprintf('等效质量 m = %.4f kg\n', m2);
fprintf('线阻尼系数 b = %.4f N·s/m\n', b2);
fprintf('==============================================================\n\n');

%% ===================== 双方法结果对比与验证 =====================
a_fit = acc_model(params_fit, t);
SS_res = sum((a - a_fit).^2);
SS_tot = sum((a - mean(a)).^2);
R2 = 1 - SS_res / SS_tot;

fprintf('==================== 结果对比 =====================\n');
fprintf('拟合优度 R² = %.4f\n', R2);
fprintf('固有频率相对偏差：%.2f %%\n', abs(fn1 - fn2)/fn1 * 100);
fprintf('阻尼比相对偏差：%.2f %%\n', abs(zeta1 - zeta2)/zeta1 * 100);
fprintf('等效质量相对偏差：%.2f %%\n', abs(m1 - m2)/m1 * 100);
fprintf('线阻尼系数相对偏差：%.2f %%\n', abs(b1 - b2)/b1 * 100);
fprintf('===================================================\n');
fprintf('结论：偏差小于10%说明两种方法结果一致，参数可信\n\n');

%% ===================== 绘图验证 =====================
% 图1：峰值识别结果
figure('Color','w','Position',[100,100,900,400]);
plot(t, a, 'b-','LineWidth',1.2, 'DisplayName','原始加速度（无滤波）');
hold on;
plot(t_peaks, pks, 'ro','MarkerFaceColor','r','MarkerSize',6, 'DisplayName','识别有效峰值');
xlabel('时间 t (s)'); ylabel('加速度 a (m/s^2)');
title('对数衰减法：峰值识别结果');
legend('Location','best');
grid on;

% 图2：拟合效果对比
figure('Color','w','Position',[100,100,900,400]);
plot(t, a, 'b-','LineWidth',1.2, 'DisplayName','原始加速度（无滤波）');
hold on;
plot(t, a_fit, 'r--','LineWidth',1.5, 'DisplayName','最小二乘拟合曲线');
xlabel('时间 t (s)'); ylabel('加速度 a (m/s^2)');
title(['非线性最小二乘拟合效果（R² = ',num2str(R2,'%.4f'),'）']);
legend('Location','best');
grid on;