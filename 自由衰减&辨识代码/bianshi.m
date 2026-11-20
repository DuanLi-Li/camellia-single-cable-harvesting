%% ===================== 枝条等效参数辨识（非线性最小二乘） =====================
% 模型：m*ddx + c*dx + k*x = 0
% 方法：对自由衰减加速度记录做非线性最小二乘拟合，输出等效质量 m、等效阻尼 c、
%       固有频率 fn、阻尼比 zeta 与拟合优度 R^2。
% 窗口：以记录中的最高峰为起点，向后取足够长的数据窗（默认 2 s，或至记录末尾）；
%       峰值仅用于确定窗口与初值，阻尼比初值取 2%，最终参数全部由最小二乘拟合得到。
% 预处理：仅去除直流零点，不滤波。
%
% 用法：修改下面的 dataFile 与 k，运行即可。
%   近根部   k = 8667.10，文件 sdd_928d_1~3.mat
%   中间位置 k = 3917.00，文件 sdd_928z_1~3.mat
%   近自由端 k = 2290.35，文件 sdd_928g_1~3.mat
clear; clc; close all;

%% ===================== 1. 数据与已知刚度 =====================
dataFile   = 'sdd_928g_3.mat';
k          = 2290.35;      % 该夹持位置的等效弯曲刚度 (N/m)
win_length = 2.5;          % 拟合窗口长度 (s)，或至记录末尾

S = load(dataFile);
fn_field = fieldnames(S);
st = S.(fn_field{1});
t_all = st.X.Data(:)';
a_all = st.Y.Data(:)';

fprintf('==================== 输入 ====================\n');
fprintf('数据文件：%s\n', dataFile);
fprintf('线刚度 k = %.2f N/m\n', k);
fprintf('==============================================\n\n');

%% ===================== 2. 以最高峰为起点截取窗口 =====================
fs   = 1/(t_all(2)-t_all(1));
a_dc = a_all - mean(a_all);

[pks, locs] = findpeaks(a_dc, 'MinPeakHeight', 0.2*max(a_dc), ...
                             'MinPeakDistance', round(0.07*fs));
if numel(pks) < 5
    error('有效峰值过少，请检查数据或峰值阈值');
end
[~, imax] = max(pks);
t0 = t_all(locs(imax));                 % 最高峰所在时刻 = 窗口起点
t1 = min(t_all(end), t0 + win_length);  % 向后取足够长的窗口

w  = t_all >= t0 & t_all <= t1;
t  = t_all(w) - t0;                     % 以最高峰为零时刻
a  = a_all(w);
a  = a - mean(a);                       % 再次去直流
n  = numel(t);

fprintf('最高峰时刻 t0 = %.4f s，拟合窗口 %.2f ~ %.2f s（%d 点）\n', t0, t0, t1, n);
fprintf('采样率 fs = %.1f Hz\n\n', fs);

%% ===================== 3. 由峰值特征给出初值 =====================
t_peaks  = t_all(locs);
Td       = mean(diff(t_peaks));            % 平均峰间距
A0       = pks(imax);                      % 最高峰幅值
zeta0    = 0.02;                           % 阻尼比初值（2%）
omega_n0 = 2*pi/Td/sqrt(1-zeta0^2);        % 无阻尼固有角频率初值
phi0     = 0;                              % 起点为峰，相位初值取 0

fprintf('平均峰间距 Td = %.4f s\n', Td);
fprintf('初值：A0 = %.3f, zeta0 = %.3f, fn0 = %.3f Hz\n\n', A0, zeta0, omega_n0/(2*pi));

%% ===================== 4. 非线性最小二乘拟合 =====================
acc_model = @(p, tt) p(1).*exp(-p(2)*p(3)*tt).*cos(p(3)*sqrt(1-p(2)^2)*tt + p(4));
params0   = [A0, zeta0, omega_n0, phi0];

options = optimset('MaxIter', 30000, 'MaxFunEvals', 100000, ...
                   'TolX', 1e-12, 'TolFun', 1e-12, 'Display', 'off');

params_fit = fminsearch(@(p) sum((a - acc_model(p, t)).^2), params0, options);

A_fit   = params_fit(1);
zeta    = abs(params_fit(2));
omega_n = params_fit(3);

m = k / omega_n^2;
c = 2*zeta*sqrt(m*k);

a_fit  = acc_model(params_fit, t);
SS_res = sum((a - a_fit).^2);
SS_tot = sum((a - mean(a)).^2);
R2     = 1 - SS_res/SS_tot;

fprintf('==================== 辨识结果 ====================\n');
fprintf('固有频率     fn = %.3f Hz\n', omega_n/(2*pi));
fprintf('阻尼比      zeta = %.4f (%.2f%%)\n', zeta, 100*zeta);
fprintf('等效质量       m = %.4f kg\n', m);
fprintf('等效阻尼       c = %.4f N·s/m\n', c);
fprintf('拟合优度      R^2 = %.2f%%\n', 100*R2);
fprintf('==================================================\n');

%% ===================== 5. 绘图 =====================
figure('Color','w','Position',[100,100,900,400]);
plot(t, a, 'b-', 'LineWidth', 1.2, 'DisplayName', '实测加速度');
hold on;
plot(t, a_fit, 'r--', 'LineWidth', 1.4, 'DisplayName', '最小二乘拟合');
xlabel('相对最高峰的时间 t (s)'); ylabel('加速度 a (m/s^2)');
title(sprintf('非线性最小二乘拟合（R^2 = %.2f%%）', 100*R2));
legend('Location', 'best'); grid on;

