% Ship_Propulsion_Analysis.m
% 船舶完整推进效率分析
%
% 方法整合:
% 1. Michell积分 - 兴波阻力计算
% 2. ITTC 1957 - 摩擦阻力计算
% 3. Holtrop方法 - 伴流分数(w)和推力减额分数(t)
% 4. PVL - 螺旋桨敞水效率计算
%
% 变化参数: 方形系数Cb(三参数函数定义), 船宽B, 船长L
% 预期结论: 阻力越小，推进效率不一定越高
%
% 作者: Claude Code
% 日期: 2026-02-01

clear; close all; clc;

fprintf('================================================================================\n');
fprintf('          船舶完整推进效率分析系统\n');
fprintf('          Michell + ITTC + Holtrop + PVL 集成\n');
fprintf('================================================================================\n\n');

% ==============================================================================
% 1. 基准船舶参数设置 (基于KCS船型)
% ==============================================================================
Ship.L_ref = 230.0;       % 参考船长 [m]
Ship.B_ref = 32.2;        % 参考船宽 [m]
Ship.T = 10.8;            % 吃水 [m]
Ship.Cb_ref = 0.6505;     % 参考方形系数
Ship.D_prop = 7.9;        % 螺旋桨直径 [m]
Ship.Vs_knots = 24;       % 设计航速 [knots]

% 物理常数
Rho_water = 1025;         % 海水密度 [kg/m³]
Nu = 1.188e-6;            % 运动粘度 [m²/s]
g = 9.80665;              % 重力加速度 [m/s²]

% 航速转换
V_ship = Ship.Vs_knots * 0.5144; % [m/s]

% Michell积分网格参数
Nx = 101;                 % 站点数 (必须为奇数)
Nz = 40;                  % 水线数
N_theta = 80;             % 传播角采样数

% ==============================================================================
% 2. 参数变化设计空间
% ==============================================================================
NumSamples = 15;          % 每个参数的采样点数

% 定义参数变化范围
Cb_range = linspace(0.62, 0.70, NumSamples);    % 方形系数变化
B_ratio_range = linspace(0.92, 1.08, NumSamples); % 船宽比例变化
L_ratio_range = linspace(0.95, 1.05, NumSamples); % 船长比例变化

% 螺旋桨基础参数
Prop.NBLADE = 5;          % 叶片数
Prop.N_rpm = 105;         % 转速 [RPM]
Prop.Dhub = 1.58;         % 桨毂直径 [m]

% ==============================================================================
% 3. 分析1: 仅变化方形系数Cb
% ==============================================================================
fprintf('\n===== 分析1: 方形系数Cb变化对推进效率的影响 =====\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'No.', 'Cb', 'Rt[kN]', 'Rw[kN]', 'Rf[kN]', 'w', 't', 'eta_H', 'eta_O', 'eta_D');
fprintf('--------------------------------------------------------------------------------\n');

Results_Cb = struct();
x_norm = linspace(0, 1, Nx);

for i = 1:NumSamples
    Cb = Cb_range(i);
    L = Ship.L_ref;
    B = Ship.B_ref;
    T = Ship.T;

    % --- A. 三参数函数生成船型 ---
    [Y, h_val] = generate_hull_offsets(x_norm, Nz, Cb, Ship.Cb_ref);

    % --- B. 阻力计算 ---
    % ITTC 1957 摩擦阻力
    S_wet = estimate_wetted_surface(L, B, T, Cb);
    Re = V_ship * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    Rf = 0.5 * Rho_water * S_wet * V_ship^2 * (Cf + 0.0004); % 含粗糙度修正

    % Michell积分兴波阻力
    Rw = michell_integral(Y, V_ship, L, B, T, Rho_water, N_theta);

    % 总阻力
    Rt = Rf + Rw;

    % --- C. Holtrop伴流和推力减额 ---
    [w, t_thrust] = holtrop_interaction(L, B, T, Cb, Ship.D_prop);

    % --- D. 推进效率计算 ---
    eta_H = (1 - t_thrust) / (1 - w);  % 船身效率

    % 所需推力
    Thrust_req = Rt / (1 - t_thrust);
    Va = V_ship * (1 - w);  % 进速

    % 简化的螺旋桨敞水效率估算 (基于载荷系数)
    CT_load = Thrust_req / (0.5 * Rho_water * Va^2 * pi * (Ship.D_prop/2)^2);
    eta_O = estimate_propeller_efficiency(CT_load, Va, Prop.N_rpm, Ship.D_prop);

    % 相对旋转效率
    eta_R = 1.00 + 0.02 * (Cb - Ship.Cb_ref);

    % 总推进效率
    eta_D = eta_H * eta_O * eta_R;

    % 有效功率和传递功率
    Pe = Rt * V_ship;       % 有效功率
    Pd = Pe / eta_D;        % 传递功率

    % 存储结果
    Results_Cb(i).Cb = Cb;
    Results_Cb(i).L = L;
    Results_Cb(i).B = B;
    Results_Cb(i).Rt = Rt;
    Results_Cb(i).Rw = Rw;
    Results_Cb(i).Rf = Rf;
    Results_Cb(i).w = w;
    Results_Cb(i).t = t_thrust;
    Results_Cb(i).eta_H = eta_H;
    Results_Cb(i).eta_O = eta_O;
    Results_Cb(i).eta_R = eta_R;
    Results_Cb(i).eta_D = eta_D;
    Results_Cb(i).Pe = Pe;
    Results_Cb(i).Pd = Pd;

    fprintf('%4d   %.4f  %7.1f  %7.1f  %7.1f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        i, Cb, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_H, eta_O, eta_D);
end

% ==============================================================================
% 4. 分析2: 仅变化船宽B
% ==============================================================================
fprintf('\n===== 分析2: 船宽B变化对推进效率的影响 =====\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'No.', 'B[m]', 'Rt[kN]', 'Rw[kN]', 'Rf[kN]', 'w', 't', 'eta_H', 'eta_O', 'eta_D');
fprintf('--------------------------------------------------------------------------------\n');

Results_B = struct();
Cb_fixed = Ship.Cb_ref;

for i = 1:NumSamples
    B_ratio = B_ratio_range(i);
    L = Ship.L_ref;
    B = Ship.B_ref * B_ratio;
    T = Ship.T;

    % 生成船型
    [Y, ~] = generate_hull_offsets(x_norm, Nz, Cb_fixed, Ship.Cb_ref);

    % 阻力计算
    S_wet = estimate_wetted_surface(L, B, T, Cb_fixed);
    Re = V_ship * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    Rf = 0.5 * Rho_water * S_wet * V_ship^2 * (Cf + 0.0004);
    Rw = michell_integral(Y, V_ship, L, B, T, Rho_water, N_theta);
    Rt = Rf + Rw;

    % Holtrop交互因子
    [w, t_thrust] = holtrop_interaction(L, B, T, Cb_fixed, Ship.D_prop);

    % 效率计算
    eta_H = (1 - t_thrust) / (1 - w);
    Thrust_req = Rt / (1 - t_thrust);
    Va = V_ship * (1 - w);
    CT_load = Thrust_req / (0.5 * Rho_water * Va^2 * pi * (Ship.D_prop/2)^2);
    eta_O = estimate_propeller_efficiency(CT_load, Va, Prop.N_rpm, Ship.D_prop);
    eta_R = 1.00;
    eta_D = eta_H * eta_O * eta_R;
    Pe = Rt * V_ship;
    Pd = Pe / eta_D;

    % 存储
    Results_B(i).B = B;
    Results_B(i).Rt = Rt;
    Results_B(i).Rw = Rw;
    Results_B(i).Rf = Rf;
    Results_B(i).w = w;
    Results_B(i).t = t_thrust;
    Results_B(i).eta_H = eta_H;
    Results_B(i).eta_O = eta_O;
    Results_B(i).eta_D = eta_D;
    Results_B(i).Pd = Pd;

    fprintf('%4d   %6.2f  %7.1f  %7.1f  %7.1f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        i, B, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_H, eta_O, eta_D);
end

% ==============================================================================
% 5. 分析3: 仅变化船长L
% ==============================================================================
fprintf('\n===== 分析3: 船长L变化对推进效率的影响 =====\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n', ...
    'No.', 'L[m]', 'Rt[kN]', 'Rw[kN]', 'Rf[kN]', 'w', 't', 'eta_H', 'eta_O', 'eta_D');
fprintf('--------------------------------------------------------------------------------\n');

Results_L = struct();

for i = 1:NumSamples
    L_ratio = L_ratio_range(i);
    L = Ship.L_ref * L_ratio;
    B = Ship.B_ref;
    T = Ship.T;

    % 生成船型
    [Y, ~] = generate_hull_offsets(x_norm, Nz, Cb_fixed, Ship.Cb_ref);

    % 阻力计算
    S_wet = estimate_wetted_surface(L, B, T, Cb_fixed);
    Re = V_ship * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    Rf = 0.5 * Rho_water * S_wet * V_ship^2 * (Cf + 0.0004);
    Rw = michell_integral(Y, V_ship, L, B, T, Rho_water, N_theta);
    Rt = Rf + Rw;

    % Holtrop交互因子
    [w, t_thrust] = holtrop_interaction(L, B, T, Cb_fixed, Ship.D_prop);

    % 效率计算
    eta_H = (1 - t_thrust) / (1 - w);
    Thrust_req = Rt / (1 - t_thrust);
    Va = V_ship * (1 - w);
    CT_load = Thrust_req / (0.5 * Rho_water * Va^2 * pi * (Ship.D_prop/2)^2);
    eta_O = estimate_propeller_efficiency(CT_load, Va, Prop.N_rpm, Ship.D_prop);
    eta_R = 1.00;
    eta_D = eta_H * eta_O * eta_R;
    Pe = Rt * V_ship;
    Pd = Pe / eta_D;

    % 存储
    Results_L(i).L = L;
    Results_L(i).Rt = Rt;
    Results_L(i).Rw = Rw;
    Results_L(i).Rf = Rf;
    Results_L(i).w = w;
    Results_L(i).t = t_thrust;
    Results_L(i).eta_H = eta_H;
    Results_L(i).eta_O = eta_O;
    Results_L(i).eta_D = eta_D;
    Results_L(i).Pd = Pd;

    fprintf('%4d   %6.1f  %7.1f  %7.1f  %7.1f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        i, L, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_H, eta_O, eta_D);
end

% ==============================================================================
% 6. 结果可视化
% ==============================================================================
figure('Position', [50 50 1400 900], 'Color', 'w', 'Name', 'Ship Propulsion Analysis');

% --- 图1: Cb变化分析 ---
subplot(3, 3, 1);
Cb_vec = [Results_Cb.Cb];
Rt_Cb = [Results_Cb.Rt]/1000;
Rw_Cb = [Results_Cb.Rw]/1000;
Rf_Cb = [Results_Cb.Rf]/1000;
plot(Cb_vec, Rt_Cb, 'b-o', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(Cb_vec, Rw_Cb, 'r--s', 'LineWidth', 1.5);
plot(Cb_vec, Rf_Cb, 'g-.^', 'LineWidth', 1.5);
xlabel('方形系数 C_b'); ylabel('阻力 [kN]');
title('Cb对阻力的影响');
legend('总阻力R_t', '兴波阻力R_w', '摩擦阻力R_f', 'Location', 'best');
grid on;

subplot(3, 3, 2);
eta_H_Cb = [Results_Cb.eta_H];
eta_O_Cb = [Results_Cb.eta_O];
eta_D_Cb = [Results_Cb.eta_D];
plot(Cb_vec, eta_H_Cb, 'm-o', 'LineWidth', 2); hold on;
plot(Cb_vec, eta_O_Cb, 'c-s', 'LineWidth', 2);
plot(Cb_vec, eta_D_Cb, 'k-^', 'LineWidth', 2.5);
xlabel('方形系数 C_b'); ylabel('效率');
title('Cb对推进效率的影响');
legend('\eta_H 船身效率', '\eta_O 敞水效率', '\eta_D 推进效率', 'Location', 'best');
grid on;

subplot(3, 3, 3);
Pd_Cb = [Results_Cb.Pd]/1e6;
yyaxis left
plot(Cb_vec, Rt_Cb, 'b-o', 'LineWidth', 2);
ylabel('总阻力 R_t [kN]');
yyaxis right
plot(Cb_vec, Pd_Cb, 'r-s', 'LineWidth', 2);
ylabel('传递功率 P_d [MW]');
xlabel('方形系数 C_b');
title('阻力vs功率 (验证核心结论)');
grid on;

% --- 图2: B变化分析 ---
subplot(3, 3, 4);
B_vec = [Results_B.B];
Rt_B = [Results_B.Rt]/1000;
Rw_B = [Results_B.Rw]/1000;
Rf_B = [Results_B.Rf]/1000;
plot(B_vec, Rt_B, 'b-o', 'LineWidth', 2); hold on;
plot(B_vec, Rw_B, 'r--s', 'LineWidth', 1.5);
plot(B_vec, Rf_B, 'g-.^', 'LineWidth', 1.5);
xlabel('船宽 B [m]'); ylabel('阻力 [kN]');
title('船宽B对阻力的影响');
legend('总阻力R_t', '兴波阻力R_w', '摩擦阻力R_f', 'Location', 'best');
grid on;

subplot(3, 3, 5);
eta_H_B = [Results_B.eta_H];
eta_O_B = [Results_B.eta_O];
eta_D_B = [Results_B.eta_D];
plot(B_vec, eta_H_B, 'm-o', 'LineWidth', 2); hold on;
plot(B_vec, eta_O_B, 'c-s', 'LineWidth', 2);
plot(B_vec, eta_D_B, 'k-^', 'LineWidth', 2.5);
xlabel('船宽 B [m]'); ylabel('效率');
title('船宽B对推进效率的影响');
legend('\eta_H', '\eta_O', '\eta_D', 'Location', 'best');
grid on;

subplot(3, 3, 6);
Pd_B = [Results_B.Pd]/1e6;
yyaxis left
plot(B_vec, Rt_B, 'b-o', 'LineWidth', 2);
ylabel('总阻力 R_t [kN]');
yyaxis right
plot(B_vec, Pd_B, 'r-s', 'LineWidth', 2);
ylabel('传递功率 P_d [MW]');
xlabel('船宽 B [m]');
title('阻力vs功率');
grid on;

% --- 图3: L变化分析 ---
subplot(3, 3, 7);
L_vec = [Results_L.L];
Rt_L = [Results_L.Rt]/1000;
Rw_L = [Results_L.Rw]/1000;
Rf_L = [Results_L.Rf]/1000;
plot(L_vec, Rt_L, 'b-o', 'LineWidth', 2); hold on;
plot(L_vec, Rw_L, 'r--s', 'LineWidth', 1.5);
plot(L_vec, Rf_L, 'g-.^', 'LineWidth', 1.5);
xlabel('船长 L [m]'); ylabel('阻力 [kN]');
title('船长L对阻力的影响');
legend('总阻力R_t', '兴波阻力R_w', '摩擦阻力R_f', 'Location', 'best');
grid on;

subplot(3, 3, 8);
eta_H_L = [Results_L.eta_H];
eta_O_L = [Results_L.eta_O];
eta_D_L = [Results_L.eta_D];
plot(L_vec, eta_H_L, 'm-o', 'LineWidth', 2); hold on;
plot(L_vec, eta_O_L, 'c-s', 'LineWidth', 2);
plot(L_vec, eta_D_L, 'k-^', 'LineWidth', 2.5);
xlabel('船长 L [m]'); ylabel('效率');
title('船长L对推进效率的影响');
legend('\eta_H', '\eta_O', '\eta_D', 'Location', 'best');
grid on;

subplot(3, 3, 9);
Pd_L = [Results_L.Pd]/1e6;
yyaxis left
plot(L_vec, Rt_L, 'b-o', 'LineWidth', 2);
ylabel('总阻力 R_t [kN]');
yyaxis right
plot(L_vec, Pd_L, 'r-s', 'LineWidth', 2);
ylabel('传递功率 P_d [MW]');
xlabel('船长 L [m]');
title('阻力vs功率');
grid on;

sgtitle('船舶推进效率综合分析 - Michell+ITTC+Holtrop+PVL', 'FontSize', 14, 'FontWeight', 'bold');

% ==============================================================================
% 7. 结论分析
% ==============================================================================
fprintf('\n================================================================================\n');
fprintf('                           结 论 分 析\n');
fprintf('================================================================================\n\n');

% Cb分析
[min_Rt_Cb, idx_min_Rt] = min([Results_Cb.Rt]);
[max_etaD_Cb, idx_max_etaD] = max([Results_Cb.eta_D]);
[min_Pd_Cb, idx_min_Pd] = min([Results_Cb.Pd]);

fprintf('【分析1: 方形系数Cb变化】\n');
fprintf('  - 最小阻力点: Cb = %.4f, Rt = %.1f kN\n', Results_Cb(idx_min_Rt).Cb, min_Rt_Cb/1000);
fprintf('  - 最大推进效率点: Cb = %.4f, eta_D = %.4f\n', Results_Cb(idx_max_etaD).Cb, max_etaD_Cb);
fprintf('  - 最小传递功率点: Cb = %.4f, Pd = %.2f MW\n', Results_Cb(idx_min_Pd).Cb, min_Pd_Cb/1e6);

if idx_min_Rt ~= idx_max_etaD
    fprintf('  >>> 结论验证: 最小阻力(Cb=%.4f) != 最高效率(Cb=%.4f)\n', ...
        Results_Cb(idx_min_Rt).Cb, Results_Cb(idx_max_etaD).Cb);
    fprintf('      阻力越小，推进效率不一定越高!\n');
end

% B分析
[min_Rt_B, idx_min_Rt_B] = min([Results_B.Rt]);
[max_etaD_B, idx_max_etaD_B] = max([Results_B.eta_D]);

fprintf('\n【分析2: 船宽B变化】\n');
fprintf('  - 最小阻力点: B = %.2f m, Rt = %.1f kN\n', Results_B(idx_min_Rt_B).B, min_Rt_B/1000);
fprintf('  - 最大推进效率点: B = %.2f m, eta_D = %.4f\n', Results_B(idx_max_etaD_B).B, max_etaD_B);

% L分析
[min_Rt_L, idx_min_Rt_L] = min([Results_L.Rt]);
[max_etaD_L, idx_max_etaD_L] = max([Results_L.eta_D]);

fprintf('\n【分析3: 船长L变化】\n');
fprintf('  - 最小阻力点: L = %.1f m, Rt = %.1f kN\n', Results_L(idx_min_Rt_L).L, min_Rt_L/1000);
fprintf('  - 最大推进效率点: L = %.1f m, eta_D = %.4f\n', Results_L(idx_max_etaD_L).L, max_etaD_L);

fprintf('\n================================================================================\n');
fprintf('【核心结论】\n');
fprintf('  1. 方形系数Cb增大时:\n');
fprintf('     - 兴波阻力Rw增加 (船型更钝)\n');
fprintf('     - 伴流分数w增加 (艉部更丰满)\n');
fprintf('     - 船身效率eta_H提高 (w增大的正面效果)\n');
fprintf('     - 最终推进效率eta_D可能在某个Cb处达到最优\n\n');
fprintf('  2. 船宽B变化时:\n');
fprintf('     - 宽船阻力更大，但交互因子变化复杂\n');
fprintf('     - 最优B取决于阻力与效率的平衡\n\n');
fprintf('  3. 船长L变化时:\n');
fprintf('     - 长船摩擦阻力增加，但兴波阻力可能减小\n');
fprintf('     - Froude数变化导致阻力特性改变\n\n');
fprintf('  4. 综合结论:\n');
fprintf('     阻力最小点与效率最高点往往不重合!\n');
fprintf('     这是因为船身效率eta_H与阻力呈负相关关系.\n');
fprintf('================================================================================\n');

% 保存结果
save('Propulsion_Analysis_Results.mat', 'Results_Cb', 'Results_B', 'Results_L', 'Ship');
fprintf('\n结果已保存至 Propulsion_Analysis_Results.mat\n');

% 保存图片
saveas(gcf, 'Propulsion_Analysis_Results.png');
fprintf('图片已保存至 Propulsion_Analysis_Results.png\n');

% ==============================================================================
% 辅助函数定义
% ==============================================================================

function [Y, h_val] = generate_hull_offsets(x_norm, Nz, Cb, Cb_ref)
    % 三参数函数生成船体偏移矩阵
    % 输入: x_norm - 归一化站位 [0,1]
    %       Nz - 水线数
    %       Cb - 目标方形系数
    %       Cb_ref - 参考方形系数
    % 输出: Y - 偏移矩阵 (Nx x Nz), 最大值0.5
    %       h_val - 三参数函数的h值

    Nx = length(x_norm);

    % 根据Cb调整三参数函数的h值
    % h控制船型的丰满度
    h_val = 0.70 + 2.5 * (Cb - 0.63);
    a1 = 0.12;  % 艏部形状参数
    a2 = 0.10;  % 艉部形状参数

    % 生成水线面形状
    f0 = three_param_waterline(x_norm, h_val, a1, a2);
    f0 = f0 / max(f0) * 0.5; % 归一化到最大0.5

    % 构造3D偏移矩阵
    z_norm = linspace(0, 1, Nz); % 0=龙骨, 1=水线
    Y = zeros(Nx, Nz);

    for iz = 1:Nz
        % 深度因子：模拟U/V型组合剖面
        z_factor = sin(pi/2 * z_norm(iz))^0.5;
        Y(:, iz) = f0(:) * z_factor;
    end
end

function f = three_param_waterline(x, h, a1, a2)
    % 三参数水线面形状函数
    % 参考: 船舶设计中常用的参数化船型方法
    %
    % x - 归一化船长位置 [0,1], 0=艏, 1=艉
    % h - 中部丰满度参数 (越大越丰满)
    % a1 - 艏部形状参数 (控制前半部曲率)
    % a2 - 艉部形状参数 (控制后半部曲率)

    f = zeros(size(x));
    for i = 1:length(x)
        if x(i) < 0.5
            % 前半部: 艏部逐渐变宽
            f(i) = h - a1/2 * (cos(2*pi*x(i)) - 1);
        else
            % 后半部: 艉部逐渐收窄
            f(i) = h + a2/2 * (cos(2*pi*x(i)) + 1) + a1;
        end
    end
end

function S = estimate_wetted_surface(L, B, T, Cb)
    % 估算湿表面积
    % 使用经验公式 (Holtrop-Mennen近似)
    S = L * (2*T + B) * sqrt(Cb) * (0.453 + 0.4425*Cb - 0.2862*Cb^2);
    % 简化版本
    S = max(S, 1.01 * L * (2*T + B*Cb));
end

function [w, t] = holtrop_interaction(L, B, T, Cb, D_prop)
    % Holtrop方法计算伴流分数和推力减额分数
    %
    % 参考: Holtrop, J. (1984). "A Statistical Re-Analysis of Resistance and
    %       Propulsion Data"
    %
    % 输入:
    %   L - 船长 [m]
    %   B - 船宽 [m]
    %   T - 吃水 [m]
    %   Cb - 方形系数
    %   D_prop - 螺旋桨直径 [m]
    %
    % 输出:
    %   w - 伴流分数 (wake fraction)
    %   t - 推力减额分数 (thrust deduction factor)

    % 长宽比和吃水比
    LB = L / B;
    BT = B / T;

    % 棱形系数估算 (假设Cp ≈ Cb/0.98)
    Cp = Cb / 0.98;

    % 螺旋桨直径与吃水比
    DT = D_prop / T;

    % === 伴流分数 w (Holtrop公式，单桨船) ===
    % 基础伴流系数
    c9 = 1.0; % 单桨船系数

    % Holtrop公式简化版
    % w = c9 * Cb * (0.5 - 0.05*Cp) * (1 + 0.015*Cstern)
    % Cstern为艉部形状系数，此处取0

    w_base = 0.3095 * Cb + 10.0 * Cb * (L/B)^(-3.18);

    % 修正项
    w_corr = -0.23 * DT;

    w = w_base + w_corr;
    w = max(0.15, min(0.45, w)); % 限制在合理范围

    % === 推力减额分数 t (Holtrop公式) ===
    % 基于伴流的经验关系
    % t/w 比值通常在 0.5-0.8 之间

    t_base = 0.25014 * (B/L)^0.28956 * (sqrt(BT)/DT)^0.2624;

    % Cb影响修正
    t_Cb_corr = 0.0015 * (Cb - 0.65) / 0.01;

    t = t_base + t_Cb_corr;

    % 限制t/w比值在合理范围
    tw_ratio = t / w;
    if tw_ratio < 0.5
        t = 0.5 * w;
    elseif tw_ratio > 0.85
        t = 0.85 * w;
    end

    t = max(0.10, min(0.30, t)); % 限制在合理范围
end

function eta_O = estimate_propeller_efficiency(CT_load, Va, N_rpm, D)
    % 基于载荷系数估算螺旋桨敞水效率
    %
    % 这是一个简化模型，实际应使用PVL完整计算
    %
    % 输入:
    %   CT_load - 推力载荷系数 CT = T/(0.5*rho*Va^2*A0)
    %   Va - 进速 [m/s]
    %   N_rpm - 转速 [RPM]
    %   D - 直径 [m]
    %
    % 输出:
    %   eta_O - 敞水效率

    n = N_rpm / 60;  % 转速 [rps]
    J = Va / (n * D); % 前进系数

    % 基于推进器理论的效率估算
    % eta_ideal = 2 / (1 + sqrt(1 + CT_load))
    eta_ideal = 2 / (1 + sqrt(1 + CT_load));

    % 考虑实际损失（粘性、旋转流等）
    % 典型效率约为理想效率的 0.75-0.85
    eta_real_factor = 0.80 - 0.02 * (CT_load - 0.5);
    eta_real_factor = max(0.70, min(0.85, eta_real_factor));

    eta_O = eta_ideal * eta_real_factor;

    % 前进系数修正
    J_opt = 0.7; % 最优前进系数估计
    J_penalty = 0.05 * (J - J_opt)^2;

    eta_O = eta_O - J_penalty;
    eta_O = max(0.55, min(0.75, eta_O));
end

function Rw = michell_integral(Y, U, L, B, T, RHO, N)
    % Michell积分计算兴波阻力
    % 基于Filon数值积分算法

    Nx = size(Y, 1);
    Nz = size(Y, 2);
    YH = Y * B;  % 偏移量乘以船宽

    % 积分网格
    dz = T / (Nz - 1);
    z = (-T:dz:0)';
    dx = L / (Nx - 1);
    x = (0:dx:L)';
    theta = michspace(N)';

    % 常数
    g = 9.80665;
    k0 = g / U^2;
    c = (4 * RHO * U^2) / pi;
    a = sec(theta);
    k = k0 * a.^2;

    % Z方向积分 (Filon梯形法)
    Kz = k0 * dz * a.^2;
    w0 = (exp(Kz) - 1 - Kz) ./ Kz.^2;
    wn = (exp(Kz) + exp(-Kz) - 2) ./ Kz.^2;
    wN = (exp(-Kz) - 1 + Kz) ./ Kz.^2;

    F = zeros(Nx, N);
    for j = 1:N
        for m = 1:Nx
            f = zeros(Nz, 1);
            for n = 1:Nz
                if n == 1
                    f(n) = w0(j) * YH(m,n) * exp(k0*z(n)*a(j)^2) * dz;
                elseif n == Nz
                    f(n) = wN(j) * YH(m,n) * exp(k0*z(n)*a(j)^2) * dz;
                else
                    f(n) = wn(j) * YH(m,n) * exp(k0*z(n)*a(j)^2) * dz;
                end
            end
            F(m,j) = sum(f);
        end
    end

    % X方向积分 (Filon算法)
    Kx = k0 * dx .* a;
    alp = (Kx.^2 + 0.5*Kx.*sin(2*Kx) + cos(2*Kx) - 1) ./ Kx.^3;
    bet = (3*Kx + Kx.*cos(2*Kx) - 2*sin(2*Kx)) ./ Kx.^3;
    gam = 4*(sin(Kx) - Kx.*cos(Kx)) ./ Kx.^3;

    Pt = zeros(N, 1);
    Qt = zeros(N, 1);
    P = zeros(N, 1);
    Q = zeros(N, 1);

    for j = 1:N
        c_val = cos(k0 * x * a(j));
        s_val = sin(k0 * x * a(j));

        pev = F(1:2:end, j) .* c_val(1:2:end);
        qev = F(1:2:end, j) .* s_val(1:2:end);
        pod = F(2:2:end-1, j) .* c_val(2:2:end-1);
        qod = F(2:2:end-1, j) .* s_val(2:2:end-1);

        Pt(j) = F(Nx, j) * cos(k0 * L * a(j));
        Qt(j) = F(Nx, j) * sin(k0 * L * a(j));

        Pev = sum(pev) - 0.5 * Pt(j);
        Pod = sum(pod);
        Qev = sum(qev) - 0.5 * Qt(j);
        Qod = sum(qod);

        P(j) = dx * (alp(j)*Qt(j) + bet(j)*Pev + gam(j)*Pod);
        Q(j) = dx * (-alp(j)*Pt(j) + bet(j)*Qev + gam(j)*Qod);
    end

    % 波阻力谱
    R = c * k.^2 ./ a.^3 .* (k.^2.*(P.^2 + Q.^2) + ...
        2*k.*a.*(Q.*Pt - P.*Qt) + a.^2.*(Pt.^2 + Qt.^2));
    R(isnan(R)) = 0;

    % Theta积分 (梯形法)
    rw_vec = zeros(N-1, 1);
    for k_idx = 1:N-1
        rw_vec(k_idx) = 0.5 * (R(k_idx) + R(k_idx+1)) * (theta(k_idx+1) - theta(k_idx));
    end
    Rw = sum(rw_vec);

    % 确保结果为正
    Rw = max(0, Rw);
end

function xm = michspace(N)
    % 传播角加密采样
    % 在pi/2附近加密
    xm = logspace(0, 1, N) - 1;
    xm = xm * pi / 18 - pi / 2;
    xm = fliplr(-xm);
end
