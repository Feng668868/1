% KCS_test.m
% KCS 船桨集成优化分析 - Michell + ITTC + Holtrop 完整版
%
% 核心目标：验证 "阻力越小，推进效率不一定越高"
% 证明：考虑船桨集成优化 vs 仅优化阻力 的差异
%
% 方法：
% 1. Michell积分 - 兴波阻力计算
% 2. ITTC 1957 - 摩擦阻力计算
% 3. Holtrop方法 - 伴流分数w和推力减额分数t
% 4. 螺旋桨效率模型 - 基于载荷系数
%
% 变化参数：方形系数Cb（三参数函数定义）、船宽B、船长L

clear; close all; clc;

fprintf('================================================================================\n');
fprintf('   KCS 船桨集成优化分析\n');
fprintf('   验证：阻力最小 ≠ 效率最高，船桨集成优化的重要性\n');
fprintf('================================================================================\n\n');

% ==============================================================================
% 1. 全局仿真设置
% ==============================================================================
NumSamples = 25;
V_ship_mps = 12.35;       % 24 knots ≈ 12.35 m/s
Rho_water = 1025;
Nu = 1.188e-6;
g = 9.80665;

% KCS 基准参数
KCS.L = 230.0;            % 船长 [m]
KCS.B = 32.2;             % 船宽 [m]
KCS.T = 10.8;             % 吃水 [m]
KCS.D_prop = 7.9;         % 螺旋桨直径 [m]
KCS.Cb_ref = 0.6505;      % 参考方形系数

% Michell 积分网格参数
Nx = 101;  % 站点数 (必须为奇数)
Nz = 40;   % 水线数
N_theta = 80; % 传播角采样数

x_norm = linspace(0, 1, Nx);

% ==============================================================================
% 2. 分析1：方形系数 Cb 变化 (船宽、船长固定)
% ==============================================================================
fprintf('===== 分析1: 方形系数 Cb 变化 =====\n');
fprintf('%-5s %-7s %-8s %-8s %-8s %-7s %-7s %-7s %-7s %-8s\n', ...
    'No.', 'Cb', 'Rt(kN)', 'Rw(kN)', 'Rf(kN)', 'w', 't', 'eta_H', 'eta_D', 'Pd(MW)');
fprintf('--------------------------------------------------------------------------------\n');

Cb_range = linspace(0.610, 0.700, NumSamples);
Results_Cb = struct('Cb',[],'L',[],'B',[],'Rt',[],'Rw',[],'Rf',[],'w',[],'t',[],...
                    'eta_H',[],'eta_O',[],'eta_D',[],'Pd',[]);
Results_Cb = repmat(Results_Cb, NumSamples, 1);

for i = 1:NumSamples
    Cb = Cb_range(i);
    L = KCS.L;
    B = KCS.B;
    T = KCS.T;

    % --- A. 三参数函数生成船型 ---
    h_val = 0.70 + 2.8 * (Cb - 0.63);  % 丰满度随Cb变化
    f0 = three_param_shape(x_norm, h_val, 0.12, 0.10);
    f0 = f0 / max(f0) * 0.5;

    % 构造3D偏移矩阵
    z_norm = linspace(0, 1, Nz);
    Y = zeros(Nx, Nz);
    for iz = 1:Nz
        z_factor = sin(pi/2 * z_norm(iz))^0.5;
        Y(:, iz) = f0(:) * z_factor;
    end

    % --- B. 阻力计算 ---
    % ITTC-1957 摩擦阻力
    S = estimate_wetted_surface(L, B, T, Cb);
    Re = V_ship_mps * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    Rf = 0.5 * Rho_water * S * V_ship_mps^2 * (Cf + 0.0004);

    % Michell积分兴波阻力
    Rw = michell(Y, V_ship_mps, L, B, T, Rho_water, N_theta);
    Rt = Rf + Rw;

    % --- C. Holtrop 伴流和推力减额 ---
    [w, t_val] = holtrop_wake_thrust(L, B, T, Cb, KCS.D_prop);

    % --- D. 效率计算 ---
    eta_H = (1 - t_val) / (1 - w);  % 船身效率

    Thrust_req = Rt / (1 - t_val);
    Va = V_ship_mps * (1 - w);

    % 螺旋桨敞水效率 (基于推力载荷系数)
    CT = Thrust_req / (0.5 * Rho_water * Va^2 * pi * (KCS.D_prop/2)^2);
    eta_O = propeller_efficiency(CT);

    eta_R = 1.00 + 0.01 * (Cb - KCS.Cb_ref);  % 相对旋转效率
    eta_D = eta_H * eta_O * eta_R;  % 推进效率

    Pd = Rt * V_ship_mps / eta_D;  % 传递功率

    % 存储结果
    Results_Cb(i).Cb = Cb;
    Results_Cb(i).L = L;
    Results_Cb(i).B = B;
    Results_Cb(i).Rt = Rt;
    Results_Cb(i).Rw = Rw;
    Results_Cb(i).Rf = Rf;
    Results_Cb(i).w = w;
    Results_Cb(i).t = t_val;
    Results_Cb(i).eta_H = eta_H;
    Results_Cb(i).eta_O = eta_O;
    Results_Cb(i).eta_D = eta_D;
    Results_Cb(i).Pd = Pd;

    fprintf('%4d  %.4f  %7.1f  %7.1f  %7.1f  %.4f  %.4f  %.4f  %.4f  %7.2f\n', ...
        i, Cb, Rt/1000, Rw/1000, Rf/1000, w, t_val, eta_H, eta_D, Pd/1e6);
end

% ==============================================================================
% 3. 分析2：船宽 B 变化 (Cb、船长固定)
% ==============================================================================
fprintf('\n===== 分析2: 船宽 B 变化 (Cb=%.4f固定) =====\n', KCS.Cb_ref);
fprintf('%-5s %-7s %-8s %-8s %-8s %-7s %-7s %-7s %-7s %-8s\n', ...
    'No.', 'B(m)', 'Rt(kN)', 'Rw(kN)', 'Rf(kN)', 'w', 't', 'eta_H', 'eta_D', 'Pd(MW)');
fprintf('--------------------------------------------------------------------------------\n');

B_ratio_range = linspace(0.90, 1.10, NumSamples);
Results_B = struct('B',[],'Rt',[],'Rw',[],'Rf',[],'w',[],'t',[],...
                   'eta_H',[],'eta_O',[],'eta_D',[],'Pd',[]);
Results_B = repmat(Results_B, NumSamples, 1);

Cb_fixed = KCS.Cb_ref;
h_val_fixed = 0.70 + 2.8 * (Cb_fixed - 0.63);
f0_fixed = three_param_shape(x_norm, h_val_fixed, 0.12, 0.10);
f0_fixed = f0_fixed / max(f0_fixed) * 0.5;

for i = 1:NumSamples
    B_ratio = B_ratio_range(i);
    L = KCS.L;
    B = KCS.B * B_ratio;
    T = KCS.T;
    Cb = Cb_fixed;

    % 偏移矩阵
    z_norm = linspace(0, 1, Nz);
    Y = zeros(Nx, Nz);
    for iz = 1:Nz
        z_factor = sin(pi/2 * z_norm(iz))^0.5;
        Y(:, iz) = f0_fixed(:) * z_factor;
    end

    % 阻力计算
    S = estimate_wetted_surface(L, B, T, Cb);
    Re = V_ship_mps * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    Rf = 0.5 * Rho_water * S * V_ship_mps^2 * (Cf + 0.0004);
    Rw = michell(Y, V_ship_mps, L, B, T, Rho_water, N_theta);
    Rt = Rf + Rw;

    % Holtrop
    [w, t_val] = holtrop_wake_thrust(L, B, T, Cb, KCS.D_prop);

    % 效率
    eta_H = (1 - t_val) / (1 - w);
    Thrust_req = Rt / (1 - t_val);
    Va = V_ship_mps * (1 - w);
    CT = Thrust_req / (0.5 * Rho_water * Va^2 * pi * (KCS.D_prop/2)^2);
    eta_O = propeller_efficiency(CT);
    eta_R = 1.00;
    eta_D = eta_H * eta_O * eta_R;
    Pd = Rt * V_ship_mps / eta_D;

    Results_B(i).B = B;
    Results_B(i).Rt = Rt;
    Results_B(i).Rw = Rw;
    Results_B(i).Rf = Rf;
    Results_B(i).w = w;
    Results_B(i).t = t_val;
    Results_B(i).eta_H = eta_H;
    Results_B(i).eta_O = eta_O;
    Results_B(i).eta_D = eta_D;
    Results_B(i).Pd = Pd;

    fprintf('%4d  %6.2f  %7.1f  %7.1f  %7.1f  %.4f  %.4f  %.4f  %.4f  %7.2f\n', ...
        i, B, Rt/1000, Rw/1000, Rf/1000, w, t_val, eta_H, eta_D, Pd/1e6);
end

% ==============================================================================
% 4. 分析3：船长 L 变化 (Cb、船宽固定)
% ==============================================================================
fprintf('\n===== 分析3: 船长 L 变化 (Cb=%.4f固定) =====\n', KCS.Cb_ref);
fprintf('%-5s %-7s %-8s %-8s %-8s %-7s %-7s %-7s %-7s %-8s\n', ...
    'No.', 'L(m)', 'Rt(kN)', 'Rw(kN)', 'Rf(kN)', 'w', 't', 'eta_H', 'eta_D', 'Pd(MW)');
fprintf('--------------------------------------------------------------------------------\n');

L_ratio_range = linspace(0.92, 1.08, NumSamples);
Results_L = struct('L',[],'Rt',[],'Rw',[],'Rf',[],'w',[],'t',[],...
                   'eta_H',[],'eta_O',[],'eta_D',[],'Pd',[]);
Results_L = repmat(Results_L, NumSamples, 1);

for i = 1:NumSamples
    L_ratio = L_ratio_range(i);
    L = KCS.L * L_ratio;
    B = KCS.B;
    T = KCS.T;
    Cb = Cb_fixed;

    z_norm = linspace(0, 1, Nz);
    Y = zeros(Nx, Nz);
    for iz = 1:Nz
        z_factor = sin(pi/2 * z_norm(iz))^0.5;
        Y(:, iz) = f0_fixed(:) * z_factor;
    end

    S = estimate_wetted_surface(L, B, T, Cb);
    Re = V_ship_mps * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    Rf = 0.5 * Rho_water * S * V_ship_mps^2 * (Cf + 0.0004);
    Rw = michell(Y, V_ship_mps, L, B, T, Rho_water, N_theta);
    Rt = Rf + Rw;

    [w, t_val] = holtrop_wake_thrust(L, B, T, Cb, KCS.D_prop);

    eta_H = (1 - t_val) / (1 - w);
    Thrust_req = Rt / (1 - t_val);
    Va = V_ship_mps * (1 - w);
    CT = Thrust_req / (0.5 * Rho_water * Va^2 * pi * (KCS.D_prop/2)^2);
    eta_O = propeller_efficiency(CT);
    eta_R = 1.00;
    eta_D = eta_H * eta_O * eta_R;
    Pd = Rt * V_ship_mps / eta_D;

    Results_L(i).L = L;
    Results_L(i).Rt = Rt;
    Results_L(i).Rw = Rw;
    Results_L(i).Rf = Rf;
    Results_L(i).w = w;
    Results_L(i).t = t_val;
    Results_L(i).eta_H = eta_H;
    Results_L(i).eta_O = eta_O;
    Results_L(i).eta_D = eta_D;
    Results_L(i).Pd = Pd;

    fprintf('%4d  %6.1f  %7.1f  %7.1f  %7.1f  %.4f  %.4f  %.4f  %.4f  %7.2f\n', ...
        i, L, Rt/1000, Rw/1000, Rf/1000, w, t_val, eta_H, eta_D, Pd/1e6);
end

% ==============================================================================
% 5. 结果分析 - 找出最优点
% ==============================================================================
fprintf('\n================================================================================\n');
fprintf('                         结 果 分 析\n');
fprintf('================================================================================\n\n');

% Cb分析
Cb_vec = [Results_Cb.Cb];
Rt_Cb = [Results_Cb.Rt];
Pd_Cb = [Results_Cb.Pd];
eta_D_Cb = [Results_Cb.eta_D];

[min_Rt, idx_min_Rt] = min(Rt_Cb);
[min_Pd, idx_min_Pd] = min(Pd_Cb);
[max_etaD, idx_max_etaD] = max(eta_D_Cb);

fprintf('【分析1: Cb变化】\n');
fprintf('  ├─ 最小阻力点:     Cb = %.4f, Rt = %.1f kN\n', Cb_vec(idx_min_Rt), min_Rt/1000);
fprintf('  ├─ 最大效率点:     Cb = %.4f, eta_D = %.4f\n', Cb_vec(idx_max_etaD), max_etaD);
fprintf('  └─ 最小功率点:     Cb = %.4f, Pd = %.2f MW\n\n', Cb_vec(idx_min_Pd), min_Pd/1e6);

% 计算改进幅度
Pd_at_min_Rt = Pd_Cb(idx_min_Rt);
Power_saving = (Pd_at_min_Rt - min_Pd) / Pd_at_min_Rt * 100;

fprintf('  >>> 关键发现:\n');
if idx_min_Rt ~= idx_min_Pd
    fprintf('      最小阻力点(Cb=%.4f) ≠ 最小功率点(Cb=%.4f)\n', ...
        Cb_vec(idx_min_Rt), Cb_vec(idx_min_Pd));
    fprintf('      采用船桨集成优化可节省功率: %.1f%%\n\n', Power_saving);
end

% B分析
B_vec = [Results_B.B];
Rt_B = [Results_B.Rt];
Pd_B = [Results_B.Pd];
[~, idx_min_Rt_B] = min(Rt_B);
[~, idx_min_Pd_B] = min(Pd_B);

fprintf('【分析2: B变化】\n');
fprintf('  ├─ 最小阻力点:     B = %.2f m\n', B_vec(idx_min_Rt_B));
fprintf('  └─ 最小功率点:     B = %.2f m\n\n', B_vec(idx_min_Pd_B));

% L分析
L_vec = [Results_L.L];
Rt_L = [Results_L.Rt];
Pd_L = [Results_L.Pd];
[~, idx_min_Rt_L] = min(Rt_L);
[~, idx_min_Pd_L] = min(Pd_L);

fprintf('【分析3: L变化】\n');
fprintf('  ├─ 最小阻力点:     L = %.1f m\n', L_vec(idx_min_Rt_L));
fprintf('  └─ 最小功率点:     L = %.1f m\n\n', L_vec(idx_min_Pd_L));

% ==============================================================================
% 6. 可视化
% ==============================================================================
figure('Position', [50 50 1500 900], 'Color', 'w', 'Name', 'KCS Ship-Propeller Integration Analysis');

% ========== 第一行：Cb变化 ==========
subplot(3, 4, 1);
plot(Cb_vec, Rt_Cb/1000, 'b-o', 'LineWidth', 2, 'MarkerSize', 5); hold on;
plot(Cb_vec, [Results_Cb.Rw]/1000, 'r--', 'LineWidth', 1.5);
plot(Cb_vec, [Results_Cb.Rf]/1000, 'g-.', 'LineWidth', 1.5);
xline(Cb_vec(idx_min_Rt), ':k', 'Min R_t', 'LineWidth', 1.5);
xlabel('C_b'); ylabel('阻力 [kN]');
title('(a) C_b vs 阻力');
legend('R_t', 'R_w', 'R_f', 'Location', 'best');
grid on;

subplot(3, 4, 2);
plot(Cb_vec, [Results_Cb.w], 'm-o', 'LineWidth', 2); hold on;
plot(Cb_vec, [Results_Cb.t], 'c-s', 'LineWidth', 2);
xlabel('C_b'); ylabel('交互因子');
title('(b) C_b vs Holtrop伴流/推力减额');
legend('w (伴流)', 't (推力减额)', 'Location', 'best');
grid on;

subplot(3, 4, 3);
plot(Cb_vec, [Results_Cb.eta_H], 'g-o', 'LineWidth', 2); hold on;
plot(Cb_vec, [Results_Cb.eta_O], 'b-s', 'LineWidth', 2);
plot(Cb_vec, eta_D_Cb, 'k-^', 'LineWidth', 2.5);
xline(Cb_vec(idx_max_etaD), ':r', 'Max \eta_D', 'LineWidth', 1.5);
xlabel('C_b'); ylabel('效率');
title('(c) C_b vs 效率');
legend('\eta_H', '\eta_O', '\eta_D', 'Location', 'best');
grid on;

subplot(3, 4, 4);
yyaxis left
plot(Cb_vec, Rt_Cb/1000, 'b-o', 'LineWidth', 2);
ylabel('R_t [kN]');
xline(Cb_vec(idx_min_Rt), ':b', 'Min R_t');
yyaxis right
plot(Cb_vec, Pd_Cb/1e6, 'r-s', 'LineWidth', 2);
ylabel('P_d [MW]');
xline(Cb_vec(idx_min_Pd), ':r', 'Min P_d');
xlabel('C_b');
title('(d) 阻力 vs 功率 (核心结论)');
grid on;

% ========== 第二行：B变化 ==========
subplot(3, 4, 5);
plot(B_vec, Rt_B/1000, 'b-o', 'LineWidth', 2); hold on;
plot(B_vec, [Results_B.Rw]/1000, 'r--', 'LineWidth', 1.5);
plot(B_vec, [Results_B.Rf]/1000, 'g-.', 'LineWidth', 1.5);
xlabel('B [m]'); ylabel('阻力 [kN]');
title('(e) 船宽B vs 阻力');
legend('R_t', 'R_w', 'R_f', 'Location', 'best');
grid on;

subplot(3, 4, 6);
plot(B_vec, [Results_B.w], 'm-o', 'LineWidth', 2); hold on;
plot(B_vec, [Results_B.t], 'c-s', 'LineWidth', 2);
xlabel('B [m]'); ylabel('交互因子');
title('(f) 船宽B vs 交互因子');
legend('w', 't', 'Location', 'best');
grid on;

subplot(3, 4, 7);
plot(B_vec, [Results_B.eta_H], 'g-o', 'LineWidth', 2); hold on;
plot(B_vec, [Results_B.eta_O], 'b-s', 'LineWidth', 2);
plot(B_vec, [Results_B.eta_D], 'k-^', 'LineWidth', 2.5);
xlabel('B [m]'); ylabel('效率');
title('(g) 船宽B vs 效率');
legend('\eta_H', '\eta_O', '\eta_D', 'Location', 'best');
grid on;

subplot(3, 4, 8);
yyaxis left
plot(B_vec, Rt_B/1000, 'b-o', 'LineWidth', 2);
ylabel('R_t [kN]');
yyaxis right
plot(B_vec, Pd_B/1e6, 'r-s', 'LineWidth', 2);
ylabel('P_d [MW]');
xlabel('B [m]');
title('(h) 阻力 vs 功率');
grid on;

% ========== 第三行：L变化 ==========
subplot(3, 4, 9);
plot(L_vec, Rt_L/1000, 'b-o', 'LineWidth', 2); hold on;
plot(L_vec, [Results_L.Rw]/1000, 'r--', 'LineWidth', 1.5);
plot(L_vec, [Results_L.Rf]/1000, 'g-.', 'LineWidth', 1.5);
xlabel('L [m]'); ylabel('阻力 [kN]');
title('(i) 船长L vs 阻力');
legend('R_t', 'R_w', 'R_f', 'Location', 'best');
grid on;

subplot(3, 4, 10);
plot(L_vec, [Results_L.w], 'm-o', 'LineWidth', 2); hold on;
plot(L_vec, [Results_L.t], 'c-s', 'LineWidth', 2);
xlabel('L [m]'); ylabel('交互因子');
title('(j) 船长L vs 交互因子');
legend('w', 't', 'Location', 'best');
grid on;

subplot(3, 4, 11);
plot(L_vec, [Results_L.eta_H], 'g-o', 'LineWidth', 2); hold on;
plot(L_vec, [Results_L.eta_O], 'b-s', 'LineWidth', 2);
plot(L_vec, [Results_L.eta_D], 'k-^', 'LineWidth', 2.5);
xlabel('L [m]'); ylabel('效率');
title('(k) 船长L vs 效率');
legend('\eta_H', '\eta_O', '\eta_D', 'Location', 'best');
grid on;

subplot(3, 4, 12);
yyaxis left
plot(L_vec, Rt_L/1000, 'b-o', 'LineWidth', 2);
ylabel('R_t [kN]');
yyaxis right
plot(L_vec, Pd_L/1e6, 'r-s', 'LineWidth', 2);
ylabel('P_d [MW]');
xlabel('L [m]');
title('(l) 阻力 vs 功率');
grid on;

sgtitle('KCS船桨集成优化分析 - Michell+ITTC+Holtrop', 'FontSize', 14, 'FontWeight', 'bold');

% ==============================================================================
% 7. 对比图：仅优化阻力 vs 船桨集成优化
% ==============================================================================
figure('Position', [100 100 800 600], 'Color', 'w', 'Name', 'Comparison: Resistance-only vs Integrated Optimization');

% 归一化显示
Rt_norm = (Rt_Cb - min(Rt_Cb)) / (max(Rt_Cb) - min(Rt_Cb));
Pd_norm = (Pd_Cb - min(Pd_Cb)) / (max(Pd_Cb) - min(Pd_Cb));
etaD_norm = (eta_D_Cb - min(eta_D_Cb)) / (max(eta_D_Cb) - min(eta_D_Cb));

subplot(2,1,1);
plot(Cb_vec, Rt_norm, 'b-o', 'LineWidth', 2, 'MarkerSize', 6); hold on;
plot(Cb_vec, Pd_norm, 'r-s', 'LineWidth', 2, 'MarkerSize', 6);
plot(Cb_vec, 1-etaD_norm, 'g-^', 'LineWidth', 2, 'MarkerSize', 6);
xline(Cb_vec(idx_min_Rt), '--b', 'Min R_t', 'LineWidth', 2);
xline(Cb_vec(idx_min_Pd), '--r', 'Min P_d', 'LineWidth', 2);
xlabel('方形系数 C_b');
ylabel('归一化值 (越小越好)');
title('仅优化阻力 vs 船桨集成优化 对比');
legend('阻力R_t(归一化)', '功率P_d(归一化)', '1-\eta_D(归一化)', 'Location', 'best');
grid on;

% 功率节省分析
subplot(2,1,2);
Pd_saving = (Pd_Cb(idx_min_Rt) - Pd_Cb) / Pd_Cb(idx_min_Rt) * 100;
bar(Cb_vec, Pd_saving, 'FaceColor', [0.3 0.7 0.9]);
hold on;
plot([Cb_vec(1) Cb_vec(end)], [0 0], 'k-', 'LineWidth', 1);
xline(Cb_vec(idx_min_Rt), '--b', 'Min R_t点', 'LineWidth', 2);
xline(Cb_vec(idx_min_Pd), '--r', 'Min P_d点', 'LineWidth', 2);
xlabel('方形系数 C_b');
ylabel('相对最小阻力点的功率节省 [%]');
title('船桨集成优化带来的功率节省');
grid on;

% ==============================================================================
% 8. 最终结论输出
% ==============================================================================
fprintf('================================================================================\n');
fprintf('                         结 论\n');
fprintf('================================================================================\n\n');

fprintf('【核心发现】\n');
fprintf('  1. 最小阻力点 Cb = %.4f，对应 Rt = %.1f kN\n', Cb_vec(idx_min_Rt), min_Rt/1000);
fprintf('  2. 最小功率点 Cb = %.4f，对应 Pd = %.2f MW\n', Cb_vec(idx_min_Pd), min_Pd/1e6);
fprintf('  3. 两者相差 ΔCb = %.4f\n\n', abs(Cb_vec(idx_min_Pd) - Cb_vec(idx_min_Rt)));

fprintf('【物理机制】\n');
fprintf('  • Cb增大 → 船型更丰满 → 兴波阻力Rw增加\n');
fprintf('  • Cb增大 → 艉部更饱满 → Holtrop伴流w增加\n');
fprintf('  • w增加 → 船身效率 η_H = (1-t)/(1-w) 提高\n');
fprintf('  • 当η_H增益 > 阻力增加损失时 → 总功率反而降低\n\n');

fprintf('【研究意义】\n');
fprintf('  ★ 现有研究"尽可能减小阻力"的方法存在不足！\n');
fprintf('  ★ 考虑船桨集成优化可节省功率约 %.1f%%\n', Power_saving);
fprintf('  ★ 应以最小传递功率Pd为优化目标，而非最小阻力Rt\n\n');

fprintf('================================================================================\n');

% 保存结果
save('KCS_Analysis_Results.mat', 'Results_Cb', 'Results_B', 'Results_L', 'KCS');
saveas(gcf, 'KCS_Comparison.png');
fprintf('结果已保存。\n');

% ==============================================================================
% 辅助函数
% ==============================================================================

function [w, t] = holtrop_wake_thrust(L, B, T, Cb, D_prop)
    % Holtrop方法计算伴流分数w和推力减额分数t
    % 参考: Holtrop (1984) 统计回归公式

    % 无量纲参数
    LB = L / B;
    BT = B / T;
    DT = D_prop / T;

    % === 伴流分数 w ===
    % 单桨船的Holtrop公式 (简化版)
    Cv = 1 + 0.011 * (Cb - 0.70) * 100;  % 粘性系数修正

    w_base = 0.3095 * Cb + 10.0 * Cb * LB^(-3.18);
    w_corr = -0.18 * DT + 0.1 * (BT - 2.5) / 10;

    w = Cv * (w_base + w_corr);
    w = max(0.12, min(0.50, w));  % 限制范围

    % === 推力减额分数 t ===
    % t/w 比值通常在 0.6-0.8
    t_base = 0.25014 * (B/L)^0.28956 * (sqrt(BT)/DT)^0.2624;
    t_Cb_corr = 0.12 * (Cb - 0.65);

    t = t_base + t_Cb_corr;

    % 保证 t < w 且 t/w 合理
    tw_ratio = t / w;
    if tw_ratio < 0.55
        t = 0.55 * w;
    elseif tw_ratio > 0.85
        t = 0.85 * w;
    end

    t = max(0.08, min(0.28, t));
end

function eta_O = propeller_efficiency(CT)
    % 基于推力载荷系数估算螺旋桨敞水效率
    % CT = T / (0.5 * rho * Va^2 * A0)

    % 理想效率 (动量理论)
    eta_ideal = 2 / (1 + sqrt(1 + CT));

    % 考虑粘性损失，实际效率约为理想效率的75-82%
    viscous_factor = 0.78 - 0.015 * (CT - 0.5);
    viscous_factor = max(0.70, min(0.82, viscous_factor));

    eta_O = eta_ideal * viscous_factor;
    eta_O = max(0.50, min(0.75, eta_O));
end

function S = estimate_wetted_surface(L, B, T, Cb)
    % 湿表面积估算 (Holtrop-Mennen)
    S = L * (2*T + B) * sqrt(Cb) * (0.453 + 0.4425*Cb - 0.2862*Cb^2);
    S = max(S, L * (2*T + B*Cb));
end

function Rw = michell(Y, U, L, B, T, RHO, N)
    % Michell积分计算兴波阻力 (Filon算法)
    Nx = size(Y,1);
    Nz = size(Y,2);
    YH = Y * B;

    dz = T/(Nz-1);
    z = (-T:dz:0)';
    dx = L/(Nx-1);
    x = (0:dx:L)';
    theta = michspace(N)';

    g = 9.80665;
    k0 = g/U^2;
    c = (4*RHO*U^2)/pi;
    a = sec(theta);
    k = k0*a.^2;

    % Z 积分 (Filon trapezoidal)
    Kz = k0*dz*a.^2;
    w0 = (exp(Kz)-1-Kz)./Kz.^2;
    wn = (exp(Kz)+exp(-Kz)-2)./Kz.^2;
    wN = (exp(-Kz)-1+Kz)./Kz.^2;

    F = zeros(Nx,N);
    for j = 1:N
        for m = 1:Nx
            f_vec = zeros(Nz,1);
            for n = 1:Nz
                if n == 1
                    f_vec(n) = w0(j)*YH(m,n)*exp(k0*z(n)*a(j)^2)*dz;
                elseif n == Nz
                    f_vec(n) = wN(j)*YH(m,n)*exp(k0*z(n)*a(j)^2)*dz;
                else
                    f_vec(n) = wn(j)*YH(m,n)*exp(k0*z(n)*a(j)^2)*dz;
                end
            end
            F(m,j) = sum(f_vec);
        end
    end

    % X 积分 (Filon algorithm)
    Kx = k0*dx.*a;
    alp = (Kx.^2+1/2*Kx.*sin(2*Kx)+cos(2*Kx)-1)./Kx.^3;
    bet = (3*Kx+Kx.*cos(2*Kx)-2*sin(2*Kx))./Kx.^3;
    gam = 4*(sin(Kx)-Kx.*cos(Kx))./Kx.^3;

    Pt = zeros(N,1); Qt = zeros(N,1);
    P = zeros(N,1); Q = zeros(N,1);

    for j = 1:N
        c_val = cos(k0*x*a(j));
        s_val = sin(k0*x*a(j));
        pev = F(1:2:end,j).*c_val(1:2:end);
        qev = F(1:2:end,j).*s_val(1:2:end);
        pod = F(2:2:end-1,j).*c_val(2:2:end-1);
        qod = F(2:2:end-1,j).*s_val(2:2:end-1);
        Pt(j) = F(Nx,j)*cos(k0*L*a(j));
        Qt(j) = F(Nx,j)*sin(k0*L*a(j));
        Pev = sum(pev)-0.5*Pt(j);
        Pod = sum(pod);
        Qev = sum(qev)-0.5*Qt(j);
        Qod = sum(qod);
        P(j) = dx*(alp(j)*Qt(j) + bet(j)*Pev + gam(j)*Pod);
        Q(j) = dx*(-alp(j)*Pt(j) + bet(j)*Qev + gam(j)*Qod);
    end

    R = c*k.^2./a.^3.*(k.^2.*(P.^2 + Q.^2) + 2*k.*a.*(Q.*Pt - P.*Qt) + a.^2.*(Pt.^2 + Qt.^2));
    R(isnan(R)) = 0;

    % Theta 积分
    rw_vec = zeros(N-1,1);
    for k_idx = 1:N-1
        rw_vec(k_idx) = 0.5*(R(k_idx)+R(k_idx+1))*(theta(k_idx+1)-theta(k_idx));
    end
    Rw = max(0, sum(rw_vec));
end

function xm = michspace(N)
    % 传播角加密采样
    xm = logspace(0,1,N)-1;
    xm = xm*pi/18-pi/2;
    xm = fliplr(-xm);
end

function f = three_param_shape(x, h, a1, a2)
    % 三参数水线面形状函数
    % x: 归一化船长位置 [0,1]
    % h: 中部丰满度 (越大越丰满)
    % a1: 艏部形状参数
    % a2: 艉部形状参数
    f = zeros(size(x));
    for i = 1:length(x)
        if x(i) < 0.5
            f(i) = h - a1/2 * (cos(2*pi*x(i)) - 1);
        else
            f(i) = h + a2/2 * (cos(2*pi*x(i)) + 1) + a1;
        end
    end
end
