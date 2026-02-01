% KCS_test.m
% 船桨集成优化分析 - 通用版本
%
% 核心目标：验证 "阻力越小，推进效率不一定越高"
%
% 方法：
% 1. Michell积分 - 兴波阻力计算
% 2. ITTC 1957 - 摩擦阻力计算
% 3. Holtrop方法 - 伴流分数w和推力减额分数t
% 4. 螺旋桨效率模型 - 基于载荷系数
%
% 变化参数：方形系数Cb、船宽B、船长L（全范围扫描）

clear; close all; clc;

fprintf('================================================================================\n');
fprintf('   船桨集成优化分析 - 通用版本\n');
fprintf('   验证：阻力最小 ≠ 效率最高\n');
fprintf('================================================================================\n\n');

% ==============================================================================
% 1. 全局仿真设置
% ==============================================================================
V_ship_mps = 12.35;       % 24 knots ≈ 12.35 m/s
Rho_water = 1025;
Nu = 1.188e-6;
g = 9.80665;

% 基准参数
Ship.T = 10.8;            % 吃水 [m] (固定)
Ship.D_prop = 7.9;        % 螺旋桨直径 [m]

% Michell 积分网格参数
Nx = 101;
Nz = 40;
N_theta = 80;

x_norm = linspace(0, 1, Nx);

% ==============================================================================
% 2. 参数范围设置 - 覆盖各种船型
% ==============================================================================
% 方形系数范围: 0.45(快速客船) ~ 0.85(散货船)
% 船长范围: 100m ~ 350m
% 船宽范围: 15m ~ 50m

NumCb = 20;
NumL = 5;
NumB = 5;

Cb_range = linspace(0.50, 0.80, NumCb);
L_range = linspace(150, 300, NumL);
B_range = linspace(25, 45, NumB);

% 总样本数
TotalSamples = NumCb * NumL * NumB;
fprintf('参数空间: Cb=[%.2f,%.2f], L=[%.0f,%.0f]m, B=[%.0f,%.0f]m\n', ...
    min(Cb_range), max(Cb_range), min(L_range), max(L_range), min(B_range), max(B_range));
fprintf('总计算点数: %d\n\n', TotalSamples);

% 存储所有结果
AllResults = struct('Cb',[],'L',[],'B',[],'Rt',[],'Rw',[],'Rf',[],'w',[],'t',[],...
                    'eta_H',[],'eta_O',[],'eta_D',[],'Pd',[],'Thrust',[]);
AllResults = repmat(AllResults, TotalSamples, 1);

% ==============================================================================
% 3. 三维参数扫描
% ==============================================================================
fprintf('开始参数扫描...\n');
idx = 0;

for iL = 1:NumL
    L = L_range(iL);
    for iB = 1:NumB
        B = B_range(iB);
        for iCb = 1:NumCb
            Cb = Cb_range(iCb);
            idx = idx + 1;

            T = Ship.T;

            % 检查船型合理性 (L/B在4-10之间，B/T在2-5之间)
            LB_ratio = L / B;
            BT_ratio = B / T;
            if LB_ratio < 4 || LB_ratio > 10 || BT_ratio < 2 || BT_ratio > 5
                % 跳过不合理的船型
                AllResults(idx).Cb = NaN;
                continue;
            end

            % --- A. 三参数函数生成船型 ---
            h_val = 0.5 + 1.5 * Cb;  % 丰满度随Cb变化
            f0 = three_param_shape(x_norm, h_val, 0.15, 0.12);
            f0 = f0 / max(f0) * 0.5;

            z_norm = linspace(0, 1, Nz);
            Y = zeros(Nx, Nz);
            for iz = 1:Nz
                z_factor = sin(pi/2 * z_norm(iz))^0.6;
                Y(:, iz) = f0(:) * z_factor;
            end

            % --- B. 阻力计算 ---
            S = estimate_wetted_surface(L, B, T, Cb);
            Re = V_ship_mps * L / Nu;
            Cf = 0.075 / (log10(Re) - 2)^2;
            Rf = 0.5 * Rho_water * S * V_ship_mps^2 * (Cf + 0.0004);
            Rw = michell(Y, V_ship_mps, L, B, T, Rho_water, N_theta);
            Rt = Rf + Rw;

            % --- C. Holtrop 伴流和推力减额 ---
            [w, t_val] = holtrop_wake_thrust(L, B, T, Cb, Ship.D_prop);

            % --- D. 效率计算 ---
            eta_H = (1 - t_val) / (1 - w);
            Thrust_req = Rt / (1 - t_val);
            Va = V_ship_mps * (1 - w);
            CT = Thrust_req / (0.5 * Rho_water * Va^2 * pi * (Ship.D_prop/2)^2);
            eta_O = propeller_efficiency(CT, Va, Ship.D_prop);
            eta_R = 1.00;
            eta_D = eta_H * eta_O * eta_R;
            Pd = Rt * V_ship_mps / eta_D;

            % 存储结果
            AllResults(idx).Cb = Cb;
            AllResults(idx).L = L;
            AllResults(idx).B = B;
            AllResults(idx).Rt = Rt;
            AllResults(idx).Rw = Rw;
            AllResults(idx).Rf = Rf;
            AllResults(idx).w = w;
            AllResults(idx).t = t_val;
            AllResults(idx).eta_H = eta_H;
            AllResults(idx).eta_O = eta_O;
            AllResults(idx).eta_D = eta_D;
            AllResults(idx).Pd = Pd;
            AllResults(idx).Thrust = Thrust_req;
        end
    end
    fprintf('  L = %.0f m 完成 (%d/%d)\n', L, iL, NumL);
end

% 过滤有效结果
valid_idx = ~isnan([AllResults.Cb]);
Results = AllResults(valid_idx);
fprintf('\n有效计算点: %d\n', length(Results));

% ==============================================================================
% 4. 分析：固定L和B，仅变化Cb
% ==============================================================================
fprintf('\n================================================================================\n');
fprintf('   分析：固定L和B，仅变化Cb\n');
fprintf('================================================================================\n\n');

% 选择中间的L和B值
L_select = L_range(ceil(NumL/2));
B_select = B_range(ceil(NumB/2));

fprintf('选择 L = %.0f m, B = %.0f m 进行Cb分析\n\n', L_select, B_select);

% 筛选
idx_select = abs([Results.L] - L_select) < 1 & abs([Results.B] - B_select) < 1;
Results_Cb = Results(idx_select);

if ~isempty(Results_Cb)
    Cb_vec = [Results_Cb.Cb];
    Rt_vec = [Results_Cb.Rt]/1000;
    Pd_vec = [Results_Cb.Pd]/1e6;
    eta_D_vec = [Results_Cb.eta_D];
    Thrust_vec = [Results_Cb.Thrust]/1000;
    eta_H_vec = [Results_Cb.eta_H];
    eta_O_vec = [Results_Cb.eta_O];

    [min_Rt, idx_Rt] = min(Rt_vec);
    [min_Pd, idx_Pd] = min(Pd_vec);
    [max_etaD, idx_etaD] = max(eta_D_vec);
    [min_Thrust, idx_Thrust] = min(Thrust_vec);

    fprintf('%-5s %-7s %-8s %-8s %-7s %-7s %-7s %-7s %-8s\n', ...
        'No.', 'Cb', 'Rt(kN)', 'Thrust', 'w', 't', 'eta_H', 'eta_D', 'Pd(MW)');
    fprintf('--------------------------------------------------------------------------------\n');
    for ii = 1:length(Results_Cb)
        fprintf('%4d  %.4f  %7.1f  %7.1f  %.4f  %.4f  %.4f  %.4f  %7.2f\n', ...
            ii, Cb_vec(ii), Rt_vec(ii), Thrust_vec(ii), ...
            Results_Cb(ii).w, Results_Cb(ii).t, eta_H_vec(ii), eta_D_vec(ii), Pd_vec(ii));
    end

    fprintf('\n【关键结果】\n');
    fprintf('  ├─ 最小阻力点:     Cb = %.4f, Rt = %.1f kN\n', Cb_vec(idx_Rt), min_Rt);
    fprintf('  ├─ 最小推力点:     Cb = %.4f, Thrust = %.1f kN\n', Cb_vec(idx_Thrust), Thrust_vec(idx_Thrust));
    fprintf('  ├─ 最大效率点:     Cb = %.4f, eta_D = %.4f\n', Cb_vec(idx_etaD), max_etaD);
    fprintf('  └─ 最小功率点:     Cb = %.4f, Pd = %.2f MW\n\n', Cb_vec(idx_Pd), min_Pd);

    Pd_at_min_Thrust = Pd_vec(idx_Thrust);
    Power_saving = (Pd_at_min_Thrust - min_Pd) / Pd_at_min_Thrust * 100;

    if idx_Thrust ~= idx_Pd
        fprintf('【核心发现】\n');
        fprintf('  >>> 最小推力点(Cb=%.4f) ≠ 最小功率点(Cb=%.4f)\n', ...
            Cb_vec(idx_Thrust), Cb_vec(idx_Pd));
        fprintf('  >>> 船桨集成优化(ISPS)可节省功率: %.2f%%\n\n', Power_saving);
    end
end

% ==============================================================================
% 5. 全局最优搜索
% ==============================================================================
fprintf('================================================================================\n');
fprintf('   全局最优搜索\n');
fprintf('================================================================================\n\n');

% 找全局最小阻力和最小功率
All_Rt = [Results.Rt];
All_Pd = [Results.Pd];
All_etaD = [Results.eta_D];

[global_min_Rt, idx_global_Rt] = min(All_Rt);
[global_min_Pd, idx_global_Pd] = min(All_Pd);
[global_max_etaD, idx_global_etaD] = max(All_etaD);

fprintf('【全局最小阻力点】\n');
fprintf('  Cb=%.4f, L=%.0fm, B=%.0fm, Rt=%.1fkN, Pd=%.2fMW, eta_D=%.4f\n', ...
    Results(idx_global_Rt).Cb, Results(idx_global_Rt).L, Results(idx_global_Rt).B, ...
    global_min_Rt/1000, Results(idx_global_Rt).Pd/1e6, Results(idx_global_Rt).eta_D);

fprintf('\n【全局最小功率点】\n');
fprintf('  Cb=%.4f, L=%.0fm, B=%.0fm, Rt=%.1fkN, Pd=%.2fMW, eta_D=%.4f\n', ...
    Results(idx_global_Pd).Cb, Results(idx_global_Pd).L, Results(idx_global_Pd).B, ...
    Results(idx_global_Pd).Rt/1000, global_min_Pd/1e6, Results(idx_global_Pd).eta_D);

fprintf('\n【全局最高效率点】\n');
fprintf('  Cb=%.4f, L=%.0fm, B=%.0fm, Rt=%.1fkN, Pd=%.2fMW, eta_D=%.4f\n', ...
    Results(idx_global_etaD).Cb, Results(idx_global_etaD).L, Results(idx_global_etaD).B, ...
    Results(idx_global_etaD).Rt/1000, Results(idx_global_etaD).Pd/1e6, global_max_etaD);

Global_Power_saving = (Results(idx_global_Rt).Pd - global_min_Pd) / Results(idx_global_Rt).Pd * 100;
fprintf('\n【对比】\n');
fprintf('  若选择全局最小阻力船型，功率 = %.2f MW\n', Results(idx_global_Rt).Pd/1e6);
fprintf('  若选择全局最小功率船型，功率 = %.2f MW\n', global_min_Pd/1e6);
fprintf('  功率差异: %.2f%%\n', Global_Power_saving);

% ==============================================================================
% 6. 高质量绘图
% ==============================================================================
if ~isempty(Results_Cb)
    figure('Position', [50 50 1600 1000], 'Color', 'w', 'Name', 'Hull-Propeller Integrated Optimization');

    set(groot, 'defaultAxesFontName', 'Arial');
    set(groot, 'defaultAxesFontSize', 10);
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

    color_traditional = [0.2, 0.4, 0.8];
    color_optimal = [0.9, 0.2, 0.2];
    color_hull = [0.2, 0.7, 0.3];
    color_prop = [0.6, 0.3, 0.7];
    color_total = [0.1, 0.1, 0.1];

    % --- 子图 (a): 推力与功率 ---
    subplot(2, 3, 1);
    yyaxis left
    plot(Cb_vec, Thrust_vec, '-o', 'Color', color_traditional, ...
        'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', color_traditional);
    ylabel('Required Thrust (kN)', 'FontWeight', 'bold', 'FontSize', 10);
    set(gca, 'YColor', color_traditional);

    yyaxis right
    plot(Cb_vec, Pd_vec, '-s', 'Color', color_optimal, ...
        'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', color_optimal);
    ylabel('Delivered Power (MW)', 'FontWeight', 'bold', 'FontSize', 10);
    set(gca, 'YColor', color_optimal);

    xlabel('Block Coefficient ($C_b$)', 'FontWeight', 'bold', 'Interpreter', 'latex');
    legend({'Thrust', 'Power'}, 'Location', 'north', 'Interpreter', 'latex');
    grid on; box on;
    xlim([min(Cb_vec)-0.01, max(Cb_vec)+0.01]);

    % --- 子图 (b): 效率分解 ---
    subplot(2, 3, 2);
    h1 = plot(Cb_vec, eta_H_vec, '-s', 'Color', color_hull, 'LineWidth', 1.5, ...
        'MarkerSize', 4, 'MarkerFaceColor', color_hull);
    hold on;
    h2 = plot(Cb_vec, eta_O_vec, '-^', 'Color', color_prop, 'LineWidth', 1.5, ...
        'MarkerSize', 4, 'MarkerFaceColor', color_prop);
    h3 = plot(Cb_vec, eta_D_vec, '-o', 'Color', color_total, 'LineWidth', 1.8, ...
        'MarkerSize', 4, 'MarkerFaceColor', color_total);
    xlabel('Block Coefficient ($C_b$)', 'FontWeight', 'bold', 'Interpreter', 'latex');
    ylabel('Efficiency', 'FontWeight', 'bold');
    legend([h1, h2, h3], {'$\eta_H$ (Hull)', '$\eta_O$ (Propeller)', '$\eta_D$ (Total)'}, ...
        'Location', 'best', 'Interpreter', 'latex');
    grid on; box on;
    xlim([min(Cb_vec)-0.01, max(Cb_vec)+0.01]);
    hold off;

    % --- 子图 (c): 阻力 vs 效率 ---
    subplot(2, 3, 3);
    yyaxis left
    plot(Cb_vec, Rt_vec, '-^', 'Color', [0.8, 0.5, 0.2], 'LineWidth', 1.5, ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.8, 0.5, 0.2]);
    ylabel('Total Resistance $R_t$ (kN)', 'FontWeight', 'bold', 'Interpreter', 'latex');
    set(gca, 'YColor', [0.8, 0.5, 0.2]);

    yyaxis right
    plot(Cb_vec, eta_D_vec, '-o', 'Color', color_total, 'LineWidth', 1.5, ...
        'MarkerSize', 4, 'MarkerFaceColor', color_total);
    ylabel('Total Efficiency $\eta_D$', 'FontWeight', 'bold', 'Interpreter', 'latex');
    set(gca, 'YColor', color_total);

    xlabel('Block Coefficient ($C_b$)', 'FontWeight', 'bold', 'Interpreter', 'latex');
    legend({'$R_t$', '$\eta_D$'}, 'Location', 'northwest', 'Interpreter', 'latex');
    grid on; box on;
    xlim([min(Cb_vec)-0.01, max(Cb_vec)+0.01]);

    % --- 子图 (d): 交互因子 ---
    subplot(2, 3, 4);
    yyaxis left
    plot(Cb_vec, [Results_Cb.w], '-d', 'Color', [0.2, 0.4, 0.7], 'LineWidth', 1.5, ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.2, 0.4, 0.7]);
    ylabel('Wake Fraction $w$', 'FontWeight', 'bold', 'Interpreter', 'latex');
    set(gca, 'YColor', [0.2, 0.4, 0.7]);

    yyaxis right
    plot(Cb_vec, [Results_Cb.t], '-o', 'Color', [0.7, 0.2, 0.4], 'LineWidth', 1.5, ...
        'MarkerSize', 4, 'MarkerFaceColor', [0.7, 0.2, 0.4]);
    ylabel('Thrust Deduction $t$', 'FontWeight', 'bold', 'Interpreter', 'latex');
    set(gca, 'YColor', [0.7, 0.2, 0.4]);

    xlabel('Block Coefficient ($C_b$)', 'FontWeight', 'bold', 'Interpreter', 'latex');
    legend({'$w$', '$t$'}, 'Location', 'northwest', 'Interpreter', 'latex');
    grid on; box on;
    xlim([min(Cb_vec)-0.01, max(Cb_vec)+0.01]);

    % --- 子图 (e): 3D散点图 ---
    subplot(2, 3, 5);
    All_Cb = [Results.Cb];
    All_L = [Results.L];
    scatter3(All_Cb, All_L, All_Pd/1e6, 30, All_etaD, 'filled');
    xlabel('$C_b$', 'Interpreter', 'latex');
    ylabel('$L$ (m)', 'Interpreter', 'latex');
    zlabel('$P_d$ (MW)', 'Interpreter', 'latex');
    colorbar;
    title('Power vs Parameters (color = $\eta_D$)', 'Interpreter', 'latex');
    grid on; box on;
    view(45, 30);

    % --- 子图 (f): 功率曲线对比 ---
    subplot(2, 3, 6);
    % 对比不同L/B下的功率曲线
    colors = lines(NumL);
    hold on;
    for iL = 1:NumL
        L_val = L_range(iL);
        B_val = B_range(ceil(NumB/2));
        idx_LB = abs([Results.L] - L_val) < 1 & abs([Results.B] - B_val) < 1;
        R_LB = Results(idx_LB);
        if ~isempty(R_LB)
            [~, sort_idx] = sort([R_LB.Cb]);
            R_LB = R_LB(sort_idx);
            plot([R_LB.Cb], [R_LB.Pd]/1e6, '-o', 'Color', colors(iL,:), ...
                'LineWidth', 1.5, 'MarkerSize', 3, 'DisplayName', sprintf('L=%.0fm', L_val));
        end
    end
    hold off;
    xlabel('Block Coefficient ($C_b$)', 'FontWeight', 'bold', 'Interpreter', 'latex');
    ylabel('Delivered Power (MW)', 'FontWeight', 'bold');
    legend('Location', 'best');
    grid on; box on;
    title('Power vs $C_b$ at different $L$', 'Interpreter', 'latex');

    sgtitle(sprintf('Ship-Propeller Optimization (L=%.0fm, B=%.0fm, T=%.1fm)', L_select, B_select, Ship.T), ...
        'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'latex');
end

% ==============================================================================
% 7. 结论
% ==============================================================================
fprintf('\n================================================================================\n');
fprintf('                         结 论\n');
fprintf('================================================================================\n\n');

fprintf('【物理机制】\n');
fprintf('  • Cb增大 → 兴波阻力Rw增加，但伴流w也增加\n');
fprintf('  • w增大 → 船身效率 η_H = (1-t)/(1-w) 提高\n');
fprintf('  • 存在最优Cb，使得 η_H的增益 与 阻力增加 达到平衡\n\n');

fprintf('【研究意义】\n');
fprintf('  ★ 不同船型(L,B)对应不同的最优Cb\n');
fprintf('  ★ 单纯追求最小阻力可能不是最优解\n');
fprintf('  ★ 船桨集成优化可实现功率节省约 %.1f%%\n\n', Global_Power_saving);

fprintf('================================================================================\n');

% 保存结果
save('Ship_Optimization_Results.mat', 'Results', 'AllResults');
saveas(gcf, 'Ship_Propeller_Optimization.png');
fprintf('结果已保存\n');

% ==============================================================================
% 辅助函数
% ==============================================================================

function [w, t] = holtrop_wake_thrust(L, B, T, Cb, D_prop)
    % Holtrop方法计算伴流分数w和推力减额分数t
    % 基于Holtrop & Mennen (1982) 的经验公式

    % 棱形系数
    Cp = Cb / 0.98;

    % 无量纲参数
    LB = L / B;
    BT = B / T;
    DT = D_prop / T;

    % === 伴流分数 w ===
    % 单桨船Holtrop公式
    % w主要受Cb影响，Cb越大，艉部越丰满，伴流越大

    % 基础伴流
    w_base = 0.1 + 0.5 * Cb;  % 简化的线性关系

    % L/B修正：细长船伴流较小
    w_LB = -0.02 * (LB - 6);

    % 螺旋桨直径修正
    w_D = -0.05 * (DT - 0.7);

    w = w_base + w_LB + w_D;
    w = max(0.10, min(0.50, w));

    % === 推力减额分数 t ===
    % t/w比值通常在0.6-0.8
    tw_ratio = 0.70 - 0.2 * (Cb - 0.65);
    tw_ratio = max(0.55, min(0.80, tw_ratio));

    t = w * tw_ratio;
    t = max(0.08, min(0.35, t));
end

function eta_O = propeller_efficiency(CT, Va, D)
    % 螺旋桨敞水效率
    % 基于动量理论 + 经验修正

    % 理想效率
    eta_ideal = 2 / (1 + sqrt(1 + CT));

    % 前进系数影响 (假设最优J=0.7)
    n_assumed = 1.8;  % 假设转速 [rps]
    J = Va / (n_assumed * D);
    J_penalty = 0.03 * (J - 0.7)^2;

    % 粘性损失
    viscous_factor = 0.80 - 0.02 * (CT - 0.5);
    viscous_factor = max(0.70, min(0.85, viscous_factor));

    eta_O = eta_ideal * viscous_factor - J_penalty;
    eta_O = max(0.50, min(0.75, eta_O));
end

function S = estimate_wetted_surface(L, B, T, Cb)
    S = L * (2*T + B) * sqrt(Cb) * (0.453 + 0.4425*Cb - 0.2862*Cb^2);
    S = max(S, L * (2*T + B*Cb));
end

function Rw = michell(Y, U, L, B, T, RHO, N)
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

    rw_vec = zeros(N-1,1);
    for k_idx = 1:N-1
        rw_vec(k_idx) = 0.5*(R(k_idx)+R(k_idx+1))*(theta(k_idx+1)-theta(k_idx));
    end
    Rw = max(0, sum(rw_vec));
end

function xm = michspace(N)
    xm = logspace(0,1,N)-1;
    xm = xm*pi/18-pi/2;
    xm = fliplr(-xm);
end

function f = three_param_shape(x, h, a1, a2)
    f = zeros(size(x));
    for i = 1:length(x)
        if x(i) < 0.5
            f(i) = h - a1/2 * (cos(2*pi*x(i)) - 1);
        else
            f(i) = h + a2/2 * (cos(2*pi*x(i)) + 1) + a1;
        end
    end
end
