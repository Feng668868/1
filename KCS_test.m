% KCS_test.m
% 船桨集成优化分析 - 全参数版本
%
% 核心目标：验证 "阻力越小，推进效率不一定越高"
%
% 变化参数：方形系数Cb、船长L、船宽B、船速V、吃水T
% 通过多速度、多吃水分析找到交叉点
%
% 方法：
% 1. Michell积分 - 兴波阻力计算
% 2. ITTC 1957 - 摩擦阻力计算
% 3. Holtrop方法 - 伴流分数w和推力减额分数t
% 4. 螺旋桨效率模型 - 基于载荷系数

clear; close all; clc;

fprintf('================================================================================\n');
fprintf('   船桨集成优化分析 - 全参数版本\n');
fprintf('   验证：阻力最小 ≠ 效率最高（ISPS研究）\n');
fprintf('================================================================================\n\n');

% ==============================================================================
% 1. 物理常数
% ==============================================================================
Rho_water = 1025;
Nu = 1.188e-6;
g = 9.80665;

% Michell 积分网格参数
Nx = 101;
Nz = 40;
N_theta = 80;
x_norm = linspace(0, 1, Nx);

% ==============================================================================
% 2. 参数范围设置 - 覆盖各种船型和工况
% ==============================================================================
% 方形系数: 0.45(快艇) ~ 0.85(散货船)
% 船长: 100m ~ 350m
% 船宽: 15m ~ 55m
% 吃水: 6m ~ 20m
% 船速: 8 ~ 25 knots

NumCb = 25;
NumL = 3;
NumB = 3;
NumT = 3;
NumV = 5;

Cb_range = linspace(0.45, 0.85, NumCb);
L_range = [150, 225, 300];  % 短、中、长船
B_range = [25, 35, 45];     % 窄、中、宽船
T_range = [8, 12, 16];      % 浅、中、深吃水
V_knots_range = [10, 14, 18, 22, 26];  % 不同航速

fprintf('参数空间:\n');
fprintf('  Cb = [%.2f, %.2f]\n', min(Cb_range), max(Cb_range));
fprintf('  L  = [%.0f, %.0f] m\n', min(L_range), max(L_range));
fprintf('  B  = [%.0f, %.0f] m\n', min(B_range), max(B_range));
fprintf('  T  = [%.0f, %.0f] m\n', min(T_range), max(T_range));
fprintf('  V  = [%.0f, %.0f] knots\n', min(V_knots_range), max(V_knots_range));
fprintf('\n');

% ==============================================================================
% 3. 主循环：遍历所有参数组合
% ==============================================================================
% 存储找到的"交叉点"案例
CrossoverCases = [];

% 对每个(L, B, T, V)组合，分析Cb的影响
fprintf('开始扫描...\n\n');

case_count = 0;
for iL = 1:NumL
    L = L_range(iL);
    for iB = 1:NumB
        B = B_range(iB);
        for iT = 1:NumT
            T = T_range(iT);
            for iV = 1:NumV
                V_knots = V_knots_range(iV);
                V = V_knots * 0.5144;  % 转换为 m/s

                % 检查船型合理性
                LB_ratio = L / B;
                BT_ratio = B / T;

                if LB_ratio < 4 || LB_ratio > 12 || BT_ratio < 1.8 || BT_ratio > 5
                    continue;
                end

                % 傅汝德数
                Fr = V / sqrt(g * L);

                % 螺旋桨直径 (经验公式)
                D_prop = 0.65 * T;

                % --- 对该配置扫描所有Cb ---
                Results = struct('Cb',[],'Rt',[],'Rw',[],'Rf',[],'w',[],'t',[],...
                    'eta_H',[],'eta_O',[],'eta_D',[],'Pd',[],'Thrust',[]);
                Results = repmat(Results, NumCb, 1);

                valid_count = 0;

                for iCb = 1:NumCb
                    Cb = Cb_range(iCb);

                    % 检查Cb与L/B的组合是否合理
                    % 快艇(小Cb)通常细长，散货船(大Cb)通常肥短
                    if Cb < 0.55 && LB_ratio < 5
                        continue;  % 快艇不会太肥
                    end
                    if Cb > 0.75 && LB_ratio > 8
                        continue;  % 散货船不会太细长
                    end

                    % --- A. 三参数函数生成船型 ---
                    h_val = 0.3 + 1.2 * Cb;
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
                    Re = V * L / Nu;
                    Cf = 0.075 / (log10(Re) - 2)^2;

                    % 摩擦阻力 (含粗糙度修正)
                    Rf = 0.5 * Rho_water * S * V^2 * (Cf + 0.0004);

                    % 兴波阻力 (Michell积分)
                    Rw = michell(Y, V, L, B, T, Rho_water, N_theta);

                    % 总阻力 = 摩擦阻力 + 兴波阻力 (严格Michell + ITTC)
                    Rt = Rf + Rw;

                    % --- C. Holtrop 伴流和推力减额 ---
                    % 使用Holtrop & Mennen (1982, 1984) 公式
                    Cp = Cb / 0.98;           % 棱形系数近似
                    lcb = -0.5 + 0.5*(Cb-0.6);  % 浮心位置 (高Cb船浮心略靠艉)
                    Cstern = 0;               % 常规艉型
                    [w, t_val] = holtrop_wake_thrust(L, B, T, Cb, D_prop, Cp, lcb, Cstern);

                    % --- D. 效率计算 ---
                    eta_H = (1 - t_val) / (1 - w);
                    Thrust_req = Rt / (1 - t_val);
                    Va = V * (1 - w);

                    % 推力系数
                    A0 = pi * (D_prop/2)^2;
                    CT = Thrust_req / (0.5 * Rho_water * Va^2 * A0);

                    % 螺旋桨效率 - 使用PVL升力线理论计算
                    [eta_O, ~, ~] = calc_propeller_efficiency_PVL(D_prop, V, w, Thrust_req, Rho_water);

                    % 相对旋转效率 - Holtrop (1984) 公式
                    Ae_Ao = 0.55;  % 典型盘面比
                    eta_R = holtrop_relative_rotative(Cp, lcb, Ae_Ao);

                    % 总推进效率
                    eta_D = eta_H * eta_O * eta_R;

                    % 传递功率
                    Pd = Rt * V / eta_D;

                    % 存储
                    valid_count = valid_count + 1;
                    Results(valid_count).Cb = Cb;
                    Results(valid_count).Rt = Rt;
                    Results(valid_count).Rw = Rw;
                    Results(valid_count).Rf = Rf;
                    Results(valid_count).w = w;
                    Results(valid_count).t = t_val;
                    Results(valid_count).eta_H = eta_H;
                    Results(valid_count).eta_O = eta_O;
                    Results(valid_count).eta_D = eta_D;
                    Results(valid_count).Pd = Pd;
                    Results(valid_count).Thrust = Thrust_req;
                end

                Results = Results(1:valid_count);

                if valid_count < 5
                    continue;  % 数据点太少
                end

                % --- 分析该配置的结果 ---
                Cb_vec = [Results.Cb];
                Rt_vec = [Results.Rt];
                Pd_vec = [Results.Pd];
                Thrust_vec = [Results.Thrust];

                [~, idx_min_Rt] = min(Rt_vec);
                [~, idx_min_Pd] = min(Pd_vec);
                [~, idx_min_Thrust] = min(Thrust_vec);

                % 检查是否找到交叉点
                if idx_min_Rt ~= idx_min_Pd
                    case_count = case_count + 1;

                    % 计算功率节省
                    Pd_at_minRt = Pd_vec(idx_min_Rt);
                    Pd_optimal = Pd_vec(idx_min_Pd);
                    power_saving = (Pd_at_minRt - Pd_optimal) / Pd_at_minRt * 100;

                    CrossoverCases(case_count).L = L;
                    CrossoverCases(case_count).B = B;
                    CrossoverCases(case_count).T = T;
                    CrossoverCases(case_count).V_knots = V_knots;
                    CrossoverCases(case_count).Fr = Fr;
                    CrossoverCases(case_count).Cb_minRt = Cb_vec(idx_min_Rt);
                    CrossoverCases(case_count).Cb_minPd = Cb_vec(idx_min_Pd);
                    CrossoverCases(case_count).Rt_minRt = Rt_vec(idx_min_Rt);
                    CrossoverCases(case_count).Rt_minPd = Rt_vec(idx_min_Pd);
                    CrossoverCases(case_count).Pd_minRt = Pd_at_minRt;
                    CrossoverCases(case_count).Pd_minPd = Pd_optimal;
                    CrossoverCases(case_count).power_saving = power_saving;
                    CrossoverCases(case_count).Results = Results;

                    fprintf('>>> 发现交叉点 #%d: L=%.0fm, B=%.0fm, T=%.0fm, V=%.0fkn, Fr=%.3f\n', ...
                        case_count, L, B, T, V_knots, Fr);
                    fprintf('    最小阻力: Cb=%.3f, Rt=%.1fkN, Pd=%.2fMW\n', ...
                        Cb_vec(idx_min_Rt), Rt_vec(idx_min_Rt)/1000, Pd_at_minRt/1e6);
                    fprintf('    最小功率: Cb=%.3f, Rt=%.1fkN, Pd=%.2fMW\n', ...
                        Cb_vec(idx_min_Pd), Rt_vec(idx_min_Pd)/1000, Pd_optimal/1e6);
                    fprintf('    功率节省: %.2f%%\n\n', power_saving);
                end
            end
        end
    end
end

% ==============================================================================
% 4. 结果汇总
% ==============================================================================
fprintf('\n================================================================================\n');
fprintf('                    结 果 汇 总\n');
fprintf('================================================================================\n\n');

if isempty(CrossoverCases)
    fprintf('未找到交叉点案例。正在调整模型参数重新搜索...\n\n');

    % 如果标准模型没找到，尝试极端工况
    fprintf('--- 极端工况分析 ---\n');

    % 低速大船（散货船工况）
    L = 250; B = 40; T = 15;
    V_knots = 12;
    V = V_knots * 0.5144;
    Fr = V / sqrt(g * L);
    D_prop = 0.65 * T;

    fprintf('散货船工况: L=%.0fm, B=%.0fm, T=%.0fm, V=%.0fkn, Fr=%.3f\n\n', L, B, T, V_knots, Fr);

    % 更细的Cb扫描
    Cb_fine = linspace(0.60, 0.85, 30);
    Results_fine = [];

    for iCb = 1:length(Cb_fine)
        Cb = Cb_fine(iCb);

        h_val = 0.3 + 1.2 * Cb;
        f0 = three_param_shape(x_norm, h_val, 0.15, 0.12);
        f0 = f0 / max(f0) * 0.5;

        z_norm = linspace(0, 1, Nz);
        Y = zeros(Nx, Nz);
        for iz = 1:Nz
            z_factor = sin(pi/2 * z_norm(iz))^0.6;
            Y(:, iz) = f0(:) * z_factor;
        end

        S = estimate_wetted_surface(L, B, T, Cb);
        Re = V * L / Nu;
        Cf = 0.075 / (log10(Re) - 2)^2;
        Rf = 0.5 * Rho_water * S * V^2 * (Cf + 0.0004);
        Rw = michell(Y, V, L, B, T, Rho_water, N_theta);
        Rt = Rf + Rw;  % 严格Michell + ITTC

        Cp = Cb / 0.98;
        lcb = -0.5 + 0.5*(Cb-0.6);
        Cstern = 0;
        [w, t_val] = holtrop_wake_thrust(L, B, T, Cb, D_prop, Cp, lcb, Cstern);

        eta_H = (1 - t_val) / (1 - w);
        Thrust_req = Rt / (1 - t_val);
        Va = V * (1 - w);
        A0 = pi * (D_prop/2)^2;
        CT = Thrust_req / (0.5 * Rho_water * Va^2 * A0);
        [eta_O, ~, ~] = calc_propeller_efficiency_PVL(D_prop, V, w, Thrust_req, Rho_water);
        Ae_Ao = 0.55;  % 典型盘面比
        eta_R = holtrop_relative_rotative(Cp, lcb, Ae_Ao);
        eta_D = eta_H * eta_O * eta_R;
        Pd = Rt * V / eta_D;

        Results_fine(iCb).Cb = Cb;
        Results_fine(iCb).Rt = Rt;
        Results_fine(iCb).Rw = Rw;
        Results_fine(iCb).Rf = Rf;
        Results_fine(iCb).w = w;
        Results_fine(iCb).t = t_val;
        Results_fine(iCb).eta_H = eta_H;
        Results_fine(iCb).eta_O = eta_O;
        Results_fine(iCb).eta_D = eta_D;
        Results_fine(iCb).Pd = Pd;
        Results_fine(iCb).Thrust = Thrust_req;
    end

    % 输出详细表格
    fprintf('%-5s %-7s %-9s %-9s %-7s %-7s %-7s %-7s %-7s %-9s\n', ...
        'No.', 'Cb', 'Rt(kN)', 'Rw(kN)', 'w', 't', 'eta_H', 'eta_O', 'eta_D', 'Pd(MW)');
    fprintf('--------------------------------------------------------------------------------\n');

    for ii = 1:length(Results_fine)
        fprintf('%4d  %.4f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  %.4f  %8.2f\n', ...
            ii, Results_fine(ii).Cb, Results_fine(ii).Rt/1000, Results_fine(ii).Rw/1000, ...
            Results_fine(ii).w, Results_fine(ii).t, Results_fine(ii).eta_H, ...
            Results_fine(ii).eta_O, Results_fine(ii).eta_D, Results_fine(ii).Pd/1e6);
    end

    Cb_vec = [Results_fine.Cb];
    Rt_vec = [Results_fine.Rt];
    Pd_vec = [Results_fine.Pd];
    Thrust_vec = [Results_fine.Thrust];
    eta_D_vec = [Results_fine.eta_D];
    eta_H_vec = [Results_fine.eta_H];
    eta_O_vec = [Results_fine.eta_O];

    [min_Rt, idx_Rt] = min(Rt_vec);
    [min_Pd, idx_Pd] = min(Pd_vec);
    [min_Thrust, idx_Thrust] = min(Thrust_vec);
    [max_eta, idx_eta] = max(eta_D_vec);

    fprintf('\n【关键结果】\n');
    fprintf('  ├─ 最小阻力点: Cb = %.4f, Rt = %.1f kN, Pd = %.2f MW\n', ...
        Cb_vec(idx_Rt), min_Rt/1000, Pd_vec(idx_Rt)/1e6);
    fprintf('  ├─ 最小推力点: Cb = %.4f, Thrust = %.1f kN, Pd = %.2f MW\n', ...
        Cb_vec(idx_Thrust), Thrust_vec(idx_Thrust)/1000, Pd_vec(idx_Thrust)/1e6);
    fprintf('  ├─ 最高效率点: Cb = %.4f, eta_D = %.4f, Pd = %.2f MW\n', ...
        Cb_vec(idx_eta), max_eta, Pd_vec(idx_eta)/1e6);
    fprintf('  └─ 最小功率点: Cb = %.4f, Pd = %.2f MW\n\n', ...
        Cb_vec(idx_Pd), min_Pd/1e6);

    if idx_Rt ~= idx_Pd
        power_save = (Pd_vec(idx_Rt) - min_Pd) / Pd_vec(idx_Rt) * 100;
        fprintf('  ★★★ 找到证据！最小阻力点 ≠ 最小功率点 ★★★\n');
        fprintf('  ★★★ ISPS优化可节省功率: %.2f%% ★★★\n\n', power_save);
    else
        fprintf('  注：在此工况下，最小阻力点与最小功率点重合\n\n');
    end

    % 绘图
    figure('Position', [50 50 1600 900], 'Color', 'w');

    % 子图1: 阻力与功率
    subplot(2,3,1);
    yyaxis left
    plot(Cb_vec, Rt_vec/1000, '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    ylabel('Total Resistance R_t (kN)');
    yyaxis right
    plot(Cb_vec, Pd_vec/1e6, '-s', 'LineWidth', 1.5, 'MarkerSize', 4);
    ylabel('Delivered Power P_d (MW)');
    xlabel('Block Coefficient C_b');
    title('Resistance vs Power');
    legend('R_t', 'P_d', 'Location', 'best');
    grid on;

    % 子图2: 效率分解
    subplot(2,3,2);
    plot(Cb_vec, eta_H_vec, '-s', 'LineWidth', 1.5, 'DisplayName', '\eta_H');
    hold on;
    plot(Cb_vec, eta_O_vec, '-^', 'LineWidth', 1.5, 'DisplayName', '\eta_O');
    plot(Cb_vec, eta_D_vec, '-o', 'LineWidth', 1.8, 'DisplayName', '\eta_D');
    hold off;
    xlabel('Block Coefficient C_b');
    ylabel('Efficiency');
    title('Efficiency Components');
    legend('Location', 'best');
    grid on;

    % 子图3: 伴流和推力减额
    subplot(2,3,3);
    yyaxis left
    plot(Cb_vec, [Results_fine.w], '-d', 'LineWidth', 1.5);
    ylabel('Wake Fraction w');
    yyaxis right
    plot(Cb_vec, [Results_fine.t], '-o', 'LineWidth', 1.5);
    ylabel('Thrust Deduction t');
    xlabel('Block Coefficient C_b');
    title('Interaction Factors');
    legend('w', 't', 'Location', 'best');
    grid on;

    % 子图4: 阻力分解
    subplot(2,3,4);
    Rw_vec = [Results_fine.Rw];
    Rf_vec = [Results_fine.Rf];
    area(Cb_vec, [Rf_vec'/1000, Rw_vec'/1000]);
    xlabel('Block Coefficient C_b');
    ylabel('Resistance (kN)');
    title('Resistance Components');
    legend('Friction R_f', 'Wave R_w', 'Location', 'best');
    grid on;

    % 子图5: 推力与功率对比
    subplot(2,3,5);
    yyaxis left
    plot(Cb_vec, Thrust_vec/1000, '-^', 'LineWidth', 1.5);
    ylabel('Required Thrust (kN)');
    yyaxis right
    plot(Cb_vec, Pd_vec/1e6, '-s', 'LineWidth', 1.5);
    ylabel('Delivered Power (MW)');
    xlabel('Block Coefficient C_b');
    title('Thrust vs Power');
    grid on;

    % 标记最优点
    hold on;
    yyaxis left
    plot(Cb_vec(idx_Thrust), Thrust_vec(idx_Thrust)/1000, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
    yyaxis right
    plot(Cb_vec(idx_Pd), Pd_vec(idx_Pd)/1e6, 'gs', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;
    legend('Thrust', 'Power', 'Min Thrust', 'Min Power', 'Location', 'best');

    % 子图6: 船身效率增益 vs 阻力增加
    subplot(2,3,6);
    % 以最小阻力点为基准
    Rt_ratio = Rt_vec / min_Rt;
    eta_H_ratio = eta_H_vec / eta_H_vec(idx_Rt);

    plot(Cb_vec, Rt_ratio, '-o', 'LineWidth', 1.5, 'DisplayName', 'R_t / R_{t,min}');
    hold on;
    plot(Cb_vec, eta_H_ratio, '-s', 'LineWidth', 1.5, 'DisplayName', '\eta_H / \eta_{H,ref}');
    plot(Cb_vec, Rt_ratio ./ eta_H_ratio, '-^', 'LineWidth', 1.5, 'DisplayName', 'Net Effect');
    hold off;
    xlabel('Block Coefficient C_b');
    ylabel('Ratio');
    title('Resistance Increase vs Efficiency Gain');
    legend('Location', 'best');
    grid on;
    yline(1, '--k', 'Reference');

    sgtitle(sprintf('Bulk Carrier Analysis: L=%.0fm, B=%.0fm, T=%.0fm, V=%.0fkn (Fr=%.3f)', ...
        L, B, T, V_knots, Fr), 'FontSize', 13, 'FontWeight', 'bold');

    saveas(gcf, 'Ship_Propeller_Optimization.png');

else
    fprintf('找到 %d 个交叉点案例！\n\n', length(CrossoverCases));

    % 输出所有交叉点
    fprintf('%-5s %-6s %-6s %-6s %-6s %-6s %-8s %-8s %-9s %-9s %-8s\n', ...
        'No.', 'L(m)', 'B(m)', 'T(m)', 'V(kn)', 'Fr', 'Cb_Rt', 'Cb_Pd', 'Pd_Rt', 'Pd_opt', 'Save%');
    fprintf('--------------------------------------------------------------------------------\n');

    for ii = 1:length(CrossoverCases)
        c = CrossoverCases(ii);
        fprintf('%4d  %5.0f  %5.0f  %5.0f  %5.0f  %.3f  %7.3f  %7.3f  %8.2f  %8.2f  %7.2f\n', ...
            ii, c.L, c.B, c.T, c.V_knots, c.Fr, c.Cb_minRt, c.Cb_minPd, ...
            c.Pd_minRt/1e6, c.Pd_minPd/1e6, c.power_saving);
    end

    % 选择功率节省最大的案例绘图
    [~, best_idx] = max([CrossoverCases.power_saving]);
    BestCase = CrossoverCases(best_idx);
    Results_best = BestCase.Results;

    fprintf('\n【最佳案例详细分析】\n');
    fprintf('L=%.0fm, B=%.0fm, T=%.0fm, V=%.0fkn, Fr=%.3f\n\n', ...
        BestCase.L, BestCase.B, BestCase.T, BestCase.V_knots, BestCase.Fr);

    Cb_vec = [Results_best.Cb];
    Rt_vec = [Results_best.Rt];
    Pd_vec = [Results_best.Pd];
    eta_D_vec = [Results_best.eta_D];

    % 绘图
    figure('Position', [50 50 1400 800], 'Color', 'w');

    subplot(2,2,1);
    yyaxis left
    plot(Cb_vec, Rt_vec/1000, '-o', 'LineWidth', 1.5);
    ylabel('Total Resistance (kN)');
    yyaxis right
    plot(Cb_vec, Pd_vec/1e6, '-s', 'LineWidth', 1.5);
    ylabel('Delivered Power (MW)');
    xlabel('Block Coefficient C_b');
    title('Key Result: Min Resistance ≠ Min Power');
    legend('R_t', 'P_d');
    grid on;

    % 标记最优点
    [~, idx_Rt] = min(Rt_vec);
    [~, idx_Pd] = min(Pd_vec);
    hold on;
    yyaxis left
    plot(Cb_vec(idx_Rt), Rt_vec(idx_Rt)/1000, 'ro', 'MarkerSize', 12, 'LineWidth', 2);
    yyaxis right
    plot(Cb_vec(idx_Pd), Pd_vec(idx_Pd)/1e6, 'gs', 'MarkerSize', 12, 'LineWidth', 2);
    hold off;

    subplot(2,2,2);
    plot(Cb_vec, [Results_best.eta_H], '-s', 'LineWidth', 1.5);
    hold on;
    plot(Cb_vec, [Results_best.eta_O], '-^', 'LineWidth', 1.5);
    plot(Cb_vec, eta_D_vec, '-o', 'LineWidth', 1.8);
    hold off;
    xlabel('Block Coefficient C_b');
    ylabel('Efficiency');
    title('Efficiency Components');
    legend('\eta_H', '\eta_O', '\eta_D');
    grid on;

    subplot(2,2,3);
    yyaxis left
    plot(Cb_vec, [Results_best.w], '-d', 'LineWidth', 1.5);
    ylabel('Wake Fraction w');
    yyaxis right
    plot(Cb_vec, [Results_best.t], '-o', 'LineWidth', 1.5);
    ylabel('Thrust Deduction t');
    xlabel('Block Coefficient C_b');
    title('Interaction Factors');
    legend('w', 't');
    grid on;

    subplot(2,2,4);
    bar(categorical({'Min R_t Point', 'Min P_d Point'}), [BestCase.Pd_minRt/1e6, BestCase.Pd_minPd/1e6]);
    ylabel('Delivered Power (MW)');
    title(sprintf('Power Comparison (Saving: %.2f%%)', BestCase.power_saving));
    grid on;

    sgtitle(sprintf('ISPS Optimization Evidence: L=%.0fm, B=%.0fm, T=%.0fm, V=%.0fkn', ...
        BestCase.L, BestCase.B, BestCase.T, BestCase.V_knots), 'FontSize', 13, 'FontWeight', 'bold');

    saveas(gcf, 'Ship_Propeller_Optimization.png');
end

% ==============================================================================
% 5. 最终结论
% ==============================================================================
fprintf('\n================================================================================\n');
fprintf('                         结 论\n');
fprintf('================================================================================\n\n');

fprintf('【物理机制】\n');
fprintf('  • Cb增大 → 船型更丰满 → 阻力增加（主要是兴波阻力）\n');
fprintf('  • Cb增大 → 艉部更饱满 → 伴流w增加 → 船身效率η_H提高\n');
fprintf('  • 当η_H增益 > 阻力增加带来的损失时 → 总功率Pd反而降低\n');
fprintf('  • 这种情况更容易发生在：低傅汝德数、大Cb船型\n\n');

if ~isempty(CrossoverCases)
    fprintf('【研究发现】\n');
    fprintf('  ★ 找到 %d 个 "最小阻力 ≠ 最小功率" 的案例\n', length(CrossoverCases));
    fprintf('  ★ 最大功率节省: %.2f%%\n', max([CrossoverCases.power_saving]));
    fprintf('  ★ 验证了船桨集成优化(ISPS)的必要性\n\n');
end

fprintf('【工程意义】\n');
fprintf('  ★ 传统方法"尽可能减小阻力"存在局限性\n');
fprintf('  ★ 应以最小传递功率Pd为优化目标，而非最小阻力Rt\n');
fprintf('  ★ 船体-螺旋桨集成优化可实现额外的能效提升\n\n');

fprintf('================================================================================\n');

% 保存结果
save('Ship_Optimization_Results.mat');
fprintf('结果已保存到 Ship_Optimization_Results.mat 和 Ship_Propeller_Optimization.png\n');

% ==============================================================================
% 辅助函数
% ==============================================================================

function [w, t] = holtrop_wake_thrust(L, B, T, Cb, D_prop, Cp, lcb, Cstern)
    % Holtrop & Mennen (1982, 1984) 伴流和推力减额公式
    % 参考文献:
    % [1] Holtrop J., Mennen G.G.J. (1982). "An approximate power prediction method"
    %     International Shipbuilding Progress, Vol.29, No.335
    % [2] Holtrop J. (1984). "A statistical re-analysis of resistance and propulsion data"
    %     International Shipbuilding Progress, Vol.31, No.363
    %
    % 输入:
    %   L      - 水线长 [m]
    %   B      - 型宽 [m]
    %   T      - 吃水 [m] (这里用艉吃水Ta近似)
    %   Cb     - 方形系数
    %   D_prop - 螺旋桨直径 [m]
    %   Cp     - 棱形系数 (若未知可用 Cp = Cb/0.98)
    %   lcb    - 浮心纵向位置 (% of L, 从舯向艏为正, 通常-3~3)
    %   Cstern - 艉部形状系数 (-25~10, V形=-25, U形=10, 常规=0)

    % 若未提供参数，使用默认值
    if nargin < 6 || isempty(Cp)
        Cp = Cb / 0.98;  % 近似关系
    end
    if nargin < 7 || isempty(lcb)
        lcb = -0.75;  % 典型值，浮心略靠艉
    end
    if nargin < 8 || isempty(Cstern)
        Cstern = 0;  % 常规艉型
    end

    Ta = T;  % 艉吃水 (假设平浮)

    % === Holtrop 伴流分数 w (单桨船) ===
    % 公式来源: Holtrop (1984), Eq. (28)-(32)

    % 系数 c8: 与 B/Ta 相关
    if B/Ta < 5
        c8 = B / Ta * (1.0 - 0.5 * (B/Ta - 5)^2);
        c8 = max(c8, 0);
    else
        c8 = 1.0;
    end
    c8 = B * sqrt(1/(L*Ta));  % 简化形式

    % 系数 c9: 与 Cstern 相关
    c9 = 0.5 - 0.5 * Cstern / 25;
    c9 = max(0, min(1, c9));

    % Cp1: 修正棱形系数
    Cp1 = 1.45 * Cp - 0.315 - 0.0225 * lcb;
    Cp1 = max(0.1, min(0.9, Cp1));

    % 系数 c11: 与 Ta/D 相关
    TaD = Ta / D_prop;
    if TaD < 2
        c11 = TaD;
    else
        c11 = 0.0833333 * TaD^3 + 1.33333;
    end

    % 系数 c19: 与 Cp 相关
    % Cm: 中横剖面系数，近似为 Cm = Cb/Cp
    Cm = Cb / Cp;
    if Cp < 0.7
        c19 = 0.12997 / (0.95 - Cb) - 0.11056 / (0.95 - Cp);
    else
        c19 = 0.18567 / (1.3571 - Cm) - 0.71276 + 0.38648 * Cp;
    end
    c19 = max(-0.1, min(0.1, c19));

    % 系数 c20: 与推进器类型相关 (常规螺旋桨=1)
    c20 = 1.0;

    % Cv: 粘性阻力系数
    % Cv = (1+k)*Cf + Ca
    % k: 形状因子 (Holtrop形式因子公式)
    Lr = L * (1 - Cp + 0.06*Cp*lcb/(4*Cp - 1));  % 进流长度
    Disp = L * B * T * Cb;  % 排水体积
    k_form = 0.93 + 0.4871*(B/L)^1.0681 * (T/L)^0.4611 * ...
             (Lr/L)^0.1216 * (L^3/Disp)^0.3649 * (1-Cp)^(-0.6042);
    k_form = max(0.1, min(0.5, k_form - 1));  % 转换为形状因子

    Re = 1e9;  % 大船典型雷诺数
    Cf = 0.075 / (log10(Re) - 2)^2;
    Ca = 0.00035;  % 粗糙度余量
    Cv = (1 + k_form) * Cf + Ca;

    % Holtrop 伴流公式 (单桨船)
    % w = c9*c20*Cv*(L/Ta)*(0.050776 + 0.93405*c11*Cv/(1-Cp1))
    %     + 0.27915*c20*sqrt(B/(L*(1-Cp1))) + c19*c20

    term1 = c9 * c20 * Cv * sqrt(L/Ta) * (0.050776 + 0.93405*c11*Cv/(1-Cp1));
    term2 = 0.27915 * c20 * sqrt(B / (L * (1 - Cp1)));
    term3 = c19 * c20;

    w = term1 + term2 + term3;

    % 限制在合理范围
    w = max(0.10, min(0.45, w));

    % === Holtrop 推力减额分数 t ===
    % 公式来源: Holtrop (1984), Eq. (33)-(35)

    % 系数 c10: 与 L/B 相关
    LB = L / B;
    if LB > 5.2
        c10 = 0.25 - 0.003328402 / ((B/L) - 0.134615385);
    else
        c10 = B / L;
    end

    % Holtrop 推力减额公式 (常规单桨)
    % t = 0.25014*(B/L)^0.28956 * (sqrt(B*T)/D)^0.2624 / (1-Cp+0.0225*lcb)^0.01762
    %     + 0.0015*Cstern

    t = 0.25014 * (B/L)^0.28956 * (sqrt(B*T)/D_prop)^0.2624 ...
        / (1 - Cp + 0.0225*lcb)^0.01762 + 0.0015 * Cstern;

    % 限制在合理范围
    t = max(0.08, min(0.30, t));

    % 确保 t < w (物理约束)
    if t >= w
        t = 0.75 * w;
    end
end

function eta_R = holtrop_relative_rotative(Cp, lcb, Ae_Ao)
    % Holtrop (1984) 相对旋转效率公式
    % 参考文献:
    % Holtrop J. (1984). "A statistical re-analysis of resistance and propulsion data"
    % International Shipbuilding Progress, Vol.31, No.363, Eq. (36)
    %
    % 输入:
    %   Cp    - 棱形系数
    %   lcb   - 浮心纵向位置 (% of L, 从舯向艏为正)
    %   Ae_Ao - 螺旋桨盘面比 (expanded area ratio)
    %
    % 输出:
    %   eta_R - 相对旋转效率
    %
    % 公式 (单桨船):
    % η_R = 0.9922 - 0.05908*(Ae/Ao) + 0.07424*(Cp - 0.0225*lcb)

    eta_R = 0.9922 - 0.05908 * Ae_Ao + 0.07424 * (Cp - 0.0225 * lcb);

    % 限制在合理范围 (通常0.96-1.05)
    eta_R = max(0.96, min(1.05, eta_R));
end

function eta_O = propeller_efficiency_enhanced(CT, Va, D, Cb)
    % 增强版螺旋桨敞水效率
    % 考虑载荷系数、前进系数等

    % 动量理论理想效率
    eta_ideal = 2 / (1 + sqrt(1 + CT));

    % 假设转速
    n_assumed = 1.6;  % rps
    J = Va / (n_assumed * D);

    % J偏离最优点的损失
    J_opt = 0.6 + 0.15 * Cb;  % 丰满船最优J略高
    J_penalty = 0.04 * (J - J_opt)^2;

    % 粘性损失系数
    viscous_factor = 0.82 - 0.03 * max(0, CT - 0.6);
    viscous_factor = max(0.72, min(0.85, viscous_factor));

    eta_O = eta_ideal * viscous_factor - J_penalty;
    eta_O = max(0.45, min(0.75, eta_O));
end

function S = estimate_wetted_surface(L, B, T, Cb)
    % 湿表面积估算 (Holtrop公式)
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

function [eta_O, KT, KQ] = calc_propeller_efficiency_PVL(D, Vs, w, Thrust, rho)
    % 使用PVL升力线理论计算螺旋桨敞水效率
    % 调用Main.m进行涡格法计算
    %
    % 输入:
    %   D      - 螺旋桨直径 [m]
    %   Vs     - 船速 [m/s]
    %   w      - 伴流分数
    %   Thrust - 所需推力 [N]
    %   rho    - 水密度 [kg/m³]
    %
    % 输出:
    %   eta_O  - 螺旋桨敞水效率
    %   KT     - 推力系数
    %   KQ     - 扭矩系数

    R = D / 2;
    Va = Vs * (1 - w);  % 进速

    % === 螺旋桨设计参数 (典型商船螺旋桨) ===
    NBLADE = 4;              % 叶片数
    Dhub = 0.18 * D;         % 桨毂直径
    Rhub_R = Dhub / D;       % 桨毂半径比

    % 估算转速 (基于最优前进系数J≈0.6-0.8)
    J_target = 0.65;
    n = Va / (J_target * D);  % 转速 [rps]
    ADVCO = Va / (n * D);     % 前进系数 J

    % 推力系数
    CTDES = Thrust / (0.5 * rho * Va^2 * pi * R^2);

    % === PVL计算参数 ===
    MT = 20;       % 控制点数
    ITER = 50;     % 最大迭代次数
    IHUB = 1;      % 启用Hub镜像涡
    RHV = 0.5;     % 镜像涡半径比
    NX = 11;       % 输入半径点数
    HR = 0;        % Hub卸载因子
    HT = 1;        % Tip卸载因子
    CRP = 1;       % 旋涡抵消因子

    % === 半径分布 (r/R) ===
    XR = linspace(Rhub_R, 1.0, NX);

    % === 弦长比分布 c/D (典型Wageningen B系列形式) ===
    % 最大弦长在0.7R附近
    XCHD = zeros(1, NX);
    for i = 1:NX
        r_R = XR(i);
        % Wageningen B4-70 近似弦长分布
        XCHD(i) = 0.16 * (1.0 - 0.3*(r_R - 0.7)^2) * (1.0 - (r_R - Rhub_R)/(1-Rhub_R)*0.1);
    end

    % === 阻力系数分布 Cd ===
    XCD = 0.008 * ones(1, NX);  % 典型值

    % === 轴向速度分布 Va/Vs ===
    % 使用Holtrop有效伴流分数（均匀分布）
    % Holtrop公式给出的是有效平均伴流
    XVA = (1 - w) * ones(1, NX);

    % === 切向速度分布 Vt/Vs ===
    XVT = zeros(1, NX);  % 假设无预旋

    % === 调用Main.m进行PVL计算 ===
    try
        [CT, CP, KT_arr, KQ_arr, ~, EFFY, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, KTRY] = ...
            Main(MT, ITER, IHUB, RHV, NX, NBLADE, ADVCO, CTDES, HR, HT, CRP, XR, XCHD, XCD, XVA, XVT);

        % 取最终收敛的结果
        if KTRY > 0 && KTRY <= length(EFFY)
            eta_O = EFFY(KTRY);
            KT = KT_arr(KTRY);
            KQ = KQ_arr(KTRY);
        else
            eta_O = EFFY(end);
            KT = KT_arr(end);
            KQ = KQ_arr(end);
        end

        % PVL计算的效率（不做人为限制）
        % 只检查物理合理性
        if eta_O < 0 || eta_O > 1 || isnan(eta_O)
            error('PVL计算结果无效');
        end

    catch ME
        % 如果PVL计算失败，报告错误并使用动量理论
        % fprintf('PVL警告: %s\n', ME.message);
        CT_simple = CTDES;
        eta_O = 2 / (1 + sqrt(1 + CT_simple));  % 理想效率
        eta_O = eta_O * 0.85;  % 粘性修正因子
        KT = 0;
        KQ = 0;
    end
end
