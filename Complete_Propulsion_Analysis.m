% Complete_Propulsion_Analysis.m
% 船舶完整推进效率分析 - 参数化研究
%
% 本程序直接调用:
% 1. michell.m      - Michell积分计算兴波阻力
% 2. ITTC 1957公式  - 摩擦阻力计算
% 3. PVL.m          - 螺旋桨敞水效率计算(升力线理论) [必须使用]
% 4. Holtrop_Propulsion_Calculation_1982.m - 伴流分数(w)、推力减额分数(t)、相对旋转效率(eta_R)
%
% 可调参数: 船长L, 船宽B, 方形系数Cb, 船速Vs, 吃水T
%
% 核心结论: 阻力越小，推进效率不一定越高
%
% 作者: Claude Code
% 日期: 2026-02-01

clear; close all; clc;

fprintf('================================================================================\n');
fprintf('          船舶完整推进效率分析系统\n');
fprintf('          Michell积分 + ITTC + PVL + Holtrop (1982) 完全集成\n');
fprintf('================================================================================\n\n');

% ==============================================================================
% 1. 基准船舶参数设置 (基于KCS集装箱船)
% ==============================================================================
fprintf('>>> 设置基准船舶参数...\n');

Ship.L_ref = 230.0;       % 参考船长 [m]
Ship.B_ref = 32.2;        % 参考船宽 [m]
Ship.T_ref = 10.8;        % 参考吃水 [m]
Ship.Cb_ref = 0.6505;     % 参考方形系数
Ship.Vs_knots_ref = 24;   % 参考航速 [knots]
Ship.D_prop = 7.9;        % 螺旋桨直径 [m]
Ship.AE_AO = 0.80;        % 展开面积比 (典型值)
Ship.Cstern = 0;          % 艉部形状系数 (正常艉部)
Ship.LCB_percent = -2.0;  % 纵向浮心位置 [%L, 艉向为负]

% 物理常数
Rho_water = 1025;         % 海水密度 [kg/m³]
Nu = 1.188e-6;            % 运动粘度 [m²/s] (15°C)
g = 9.80665;              % 重力加速度 [m/s²]

% Michell积分网格参数
Nx = 101;                 % 站点数 (必须为奇数)
Nz = 40;                  % 水线数
N_theta = 80;             % 传播角采样数

% ==============================================================================
% 2. PVL螺旋桨参数配置
% ==============================================================================
fprintf('>>> 配置PVL螺旋桨参数...\n');

% 单桨设计参数
Single_def1 = zeros(1, 3);
Single_def1(1) = 5;           % NBLADE - 叶片数
Single_def1(2) = 105;         % N - 转速 [RPM]
Single_def1(3) = Ship.D_prop; % D - 直径 [m]

Single_def2 = zeros(1, 4);
Single_def2(1) = 8.0;         % H - 轴心深度 [m]
Single_def2(2) = 0.1;         % dV - 速度扰动 [m/s]
Single_def2(3) = 1.0;         % AlphaI - 理想攻角 [deg]
Single_def2(4) = 60;          % NP - 剖面点数

% 公共参数 (将根据每个工况动态更新)
Common_Def_base = zeros(1, 11);
Common_Def_base(1) = 0;            % 占位
Common_Def_base(3) = 1.58;         % Dhub - 桨毂直径 [m]
Common_Def_base(4) = 20;           % MT - 控制点数
Common_Def_base(5) = 20;           % ITER - 最大迭代次数
Common_Def_base(6) = 0.15;         % RHV - hub镜像涡半径比
Common_Def_base(7) = 10;           % NX - 径向节点数
Common_Def_base(8) = 0;            % HR - hub卸载因子
Common_Def_base(9) = 0;            % HT - tip卸载因子
Common_Def_base(10) = 1;           % CRP - 旋涡抵消因子
Common_Def_base(11) = Rho_water;   % rho - 水密度

% 螺旋桨剖面参数
Rhub = Common_Def_base(3) / 2;
R_prop = Ship.D_prop / 2;
XR0 = [Rhub/R_prop, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0];

% 标准弦长分布 c/D
XCHD_def = [0.174, 0.202, 0.229, 0.257, 0.277, 0.280, 0.267, 0.231, 0.196, 0.0];

% 阻力系数分布
XCD_def = [0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008];

% 最大弯度比 f0/c
f0oc_def = [0.0315, 0.0220, 0.0180, 0.0155, 0.0132, 0.0115, 0.0100, 0.0085, 0.0075, 0.0];

% 最大厚度比 t0/c
t0oc_def = [0.2200, 0.1650, 0.1300, 0.1020, 0.0800, 0.0650, 0.0500, 0.0380, 0.0320, 0.003];

% 掠角分布 [deg]
skew_def = [0, 5, 10, 15, 20, 25, 30, 35, 38, 40];

% 倾角分布 Xs/D
rake_def = [0, 0.01, 0.02, 0.03, 0.035, 0.04, 0.045, 0.05, 0.052, 0.055];

% 平均线和厚度类型
Mean = 1;    % NACA a=0.8
Thick = 1;   % NACA 65A010

% ==============================================================================
% 3. 参数变化设计空间
% ==============================================================================
NumSamples = 8;          % 每个参数的采样点数

% 定义各参数的变化范围
Cb_range = linspace(0.60, 0.72, NumSamples);      % 方形系数变化
B_range = linspace(30.0, 35.0, NumSamples);       % 船宽变化 [m]
L_range = linspace(220, 240, NumSamples);         % 船长变化 [m]
Vs_range = linspace(20, 28, NumSamples);          % 船速变化 [knots]
T_range = linspace(9.5, 12.0, NumSamples);        % 吃水变化 [m]

% 基准x坐标(用于船型生成)
x_norm = linspace(0, 1, Nx);

% 有效结果计数
valid_count_Cb = 0;
valid_count_B = 0;
valid_count_L = 0;
valid_count_Vs = 0;
valid_count_T = 0;

% ==============================================================================
% 4. 分析1: 方形系数Cb变化
% ==============================================================================
fprintf('\n===== 分析1: 方形系数Cb变化对推进效率的影响 =====\n');
fprintf('%-5s %-7s %-9s %-9s %-9s %-7s %-7s %-7s %-7s %-7s %-7s\n', ...
    'No.', 'Cb', 'Rt[kN]', 'Rw[kN]', 'Rf[kN]', 'w', 't', 'eta_R', 'eta_H', 'eta_O', 'eta_D');
fprintf('---------------------------------------------------------------------------------------------\n');

Results_Cb = [];
L = Ship.L_ref;
B = Ship.B_ref;
T = Ship.T_ref;
Vs_knots = Ship.Vs_knots_ref;
V_ship = Vs_knots * 0.5144;

for i = 1:NumSamples
    Cb = Cb_range(i);

    % --- A. 生成船型偏移矩阵 ---
    Y = generate_hull_offsets(x_norm, Nz, Cb);

    % --- B. 阻力计算 ---
    % ITTC 1957 摩擦阻力
    S_wet = estimate_wetted_surface(L, B, T, Cb);
    Re = V_ship * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    K1 = 1 + calculate_form_factor(L, B, T, Cb);
    Rf = 0.5 * Rho_water * S_wet * V_ship^2 * Cf * K1;

    % Michell积分兴波阻力 (直接调用michell.m)
    Rw = michell(Y, V_ship, L, B, T, Rho_water, N_theta);
    Rw = max(0, Rw);

    Rt = Rf + Rw;

    % --- C. Holtrop推进因子 (直接调用) ---
    Cp = Cb / 0.98;
    [w, t_thrust, eta_R] = Holtrop_Propulsion_Calculation_1982(...
        L, B, T, Ship.D_prop, S_wet, Cb, Cp, Ship.LCB_percent, ...
        Ship.Cstern, Ship.AE_AO, Cf, K1);

    w = max(0.15, min(0.45, w));
    t_thrust = max(0.10, min(0.30, t_thrust));
    eta_R = max(0.95, min(1.05, eta_R));

    % --- D. 推进效率计算 ---
    eta_H = (1 - t_thrust) / (1 - w);
    Thrust_req = Rt / (1 - t_thrust);
    Va = V_ship * (1 - w);

    % 更新PVL参数
    Common_Def = Common_Def_base;
    Common_Def(2) = Va;

    % 构建伴流场速度分布
    XVA_def = ones(size(XR0)) * (1 - w) + 0.08*(1-XR0);
    XVT_def = zeros(size(XR0)) + 0.01 * (1-XR0);

    % 调用PVL计算螺旋桨敞水效率 (必须使用PVL.m)
    eta_O = call_PVL_for_efficiency(Common_Def, XR0, XCHD_def, XCD_def, ...
        XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
        Single_def1, Single_def2, Mean, Thick, Thrust_req);

    % 检查PVL计算是否有效
    if isnan(eta_O)
        fprintf('%4d  %.4f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  [PVL无效-跳过]\n', ...
            i, Cb, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H);
        continue;  % 跳过该工况
    end

    eta_D = eta_H * eta_O * eta_R;
    Pe = Rt * V_ship;
    Pd = Pe / eta_D;

    % 存储有效结果
    valid_count_Cb = valid_count_Cb + 1;
    Results_Cb(valid_count_Cb).Cb = Cb;
    Results_Cb(valid_count_Cb).Rt = Rt;
    Results_Cb(valid_count_Cb).Rw = Rw;
    Results_Cb(valid_count_Cb).Rf = Rf;
    Results_Cb(valid_count_Cb).w = w;
    Results_Cb(valid_count_Cb).t = t_thrust;
    Results_Cb(valid_count_Cb).eta_R = eta_R;
    Results_Cb(valid_count_Cb).eta_H = eta_H;
    Results_Cb(valid_count_Cb).eta_O = eta_O;
    Results_Cb(valid_count_Cb).eta_D = eta_D;
    Results_Cb(valid_count_Cb).Pe = Pe;
    Results_Cb(valid_count_Cb).Pd = Pd;

    fprintf('%4d  %.4f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        i, Cb, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H, eta_O, eta_D);

    close all hidden;
end

% ==============================================================================
% 5. 分析2: 船宽B变化
% ==============================================================================
fprintf('\n===== 分析2: 船宽B变化对推进效率的影响 =====\n');
fprintf('%-5s %-7s %-9s %-9s %-9s %-7s %-7s %-7s %-7s %-7s %-7s\n', ...
    'No.', 'B[m]', 'Rt[kN]', 'Rw[kN]', 'Rf[kN]', 'w', 't', 'eta_R', 'eta_H', 'eta_O', 'eta_D');
fprintf('---------------------------------------------------------------------------------------------\n');

Results_B = [];
Cb = Ship.Cb_ref;
L = Ship.L_ref;
T = Ship.T_ref;
V_ship = Ship.Vs_knots_ref * 0.5144;

for i = 1:NumSamples
    B = B_range(i);

    Y = generate_hull_offsets(x_norm, Nz, Cb);
    S_wet = estimate_wetted_surface(L, B, T, Cb);
    Re = V_ship * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    K1 = 1 + calculate_form_factor(L, B, T, Cb);
    Rf = 0.5 * Rho_water * S_wet * V_ship^2 * Cf * K1;
    Rw = michell(Y, V_ship, L, B, T, Rho_water, N_theta);
    Rw = max(0, Rw);
    Rt = Rf + Rw;

    Cp = Cb / 0.98;
    [w, t_thrust, eta_R] = Holtrop_Propulsion_Calculation_1982(...
        L, B, T, Ship.D_prop, S_wet, Cb, Cp, Ship.LCB_percent, Ship.Cstern, Ship.AE_AO, Cf, K1);
    w = max(0.15, min(0.45, w));
    t_thrust = max(0.10, min(0.30, t_thrust));
    eta_R = max(0.95, min(1.05, eta_R));

    eta_H = (1 - t_thrust) / (1 - w);
    Thrust_req = Rt / (1 - t_thrust);
    Va = V_ship * (1 - w);

    Common_Def = Common_Def_base;
    Common_Def(2) = Va;
    XVA_def = ones(size(XR0)) * (1 - w) + 0.08*(1-XR0);
    XVT_def = zeros(size(XR0)) + 0.01 * (1-XR0);

    eta_O = call_PVL_for_efficiency(Common_Def, XR0, XCHD_def, XCD_def, ...
        XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
        Single_def1, Single_def2, Mean, Thick, Thrust_req);

    if isnan(eta_O)
        fprintf('%4d  %6.2f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  [PVL无效-跳过]\n', ...
            i, B, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H);
        continue;
    end

    eta_D = eta_H * eta_O * eta_R;
    Pd = Rt * V_ship / eta_D;

    valid_count_B = valid_count_B + 1;
    Results_B(valid_count_B).B = B;
    Results_B(valid_count_B).Rt = Rt;
    Results_B(valid_count_B).Rw = Rw;
    Results_B(valid_count_B).Rf = Rf;
    Results_B(valid_count_B).w = w;
    Results_B(valid_count_B).t = t_thrust;
    Results_B(valid_count_B).eta_R = eta_R;
    Results_B(valid_count_B).eta_H = eta_H;
    Results_B(valid_count_B).eta_O = eta_O;
    Results_B(valid_count_B).eta_D = eta_D;
    Results_B(valid_count_B).Pd = Pd;

    fprintf('%4d  %6.2f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        i, B, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H, eta_O, eta_D);

    close all hidden;
end

% ==============================================================================
% 6. 分析3: 船长L变化
% ==============================================================================
fprintf('\n===== 分析3: 船长L变化对推进效率的影响 =====\n');
fprintf('%-5s %-7s %-9s %-9s %-9s %-7s %-7s %-7s %-7s %-7s %-7s\n', ...
    'No.', 'L[m]', 'Rt[kN]', 'Rw[kN]', 'Rf[kN]', 'w', 't', 'eta_R', 'eta_H', 'eta_O', 'eta_D');
fprintf('---------------------------------------------------------------------------------------------\n');

Results_L = [];
Cb = Ship.Cb_ref;
B = Ship.B_ref;
T = Ship.T_ref;
V_ship = Ship.Vs_knots_ref * 0.5144;

for i = 1:NumSamples
    L = L_range(i);

    Y = generate_hull_offsets(x_norm, Nz, Cb);
    S_wet = estimate_wetted_surface(L, B, T, Cb);
    Re = V_ship * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    K1 = 1 + calculate_form_factor(L, B, T, Cb);
    Rf = 0.5 * Rho_water * S_wet * V_ship^2 * Cf * K1;
    Rw = michell(Y, V_ship, L, B, T, Rho_water, N_theta);
    Rw = max(0, Rw);
    Rt = Rf + Rw;

    Cp = Cb / 0.98;
    [w, t_thrust, eta_R] = Holtrop_Propulsion_Calculation_1982(...
        L, B, T, Ship.D_prop, S_wet, Cb, Cp, Ship.LCB_percent, Ship.Cstern, Ship.AE_AO, Cf, K1);
    w = max(0.15, min(0.45, w));
    t_thrust = max(0.10, min(0.30, t_thrust));
    eta_R = max(0.95, min(1.05, eta_R));

    eta_H = (1 - t_thrust) / (1 - w);
    Thrust_req = Rt / (1 - t_thrust);
    Va = V_ship * (1 - w);

    Common_Def = Common_Def_base;
    Common_Def(2) = Va;
    XVA_def = ones(size(XR0)) * (1 - w) + 0.08*(1-XR0);
    XVT_def = zeros(size(XR0)) + 0.01 * (1-XR0);

    eta_O = call_PVL_for_efficiency(Common_Def, XR0, XCHD_def, XCD_def, ...
        XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
        Single_def1, Single_def2, Mean, Thick, Thrust_req);

    if isnan(eta_O)
        fprintf('%4d  %6.1f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  [PVL无效-跳过]\n', ...
            i, L, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H);
        continue;
    end

    eta_D = eta_H * eta_O * eta_R;
    Pd = Rt * V_ship / eta_D;

    valid_count_L = valid_count_L + 1;
    Results_L(valid_count_L).L = L;
    Results_L(valid_count_L).Rt = Rt;
    Results_L(valid_count_L).Rw = Rw;
    Results_L(valid_count_L).Rf = Rf;
    Results_L(valid_count_L).w = w;
    Results_L(valid_count_L).t = t_thrust;
    Results_L(valid_count_L).eta_R = eta_R;
    Results_L(valid_count_L).eta_H = eta_H;
    Results_L(valid_count_L).eta_O = eta_O;
    Results_L(valid_count_L).eta_D = eta_D;
    Results_L(valid_count_L).Pd = Pd;

    fprintf('%4d  %6.1f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        i, L, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H, eta_O, eta_D);

    close all hidden;
end

% ==============================================================================
% 7. 分析4: 船速Vs变化
% ==============================================================================
fprintf('\n===== 分析4: 船速Vs变化对推进效率的影响 =====\n');
fprintf('%-5s %-8s %-9s %-9s %-9s %-7s %-7s %-7s %-7s %-7s %-7s\n', ...
    'No.', 'Vs[kn]', 'Rt[kN]', 'Rw[kN]', 'Rf[kN]', 'w', 't', 'eta_R', 'eta_H', 'eta_O', 'eta_D');
fprintf('---------------------------------------------------------------------------------------------\n');

Results_Vs = [];
Cb = Ship.Cb_ref;
L = Ship.L_ref;
B = Ship.B_ref;
T = Ship.T_ref;

for i = 1:NumSamples
    Vs_knots = Vs_range(i);
    V_ship = Vs_knots * 0.5144;

    Y = generate_hull_offsets(x_norm, Nz, Cb);
    S_wet = estimate_wetted_surface(L, B, T, Cb);
    Re = V_ship * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    K1 = 1 + calculate_form_factor(L, B, T, Cb);
    Rf = 0.5 * Rho_water * S_wet * V_ship^2 * Cf * K1;
    Rw = michell(Y, V_ship, L, B, T, Rho_water, N_theta);
    Rw = max(0, Rw);
    Rt = Rf + Rw;

    Cp = Cb / 0.98;
    [w, t_thrust, eta_R] = Holtrop_Propulsion_Calculation_1982(...
        L, B, T, Ship.D_prop, S_wet, Cb, Cp, Ship.LCB_percent, Ship.Cstern, Ship.AE_AO, Cf, K1);
    w = max(0.15, min(0.45, w));
    t_thrust = max(0.10, min(0.30, t_thrust));
    eta_R = max(0.95, min(1.05, eta_R));

    eta_H = (1 - t_thrust) / (1 - w);
    Thrust_req = Rt / (1 - t_thrust);
    Va = V_ship * (1 - w);

    Common_Def = Common_Def_base;
    Common_Def(2) = Va;
    XVA_def = ones(size(XR0)) * (1 - w) + 0.08*(1-XR0);
    XVT_def = zeros(size(XR0)) + 0.01 * (1-XR0);

    eta_O = call_PVL_for_efficiency(Common_Def, XR0, XCHD_def, XCD_def, ...
        XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
        Single_def1, Single_def2, Mean, Thick, Thrust_req);

    if isnan(eta_O)
        fprintf('%4d  %7.2f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  [PVL无效-跳过]\n', ...
            i, Vs_knots, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H);
        continue;
    end

    eta_D = eta_H * eta_O * eta_R;
    Pe = Rt * V_ship;
    Pd = Pe / eta_D;

    valid_count_Vs = valid_count_Vs + 1;
    Results_Vs(valid_count_Vs).Vs = Vs_knots;
    Results_Vs(valid_count_Vs).V_ship = V_ship;
    Results_Vs(valid_count_Vs).Rt = Rt;
    Results_Vs(valid_count_Vs).Rw = Rw;
    Results_Vs(valid_count_Vs).Rf = Rf;
    Results_Vs(valid_count_Vs).w = w;
    Results_Vs(valid_count_Vs).t = t_thrust;
    Results_Vs(valid_count_Vs).eta_R = eta_R;
    Results_Vs(valid_count_Vs).eta_H = eta_H;
    Results_Vs(valid_count_Vs).eta_O = eta_O;
    Results_Vs(valid_count_Vs).eta_D = eta_D;
    Results_Vs(valid_count_Vs).Pe = Pe;
    Results_Vs(valid_count_Vs).Pd = Pd;

    fprintf('%4d  %7.2f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        i, Vs_knots, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H, eta_O, eta_D);

    close all hidden;
end

% ==============================================================================
% 8. 分析5: 吃水T变化
% ==============================================================================
fprintf('\n===== 分析5: 吃水T变化对推进效率的影响 =====\n');
fprintf('%-5s %-7s %-9s %-9s %-9s %-7s %-7s %-7s %-7s %-7s %-7s\n', ...
    'No.', 'T[m]', 'Rt[kN]', 'Rw[kN]', 'Rf[kN]', 'w', 't', 'eta_R', 'eta_H', 'eta_O', 'eta_D');
fprintf('---------------------------------------------------------------------------------------------\n');

Results_T = [];
Cb = Ship.Cb_ref;
L = Ship.L_ref;
B = Ship.B_ref;
Vs_knots = Ship.Vs_knots_ref;
V_ship = Vs_knots * 0.5144;

for i = 1:NumSamples
    T = T_range(i);

    Y = generate_hull_offsets(x_norm, Nz, Cb);
    S_wet = estimate_wetted_surface(L, B, T, Cb);
    Re = V_ship * L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    K1 = 1 + calculate_form_factor(L, B, T, Cb);
    Rf = 0.5 * Rho_water * S_wet * V_ship^2 * Cf * K1;
    Rw = michell(Y, V_ship, L, B, T, Rho_water, N_theta);
    Rw = max(0, Rw);
    Rt = Rf + Rw;

    Cp = Cb / 0.98;
    [w, t_thrust, eta_R] = Holtrop_Propulsion_Calculation_1982(...
        L, B, T, Ship.D_prop, S_wet, Cb, Cp, Ship.LCB_percent, Ship.Cstern, Ship.AE_AO, Cf, K1);
    w = max(0.15, min(0.45, w));
    t_thrust = max(0.10, min(0.30, t_thrust));
    eta_R = max(0.95, min(1.05, eta_R));

    eta_H = (1 - t_thrust) / (1 - w);
    Thrust_req = Rt / (1 - t_thrust);
    Va = V_ship * (1 - w);

    Common_Def = Common_Def_base;
    Common_Def(2) = Va;
    XVA_def = ones(size(XR0)) * (1 - w) + 0.08*(1-XR0);
    XVT_def = zeros(size(XR0)) + 0.01 * (1-XR0);

    eta_O = call_PVL_for_efficiency(Common_Def, XR0, XCHD_def, XCD_def, ...
        XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
        Single_def1, Single_def2, Mean, Thick, Thrust_req);

    if isnan(eta_O)
        fprintf('%4d  %6.2f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  [PVL无效-跳过]\n', ...
            i, T, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H);
        continue;
    end

    eta_D = eta_H * eta_O * eta_R;
    Pd = Rt * V_ship / eta_D;

    valid_count_T = valid_count_T + 1;
    Results_T(valid_count_T).T = T;
    Results_T(valid_count_T).Rt = Rt;
    Results_T(valid_count_T).Rw = Rw;
    Results_T(valid_count_T).Rf = Rf;
    Results_T(valid_count_T).w = w;
    Results_T(valid_count_T).t = t_thrust;
    Results_T(valid_count_T).eta_R = eta_R;
    Results_T(valid_count_T).eta_H = eta_H;
    Results_T(valid_count_T).eta_O = eta_O;
    Results_T(valid_count_T).eta_D = eta_D;
    Results_T(valid_count_T).Pd = Pd;

    fprintf('%4d  %6.2f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        i, T, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H, eta_O, eta_D);

    close all hidden;
end

% ==============================================================================
% 9. 核心结论可视化: 阻力越小，推进效率不一定越高
% ==============================================================================
fprintf('\n>>> 生成核心结论证明图表...\n');

figure('Position', [50 50 1200 800], 'Color', 'w', 'Name', 'Core Conclusion: Resistance vs Efficiency');

% --- 子图1: Cb变化 - 阻力与效率对比 ---
if valid_count_Cb >= 2
    subplot(2, 3, 1);
    Cb_vec = [Results_Cb.Cb];
    Rt_Cb = [Results_Cb.Rt]/1000;
    eta_D_Cb = [Results_Cb.eta_D];

    yyaxis left;
    plot(Cb_vec, Rt_Cb, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    ylabel('总阻力 R_t [kN]', 'FontSize', 11);

    yyaxis right;
    plot(Cb_vec, eta_D_Cb, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    ylabel('推进效率 \eta_D', 'FontSize', 11);

    xlabel('方形系数 C_b', 'FontSize', 11);
    title('C_b变化: 阻力 vs 效率', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;

    % 标注最小阻力点和最高效率点
    [~, idx_min_Rt] = min(Rt_Cb);
    [~, idx_max_eta] = max(eta_D_Cb);
    hold on;
    yyaxis left;
    plot(Cb_vec(idx_min_Rt), Rt_Cb(idx_min_Rt), 'bp', 'MarkerSize', 15, 'MarkerFaceColor', 'c');
    yyaxis right;
    plot(Cb_vec(idx_max_eta), eta_D_Cb(idx_max_eta), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'm');
    legend('R_t', '\eta_D', '最小阻力', '最高效率', 'Location', 'best');
end

% --- 子图2: B变化 - 阻力与效率对比 ---
if valid_count_B >= 2
    subplot(2, 3, 2);
    B_vec = [Results_B.B];
    Rt_B = [Results_B.Rt]/1000;
    eta_D_B = [Results_B.eta_D];

    yyaxis left;
    plot(B_vec, Rt_B, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    ylabel('R_t [kN]', 'FontSize', 11);

    yyaxis right;
    plot(B_vec, eta_D_B, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    ylabel('\eta_D', 'FontSize', 11);

    xlabel('船宽 B [m]', 'FontSize', 11);
    title('B变化: 阻力 vs 效率', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;

    [~, idx_min_Rt] = min(Rt_B);
    [~, idx_max_eta] = max(eta_D_B);
    hold on;
    yyaxis left;
    plot(B_vec(idx_min_Rt), Rt_B(idx_min_Rt), 'bp', 'MarkerSize', 15, 'MarkerFaceColor', 'c');
    yyaxis right;
    plot(B_vec(idx_max_eta), eta_D_B(idx_max_eta), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'm');
end

% --- 子图3: L变化 - 阻力与效率对比 ---
if valid_count_L >= 2
    subplot(2, 3, 3);
    L_vec = [Results_L.L];
    Rt_L = [Results_L.Rt]/1000;
    eta_D_L = [Results_L.eta_D];

    yyaxis left;
    plot(L_vec, Rt_L, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    ylabel('R_t [kN]', 'FontSize', 11);

    yyaxis right;
    plot(L_vec, eta_D_L, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    ylabel('\eta_D', 'FontSize', 11);

    xlabel('船长 L [m]', 'FontSize', 11);
    title('L变化: 阻力 vs 效率', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;

    [~, idx_min_Rt] = min(Rt_L);
    [~, idx_max_eta] = max(eta_D_L);
    hold on;
    yyaxis left;
    plot(L_vec(idx_min_Rt), Rt_L(idx_min_Rt), 'bp', 'MarkerSize', 15, 'MarkerFaceColor', 'c');
    yyaxis right;
    plot(L_vec(idx_max_eta), eta_D_L(idx_max_eta), 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'm');
end

% --- 子图4: Vs变化 - 阻力与效率对比 ---
if valid_count_Vs >= 2
    subplot(2, 3, 4);
    Vs_vec = [Results_Vs.Vs];
    Rt_Vs = [Results_Vs.Rt]/1000;
    eta_D_Vs = [Results_Vs.eta_D];

    yyaxis left;
    plot(Vs_vec, Rt_Vs, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    ylabel('R_t [kN]', 'FontSize', 11);

    yyaxis right;
    plot(Vs_vec, eta_D_Vs, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    ylabel('\eta_D', 'FontSize', 11);

    xlabel('船速 V_s [knots]', 'FontSize', 11);
    title('V_s变化: 阻力 vs 效率', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
end

% --- 子图5: T变化 - 阻力与效率对比 ---
if valid_count_T >= 2
    subplot(2, 3, 5);
    T_vec = [Results_T.T];
    Rt_T = [Results_T.Rt]/1000;
    eta_D_T = [Results_T.eta_D];

    yyaxis left;
    plot(T_vec, Rt_T, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    ylabel('R_t [kN]', 'FontSize', 11);

    yyaxis right;
    plot(T_vec, eta_D_T, 'r-s', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    ylabel('\eta_D', 'FontSize', 11);

    xlabel('吃水 T [m]', 'FontSize', 11);
    title('T变化: 阻力 vs 效率', 'FontSize', 12, 'FontWeight', 'bold');
    grid on;
end

% --- 子图6: 物理解释文字 ---
subplot(2, 3, 6);
axis off;
text(0.1, 0.9, '核心结论证明:', 'FontSize', 14, 'FontWeight', 'bold');
text(0.1, 0.75, '阻力越小，推进效率不一定越高', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
text(0.1, 0.55, '物理解释:', 'FontSize', 11, 'FontWeight', 'bold');
text(0.1, 0.42, '1. 船型更丰满(Cb增大) \rightarrow 阻力增加', 'FontSize', 10);
text(0.1, 0.32, '2. 但艉部更饱满 \rightarrow 伴流分数w增加', 'FontSize', 10);
text(0.1, 0.22, '3. w增加 \rightarrow 船身效率\eta_H提高', 'FontSize', 10);
text(0.1, 0.12, '4. 当\eta_H增益 > 阻力损失 \rightarrow 总效率更高', 'FontSize', 10);
text(0.1, 0.02, '5. 因此最小阻力点 \neq 最高效率点', 'FontSize', 10, 'Color', 'b');

sgtitle({'船舶推进效率分析: 阻力与效率的非单调关系', ...
    '(Michell积分 + ITTC + PVL + Holtrop 1982)'}, 'FontSize', 14, 'FontWeight', 'bold');

saveas(gcf, 'Core_Conclusion_Resistance_vs_Efficiency.png');

% ==============================================================================
% 10. 结论分析输出
% ==============================================================================
fprintf('\n================================================================================\n');
fprintf('                           结 论 分 析\n');
fprintf('================================================================================\n\n');

% Cb分析
if valid_count_Cb >= 2
    [~, idx_min_Rt_Cb] = min([Results_Cb.Rt]);
    [~, idx_max_etaD_Cb] = max([Results_Cb.eta_D]);
    [~, idx_min_Pd_Cb] = min([Results_Cb.Pd]);
    fprintf('【分析1: 方形系数Cb变化】(有效工况: %d/%d)\n', valid_count_Cb, NumSamples);
    fprintf('  - 最小阻力点: Cb = %.4f, Rt = %.1f kN\n', Results_Cb(idx_min_Rt_Cb).Cb, Results_Cb(idx_min_Rt_Cb).Rt/1000);
    fprintf('  - 最高推进效率点: Cb = %.4f, eta_D = %.4f\n', Results_Cb(idx_max_etaD_Cb).Cb, Results_Cb(idx_max_etaD_Cb).eta_D);
    fprintf('  - 最小功率点: Cb = %.4f, Pd = %.2f MW\n', Results_Cb(idx_min_Pd_Cb).Cb, Results_Cb(idx_min_Pd_Cb).Pd/1e6);
    if idx_min_Rt_Cb ~= idx_max_etaD_Cb
        fprintf('  >>> 结论验证: 最小阻力点(Cb=%.4f) 与 最高效率点(Cb=%.4f) 不同!\n', ...
            Results_Cb(idx_min_Rt_Cb).Cb, Results_Cb(idx_max_etaD_Cb).Cb);
    end
    fprintf('\n');
end

% B分析
if valid_count_B >= 2
    [~, idx_min_Rt_B] = min([Results_B.Rt]);
    [~, idx_max_etaD_B] = max([Results_B.eta_D]);
    fprintf('【分析2: 船宽B变化】(有效工况: %d/%d)\n', valid_count_B, NumSamples);
    fprintf('  - 最小阻力点: B = %.2f m\n', Results_B(idx_min_Rt_B).B);
    fprintf('  - 最高推进效率点: B = %.2f m\n', Results_B(idx_max_etaD_B).B);
    if idx_min_Rt_B ~= idx_max_etaD_B
        fprintf('  >>> 结论验证: 最小阻力点(B=%.2f) 与 最高效率点(B=%.2f) 不同!\n', ...
            Results_B(idx_min_Rt_B).B, Results_B(idx_max_etaD_B).B);
    end
    fprintf('\n');
end

% L分析
if valid_count_L >= 2
    [~, idx_min_Rt_L] = min([Results_L.Rt]);
    [~, idx_max_etaD_L] = max([Results_L.eta_D]);
    fprintf('【分析3: 船长L变化】(有效工况: %d/%d)\n', valid_count_L, NumSamples);
    fprintf('  - 最小阻力点: L = %.1f m\n', Results_L(idx_min_Rt_L).L);
    fprintf('  - 最高推进效率点: L = %.1f m\n', Results_L(idx_max_etaD_L).L);
    if idx_min_Rt_L ~= idx_max_etaD_L
        fprintf('  >>> 结论验证: 最小阻力点(L=%.1f) 与 最高效率点(L=%.1f) 不同!\n', ...
            Results_L(idx_min_Rt_L).L, Results_L(idx_max_etaD_L).L);
    end
    fprintf('\n');
end

% Vs分析
if valid_count_Vs >= 2
    [~, idx_max_etaD_Vs] = max([Results_Vs.eta_D]);
    fprintf('【分析4: 船速Vs变化】(有效工况: %d/%d)\n', valid_count_Vs, NumSamples);
    fprintf('  - 最高推进效率点: Vs = %.1f knots, eta_D = %.4f\n', ...
        Results_Vs(idx_max_etaD_Vs).Vs, Results_Vs(idx_max_etaD_Vs).eta_D);
    fprintf('\n');
end

% T分析
if valid_count_T >= 2
    [~, idx_min_Rt_T] = min([Results_T.Rt]);
    [~, idx_max_etaD_T] = max([Results_T.eta_D]);
    fprintf('【分析5: 吃水T变化】(有效工况: %d/%d)\n', valid_count_T, NumSamples);
    fprintf('  - 最小阻力点: T = %.2f m\n', Results_T(idx_min_Rt_T).T);
    fprintf('  - 最高推进效率点: T = %.2f m\n', Results_T(idx_max_etaD_T).T);
    if idx_min_Rt_T ~= idx_max_etaD_T
        fprintf('  >>> 结论验证: 最小阻力点(T=%.2f) 与 最高效率点(T=%.2f) 不同!\n', ...
            Results_T(idx_min_Rt_T).T, Results_T(idx_max_etaD_T).T);
    end
    fprintf('\n');
end

fprintf('================================================================================\n');
fprintf('【总结】\n');
fprintf('  本分析采用以下方法计算船舶完整推进效率:\n');
fprintf('  1. Michell积分 (michell.m) - 兴波阻力 Rw\n');
fprintf('  2. ITTC 1957 公式 - 摩擦阻力 Rf\n');
fprintf('  3. Holtrop 1982 公式 - 伴流分数w、推力减额分数t、相对旋转效率eta_R\n');
fprintf('  4. PVL升力线理论 (PVL.m) - 螺旋桨敞水效率 eta_O [直接调用]\n');
fprintf('\n');
fprintf('  完整推进效率公式: eta_D = eta_H * eta_O * eta_R\n');
fprintf('  其中船身效率: eta_H = (1-t)/(1-w)\n');
fprintf('\n');
fprintf('  核心结论: 阻力越小，推进效率不一定越高!\n');
fprintf('================================================================================\n');

% 保存结果
save('Complete_Propulsion_Results.mat', 'Results_Cb', 'Results_B', 'Results_L', 'Results_Vs', 'Results_T', 'Ship');
fprintf('\n结果已保存至 Complete_Propulsion_Results.mat\n');
fprintf('图表已保存至 Core_Conclusion_Resistance_vs_Efficiency.png\n');

% ==============================================================================
% 辅助函数
% ==============================================================================

function Y = generate_hull_offsets(x_norm, Nz, Cb)
    % 基于方形系数生成船体偏移矩阵(用于Michell积分)
    Nx = length(x_norm);

    h_val = 0.70 + 2.5 * (Cb - 0.63);
    a1 = 0.12;
    a2 = 0.10;

    f0 = zeros(size(x_norm));
    for ii = 1:length(x_norm)
        if x_norm(ii) < 0.5
            f0(ii) = h_val - a1/2 * (cos(2*pi*x_norm(ii)) - 1);
        else
            f0(ii) = h_val + a2/2 * (cos(2*pi*x_norm(ii)) + 1) + a1;
        end
    end
    f0 = f0 / max(f0) * 0.5;

    z_norm = linspace(0, 1, Nz);
    Y = zeros(Nx, Nz);
    for iz = 1:Nz
        z_factor = sin(pi/2 * z_norm(iz))^0.5;
        Y(:, iz) = f0(:) * z_factor;
    end
end

function S = estimate_wetted_surface(L, B, T, Cb)
    Cm = 0.98;
    Cwp = 0.75 + 0.125 * Cb;
    Abt = 0;
    S = L * (2*T + B) * sqrt(Cm) * (0.453 + 0.4425*Cb - 0.2862*Cm - 0.003467*B/T + 0.3696*Cwp) + 2.38*Abt/Cb;
    S = max(S, 1.01 * L * (2*T + B*Cb));
end

function k1 = calculate_form_factor(L, B, T, Cb)
    LR = L * (1 - Cb + 0.06*Cb/(4*Cb-1));
    c14 = 1;
    Cp = Cb / 0.98;
    k1 = 0.93 + 0.4871 * c14 * (B/LR)^1.0681 * (T/L)^0.4611 * (L/LR)^0.1216 * ...
         (L^3/(L*B*T*Cb))^0.3649 * (1-Cp)^(-0.6042);
    k1 = max(0.05, min(0.5, k1));
end

function eta_O = call_PVL_for_efficiency(Common_Def, XR0, XCHD_def, XCD_def, ...
    XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
    Single_def1, Single_def2, Mean, Thick, Thrust_req)
    % 直接调用用户的PVL.m计算螺旋桨敞水效率
    % PVL返回: [Sigma, skew, rake, EFFY, eta_final, XVA]
    % 第5个返回值 eta_final 才是最终效率

    try
        [~, ~, ~, ~, eta_final, ~] = PVL(Common_Def, XR0, XCHD_def, XCD_def, ...
            XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
            Single_def1, Single_def2, Mean, Thick, Thrust_req);

        if isnan(eta_final) || eta_final <= 0 || eta_final > 1
            eta_O = NaN;
        else
            eta_O = eta_final;
        end
    catch
        eta_O = NaN;
    end
end
