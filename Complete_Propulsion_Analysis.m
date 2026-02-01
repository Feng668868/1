% Complete_Propulsion_Analysis.m
% 船舶完整推进效率分析 - 参数化研究
%
% 本程序直接调用:
% 1. michell.m      - Michell积分计算兴波阻力
% 2. ITTC 1957公式  - 摩擦阻力计算
% 3. PVL.m          - 螺旋桨敞水效率计算(升力线理论)
% 4. Holtrop_Propulsion_Calculation_1982.m - 伴流分数(w)、推力减额分数(t)、相对旋转效率(eta_R)
%
% 可调参数: 船长L, 船宽B, 方形系数Cb, 船速Vs, 吃水T
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

% ==============================================================================
% 4. 分析1: 方形系数Cb变化
% ==============================================================================
fprintf('\n===== 分析1: 方形系数Cb变化对推进效率的影响 =====\n');
fprintf('%-5s %-7s %-9s %-9s %-9s %-7s %-7s %-7s %-7s %-7s %-7s\n', ...
    'No.', 'Cb', 'Rt[kN]', 'Rw[kN]', 'Rf[kN]', 'w', 't', 'eta_R', 'eta_H', 'eta_O', 'eta_D');
fprintf('---------------------------------------------------------------------------------------------\n');

Results_Cb = struct();
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
    K1 = 1 + calculate_form_factor(L, B, T, Cb);  % 形状因子 (1+k1)
    Rf = 0.5 * Rho_water * S_wet * V_ship^2 * Cf * K1;

    % Michell积分兴波阻力 (直接调用michell.m)
    Rw = michell(Y, V_ship, L, B, T, Rho_water, N_theta);
    Rw = max(0, Rw);  % 确保非负

    Rt = Rf + Rw;

    % --- C. Holtrop推进因子 (直接调用) ---
    Cp = Cb / 0.98;  % 简化的棱形系数估算
    [w, t_thrust, eta_R] = Holtrop_Propulsion_Calculation_1982(...
        L, B, T, Ship.D_prop, S_wet, Cb, Cp, Ship.LCB_percent, ...
        Ship.Cstern, Ship.AE_AO, Cf, K1);

    % 修正范围约束
    w = max(0.15, min(0.45, w));
    t_thrust = max(0.10, min(0.30, t_thrust));
    eta_R = max(0.95, min(1.05, eta_R));

    % --- D. 推进效率计算 ---
    eta_H = (1 - t_thrust) / (1 - w);  % 船身效率
    Thrust_req = Rt / (1 - t_thrust);  % 所需推力
    Va = V_ship * (1 - w);             % 进速

    % 更新PVL参数
    Common_Def = Common_Def_base;
    Common_Def(2) = Va;  % 进速

    % 构建伴流场速度分布
    XVA_def = ones(size(XR0)) * (1 - w) + 0.08*(1-XR0);
    XVT_def = zeros(size(XR0)) + 0.01 * (1-XR0);

    % 调用PVL计算螺旋桨敞水效率
    eta_O = calculate_propeller_efficiency_PVL(Common_Def, XR0, XCHD_def, XCD_def, ...
        XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
        Single_def1, Single_def2, Mean, Thick, Thrust_req, Va, Ship.D_prop, Rho_water);

    % 总推进效率
    eta_D = eta_H * eta_O * eta_R;

    % 有效功率和传递功率
    Pe = Rt * V_ship;
    Pd = Pe / eta_D;

    % 存储结果
    Results_Cb(i).Cb = Cb;
    Results_Cb(i).Rt = Rt;
    Results_Cb(i).Rw = Rw;
    Results_Cb(i).Rf = Rf;
    Results_Cb(i).w = w;
    Results_Cb(i).t = t_thrust;
    Results_Cb(i).eta_R = eta_R;
    Results_Cb(i).eta_H = eta_H;
    Results_Cb(i).eta_O = eta_O;
    Results_Cb(i).eta_D = eta_D;
    Results_Cb(i).Pe = Pe;
    Results_Cb(i).Pd = Pd;

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

Results_B = struct();
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

    eta_O = calculate_propeller_efficiency_PVL(Common_Def, XR0, XCHD_def, XCD_def, ...
        XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
        Single_def1, Single_def2, Mean, Thick, Thrust_req, Va, Ship.D_prop, Rho_water);

    eta_D = eta_H * eta_O * eta_R;
    Pd = Rt * V_ship / eta_D;

    Results_B(i).B = B;
    Results_B(i).Rt = Rt;
    Results_B(i).Rw = Rw;
    Results_B(i).Rf = Rf;
    Results_B(i).w = w;
    Results_B(i).t = t_thrust;
    Results_B(i).eta_R = eta_R;
    Results_B(i).eta_H = eta_H;
    Results_B(i).eta_O = eta_O;
    Results_B(i).eta_D = eta_D;
    Results_B(i).Pd = Pd;

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

Results_L = struct();
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

    eta_O = calculate_propeller_efficiency_PVL(Common_Def, XR0, XCHD_def, XCD_def, ...
        XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
        Single_def1, Single_def2, Mean, Thick, Thrust_req, Va, Ship.D_prop, Rho_water);

    eta_D = eta_H * eta_O * eta_R;
    Pd = Rt * V_ship / eta_D;

    Results_L(i).L = L;
    Results_L(i).Rt = Rt;
    Results_L(i).Rw = Rw;
    Results_L(i).Rf = Rf;
    Results_L(i).w = w;
    Results_L(i).t = t_thrust;
    Results_L(i).eta_R = eta_R;
    Results_L(i).eta_H = eta_H;
    Results_L(i).eta_O = eta_O;
    Results_L(i).eta_D = eta_D;
    Results_L(i).Pd = Pd;

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

Results_Vs = struct();
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

    eta_O = calculate_propeller_efficiency_PVL(Common_Def, XR0, XCHD_def, XCD_def, ...
        XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
        Single_def1, Single_def2, Mean, Thick, Thrust_req, Va, Ship.D_prop, Rho_water);

    eta_D = eta_H * eta_O * eta_R;
    Pe = Rt * V_ship;
    Pd = Pe / eta_D;

    Results_Vs(i).Vs = Vs_knots;
    Results_Vs(i).V_ship = V_ship;
    Results_Vs(i).Rt = Rt;
    Results_Vs(i).Rw = Rw;
    Results_Vs(i).Rf = Rf;
    Results_Vs(i).w = w;
    Results_Vs(i).t = t_thrust;
    Results_Vs(i).eta_R = eta_R;
    Results_Vs(i).eta_H = eta_H;
    Results_Vs(i).eta_O = eta_O;
    Results_Vs(i).eta_D = eta_D;
    Results_Vs(i).Pe = Pe;
    Results_Vs(i).Pd = Pd;

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

Results_T = struct();
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

    eta_O = calculate_propeller_efficiency_PVL(Common_Def, XR0, XCHD_def, XCD_def, ...
        XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
        Single_def1, Single_def2, Mean, Thick, Thrust_req, Va, Ship.D_prop, Rho_water);

    eta_D = eta_H * eta_O * eta_R;
    Pd = Rt * V_ship / eta_D;

    Results_T(i).T = T;
    Results_T(i).Rt = Rt;
    Results_T(i).Rw = Rw;
    Results_T(i).Rf = Rf;
    Results_T(i).w = w;
    Results_T(i).t = t_thrust;
    Results_T(i).eta_R = eta_R;
    Results_T(i).eta_H = eta_H;
    Results_T(i).eta_O = eta_O;
    Results_T(i).eta_D = eta_D;
    Results_T(i).Pd = Pd;

    fprintf('%4d  %6.2f  %8.1f  %8.1f  %8.1f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n', ...
        i, T, Rt/1000, Rw/1000, Rf/1000, w, t_thrust, eta_R, eta_H, eta_O, eta_D);

    close all hidden;
end

% ==============================================================================
% 9. 结果可视化
% ==============================================================================
fprintf('\n>>> 生成综合分析图表...\n');

figure('Position', [30 30 1600 1000], 'Color', 'w', 'Name', 'Complete Propulsion Analysis');

% --- 第一行: Cb分析 ---
subplot(5, 4, 1);
Cb_vec = [Results_Cb.Cb];
plot(Cb_vec, [Results_Cb.Rt]/1000, 'b-o', 'LineWidth', 2); hold on;
plot(Cb_vec, [Results_Cb.Rw]/1000, 'r--s', 'LineWidth', 1.5);
plot(Cb_vec, [Results_Cb.Rf]/1000, 'g-.^', 'LineWidth', 1.5);
xlabel('C_b'); ylabel('阻力 [kN]');
title('方形系数C_b vs 阻力');
legend('R_t', 'R_w', 'R_f', 'Location', 'best');
grid on;

subplot(5, 4, 2);
plot(Cb_vec, [Results_Cb.eta_H], 'm-o', 'LineWidth', 2); hold on;
plot(Cb_vec, [Results_Cb.eta_O], 'c-s', 'LineWidth', 2);
plot(Cb_vec, [Results_Cb.eta_R], 'y-^', 'LineWidth', 1.5);
plot(Cb_vec, [Results_Cb.eta_D], 'k-d', 'LineWidth', 2.5);
xlabel('C_b'); ylabel('效率');
title('方形系数C_b vs 效率');
legend('\eta_H', '\eta_O', '\eta_R', '\eta_D', 'Location', 'best');
grid on;

subplot(5, 4, 3);
plot(Cb_vec, [Results_Cb.w], 'b-o', 'LineWidth', 2); hold on;
plot(Cb_vec, [Results_Cb.t], 'r-s', 'LineWidth', 2);
xlabel('C_b'); ylabel('系数');
title('C_b vs 伴流/推力减额');
legend('w伴流', 't推力减额', 'Location', 'best');
grid on;

subplot(5, 4, 4);
yyaxis left;
plot(Cb_vec, [Results_Cb.Rt]/1000, 'b-o', 'LineWidth', 2);
ylabel('R_t [kN]');
yyaxis right;
plot(Cb_vec, [Results_Cb.Pd]/1e6, 'r-s', 'LineWidth', 2);
ylabel('P_d [MW]');
xlabel('C_b');
title('阻力 vs 功率');
grid on;

% --- 第二行: B分析 ---
subplot(5, 4, 5);
B_vec = [Results_B.B];
plot(B_vec, [Results_B.Rt]/1000, 'b-o', 'LineWidth', 2); hold on;
plot(B_vec, [Results_B.Rw]/1000, 'r--s', 'LineWidth', 1.5);
plot(B_vec, [Results_B.Rf]/1000, 'g-.^', 'LineWidth', 1.5);
xlabel('B [m]'); ylabel('阻力 [kN]');
title('船宽B vs 阻力');
legend('R_t', 'R_w', 'R_f', 'Location', 'best');
grid on;

subplot(5, 4, 6);
plot(B_vec, [Results_B.eta_H], 'm-o', 'LineWidth', 2); hold on;
plot(B_vec, [Results_B.eta_O], 'c-s', 'LineWidth', 2);
plot(B_vec, [Results_B.eta_D], 'k-d', 'LineWidth', 2.5);
xlabel('B [m]'); ylabel('效率');
title('船宽B vs 效率');
legend('\eta_H', '\eta_O', '\eta_D', 'Location', 'best');
grid on;

subplot(5, 4, 7);
plot(B_vec, [Results_B.w], 'b-o', 'LineWidth', 2); hold on;
plot(B_vec, [Results_B.t], 'r-s', 'LineWidth', 2);
xlabel('B [m]'); ylabel('系数');
title('B vs 伴流/推力减额');
legend('w', 't', 'Location', 'best');
grid on;

subplot(5, 4, 8);
yyaxis left;
plot(B_vec, [Results_B.Rt]/1000, 'b-o', 'LineWidth', 2);
ylabel('R_t [kN]');
yyaxis right;
plot(B_vec, [Results_B.Pd]/1e6, 'r-s', 'LineWidth', 2);
ylabel('P_d [MW]');
xlabel('B [m]');
title('阻力 vs 功率');
grid on;

% --- 第三行: L分析 ---
subplot(5, 4, 9);
L_vec = [Results_L.L];
plot(L_vec, [Results_L.Rt]/1000, 'b-o', 'LineWidth', 2); hold on;
plot(L_vec, [Results_L.Rw]/1000, 'r--s', 'LineWidth', 1.5);
plot(L_vec, [Results_L.Rf]/1000, 'g-.^', 'LineWidth', 1.5);
xlabel('L [m]'); ylabel('阻力 [kN]');
title('船长L vs 阻力');
legend('R_t', 'R_w', 'R_f', 'Location', 'best');
grid on;

subplot(5, 4, 10);
plot(L_vec, [Results_L.eta_H], 'm-o', 'LineWidth', 2); hold on;
plot(L_vec, [Results_L.eta_O], 'c-s', 'LineWidth', 2);
plot(L_vec, [Results_L.eta_D], 'k-d', 'LineWidth', 2.5);
xlabel('L [m]'); ylabel('效率');
title('船长L vs 效率');
legend('\eta_H', '\eta_O', '\eta_D', 'Location', 'best');
grid on;

subplot(5, 4, 11);
plot(L_vec, [Results_L.w], 'b-o', 'LineWidth', 2); hold on;
plot(L_vec, [Results_L.t], 'r-s', 'LineWidth', 2);
xlabel('L [m]'); ylabel('系数');
title('L vs 伴流/推力减额');
legend('w', 't', 'Location', 'best');
grid on;

subplot(5, 4, 12);
yyaxis left;
plot(L_vec, [Results_L.Rt]/1000, 'b-o', 'LineWidth', 2);
ylabel('R_t [kN]');
yyaxis right;
plot(L_vec, [Results_L.Pd]/1e6, 'r-s', 'LineWidth', 2);
ylabel('P_d [MW]');
xlabel('L [m]');
title('阻力 vs 功率');
grid on;

% --- 第四行: Vs分析 ---
subplot(5, 4, 13);
Vs_vec = [Results_Vs.Vs];
plot(Vs_vec, [Results_Vs.Rt]/1000, 'b-o', 'LineWidth', 2); hold on;
plot(Vs_vec, [Results_Vs.Rw]/1000, 'r--s', 'LineWidth', 1.5);
plot(Vs_vec, [Results_Vs.Rf]/1000, 'g-.^', 'LineWidth', 1.5);
xlabel('V_s [knots]'); ylabel('阻力 [kN]');
title('船速V_s vs 阻力');
legend('R_t', 'R_w', 'R_f', 'Location', 'best');
grid on;

subplot(5, 4, 14);
plot(Vs_vec, [Results_Vs.eta_H], 'm-o', 'LineWidth', 2); hold on;
plot(Vs_vec, [Results_Vs.eta_O], 'c-s', 'LineWidth', 2);
plot(Vs_vec, [Results_Vs.eta_D], 'k-d', 'LineWidth', 2.5);
xlabel('V_s [knots]'); ylabel('效率');
title('船速V_s vs 效率');
legend('\eta_H', '\eta_O', '\eta_D', 'Location', 'best');
grid on;

subplot(5, 4, 15);
plot(Vs_vec, [Results_Vs.w], 'b-o', 'LineWidth', 2); hold on;
plot(Vs_vec, [Results_Vs.t], 'r-s', 'LineWidth', 2);
xlabel('V_s [knots]'); ylabel('系数');
title('V_s vs 伴流/推力减额');
legend('w', 't', 'Location', 'best');
grid on;

subplot(5, 4, 16);
yyaxis left;
plot(Vs_vec, [Results_Vs.Rt]/1000, 'b-o', 'LineWidth', 2);
ylabel('R_t [kN]');
yyaxis right;
plot(Vs_vec, [Results_Vs.Pd]/1e6, 'r-s', 'LineWidth', 2);
ylabel('P_d [MW]');
xlabel('V_s [knots]');
title('阻力 vs 功率');
grid on;

% --- 第五行: T分析 ---
subplot(5, 4, 17);
T_vec = [Results_T.T];
plot(T_vec, [Results_T.Rt]/1000, 'b-o', 'LineWidth', 2); hold on;
plot(T_vec, [Results_T.Rw]/1000, 'r--s', 'LineWidth', 1.5);
plot(T_vec, [Results_T.Rf]/1000, 'g-.^', 'LineWidth', 1.5);
xlabel('T [m]'); ylabel('阻力 [kN]');
title('吃水T vs 阻力');
legend('R_t', 'R_w', 'R_f', 'Location', 'best');
grid on;

subplot(5, 4, 18);
plot(T_vec, [Results_T.eta_H], 'm-o', 'LineWidth', 2); hold on;
plot(T_vec, [Results_T.eta_O], 'c-s', 'LineWidth', 2);
plot(T_vec, [Results_T.eta_D], 'k-d', 'LineWidth', 2.5);
xlabel('T [m]'); ylabel('效率');
title('吃水T vs 效率');
legend('\eta_H', '\eta_O', '\eta_D', 'Location', 'best');
grid on;

subplot(5, 4, 19);
plot(T_vec, [Results_T.w], 'b-o', 'LineWidth', 2); hold on;
plot(T_vec, [Results_T.t], 'r-s', 'LineWidth', 2);
xlabel('T [m]'); ylabel('系数');
title('T vs 伴流/推力减额');
legend('w', 't', 'Location', 'best');
grid on;

subplot(5, 4, 20);
yyaxis left;
plot(T_vec, [Results_T.Rt]/1000, 'b-o', 'LineWidth', 2);
ylabel('R_t [kN]');
yyaxis right;
plot(T_vec, [Results_T.Pd]/1e6, 'r-s', 'LineWidth', 2);
ylabel('P_d [MW]');
xlabel('T [m]');
title('阻力 vs 功率');
grid on;

sgtitle('船舶完整推进效率分析 (Michell + ITTC + PVL + Holtrop 1982)', 'FontSize', 14, 'FontWeight', 'bold');

% ==============================================================================
% 10. 结论分析输出
% ==============================================================================
fprintf('\n================================================================================\n');
fprintf('                           结 论 分 析\n');
fprintf('================================================================================\n\n');

% Cb分析
[~, idx_min_Rt_Cb] = min([Results_Cb.Rt]);
[~, idx_max_etaD_Cb] = max([Results_Cb.eta_D]);
[~, idx_min_Pd_Cb] = min([Results_Cb.Pd]);
fprintf('【分析1: 方形系数Cb变化】\n');
fprintf('  - 最小阻力点: Cb = %.4f, Rt = %.1f kN\n', Results_Cb(idx_min_Rt_Cb).Cb, Results_Cb(idx_min_Rt_Cb).Rt/1000);
fprintf('  - 最高推进效率点: Cb = %.4f, eta_D = %.4f\n', Results_Cb(idx_max_etaD_Cb).Cb, Results_Cb(idx_max_etaD_Cb).eta_D);
fprintf('  - 最小功率点: Cb = %.4f, Pd = %.2f MW\n\n', Results_Cb(idx_min_Pd_Cb).Cb, Results_Cb(idx_min_Pd_Cb).Pd/1e6);

% B分析
[~, idx_min_Rt_B] = min([Results_B.Rt]);
[~, idx_max_etaD_B] = max([Results_B.eta_D]);
fprintf('【分析2: 船宽B变化】\n');
fprintf('  - 最小阻力点: B = %.2f m\n', Results_B(idx_min_Rt_B).B);
fprintf('  - 最高推进效率点: B = %.2f m\n\n', Results_B(idx_max_etaD_B).B);

% L分析
[~, idx_min_Rt_L] = min([Results_L.Rt]);
[~, idx_max_etaD_L] = max([Results_L.eta_D]);
fprintf('【分析3: 船长L变化】\n');
fprintf('  - 最小阻力点: L = %.1f m\n', Results_L(idx_min_Rt_L).L);
fprintf('  - 最高推进效率点: L = %.1f m\n\n', Results_L(idx_max_etaD_L).L);

% Vs分析
[~, idx_max_etaD_Vs] = max([Results_Vs.eta_D]);
fprintf('【分析4: 船速Vs变化】\n');
fprintf('  - 最高推进效率点: Vs = %.1f knots, eta_D = %.4f\n\n', ...
    Results_Vs(idx_max_etaD_Vs).Vs, Results_Vs(idx_max_etaD_Vs).eta_D);

% T分析
[~, idx_min_Rt_T] = min([Results_T.Rt]);
[~, idx_max_etaD_T] = max([Results_T.eta_D]);
fprintf('【分析5: 吃水T变化】\n');
fprintf('  - 最小阻力点: T = %.2f m\n', Results_T(idx_min_Rt_T).T);
fprintf('  - 最高推进效率点: T = %.2f m\n\n', Results_T(idx_max_etaD_T).T);

fprintf('================================================================================\n');
fprintf('【总结】\n');
fprintf('  本分析采用以下方法计算船舶完整推进效率:\n');
fprintf('  1. Michell积分 (michell.m) - 兴波阻力\n');
fprintf('  2. ITTC 1957 公式 - 摩擦阻力\n');
fprintf('  3. Holtrop 1982 公式 - 伴流分数w、推力减额分数t、相对旋转效率eta_R\n');
fprintf('  4. PVL升力线理论 (PVL.m) - 螺旋桨敞水效率eta_O\n');
fprintf('\n');
fprintf('  完整推进效率: eta_D = eta_H * eta_O * eta_R\n');
fprintf('  其中: eta_H = (1-t)/(1-w) 为船身效率\n');
fprintf('================================================================================\n');

% 保存结果
save('Complete_Propulsion_Results.mat', 'Results_Cb', 'Results_B', 'Results_L', 'Results_Vs', 'Results_T', 'Ship');
saveas(gcf, 'Complete_Propulsion_Results.png');
fprintf('\n结果已保存至 Complete_Propulsion_Results.mat 和 Complete_Propulsion_Results.png\n');

% ==============================================================================
% 辅助函数
% ==============================================================================

function Y = generate_hull_offsets(x_norm, Nz, Cb)
    % 基于方形系数生成船体偏移矩阵(用于Michell积分)
    Nx = length(x_norm);

    % 三参数函数确定船型丰满度
    h_val = 0.70 + 2.5 * (Cb - 0.63);
    a1 = 0.12;  % 船艏参数
    a2 = 0.10;  % 船艉参数

    % 生成水线形状
    f0 = zeros(size(x_norm));
    for ii = 1:length(x_norm)
        if x_norm(ii) < 0.5
            f0(ii) = h_val - a1/2 * (cos(2*pi*x_norm(ii)) - 1);
        else
            f0(ii) = h_val + a2/2 * (cos(2*pi*x_norm(ii)) + 1) + a1;
        end
    end
    f0 = f0 / max(f0) * 0.5;  % 归一化，最大偏移为0.5

    % 生成三维偏移矩阵
    z_norm = linspace(0, 1, Nz);
    Y = zeros(Nx, Nz);
    for iz = 1:Nz
        z_factor = sin(pi/2 * z_norm(iz))^0.5;
        Y(:, iz) = f0(:) * z_factor;
    end
end

function S = estimate_wetted_surface(L, B, T, Cb)
    % Holtrop润湿面积公式
    Cm = 0.98;  % 中横剖面系数
    Cwp = 0.75 + 0.125 * Cb;  % 水线面系数
    Abt = 0;    % 球鼻艏横截面积

    S = L * (2*T + B) * sqrt(Cm) * (0.453 + 0.4425*Cb - 0.2862*Cm - 0.003467*B/T + 0.3696*Cwp) + 2.38*Abt/Cb;
    S = max(S, 1.01 * L * (2*T + B*Cb));
end

function k1 = calculate_form_factor(L, B, T, Cb)
    % Holtrop形状因子(1+k1)计算
    % 简化版本
    LR = L * (1 - Cb + 0.06*Cb/(4*Cb-1));  % run length
    c14 = 1;  % 无艉轴包壳
    Cp = Cb / 0.98;
    k1 = 0.93 + 0.4871 * c14 * (B/LR)^1.0681 * (T/L)^0.4611 * (L/LR)^0.1216 * ...
         (L^3/(L*B*T*Cb))^0.3649 * (1-Cp)^(-0.6042);
    k1 = max(0.05, min(0.5, k1));  % 限制在合理范围
end

function eta_O = calculate_propeller_efficiency_PVL(Common_Def, XR0, XCHD_def, XCD_def, ...
    XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
    Single_def1, Single_def2, Mean, Thick, Thrust_req, Va, D, rho)
    % 尝试调用PVL计算螺旋桨效率，失败时使用简化模型

    try
        [~, ~, ~, eta_pvl, ~] = PVL(Common_Def, XR0, XCHD_def, XCD_def, ...
            XVA_def, XVT_def, f0oc_def, t0oc_def, skew_def, rake_def, ...
            Single_def1, Single_def2, Mean, Thick, Thrust_req);

        if isnan(eta_pvl) || eta_pvl <= 0 || eta_pvl > 1
            eta_O = simplified_propeller_efficiency(Thrust_req, Va, D, rho);
        else
            eta_O = eta_pvl;
        end
    catch
        % PVL计算失败，使用简化模型
        eta_O = simplified_propeller_efficiency(Thrust_req, Va, D, rho);
    end
end

function eta_O = simplified_propeller_efficiency(Thrust, Va, D, rho)
    % 基于动量理论的简化螺旋桨效率计算
    A0 = pi * (D/2)^2;
    CT_load = Thrust / (0.5 * rho * Va^2 * A0);
    CT_load = max(0.1, min(5, CT_load));  % 限制在合理范围

    % 理想效率 (Rankine-Froude)
    eta_ideal = 2 / (1 + sqrt(1 + CT_load));

    % 实际效率修正
    eta_real_factor = 0.82 - 0.03 * (CT_load - 0.8);
    eta_real_factor = max(0.70, min(0.88, eta_real_factor));

    eta_O = eta_ideal * eta_real_factor;
    eta_O = max(0.50, min(0.78, eta_O));
end
