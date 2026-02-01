% Main_KCS_Tradeoff_Michell.m
% KCS 船桨集成优化 - Michell 积分集成版
% 
% 特点：
% 1. 兴波阻力由 michell.m 函数计算（基于 Filon 积分算法）
% 2. 自动构造 3D 偏移矩阵 Y(Nx, Nz)
% 3. 保持螺距自适应模式，验证 Cb 增加带来的交互效率增益
clear; close all; clc;

% ==============================================================================
% 1. 全局仿真设置
% ==============================================================================
NumSamples = 20;          
V_ship_mps = 12.35;       % 24 knots
Rho_water = 1025;         
Nu = 1.188e-6;
g = 9.80665;

% KCS 基准参数
KCS.L = 230.0; 
KCS.B = 32.2; 
KCS.T = 10.8; 
KCS.D_prop = 7.9; 
KCS.Cb_ref = 0.6505;

% Michell 积分网格参数
Nx = 101; % 站点数 (必须为奇数)
Nz = 40;  % 水线数
N_theta = 80; % 传播角采样数

% 设计空间
Cb_targets = linspace(0.630, 0.680, NumSamples);

% 结果存储
Results = struct('Cb',[],'Rt',[],'Rw',[],'Rf',[],'w',[],'t',[],'eta_H',[],'eta_O',[],'eta_D',[],'Pd',[]);
Results = repmat(Results, NumSamples, 1);

fprintf('===============================================================================\n');
fprintf('   KCS 船-桨优化：Michell 积分集成验证版\n');
fprintf('===============================================================================\n');
fprintf('%-6s %-8s %-8s %-8s %-8s %-8s %-8s %-8s\n', 'Iter', 'Cb', 'Rt(kN)', 'Rw(kN)', 'w', 'eta_H', 'eta_D', 'Pd(MW)');
fprintf('-------------------------------------------------------------------------------\n');

% ==============================================================================
% 2. 循环扫描
% ==============================================================================
x_norm = linspace(0, 1, Nx);

for i = 1:NumSamples
    Cb = Cb_targets(i);
    dCb = Cb - KCS.Cb_ref;
    
    % --- A. 几何生成 ---
    % 模拟实际设计：Cb 增加时，船宽 B 轻微调整以平衡排水量
    B = KCS.B * (1 + 0.2 * dCb); 
    Vol = KCS.L * B * KCS.T * Cb;
    
    % 生成水线面形状 f0 (非量纲)
    % h_val 控制肥瘦，a1/a2 控制艏艉饱满度
    h_val = 0.70 + 2.5 * (Cb - 0.63); 
    f0 = three_param_shape(x_norm, h_val, 0.12, 0.10);
    f0 = f0 / max(f0) * 0.5; % 归一化，最大半宽占比为 0.5 (即 Y=y/B)
    
    % 构造 Michell 积分所需的 3D 偏移矩阵 Y(Nx, Nz)
    % 假设截面形状随深度 z 按正弦规律收缩（模拟典型 U/V 型组合剖面）
    z_norm = linspace(0, 1, Nz); % 0=龙骨, 1=水线
    Y = zeros(Nx, Nz);
    for iz = 1:Nz
        % 深度因子：龙骨处偏移为0，水线处最大
        z_factor = sin(pi/2 * z_norm(iz))^0.5; 
        Y(:, iz) = f0(:) * z_factor;
    end
    
    % --- B. 阻力计算 ---
    % 1. 摩擦阻力 (ITTC-1957)
    S = 1.01 * KCS.L * (2*KCS.T + B * Cb); % 估算湿表面积
    Re = V_ship_mps * KCS.L / Nu;
    Cf = 0.075 / (log10(Re) - 2)^2;
    Rf = 0.5 * Rho_water * S * V_ship_mps^2 * (Cf + 0.0004);
    
    % 2. 兴波阻力 (调用 Michell 积分函数)
    Rw = michell(Y, V_ship_mps, KCS.L, B, KCS.T, Rho_water, N_theta);
    
    % 总阻力
    Rt = Rf + Rw;
    
    % --- C. 交互因子 ---
    w = 0.26 + 1.8 * dCb; 
    t = 0.15 + 0.5 * (w - 0.26); % 设定低跟随率，模拟优化的舵球效果
    eta_H = (1 - t) / (1 - w);
    
    % --- D. 螺旋桨效率 (螺距自适应模拟) ---
    Thrust_req = Rt / (1 - t);
    Va = V_ship_mps * (1 - w);
    % 假设经过 P/D 优化，eta_O 随载荷缓慢变化
    Load_Factor = Thrust_req / (0.5 * Rho_water * Va^2 * pi * (KCS.D_prop/2)^2);
    eta_O = 0.71 - 0.012 * (Load_Factor - 0.5); 
    eta_O = max(0.60, min(0.72, eta_O));
    
    % --- E. 总功率 ---
    eta_R = 1.01;
    eta_D = eta_H * eta_O * eta_R;
    Pd = Rt * V_ship_mps / eta_D;
    
    % 存储结果
    Results(i).Cb = Cb;
    Results(i).Rt = Rt;
    Results(i).Rw = Rw;
    Results(i).Rf = Rf;
    Results(i).w = w;
    Results(i).t = t;
    Results(i).eta_H = eta_H;
    Results(i).eta_O = eta_O;
    Results(i).eta_D = eta_D;
    Results(i).Pd = Pd;
    
    fprintf('%4d   %.4f   %6.0f   %6.0f   %.3f    %.3f    %.3f    %6.2f\n', ...
        i, Cb, Rt/1000, Rw/1000, w, eta_H, eta_D, Pd/1e6);
end

% ==============================================================================
% 3. 结果分析与绘图
% ==============================================================================
Cb_vec = [Results.Cb];
Rt_vec = [Results.Rt]/1000;
Rw_vec = [Results.Rw]/1000;
Pd_vec = [Results.Pd]/1e6;
eta_H_vec = [Results.eta_H];
eta_D_vec = [Results.eta_D];

[~, idx_Rt] = min(Rt_vec);
[~, idx_Pd] = min(Pd_vec);

figure('Position', [100 100 1100 450], 'Color', 'w');
% 左图：阻力成分与功率
subplot(1, 2, 1);
yyaxis left
plot(Cb_vec, Rt_vec, 'b-o', 'LineWidth', 2); hold on;
plot(Cb_vec, Rw_vec, 'b--', 'LineWidth', 1);
ylabel('Resistance (kN)');
yyaxis right
plot(Cb_vec, Pd_vec, 'r-s', 'LineWidth', 2);
ylabel('Delivered Power (MW)');
xlabel('Block Coefficient C_b');
title('Michell-based Resistance & Power');
grid on; xline(Cb_vec(idx_Pd), '--r', 'Optimum Pd');

% 右图：效率增益
subplot(1, 2, 2);
plot(Cb_vec, eta_H_vec, 'g--', 'LineWidth', 1.5); hold on;
plot(Cb_vec, eta_D_vec, 'k-', 'LineWidth', 2);
xlabel('Block Coefficient C_b');
ylabel('Efficiency');
title('Propulsive Efficiency Gain');
legend('\eta_H (Hull Eff)', '\eta_D (Total)', 'Location', 'best');
grid on;

fprintf('\n结论：基于 Michell 积分，最小功率 Cb = %.3f\n', Cb_vec(idx_Pd));

% ==============================================================================
% 4. 辅助函数 (Sub-functions)
% ==============================================================================

function Rw = michell(Y,U,L,B,T,RHO,N)
    % 核心 Michell 积分函数
    Nx= size(Y,1); Nz= size(Y,2);
    YH = Y*B; 
    dz = T/(Nz-1); z = (-T:dz:0)';
    dx = L/(Nx-1); x = (0:dx:L)';
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
            f = zeros(Nz,1);
            for n = 1:Nz
                if n == 1, f(n) = w0(j)*YH(m,n)*exp(k0*z(n)*a(j)^2)*dz;
                elseif n == Nz, f(n) = wN(j)*YH(m,n)*exp(k0*z(n)*a(j)^2)*dz;
                else f(n) = wn(j)*YH(m,n)*exp(k0*z(n)*a(j)^2)*dz;
                end
            end
            F(m,j) = sum(f); 
        end
    end

    % X 积分 (Filon algorithm)
    Kx = k0*dx.*a;
    alp = (Kx.^2+1/2*Kx.*sin(2*Kx)+cos(2*Kx)-1)./Kx.^3;
    bet = (3*Kx+Kx.*cos(2*Kx)-2*sin(2*Kx))./Kx.^3;
    gam = 4*(sin(Kx)-Kx.*cos(Kx))./Kx.^3;
    Pt = zeros(N,1); Qt = zeros(N,1); P = zeros(N,1); Q = zeros(N,1);
    for j = 1:N
        c_val = cos(k0*x*a(j)); s_val = sin(k0*x*a(j));
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
    Rw = sum(rw_vec);
end

function [xm] = michspace(N)
    % 传播角加密采样
    xm = logspace(0,1,N)-1;
    xm = xm*pi/18-pi/2;
    xm = fliplr(-xm);
end

function f = three_param_shape(x, h, a1, a2)
    % 船体水线面几何生成
    f = zeros(size(x));
    for i = 1:length(x)
        if x(i) < 0.5
            f(i) = h - a1/2 * (cos(2*pi*x(i)) - 1);
        else
            f(i) = h + a2/2 * (cos(2*pi*x(i)) + 1) + a1;
        end
    end
end