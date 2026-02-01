function [w, t, eta_R] = Holtrop_Propulsion_Calculation_1982(L, B, T, D, S, Cb, Cp, LCB_percent, Cstern, AE_AO_ratio, Cf, K1)
%HOLTROP_PROPULSION_CALCULATION_1982 Calculate propulsion factors based on Holtrop (1982).
%   Strict implementation of:
%   J. Holtrop and G.G.J. Mennen, "An approximate power prediction method",
%   International Shipbuilding Progress, Vol. 29, July 1982.
%
%   Inputs:
%       L           - Waterline length [m]
%       B           - Moulded breadth [m]
%       T           - Moulded draught [m] (Ta)
%       D           - Propeller diameter [m]
%       S           - Wetted surface area [m^2] (CRITICAL for 1982 method)
%       Cb          - Block coefficient
%       Cp          - Prismatic coefficient
%       LCB_percent - Longitudinal Center of Buoyancy [%L]
%                     (Forward +, Aft -. e.g., 2% aft = -2.0)
%       Cstern      - Stern shape coefficient (-10 to +10)
%       AE_AO_ratio - Expanded area ratio
%       Cf          - Frictional resistance coefficient (ITTC-1957)
%       K1          - Form factor (1 + k1)

    %% 1. Parameter Preparation
    % Source[cite: 145]: Correlation Allowance (Standard)
    CA = 0.006 * (L + 100)^(-0.16) - 0.00205;
    
    % Source[cite: 156]: Viscous Resistance Coefficient
    Cv = K1 * Cf + CA;

    % Source[cite: 157]: Corrected Prismatic Coefficient Cp1
    % Formula: 1.45 Cp - 0.315 - 0.0225 lcb
    Cp1 = 1.45 * Cp - 0.315 - 0.0225 * LCB_percent;
    Cp1 = min(0.99, max(0.1, Cp1)); % Numerical safety only

    %% 2. Wake Fraction (w)
    % Source: PDF Page 3

    % --- Coefficient c8 (Dependent on S) ---
    % Source [cite: 122-126]
    if B/T < 5
        c8 = (B * S) / (L * D * T);
    else
        c8 = (S * (7 * B/T - 25)) / (L * D * (B/T - 3));
    end
    
    % --- Coefficient c9 ---
    % Source [cite: 131-135]
    if c8 < 28
        c9 = c8;
    else
        c9 = 32 - 16 / (c8 - 24);
    end
    
    % --- Coefficient c11 ---
    % Source [cite: 136-139, 155]
    if T/D < 2
        c11 = T/D;
    else
        c11 = 0.0833333 * (T/D)^3 + 1.33333;
    end
    
    % --- Main w Formula ---
    % Source 
    % Term 1: Viscous/Form
    term_w1 = c9 * Cv * (L / T) * (0.0661875 + 1.21756 * c11 * Cv / (1 - Cp1));
    % Term 2: Beam/Propeller
    term_w2 = 0.24558 * sqrt(B / (L * (1 - Cp1)));
    % Term 3: Fullness (Cp and Cb terms)
    safe_Cp = min(0.945, Cp); % Prevent singularity at 0.95
    safe_Cb = min(0.945, Cb);
    term_w3 = -0.09726 / (0.95 - safe_Cp) + 0.11434 / (0.95 - safe_Cb);
    % Term 4: Stern Shape
    term_w4 = 0.75 * Cstern * Cv + 0.002 * Cstern;
    
    w = term_w1 + term_w2 + term_w3 + term_w4;
    
    %% 3. Thrust Deduction (t)
    % Source: PDF Page 3

    % --- Coefficient c10 ---
    % Source [cite: 161-166]
    LB_ratio = L / B;
    if LB_ratio > 5.2
        c10 = B / L;
    else
        % CRITICAL: Denominator is B/L 
        BL_ratio = B / L; 
        c10 = 0.25 - 0.003328402 / (BL_ratio - 0.134615385);
    end
    
    % --- Main t Formula ---
    % Source 
    term_t1 = 0.001979 * L / (B * (1 - Cp1));
    term_t2 = 1.0585 * c10;
    term_t3 = -0.00524;
    term_t4 = -0.1418 * (D^2) / (B * T); % Source [160] says BT
    term_t5 = 0.0015 * Cstern;
    
    t = term_t1 + term_t2 + term_t3 + term_t4 + term_t5;
    
    %% 4. Relative Rotative Efficiency (eta_R)
    % Source 
    eta_R = 0.9922 - 0.05908 * AE_AO_ratio + 0.07424 * (Cp - 0.0225 * LCB_percent);

end