% =============================================== FORTRAN Subroutine WRENCH    
function [UAIF,UTIF]=WRENCH(NBLADE,TANBIW,RCW,RVW)
% =========================================================================
% 【函数名称】：WRENCH
% 【功能】：计算单个马蹄涡（控制点 RCW 受 RVW 上涡影响）在升力线理论中的
%          诱导速度，包括轴向诱导速度 UAIF 和切向诱导速度 UTIF。
% 【适用范围】：升力线理论的涡感应影响系数构建（用于 Main 函数中的 A 矩阵）
% 【参数说明】：
%     NBLADE  - 螺旋桨叶片数
%     TANBIW  - 马蹄涡的诱导角（tan(beta_i)）
%     RCW     - 控制点半径比（被感应点）
%     RVW     - 涡点半径比（感应源位置）
% =========================================================================
    if NBLADE>20    % 如果叶片数大于 20，使用无限叶片近似（解析解更快）
        if RCW>RVW  % 控制点在感应源之后（常见情况）
            UAIF = 0;  % 没有轴向诱导
            UTIF = NBLADE*(RCW-RVW)/RCW;  % 切向诱导速度解析解
        else        % 控制点在涡点之前
            UAIF = -NBLADE*(RCW-RVW)/(RVW*TANBIW);  % 有轴向诱导
            UTIF = 0;  % 切向诱导为零
        end
        return;
    end
 % 有限叶片数（Goldstein 模型）
    XG = 1/TANBIW;     % 求解诱导角反比，用于坐标转换
    ETA = RVW/RCW;     % 源点与目标点半径比
    H = XG/ETA;        % 诱导角相关几何参数
    XS = 1 + H^2;      
    TW = sqrt(XS);     % 中间变量，平方根形式
    V = 1 + XG^2;      
    W = sqrt(V);       % 另一个根号中间变量
    AE = TW - W;       
    U = exp(AE);       % 指数项用于诱导模型
    R = (((TW - 1)/H * (XG / (W - 1))) * U)^NBLADE;  % 旋涡强度函数（带叶片数）
    XX = (1/(2 * NBLADE * XG)) * ((V / XS)^0.25);    % 基础诱导系数
    Y = ((9 * XG^2) + 2) / (V^1.5) + ((3 * H^2 - 2) / (XS^1.5));  % 调整项
    Z = 1 / (24 * NBLADE) * Y;     % 用于修正指数诱导项

    if H >= XG    % 如果诱导点 H 位于主要涡区之后
        AF = 1 + 1 / (R - 1);     % 修正项
        if AF == 0
            UAIF = NaN; UTIF = NaN; return;   % 避免除 0 错误
        else
            AA = XX * (1 / (R - 1) - Z * log(AF));  % 计算诱导系数
            UAIF = 2 * NBLADE^2 * XG * H * (1 - ETA) * AA;     % 轴向诱导速度
            UTIF = NBLADE * (1 - ETA) * (1 + 2 * NBLADE * XG * AA); % 切向诱导速度
        end
    else          % 如果诱导点在涡前面
        if R > 1e-12
            RATIO = 1 / (1 / R - 1);   % 修正项
        else
            RATIO = 0;
        end
        AG = 1 + RATIO;
        if AG == 0
            UAIF = NaN; UTIF = NaN; return;   % 避免无效运算
        else
            AB = -XX * (RATIO + Z * log(AG));  % 另一种诱导系数公式
            UAIF = NBLADE * XG * (1 - 1 / ETA) * (1 - 2 * NBLADE * XG * AB);  % 轴向诱导速度
            UTIF = 2 * NBLADE^2 * XG * (1 - ETA) * AB;  % 切向诱导速度
        end
    end
