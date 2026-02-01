% ====================================================== Main Function Code
function [CT,CP,KT,KQ,WAKE,EFFY,RC,G,VAC,VTC,UASTAR,UTSTAR,TANBC,TANBIC,CDC,CD,KTRY]=...
    Main(MT,ITER,IHUB,RHV,NX,NBLADE,ADVCO,CTDES,HR,HT,CRP,XR,XCHD,XCD,XVA,XVT)
%% =========================================================================
% 输入变量解释：
% -------------------------------------------------------------------------
% MT       - 控制点数量（沿半径方向的离散点个数）
% ITER     - 最大迭代次数（用于收敛控制）
% IHUB     - 是否启用 Hub 镜像涡流（1=启用, 0=不启用）
% RHV      - 镜像涡半径/Hub 半径比
% NX       - 用户输入的半径比 XR 的节点数（例如 10 或 11）
% NBLADE   - 螺旋桨叶片数
% ADVCO    - 前进系数 J = V / (nD)
% CTDES    - 目标推力系数（用于迭代收敛判断）
% HR       - Hub 卸载因子（0=不卸载，>0 对 hub 区减载）
% HT       - Tip 卸载因子（0=不卸载，>0 对尖部减载）
% CRP      - 旋涡抵消因子（1=不抵消，0=完全抵消）
% XR       - 用户定义的半径比 r/R 分布（用于插值）
% XCHD     - 用户输入的弦长比 c/D 分布
% XCD      - 用户输入的剖面阻力系数 Cd 分布
% XVA      - 用户输入的轴向流速比 Va/Vs 分布
% XVT      - 用户输入的切向流速比 Vt/Vs 分布
% =========================================================================
% 输出变量解释：
% -------------------------------------------------------------------------
% CT       - 推力系数（每次迭代的 CT(k)）
% CP       - 功率系数（每次迭代的 CP(k)）
% KT       - 非维推力系数（最终）
% KQ       - 非维扭矩系数（最终）
% WAKE     - 平均诱导速度系数（轴向速度比）
% EFFY     - 效率 η = CT * WAKE / CP
% RC       - 控制点位置（cosine spacing 的 r/R 分布）/////////
% G        - 环量分布（Γ）
% VAC      - 控制点的轴向速度（插值得到）
% VTC      - 控制点的切向速度（插值得到）
% UASTAR   - 环量诱导的轴向速度（非维）
% UTSTAR   - 环量诱导的切向速度（非维）
% TANBC    - 来流攻角（tanβc）
% TANBIC   - 诱导攻角（tanβi），迭代过程中更新
% CDC      - 控制点的弦长比 c/D
% CD       - 控制点的阻力系数 Cd
% KTRY     - 成功收敛所需的迭代次数
% =========================================================================
    global Parametric_Flag Single_Flag
    % ---------------------------------------------- FORTRAN Function VOLWK
    % 计算“尾流诱导速度因子（WAKE）
    YW = XR.*XVA;
    YDX = trapz(XR,YW);
    WAKE = 2*YDX/(1-XR(1)^2);
    % -------------------------------------------------------- End of VOLWK
    % =========================================================================
    %【目的】构造用于升力线理论的马蹄涡离散点 RV，以及对应的速度插值与攻角初估
    %        本段采用余弦加密方式布置 RV，用于后续环量感应计算
    % =========================================================================
    XRC = 1 - sqrt(1 - XR);             % 将原始半径比 XR 做映射，用于弦长比插值（后续 RCWG 用到）
    DEL = pi / (2 * MT);                % 控制点间的余弦角度间隔，用于布置马蹄涡点 RV（0 到 π/2） % Compute cosine spaced vortex radii
    HRR = 0.5 * (XR(NX) - XR(1));       % 半径方向的一半跨度（从 hub 到叶尖的间距的一半）
    for i = 1:MT+1
        RV(i) = XR(1) + HRR * (1 - cos(2 * (i - 1) * DEL));  
        % RV 是用于环量感应计算的“马蹄涡放置点”；
        % 采用余弦 spacing，更密集分布在 hub 和 tip 附近
    end
    % 在马蹄涡位置 RV 上插值得到轴向、切向速度
    VAV = pchip(XR, XVA, RV);           % RV 位置上的轴向速度 Va/Vs（插值）
    VTV = pchip(XR, XVT, RV);           % RV 位置上的切向速度 Vt/Vs（插值）
    % 根据局部速度分布，估算马蹄涡处的来流攻角（tanβ）
    TANBV = VAV ./ ((pi .* RV ./ ADVCO) + VTV);  
    % 式中分母是：πr/J + Vt，来自螺旋桨前进几何
    % 用于估算诱导速度比 Vb/Va，即切向速度对轴向速度的比例
    VBAV = VTV .* TANBV ./ VAV;
    % 注意：VBAV 后续用于估算诱导角与流动旋涡结构

    % Cosine spaced control point radii: Evaluate c/D,Va,Vt,tanB,Cd,Vt*,tanB/Va
        % === Step：构建控制点（RC）并插值各物理量用于升力线理论 ===
    % 使用 cosine spacing 构造控制点（r/R），并在控制点处插值剖面参数与流动参数
    % 包括弦长比、阻力系数、轴向/切向速度比、来流攻角等
    for i=1:MT      
        RC(i) = XR(1) + HRR * (1 - cos((2*i - 1) * DEL));  % 控制点 r/R，使用余弦分布（cosine spacing）
        RCWG(i) = 1 - sqrt(1 - RC(i));                     % 用于弦长插值的变换变量（从圆盘投影转换）
    end
    CDC = pchip(XRC, XCHD, RCWG);                          % 插值得到控制点上的弦长比 c/D
    CD  = pchip(XR, XCD, RC);                              % 插值得到控制点上的阻力系数 Cd
    VAC = pchip(XR, XVA, RC);                              % 插值得到控制点上的轴向流速比 Va/Vs
    VTC = pchip(XR, XVT, RC);                              % 插值得到控制点上的切向流速比 Vt/Vs
    TANBC = VAC ./ (pi .* RC ./ ADVCO + VTC);              % 计算控制点的来流攻角 tanβc
    VBAC  = VTC .* TANBC ./ VAC;                           % 计算控制点的切向速度比 Vt/Va

    % First estimation of tanBi based on 90% of actuator disk efficiency
        % === 初始估算诱导攻角 tanβi，用于构造环量方程的初始解，基于90%推进盘效率模型 ===
    EDISK = 1.8/(1+sqrt(1+CTDES/WAKE^2));                      % 基于目标推力系数 CTDES 和 WAKE 推导的理论效率（近似90%推进盘效率）
    TANBXV = TANBV.*sqrt(WAKE./(VAV-VBAV))/EDISK;              % 在马蹄涡点（RV）处估算诱导攻角 tanβi（非维），用于诱导速度计算
    TANBXC = TANBC.*sqrt(WAKE./(VAC-VBAC))/EDISK;              % 在控制点（RC）处估算诱导攻角 tanβi（非维），用于控制点速度分布

    % Unload hub and tip as specified by HR and HT
       % 施加 Hub（根部）与 Tip（尖部）区域的卸载系数 HR 和 HT，用于调整诱导角分布，使载荷分布更合理。
       % 总结：“在叶片根部和叶尖区域，把诱导攻角 tanβi 稍微减小一点，减小载荷，模拟工程中真实的螺旋桨卸载设计。”
    RM = (XR(1)+XR(NX))/2;  % 计算半径中点 RM，用于区分 Hub 区（r<RM）和 Tip 区（r>=RM）

    for i=1:MT+1
        if RV(i)<RM           % 若当前马蹄涡位置在 RM 左侧（Hub 区）
            HRF=HR;           % 使用 Hub 卸载系数
        else                 
            HRF=HT;           % 使用 Tip 卸载系数
        end
        % 计算卸载量 DTANB：根据卸载因子 HRF 和位置偏移（越远越大），将诱导攻角减小一定比例
        DTANB = HRF*(TANBXV(i)-TANBV(i))*((RV(i)-RM)/(XR(1)-RM))^2;
        TANBXV(i) = TANBXV(i)-DTANB;  % 应用卸载，调整诱导攻角 tanβ
    end

    for i=1:MT
        if RC(i)<RM           % 若控制点位于 Hub 区
            HRF=HR;           % 使用 Hub 卸载系数
        else                 
            HRF=HT;           % 使用 Tip 卸载系数
        end
        % 同样地计算控制点处的卸载量，并修正诱导攻角
        DTANB = HRF*(TANBXC(i)-TANBC(i))*((RC(i)-RM)/(XR(1)-RM))^2;
        TANBXC(i) = TANBXC(i)-DTANB;  % 应用卸载到诱导攻角 tanβi
    end

    % Iterations to scale tanBi to get desired value of thrust coefficient
        % =======【步骤5】通过迭代调整 tanβi，使推力系数 CT 收敛到目标值 CTDES =======
        % 这段代码通过多次迭代调整诱导角 tanβi 的缩放比例，使计算出来的螺旋桨推力系数 CT 尽量接近用户设定的目标值 CTDES，
        % 从而实现自适应调节环量强度，保证性能符合设计要求。
    for KTRY=1:ITER     % 最多迭代 ITER 次（通常为10）
        if KTRY==1
            T(KTRY) = 1;  % 第一次迭代缩放系数 T=1，不缩放
        elseif KTRY==2
            T(KTRY) = 1 + (CTDES - CT(1)) / (5 * CTDES);  % 第二次采用线性预测，微调缩放比例
        elseif KTRY > 2
            if CT(KTRY-1) - CT(KTRY-2) == 0  % 防止分母为0，若两次CT一致则跳出
                break
            else
                % 使用 Aitken 加速法修正缩放系数 T：让 CT 尽快接近目标 CTDES
                T(KTRY) = T(KTRY-1) + ...
                    (T(KTRY-1) - T(KTRY-2)) * (CTDES - CT(KTRY-1)) / ...
                    (CT(KTRY-1) - CT(KTRY-2));
            end
        end
        TANBIV = T(KTRY) .* TANBXV;  % 将缩放系数应用到 vortex 点处的诱导角 tanβi
        TANBIC = T(KTRY) .* TANBXC;  % 将缩放系数应用到 control 点处的诱导角 tanβi

        % Compute axial and tangential horseshoe influence coefficients
                % 计算轴向与切向的马蹄涡诱导速度影响系数矩阵（用于求解环量）
        for i=1:MT      
            RCW = RC(i);  % 控制点半径比 r/R（接收诱导速度的位置）
            for j = 1:MT+1
                RVW = RV(j);  % 涡段位置（马蹄涡放置点的 r/R）
                TANBIW = TANBIV(j);  % 当前涡段的诱导攻角 tanβi
                [UAIF,UTIF]=WRENCH(NBLADE,TANBIW,RCW,RVW);  % 调用子程序计算该涡段对控制点的轴向/切向诱导速度

                if (isnan(UAIF)==1)||(isnan(UTIF)==1)  % 若结果为非法数，报错退出
                    if Single_Flag==1
                       set(Err_Single,'visible','on','enable','on','string','Error in MPVL.'); 
                    end
                    return;
                end

                UAW(j) = -UAIF / (2 * (RC(i) - RV(j)));  % 计算马蹄涡在该控制点处的单位轴向诱导速度系数
                UTIF = UTIF * CRP;  % 根据涡流抵消因子 CRP 修改切向诱导速度（1 = 不抵消）
                UTW(j) = UTIF / (2 * (RC(i) - RV(j)));  % 计算单位切向诱导速度系数
                % Induction of corresponding hub-image trailing vorticies
                           % ======== 若启用 Hub 镜像涡流（模拟轴对称影响），计算其诱导速度 ========
                if IHUB==1      
                    RVW = XR(1)^2/RV(j);                     % 镜像涡点的半径位置（关于 hub 半径对称）
                    TANBIW = TANBIV(1)*RV(1)/RVW;            % 镜像涡点处的诱导角估计（与原涡对应）
                    [UAIF,UTIF]=WRENCH(NBLADE,TANBIW,RCW,RVW);  % 调用诱导速度函数，计算诱导速度分量
                    if (isnan(UAIF)==1)||(isnan(UTIF)==1)    % 若计算结果为 NaN，视为出错
                        if Single_Flag==1
                           set(Err_Single,'visible','on','enable','on',...  % 显示错误提示信息
                               'string','Error in MPVL.');
                        end
                        return;                              % 终止运行
                    end
                    UAW(j) = UAW(j)+UAIF/(2*(RC(i)-RVW));    % 将轴向诱导速度叠加至总诱导（Hub 镜像贡献）
                    UTIF = UTIF*CRP;                         % 施加旋涡抵消因子（如设置为 0，则完全抵消）
                    UTW(j) = UTW(j)-UTIF/(2*(RC(i)-RVW));    % 将切向诱导速度叠加（方向相反，因此为减）
                end

            end
            % Final step in building influence functions
                           % 最后一步：构建马蹄涡诱导速度的影响系数矩阵（轴向 + 切向）
            for k=1:MT  
                UAHIF(i,k) = UAW(k+1)-UAW(k);            % 控制点 i 处第 k 个轴向诱导速度影响（由 UAW 差值得出）
                UTHIF(i,k) = UTW(k+1)-UTW(k);            % 控制点 i 处第 k 个切向诱导速度影响（由 UTW 差值得出）
            end

        end
        % Solve simutaneous equations for circulation strengths G(i)
                % ======= 求解升力线理论中的线性方程组 A * G = B，以获得环量 G(i) =======
        for m=1:MT  
            B(m) = VAC(m) * ((TANBIC(m)/TANBC(m)) - 1);             % 构建右端项 B：基于来流与诱导攻角比差计算所需诱导速度
            for n=1:MT
                A(m,n) = UAHIF(m,n) - UTHIF(m,n) * TANBIC(m);       % 构建系数矩阵 A：轴向和切向诱导系数组合
            end
        end
        % ======================================= FORTRAN Subroutine SIMEQN    
        NEQ = length(B);                            % 方程个数（即控制点个数）
        IERR = 1;                                   % 初始化错误标志为 1，表示尚未成功求解
        % Find |maximum| element in each row and exit if a zero row is detected
               % ======================= SIMEQN 子模块：高斯消元解线性方程组 A·G = B =======================
        % 功能：使用带主元选取的高斯消元法，稳定求解升力线理论中的 G（环量）分布
        % 若矩阵 A 是奇异矩阵（不可逆），则返回 NaN 并提示错误，防止程序崩溃
        % 变量说明：
        %   A       - 系数矩阵（MT×MT）
        %   B       - 右端项（长度 MT）
        %   G       - 解向量（环量分布）
        %   IPIVOT  - 用于行交换的行指针
        %   D       - 每行最大元素，用于数值稳定判断

        % 逐行寻找最大元素，判断是否存在零行（奇异矩阵）
        for i=1:NEQ
            IPIVOT(i) = i;                         % 初始化行指针数组
            ROWMAX = 0;                            % 初始化该行最大值
            for j=1:NEQ
                ROWMAX = max(ROWMAX, abs(A(i,j))); % 寻找每行最大绝对值
            end
            if ROWMAX==0                           % 若全为 0，说明奇异
                fprintf('Matrix is Singular-1.\n') % 打印错误信息
                G = NaN;                           % 返回空值，避免报错
                if Single_Flag==1                  % 若是单桨模式，弹出GUI提示
                    set(Err_Single,'visible','on','enable','on','string','Error in MPVL.');
                end
                return;                            % 中断函数
            end
            D(i) = ROWMAX;                         % 保存每行最大值
        end

        NM1=NEQ-1;                                  % NEQ = 方程个数；NM1 = NEQ-1
        if NM1>0                                    % 多于1个变量才进行消元
            for k=1:NM1
                j = k;
                KP1 = k+1;
                IP = IPIVOT(k);                     % 当前主行行号
                COLMAX = abs(A(IP,k))/D(IP);        % 当前列最大比值
                for i=KP1:NEQ
                    IP = IPIVOT(i);                 % 遍历其余行
                    AWIKOV = abs(A(IP,k))/D(IP);    % 当前行该列缩放值
                    if AWIKOV>COLMAX                % 找更大主元
                        COLMAX = AWIKOV;
                        j=i;
                    end
                end
                if COLMAX==0                        % 整列为0 → 奇异
                    fprintf('Matrix is Singular-2.\n')
                    G = NaN;
                    if Single_Flag==1
                        set(Err_Single,'visible','on','enable','on','string','Error in MPVL.');
                    end
                    return;
                end

                % 行交换，确保主元最大
                IPK = IPIVOT(j);
                IPIVOT(j) = IPIVOT(k);
                IPIVOT(k) = IPK;

                % 消元：将当前列以下元素置为 0
                for i=KP1:NEQ
                    IP = IPIVOT(i);
                    A(IP,k) = A(IP,k)/A(IPK,k);             % 归一化
                    RATIO = -A(IP,k);
                    for j=KP1:NEQ
                        A(IP,j) = RATIO*A(IPK,j)+A(IP,j);   % 消元
                    end
                end
            end
            if A(IP,NEQ)==0                         % 最后一项为 0 → 奇异
                fprintf('Matrix is Singular-3.\n')
                G = NaN;
                if Single_Flag==1
                    set(Err_Single,'visible','on','enable','on','string','Error in MPVL.');
                end
                return;
            end
        end

        IERR = 0;                                    % 矩阵通过了奇异性测试，可继续求解

        if NEQ==1                                    % 如果只有一个变量，直接除法
            G(1) = B(1)/A(1,1);
        else
            % 前向替代：求 G(1), G(2), ..., G(n-1)
            IP = IPIVOT(1);
            G(1) = B(IP);
            for k=2:NEQ
                IP = IPIVOT(k);
                KM1 = k-1;
                SUMM = 0;
                for j=1:KM1
                    SUMM = A(IP,j)*G(j)+SUMM;        % 累加已知项
                end
                G(k) = B(IP)-SUMM;                   % 得到中间结果
            end

            % 后向回代：回代解出 G(n)
                G(NEQ) = G(NEQ)/A(IP,NEQ);
            k = NEQ;
            for NP1MK=2:NEQ
                KP1 = k;
                k = k-1;
                IP = IPIVOT(k);
                SUMM = 0;
                for j=KP1:NEQ
                    SUMM = A(IP,j)*G(j)+SUMM;        % 累加已知项
                end
                G(k) = (G(k)-SUMM)/A(IP,k);          % 得到最终解
            end
        end
        % =================================================== End of SIMEQN    
                % ===== 判断线性方程组是否求解失败（奇异矩阵），若失败则报错并退出函数 =====
        if IERR==1  % Matrix is singular（矩阵奇异，无法求解环量 G）
            if Single_Flag==1  % 如果当前是单桨设计模式
                set(Err_Single,'visible','on','enable','on','string','Error in MPVL.'); % 在界面显示错误提示按钮
            end
            return;  % 终止当前函数，退出求解
        end
        % Evalutate the induced velocities from the circulation G(i)
               %========================= 计算由环量 G(i) 产生的诱导速度 =========================
        for p = 1:MT       % 对于每个控制点 p（沿半径方向）
            UASTAR(p) = 0;         % 初始化轴向诱导速度为 0（非维）
            UTSTAR(p) = 0;         % 初始化切向诱导速度为 0（非维）
            for q = 1:MT           % 遍历所有马蹄涡元（源点 q）
                UASTAR(p) = UASTAR(p) + G(q) * UAHIF(p,q);  % 累加 q 段马蹄涡对控制点 p 的轴向感应速度
                UTSTAR(p) = UTSTAR(p) + G(q) * UTHIF(p,q);  % 累加 q 段马蹄涡对控制点 p 的切向感应速度
            end
        end
        % ======================================= FORTRAN Subroutine FORCES
               % ============================= 计算推力、扭矩、功率等性能指标（FORTRAN 子程序 FORCES 对应部分）
        LD = 0;     % 默认认为 CD 是阻力系数（非升阻比）
        if CD>1
            LD = 1;  % 若 CD > 1，表示输入的是升阻比（L/D），需要特殊处理
        end
        CT(KTRY) = 0;        % 初始化推力系数累加器
        CQ(KTRY) = 0;        % 初始化扭矩系数累加器
        for m=1:MT
            DR = RV(m+1)-RV(m);             % 控制点之间的径向间距（Δr）
            VTSTAR = VAC(m)/TANBC(m)+UTSTAR(m); % 切向合速度（Va + Ut*）
            VASTAR = VAC(m)+UASTAR(m);      % 轴向合速度（Va + Ua*）
            VSTAR = sqrt(VTSTAR^2+VASTAR^2);% 合速度模长
            if LD==0
                DVISC = (VSTAR^2*CDC(m)*CD(m))/(2*pi); % 若 CD 是阻力系数，用公式计算粘性阻力
            else
                FKJ = VSTAR*G(m);           % 若输入为升阻比，FKJ 表示升力功
                DVISC = FKJ/CD(m);          % 粘性阻力通过升阻比换算
            end
            CT(KTRY) = CT(KTRY)+(VTSTAR*G(m)-DVISC*VASTAR/VSTAR)*DR; % 当前点的推力贡献累加
            CQ(KTRY) = CQ(KTRY)+(VASTAR*G(m)+DVISC*VTSTAR/VSTAR)*RC(m)*DR; % 当前点的扭矩贡献累加
        end
        if IHUB~=0
            CTH = .5*(log(1/RHV)+3)*(NBLADE*G(1))^2; % 若启用 hub 镜像涡流，添加 hub 区涡流影响
        else
            CTH = 0;                                % 否则不计算 hub 区涡流影响
        end
        CT(KTRY) = CT(KTRY)*4*NBLADE-CTH;           % 总推力系数，乘以螺旋桨叶片数并减去 hub 修正
        CQ(KTRY) = CQ(KTRY)*2*NBLADE;               % 总扭矩系数
        CP(KTRY) = CQ(KTRY)*2*pi/ADVCO;             % 功率系数（与扭矩系数和前进系数有关）
        KT(KTRY) = CT(KTRY)*ADVCO^2*pi/8;           % 非维推力系数，用于换算
        KQ(KTRY) = CQ(KTRY)*ADVCO^2*pi/8;           % 非维扭矩系数
        EFFY(KTRY) = CT(KTRY)*WAKE/CP(KTRY);        % 螺旋桨推进效率 η = 推力功 / 输入功
        for i=1:length(RC)
            if (isreal(TANBIC(i))==0)||(isreal(EFFY(KTRY))==0)||(EFFY(KTRY)<=0)
                TANBIC(i) = NaN;                    % 若诱导攻角无效或效率为负，则置为无效值
                EFFY(KTRY) = NaN;                   % 整体效率标记为无效
            end
        end
        % =================================================== End of FORCES
                %=============================== 迭代终止条件：判断推力是否收敛 ===============================
        if abs(CT(KTRY)-CTDES)<(5e-6)   % 如果当前迭代得到的推力系数接近目标值，说明收敛，退出循环
            break                       % 结束迭代
        end
    end                                 % KTRY=1:ITER 的 for 循环结束（最多迭代 ITER 次）

    %======================== 如果线性方程组奇异（无法求解），直接终止程序 ==========================
    if IERR==1                          % 如果迭代过程中矩阵奇异，IERR=1
       fprintf('Matrix is Singular. Run Terminated.\n');     % 打印终止信息
       return;                         % 终止函数执行
    else
        %======================= 如果是单螺旋桨设计模式，输出详细结果文件 ==========================
        % Single_Flag==1;
        % if Single_Flag==1              % 如果当前处于“单桨设计”模式
            TANBC = atand(TANBC);      % 将 tanβc 转换为角度（°）
            TANBIC = atand(TANBIC);    % 将 tanβi 转换为角度（°）
            fid = fopen('PVL_Output.txt','w');              % 打开结果输出文件
            fprintf(fid,'\t\t\t\t\t PVL_Output.txt\n');      
            fprintf(fid,'\t\t\t\t\tPVL Output Table\n');
            fprintf(fid,'Ct= % 5.4f\n' ,CT(KTRY));           % 输出最终推力系数
            fprintf(fid,'Cp= % 5.4f\n' ,CP(KTRY));           % 输出最终功率系数
            fprintf(fid,'Kt= % 5.4f\n' ,KT(KTRY));           % 输出非维推力系数
            fprintf(fid,'Kq= % 5.4f\n' ,KQ(KTRY));           % 输出非维扭矩系数
            fprintf(fid,'Va/Vs= % 5.4f\n' ,WAKE);            % 输出平均轴向诱导速度比
            fprintf(fid,'Efficiency= %5.4f\n' ,EFFY(KTRY));  % 输出效率 η
            fprintf(fid,' r/R\t  \tG\t\t  Va\t  Vt\t  Ua\t  \tUt\t \tBeta\tBetaI\t c/D\t  Cd\t\n');
            for i = 1:length(RC)                              % 遍历所有控制点
            fprintf(fid,'%5.5f  %5.6f  %5.5f  %5.4f  %5.5f  %5.5f  %5.3f  %5.3f  %5.5f  %5.5f\n',...
                RC(i),G(i),VAC(i),VTC(i),UASTAR(i),UTSTAR(i),TANBC(i),TANBIC(i),CDC(i),CD(i));
                % 输出每个控制点的各项参数：半径比、环量、速度分布、攻角、弦长比、阻力系数等
            end
            fclose(fid);               % 关闭文件
        % end
    end

% ==================================================== End of Main Function    