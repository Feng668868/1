
function [Sigma,skew,rake,EFFY,eta_final,XVA]=...
    PVL(Common_Def,XR0,XCHD_def,XCD_def,XVA_def,XVT_def,f0oc_def,t0oc_def,skew_def,rake_def,Single_def1,Single_def2,Mean,Thick,T)
% ---------------------------------------------------- Define Variables
 % ------------------ 定义参数变量（来自用户输入）-----------------------
 NBLADE = Single_def1(1);        % 螺旋桨叶片数
    N = Single_def1(2);             % 螺旋桨转速（RPM）
    D = Single_def1(3);             % 螺旋桨直径（m）

    THRUST = T;        % 目标推力（单位N）
    V =  Common_Def(2);             % 船速/前进速度（m/s）
    Dhub =  Common_Def(3);          % hub（轮毂）直径（m）
    MT =  Common_Def(4);            % 半径方向上的面板数量（用于升力线控制点数量）
    ITER =  Common_Def(5);          % 最大迭代次数（用于收敛控制）
    RHV =  Common_Def(6);           % hub镜像涡半径 / hub半径的比值

    NX =  Common_Def(7);            % 用户输入的径向分布节点数
    HR =  Common_Def(8);            % hub卸载因子（控制根部减载程度）
    HT =  Common_Def(9);            % tip卸载因子（控制叶尖减载程度）
    CRP =  Common_Def(10);          % 旋涡抵消因子（1=不抵消，0=完全抵消）
    rho =  Common_Def(11);          % 水体密度（kg/m³）

    H = Single_def2(1);             % 螺旋桨轴心处水深（shaft centerline depth）
    dV = Single_def2(2);            % 来流速度扰动量（m/s）
    AlphaI = Single_def2(3);        % 理想攻角（度）
    NP = Single_def2(4);            % 弦向剖面点数（用于剖面几何建模）

    R = D/2;                    % 螺旋桨半径
    Rhub = Dhub/2;              % hub 半径
    n = N/60;                   % 螺旋桨转速（单位：转/秒）
    ADVCO = V/(n*D);            % 前进系数 J = V / (nD)
    CTDES = THRUST/(rho*V^2*pi*R^2/2);  % 推力系数（非维）：CT = T / (0.5 * rho * V^2 * A)



% % % % %

    IHUB = 1;              % 是否启用 hub 镜像涡流（复选框）,默认勾选（1=启用，0=禁用）




    Meanline = Mean(1);        % 用户选择的翼型平均线类型(1,2)（NACA或抛物线）
    Thickness = Thick(1);       % 用户选择的翼型厚度类型(1,2,3)（65A010、椭圆、抛物线）

    XCHD0 = XCHD_def;  % 输入的弦长比 c/D（沿 r/R）
    XCD0 = XCD_def;    % 输入的阻力系数 Cd（沿 r/R）
    XVA0 = XVA_def;    % 输入的轴向速度比 Va/Vs（沿 r/R）
    XVT0 = XVT_def;    % 输入的切向速度比 Vt/Vs（沿 r/R）

    f0oc0 = f0oc_def;  % 用户输入的最大弯度比 f0/c（用于剖面）
    t0oc0 = t0oc_def;  % 用户输入的最大厚度比 t0/c
    skew0 = skew_def;  % 用户输入的掠角分布（角度）
    rake0 = rake_def;  % 用户输入的倾角分布（Xs/D）

    if Dhub/D < .15                              % 检查 hub 与总直径的比例（小于15%报错）
        set(Err_Single,'visible','on','enable','on','string','Dhub/D >= 15%');
        return;
    end
























   % ------------------------------------------------- Make PVL_Input.txt
        % 生成输入文件 PVL_Input.txt，用于记录当前设计参数和剖面分布，供用户查看或记录计算输入
    % ------------------------------------------------- 生成 XR 插值点（用于插值剖面几何参数）
    XR1 = Rhub/R:(1-Rhub/R)/(NX-3):1;                % 在 hub 到叶尖之间，均匀分布 NX-3 个点（不含两侧中点）
    half1 = (XR1(1)+XR1(2))/2;                       % 计算第一个中点，用于精细分布
    half2 = (XR1(NX-3)+XR1(NX-2))/2;                 % 计算最后一个中点
    XR = [XR1(1) half1 XR1(2:NX-3) half2 XR1(NX-2)]; % 插入两个中点，构建 NX 个半径比点 XR（最终用于插值）

    % ------------------------------------------------- 对几何参数进行插值（从 XR0 → XR）
    XCHD = pchip(XR0,XCHD0,XR);                      % 插值：弦长比 c/D 分布
    XCD = pchip(XR0,XCD0,XR);                        % 插值：阻力系数 Cd 分布
    XVA = pchip(XR0,XVA0,XR);                        % 插值：轴向速度比 Va/Vs 分布
    XVT = pchip(XR0,XVT0,XR);                        % 插值：切向速度比 Vt/Vs 分布


    % % % % % % % % % % %% 新加两行////////////////////////////////////////////////////////////////////////////////////
    % % % % % % 
    % % % % % % f0oc0 = pchip(XR0,f0oc0,XR);  % 用户输入的最大弯度比 f0/c（用于剖面）
    % % % % % % t0oc0 = pchip(XR0,t0oc0,XR);  % 用户输入的最大厚度比 t0/c


    % ------------------------------------------------- 创建并写入 MPVL_Input.txt 文件
    Flag1 = datestr(now,31);                         % 获取当前时间戳，作为记录头
    fid = fopen('PVL_Input.txt','w');               % 打开文件用于写入（w 模式会覆盖原文件）

    fprintf(fid,'%s\tPVL_Input.txt\n',Flag1);       % 写入时间戳和文件名说明
    fprintf(fid,'%.0f \t\tNumber of Vortex Panels over the Radius\n',MT);      % 写入控制点数 MT
    fprintf(fid,'%.0f \t\tMax. Iterations in Wake Alignment\n',ITER);         % 写入最大迭代次数 ITER
    fprintf(fid,'%.0f \t\tHub Image Flag: 1=YES, 0=NO\n',IHUB);               % 写入是否使用 hub 镜像涡流
    fprintf(fid,'%.1f \tHub Vortex Radius/Hub Radius\n',RHV);                 % 写入镜像涡流半径比
    fprintf(fid,'%.0f \t\tNumber of Input Radii\n',NX);                       % 写入插值点个数 NX
    fprintf(fid,'%.0f \t\tNumber of Blades\n',NBLADE);                        % 写入叶片数
    fprintf(fid,'%.3f \tAdvance Coef., J, Based on Ship Speed\n',ADVCO);      % 写入前进系数 J
    fprintf(fid,'%.3f \tDesired Thrust Coef., Ct\n',CTDES);                   % 写入目标推力系数
    fprintf(fid,'%.0f \t\tHub Unloading Factor: 0=optimum\n',HR);            % 写入 Hub 卸载因子
    fprintf(fid,'%.0f \t\tTip Unloading Factor: 1=Reduced Loading\n',HT);    % 写入 Tip 卸载因子
    fprintf(fid,'%.0f \t\tSwirl Cencellation Factor: 1=No Cancellation\n',CRP); % 写入旋涡抵消因子

    fprintf(fid,'r/R  \t  C/D  \t   XCD\t    Va/Vs  Vt/Vs\n');               % 写入表头
    for i = 1:NX
        fprintf(fid,'%6.5f  %6.5f  %6.5f  %6.2f  %6.4f\n',XR(i),XCHD(i),XCD(i),XVA(i),XVT(i)); % 写入每个 r/R 位置的剖面数据
    end
    fclose(fid);                                                          % 关闭文件

















    % ================================================== Call Main Function
   % [CT,CP,KT,KQ,WAKE,EFFY,RC,G,VAC,VTC,UASTAR,UTSTAR,TANBC,TANBIC,CDC,CD,KTRY]=...
   %  Main2(MT,ITER,IHUB,RHV,NX,NBLADE,ADVCO,CTDES,HR,HT,CRP,XR,XR0,XCHD,XCD,XVA,XVT,f0oc0,t0oc0);
  [CT,CP,KT,KQ,WAKE,EFFY,RC,G,VAC,VTC,UASTAR,UTSTAR,TANBC,TANBIC,CDC,CD,KTRY]=...
    Main(MT,ITER,IHUB,RHV,NX,NBLADE,ADVCO,CTDES,HR,HT,CRP,XR,XCHD,XCD,XVA,XVT);


    % -------------------------------------------- Create Graphical Reports
        % ---------------------------- 创建性能图形报告窗口，绘制计算结果（升力线理论解）
    Fig1_S = figure('units','normalized','position',[.01 .06 .4 .3],'name',...
        'Graphical Report','numbertitle','off');  % 创建新的图形窗口，位置左下角

    subplot(2,2,1);         % 第一个子图（2行2列中的第1个）
    plot(RC,G);             % 绘制环量分布（Γ）对 r/R 的变化
    xlabel('r/R');          % X轴标签：半径比
    ylabel('Non-Dimensional Circulation');  % Y轴标签：无量纲环量
    grid on;                % 添加网格线
    TitleString=strcat('J=',num2str(ADVCO,'%10.3f'),'; Ct=',num2str(CT(KTRY),'%10.3f'),...
        '; Kt=',num2str(KT(KTRY),'%10.3f'),'; Kq=',num2str(KQ(KTRY),'%10.3f'),...
        '; \eta=',num2str(EFFY(KTRY),'%10.3f'));  % 构造标题：包括 J、Ct、Kt、Kq、效率
    title(TitleString);     % 设置标题

    subplot(2,2,2);         % 第二个子图：速度分布
    plot(RC,VAC,'-b',RC,VTC,'--b',RC,UASTAR,'-.r',RC,UTSTAR,':r');  % 绘制轴向、切向流速及其诱导量
    xlabel('r/R');          % X轴标签
    legend('Va/Vs','Vt/Vs','Ua*/Vs','Ut*/Vs');  % 图例说明各条曲线含义
    grid on;                % 添加网格线

    subplot(2,2,3);         % 第三个子图：攻角对比
    plot(RC,TANBC,'--b',RC,TANBIC,'-r');  % 绘制原始攻角 vs 诱导攻角
    xlabel('r/R');          % X轴标签
    ylabel('Degrees');      % Y轴单位为角度（实际上是 atan 过的）
    grid on;                % 网格线
    legend('Beta','BetaI'); % 图例说明

    subplot(2,2,4);         % 第四个子图：弦长分布
    plot(RC,CDC);           % 绘制剖面弦长比 c/D 随半径变化曲线
    xlabel('r/R');          % X轴标签
    ylabel('c/D');          % Y轴标签：无量纲弦长
    grid on;                % 添加网格

    % ----------------------------------- Propeller Performance Calculation
       % -------------------- 螺旋桨性能计算部分（根据升力线结果计算各物理参数） --------------------
    w = 2*pi*n;         % w：角速度（rad/s），n 是转速（转/秒）
    for k = 1:MT
        Vstar(k) = sqrt((VAC(k)+UASTAR(k))^2 + (w*R*RC(k)+VTC(k)+UTSTAR(k))^2); % 控制点的相对流速大小
        Gamma(k) = G(k)*2*pi*R*V;                                               % 环量（真实单位），单位：m^2/s
        Cl(k) = 2*Gamma(k) / (Vstar(k)*CDC(k)*D);                               % 升力系数 Cl，基于局部速度和弦长
        dBetai(k) = atand((tand(TANBIC(k))*w*RC(k)*R+dV)/(w*RC(k)*R))...        % 迎角偏差 dBeta_i，考虑速度扰动（+dV）
                   -atand((tand(TANBIC(k))*w*RC(k)*R-dV)/(w*RC(k)*R));          % 与扰动（-dV）相比的偏差
        Sigma(k) = (101000+rho*9.81*(H-RC(k)*R)-2500)/(rho*Vstar(k)^2/2);       % 空泡数 σ，用于 cavitation 预测///////////////////////////
    end
    f0oc = pchip(XR0,f0oc0,RC).*Cl;     % 插值 f0/c（无升力弯度），并根据升力系数缩放（反映实际升力形状）
    t0oc = pchip(XR0,t0oc0,RC);         % 插值 t0/c（厚度分布），按 RC 控制点位置插值
    fid = fopen('Performance.txt','w');                 % 创建输出性能文件 Performance.txt
    fprintf(fid,'\t\t\t\t\t\tPerformance.txt\n');       % 写文件标题
    fprintf(fid,'\t\t\t\t\tPropeller Performance Table\n');
    fprintf(fid,' r/R\t\tV*\t beta\t betai\t  Gamma\t\tCl\t Sigma\tdBetai\n'); % 写表头
    for k = 1:MT
        fprintf(fid,'%.3f\t %.3f\t %.2f\t %.2f\t %.4f\t %.3f\t %.3f\t %.2f\n'... % 按行写入各控制点数据
        ,RC(k),Vstar(k),TANBC(k),TANBIC(k),Gamma(k),Cl(k),Sigma(k),dBetai(k));
    end
    fclose(fid);    % 关闭文件

    % ------------------------------------------------ Geometry Calculation
       % ------------------------------------------------ Geometry Calculation（几何参数计算与输出）
    skew = pchip(XR0,skew0,RC);                % 使用分段三次插值将掠角 skew0 从 XR0 插值到控制点 RC
    rake = pchip(XR0,rake0,RC);                % 使用分段三次插值将倾角 rake0 从 XR0 插值到 RC
    fid = fopen('Geometry.txt','w');           % 创建并打开几何信息输出文件 Geometry.txt
    fprintf(fid,'\t\t\tGeometry.txt\n');       % 写入文件标题
    fprintf(fid,'\t\tPropeller Geometry Table\n\n');    % 写入表格标题
    fprintf(fid,'Propeller Diameter = %.1f m\n',D);      % 写入螺旋桨直径
    fprintf(fid,'Number of Blades = %.0f\n',NBLADE);     % 写入叶片数
    fprintf(fid,'Propeller Speed= %.0f RPM\n',N);        % 写入转速
    fprintf(fid,'Propeller Hub Diameter = %.2f m\n',Dhub);   % 写入桨毂直径
    if Meanline==1
        fprintf(fid,'Meanline Type: NACA a=0.8\n');      % 写入平均线类型为 NACA
    elseif Meanline==2
        fprintf(fid,'Meanline Type: Parabolic\n');       % 写入平均线类型为抛物线
    end
    if Thickness==1
        fprintf(fid,'Thickness Type: NACA 65A010\n\n');  % 写入厚度分布类型为 NACA 65A010
    elseif Thickness==2
        fprintf(fid,'Thickness Type: Elliptical\n\n');   % 写入厚度分布类型为椭圆型
    elseif Thickness==3
        fprintf(fid,'Thickness Type: Parabolic\n\n');    % 写入厚度分布类型为抛物线
    end
    fprintf(fid,' r/R\t P/D\t Skew\t Xs/D\t  c/D\t  f0/c\t  t0/c\n');  % 写入列标题
    for i = 1:MT
        ThetaP(i) = TANBIC(i) + AlphaI;                      % 计算螺距角（诱导攻角 + 理想攻角）//////////////////////////
        PitchOD(i) = tand(ThetaP(i))*pi*RC(i);              % 根据螺距角计算 P/D（螺距比）///RC(i) 是该点的半径比///////////////////////////////
        fprintf(fid, '%.3f\t %.2f\t %.1f\t %.3f\t %.3f\t %.4f\t %.4f\n'...
        ,RC(i),PitchOD(i),skew(i),rake(i),CDC(i),f0oc(i),t0oc(i));    % 写入每个剖面的几何参数
    end
    fclose(fid);                                             % 关闭文件

    % --------------------------------------------------------- BASIC SHAPE
       % ---------------------------- 【基本剖面坐标初始化与几何准备阶段】----------------------------
    % 本段代码用于根据控制点计算螺旋桨各剖面的弦长、位置、厚度和弯度线坐标
    % 包括平面剖面基本形状、选定的平均线类型（NACA 或抛物线）、厚度分布类型（NACA 或椭圆/抛物线）

    c = CDC.*D;                              % 每个控制点的实际弦长 c = c/D × D
    r = RC.*R;                               % 每个控制点的实际半径位置 r = r/R × R
    theta = 0:360/NBLADE:360;               % 每片叶片的极角分布，用于后续复制剖面(圆周上的旋转角度,由叶片数决定)

    for i = 1:MT                             % 遍历每个控制点（MT 个半径分布点）
        for j = 1:NP                         % 对每个剖面，分成 NP 个点（用于构建剖面曲线）
            x1(i,j) = c(i)/2 - c(i)/(NP-1)*(j-1);   % 沿弦线的 x 坐标，从后缘到前缘（右→左）
            station(1,j) = 1/(NP-1)*(j-1);          % 归一化的剖面位置（0 ~ 1）用于插值
            z1(i,j,1) = sqrt(r(i)^2 - x1(i,j)^2);   % 计算 z 坐标（剖面在半径圆上的投影）
        end
    end
% ---------------------------------------------------- MEANLINE & THICKNESS
       % ======================== 构建剖面中弯线（Meanline）和厚度分布（Thickness） ========================
    % 此段用于根据用户选择的剖面类型（NACA或抛物线），生成每个剖面处的中弯线坐标 f(i,j)、斜率 dfdx，以及厚度分布 t(i,j)

    x = [0 .5 .75 1.25 2.5 5 7.5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 ...
         80 85 90 95 100]./100;   % 定义归一化弦长坐标（0 到 1，共 26 点）

    if Meanline==1          % 若选择 NACA a=0.8 的中弯线
        foc = [0 .287 .404 .616 1.077 1.841 2.483 3.043 3.985 4.748 5.367 5.863 6.248...
            6.528 6.709 6.79 6.77 6.644 6.405 6.037 5.514 4.771 3.683 2.435 1.163 0]./100;  % NACA a=0.8 的标准中弯线坐标（归一化）
        fscale = f0oc./max(foc);        % 缩放比例：将标准中弯线按每个剖面设定的 f0/c 进行放大
        dfdx0 = [.48535 .44925 .40359 .34104 .27718 .23868 .21050 .16892...
            .13734 .11101 .08775 .06634 .04601 .02613 .00620 -.01433 -.03611...
            -.06010 -.08790 -.12311 -.18412 -.23921 -.25583 -.24904 -.20385];  % 中弯线斜率（df/dx）预设数据
        for i = 1:MT                         % 遍历所有控制点（每个剖面）
            for j = 1:NP                     % 遍历每个剖面上的离散点（沿弦方向）
                f(i,:) = pchip(x,foc.*fscale(i).*c(i),station);       % 插值得到每个剖面的中弯线坐标 f(i,j)
                dfdx(i,:) = pchip(x(2:end),dfdx0,station);            % 插值得到每个剖面斜率 dfdx(i,j)
            end
        end
    elseif Meanline==2      % 若选择抛物线形中弯线（Parabolic）
        for i = 1:MT
            for j = 1:NP
                f(i,j) = f0oc(i)*c(i)*(1-(2*x1(i,j)/c(i))^2);         % 抛物线中弯线公式
                dfdx(i,j) = -8*f0oc(i)*x1(i,j)/c(i);                  % 中弯线斜率的导数
            end
        end
    end

    if Thickness==1         % 若选择 NACA 65A010 厚度分布
        toc_65 = [0 .765 .928 1.183 1.623 2.182 2.65 3.04 3.658 4.127 4.483 4.742 4.912...
            4.995 4.983 4.863 4.632 4.304 3.899 3.432 2.912 2.352 1.771 1.188 .604 .021]./100;  % 标准 NACA 65A010 厚度分布
        tscale = t0oc./2./max(toc_65);     % 厚度缩放因子（注意 NACA 是上下对称的，所以除以2）
        for i = 1:MT
            for j = 1:NP
                t(i,:) = pchip(x,toc_65.*tscale(i).*c(i),station);   % 插值得到每个剖面的厚度分布
            end
        end
    elseif Thickness==2     % 若选择椭圆厚度分布
        for i = 1:MT
            for j = 1:NP
                t(i,j) = t0oc(i)*c(i)*real(sqrt(1-(2*x1(i,j)/c(i))^2));   % 椭圆厚度分布公式
            end
        end
    elseif Thickness==3     % 若选择抛物线厚度分布
        for i = 1:MT
            for j = 1:NP
                t(i,j) = t0oc(i)*c(i)*(1-(2*x1(i,j)/c(i))^2);             % 抛物线厚度分布公式
            end
        end
    end

    % -------------------------------------------------- CAMBER & THICKNESS
      % ========== 剖面计算：根据厚度与弯度线计算翼型上下表面的坐标 ==========
    for i = 1:MT                    % 对每个径向控制点（MT 个站位）
        for j = 1:NP                % 对每个弦向分布点（NP 个剖面点）
            xu(i,j) = x1(i,j) + t(i,j) * sin(atan(dfdx(i,j)));     % 上表面 X 坐标：弯度线 + 厚度偏移
            xl(i,j) = x1(i,j) - t(i,j) * sin(atan(dfdx(i,j)));     % 下表面 X 坐标：弯度线 - 厚度偏移
            yu(i,j) = f(i,j) + t(i,j) * cos(atan(dfdx(i,j)));      % 上表面 Y 坐标：弯度线 + 厚度偏移
            yl(i,j) = f(i,j) - t(i,j) * cos(atan(dfdx(i,j)));      % 下表面 Y 坐标：弯度线 - 厚度偏移
            if (isreal(yu(i,j))==0)||(isreal(yl(i,j))==0)          % 如果计算结果出现虚数，说明几何出错
                % if Single_Flag==1                                  % 如果是在“单螺旋桨设计模式”
                %     set(Err_Single,'visible','on','enable','on','string','Error in MPVL.');  % 显示错误提示
                %     if ishandle(Fig1_S)~=0                         % 如果图形窗口存在
                %         close(Fig1_S);                             % 关闭报错前的图像窗口
                %     end
                %     return;                                        % 中止程序执行
                % end
            end
        end
    end

    % -------------------------------------------------- PITCH, SKEW & RAKE
       % -------------------------------------------------- PITCH, SKEW & RAKE
    % 本段功能：基于螺旋桨的螺距角（Pitch）、掠角（Skew）、倾角（Rake）参数，将每个剖面点变换到最终三维坐标系中，构造完整的螺旋桨三维几何模型（含所有叶片）
     % 螺距角：	决定每个剖面相对于旋转轴的“扭转角”
    % 掠角：改变翼型在轴向上的位置（偏斜）
    % 倾角：模拟叶片在轴向的后移
    yrake = rake.*D;  % 计算每个半径点的倾角偏移量（沿 y 方向），rake 是非维量（Xs/D），乘直径 D 得到实际偏移量///////////////////////////////////////////
    for i = 1:MT  % 遍历所有半径位置（控制点）
        for j = 1:NP  % 遍历每个剖面的弦向点
            xup(i,j,1) = xu(i,j)*cosd(ThetaP(i)) - yu(i,j)*sind(ThetaP(i));  % 上翼型绕 pitch 角旋转（x）
            xlp(i,j,1) = xl(i,j)*cosd(ThetaP(i)) - yl(i,j)*sind(ThetaP(i));  % 下翼型绕 pitch 角旋转（x）
            yup(i,j,1) = xu(i,j)*sind(ThetaP(i)) + yu(i,j)*cosd(ThetaP(i));  % 上翼型绕 pitch 角旋转（y）
            ylp(i,j,1) = xl(i,j)*sind(ThetaP(i)) + yl(i,j)*cosd(ThetaP(i));  % 下翼型绕 pitch 角旋转（y）
            z1p(i,j,1) = z1(i,j);  % pitch 变化不影响 z 坐标，原样保留

            xus(i,j,1) = xup(i,j,1)*cosd(skew(i)) - z1p(i,j,1)*sind(skew(i));  % 上翼型绕掠角 skew 旋转（x）
            xls(i,j,1) = xlp(i,j,1)*cosd(skew(i)) - z1p(i,j,1)*sind(skew(i));  % 下翼型绕掠角 skew 旋转（x）
            yus(i,j,1) = yup(i,j,1);  % 掠角 skew 不影响 y 坐标（轴向方向）
            yls(i,j,1) = ylp(i,j,1);  % 同上
            z1s(i,j,1) = x1(i,j)*sind(skew(i)) + z1p(i,j,1)*cosd(skew(i));  % 上/下翼型的 z 坐标偏移（由 skew 引起）

            xur(i,j,1) = xus(i,j,1);  % 初始化最终上翼型坐标（x）
            xlr(i,j,1) = xls(i,j,1);  % 初始化最终下翼型坐标（x）
            yur(i,j,1) = yus(i,j,1) - yrake(i);  % 上翼型坐标在 y 方向施加 rake 偏移（轴向后移）
            ylr(i,j,1) = yls(i,j,1) - yrake(i);  % 下翼型同理
            z1r(i,j,1) = z1s(i,j,1);  % 保持最终 z 坐标
            for k = 2:length(theta)-1  % 为每个叶片复制坐标并绕圆周分布
                xur(i,j,k) = xur(i,j,1)*cosd(theta(k)) - z1r(i,j,1)*sind(theta(k));  % 上翼型绕轴对称分布（x）
                xlr(i,j,k) = xlr(i,j,1)*cosd(theta(k)) - z1r(i,j,1)*sind(theta(k));  % 下翼型同理
                yur(i,j,k) = yur(i,j,1);  % y 保持不变
                ylr(i,j,k) = ylr(i,j,1);  % 同上
                z1r(i,j,k) = xur(i,j,1)*sind(theta(k)) + z1r(i,j,1)*cosd(theta(k));  % 上下翼型 z 坐标绕轴旋转（得到完整 3D 坐标）
            end
        end
    end


    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% ... 原有代码：生成 Geometry.txt 文件 ...

% ==================== 新增代码：计算螺旋桨展开面积比和盘面积 ====================
% % ==================== 计算螺旋桨展开面积比和盘面积 ====================
% % 计算盘面积 A_O
% A_O = pi * R^2;
% 
% % 实际半径和实际弦长数组
% r_actual = RC * R;       % RC为各控制点的r/R，R为螺旋桨半径
% c_actual = CDC * D;       % CDC为各控制点的c/D，D为螺旋桨直径
% 
% % 使用梯形法则计算单叶展开面积
% A_E_single = 0;
% for i = 1:length(r_actual)-1
%     dr = r_actual(i+1) - r_actual(i);           % 半径间隔
%     c_avg = (c_actual(i) + c_actual(i+1)) / 2;   % 平均弦长
%     A_E_single = A_E_single + c_avg * dr;         % 累加面积微元
% end
% 
% % 总展开面积（所有叶片）
% A_E = A_E_single * NBLADE;
% 
% % 计算展开面积比
% AE_AO_ratio = A_E / A_O;
% 
% % % 将结果写入文件
% % fprintf(fid, '\n===== 展开面积与盘面积 =====\n');
% % fprintf(fid, '盘面积 A_O = %.4f m²\n', A_O);
% % fprintf(fid, '展开面积 A_E = %.4f m²\n', A_E);
% % fprintf(fid, '展开面积比 A_E/A_O = %.4f\n', AE_AO_ratio);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    % -------------------------- Create Figure for 2D Propeller Blade Image
       % -------------------------- 绘制螺旋桨剖面的二维视图（绘制二维的螺旋桨叶片剖面图像（横截面），展示不同半径位置（r/R）的剖面外形）
    Fig2_S = figure('units','normalized','position',[0.31 .06 .4 .3],'name',...
        'Blade Image','numbertitle','off');    % 创建一个新图窗显示二维剖面图
    style=['r' 'g' 'b' 'm' 'k'];               % 设置线条颜色样式（红绿蓝紫黑）
    str_prefix = {'r/R = '};                  % 图例前缀文字（代表半径比）
    flag=1;                                   % 用于切换颜色与图例编号
    for i = 1:ceil(MT/5):MT                    % 每隔约5个控制点画一个剖面（简洁可视）
        plot(xur(i,:,1),yur(i,:,1),style(flag));     % 绘制该剖面的上表面
        str_legend(flag)=strcat(str_prefix,num2str(RC(i)));  % 生成图例文本
        hold on;                              % 保持图像，叠加多个剖面图
        flag = flag+1;                        % 切换颜色与图例索引
    end
    flag=1;                                   % 重置颜色计数器
    for i = 1:ceil(MT/5):MT                    % 同样再绘制一次下表面
        plot(xlr(i,:,1),ylr(i,:,1),style(flag));     % 绘制该剖面的下表面
        hold on;                              % 保持图像继续叠加
        flag = flag+1;                        % 颜色索引递增
    end
    legend(str_legend,'location','northwest');       % 添加图例并放在左上角
    axis equal;     grid on;                        % 设置等比例坐标轴和网格
    xlabel('X (m)'); ylabel('Y (m)');               % 添加坐标轴标签
    hold off;                                       % 释放图像叠加锁

    % ------------------------------------------- Create 3D Propeller Image
      % ---------------------------- 创建 3D 螺旋桨图像，用于显示完整的螺旋桨外形
    Fig3_S = figure('units','normalized','position',[.61 .06 .4 .3],...
           'name','Propeller Image','numbertitle','off');   % 新建一个图窗，用于绘制3D螺旋桨

    for k = 1:NBLADE
        surf(xur(:,:,k),yur(:,:,k),z1r(:,:,k));         hold on;   % 绘制第 k 个叶片的上表面
        surf(xlr(:,:,k),ylr(:,:,k),z1r(:,:,k));         hold on;   % 绘制第 k 个叶片的下表面
    end

    tick = 0:15:90;                                     % 定义圆柱 Hub 的轮廓角度步进
    [xh0,yh0,zh0] = cylinder(Rhub*sind(tick),50);       % 生成一个弯曲的底部过渡段（小圆弧）
    surf(yh0,zh0.*.3+min(ylp(1,:,1))-.3,xh0);            % 显示底部 Hub 的曲面过渡段

    [xh1,yh1,zh1] = cylinder(Rhub,50);                  % 生成中心 Hub 的圆柱部分
    surf(yh1, zh1+min(ylp(1,:,1)), xh1);                % 显示 Hub 主体圆柱表面

    hold off;       colormap gray;          grid on;       
    axis equal;   % 关闭叠加，设置灰色色图、网格和坐标比例


    xlabel('X');    ylabel('Y');            zlabel('Z');                  % 设置坐标轴标签
    % set(New_Single,'visible','on','enable','on');                        % 显示“Try Again”按钮
    % figure(Fig_Main);                                                    % 激活主界面窗口
    % toc                                                                  % 停止计时器，显示计算耗时

    
    % ------------------ 输出最终推进效率（使用 KTRY 对应值） ------------------
eta_final = EFFY(KTRY);  % 选定解对应的推进效率

fprintf('\n========== 最终推进效率 ==========\n');
fprintf('最终选定的推进效率 η = %.4f (取自第 %d 次迭代结果)\n', eta_final, KTRY);
fprintf('===================================\n\n');

  % =================================================== Do not delete