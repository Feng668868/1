% MICHELL Michell integral wave drag computation
%
%MICHELL(Y,U,L,B,T,RHO,N)computes Rw for the non-dimensional
%hull offsets Y at speed U in water of density RHO. The offsets
%are scaled bylengthL,beam Band draft T.
%
%The integration is carried out over Nz waterlines,Nz stations
% and N propagation angles.
%
%The matrix Y should be Nx by Nz ordered from bow to stern and
%keel to waterline.The maximum value ofthe offset matrixY
%must be 1/2.Nx must be odd.All values must be specified
%in metric units.
%
%Doug Read 30.7.2008

function Rw = michell(Y,U,L,B,T,RHO,N)
%定义主函数 michell，用于计算船体波阻力 Rw，输入包括
            %Y：船体非量纲偏移矩阵（Nx × Nz），从船艏到船艉，从龙骨到水线
            %U：航速
            %L, B, T：船长、型宽、吃水
            %RHO：水密度
            %N：传播角度数量


Nx= size(Y,1);                    %determine number of stations
Nz= size(Y,2);                    %determine number of waterlines
YH = Y*B;                         %scale non-dimensional offsets by the beam
%获取站点数（Nx）和水线数（Nz），并将非量纲偏移 Y 按宽度 B 还原为真实单位

if mod(Nx,2)== 0;
    warning('Nx must be odd.');      %required for x Filon algorithm
end
%如果站点数 Nx 为偶数，则警告（因为 Filon 积分算法要求奇数节点）

%------ integration variables---------
%Z方向积分准备（从龙骨到水线）
dz=T/(Nz-1); z=-T:dz:0;      z = z';
dx=L/(Nx-1); x= 0:dx:L;      x = x';
 %theta = linspace(0,pi/2,N);   theta = theta';
theta = michspace(N);        theta = theta' ;
%构造 Z（吃水方向）和 X（船长方向）上的均匀网格
%theta 是传播角度，用 michspace(N) 生成对 π/2 附近加密的对数分布


  g = 9.80665;
  k0 = g/U^2;                         % fundamental wave number;基础重力波数
  c = (4*RHO*U^2)/pi;                 % constantc ;c是波阻计算中的常系数
  a = sec(theta);                     % convenient substitution
  k = k0*a.^2;                        % dispersion relation;与角度有关的波数

   %---------- Z INTEGTRAL---------------
   %Z方向积分（Filon变换，用指数核函数计算积分）

  %--- variables for Filon trapezoidal algorithm---
  Kz = k0*dz*a.^2;
  w0 = (exp(Kz)-1-Kz)./Kz.^2;
  wn = (exp(Kz)+exp(-Kz)-2)./Kz.^2;
  wN = (exp(-Kz)-1+Kz)./Kz.^2;
  %用 Filon 算法的核函数权重（适用于指数函数）

%--- preallocate--
f = zeros(Nz,1);
F = zeros(Nx,N);
%初始化中间变量 F，储存在 X 方向每一点、每个传播角 theta(j) 下的 Z 积分结果

for j = 1:N;
  for m = 1:Nx;
    for n = 1:Nz;
     if n == 1;
       f(n) = w0(j)*YH(m,n)*exp(k0*z(n)*a(j)^2)*dz;
     elseif n == Nz;
       f(n) = wN(j)*YH(m,n)*exp(k0*z(n)*a(j)^2)*dz;
     else
       f(n) = wn(j)*YH(m,n)*exp(k0*z(n)*a(j)^2)*dz;
     end
    end
    F(m,j) = sum(f); 
  end
end
%按照指数加权积分方法，对每个 theta(j)、每个 X 截面 m，沿吃水方向 Z 积分

 %---------- X INTEGRAL----------------
 %X方向积分（使用 Filon 算法处理余弦和正弦项）

 %--- variables for Filon algorithm--
Kx = k0*dx.*a;
alp = (Kx.^2+1/2*Kx.*sin(2*Kx)+cos(2*Kx)-1)./Kx.^3;
bet = (3*Kx+Kx.*cos(2*Kx)-2*sin(2*Kx))./Kx.^3;
gam = 4*(sin(Kx)-Kx.*cos(Kx))./Kx.^3;
%Filon 积分用到的权重系数 alpha, beta, gamma（分别对应端点项、中间偶数项、奇数项）

Nev = (Nx+1)/2;      % even Filon index
Nod = (Nx-1)/2;      % odd Filon index

 %--- preallocate---
pev = zeros(Nev,1); qev = zeros(Nev,1);
pod = zeros(Nod,1); qod = zeros(Nod,1);
Pt = zeros(N,1);    Qt = zeros(N,1);
Pev = zeros(N,1);   Qev = zeros(N,1);
Pod = zeros(N,1);   Qod = zeros(N,1);
P = zeros(N,1);     Q = zeros(N,1);
%初始化积分结果变量，分成奇偶索引与终点项进行分配

for j = 1:N;
    for m = 1:Nev;
        pev(2*m-1) = F(2*m-1,j)*cos(k0*x(2*m-1)*a(j));
        qev(2*m-1) = F(2*m-1,j)*sin(k0*x(2*m-1)*a(j));
    end
    for m = 1:Nod;
        pod(2*m) = F(2*m,j)*cos(k0*x(2*m)*a(j));
        qod(2*m) = F(2*m,j)*sin(k0*x(2*m)*a(j));
    end

    Pt(j) = F(Nx,j)*cos(k0*L*a(j));  
    Qt(j) = F(Nx,j)*sin(k0*L*a(j));

    Pev(j) = sum(pev)-1/2*Pt(j);
    Pod(j) = sum(pod);
    Qev(j) = sum(qev)-1/2*Qt(j);
    Qod(j) = sum(qod);

    P(j)=dx*( alp(j)*Qt(j)+bet(j)*Pev(j)+gam(j)*Pod(j));
    Q(j)=dx*(-alp(j)*Pt(j)+bet(j)*Qev(j)+gam(j)*Qod(j));
end
%对每个角度 theta(j)，分别完成 X 方向积分，得到 P 和 Q 两个向量（与余弦、正弦部分相关）

%波阻抗谱计算与传播角积分:
 R = c*k.^2./a.^3.*...
     (k.^2.*(P.^2 + Q.^2)+...
      2*k.*a.*(Q.*Pt- P.*Qt )+...
       a.^2.*( Pt.^2 + Qt.^2 ));
  R(isnan(R)) = 0;
  %计算在每个传播角 theta 下的波阻力谱 R(θ)，最后处理 NaN（如出现除 0）

 %---------- THETA INTEGRAL------------

 rw = zeros(N-1,1);

 for k = 1:N-1;
     rw(k) = 1/2*(R(k)+R(k+1))*(theta(k+1)-theta(k));
 end 

Rw = sum(rw);
%使用复化梯形积分法对 theta 方向进行数值积分，得到总波阻力 Rw

 %------------ END---------------------

 function [xm]=michspace(N);

  % MICHSPACE log spacing for Michell integral
  %
  % MICHSPACE(N) produces log base 10 spacing over N propagation
  % angles between 0 and pi/2. Points are more closely spaced
  % near pi/2.

   xm = logspace(0,1,N)-1;   %michspace(N) 创建 N 个传播角，以 π/2 附近为密集区域的 log 分布,返回的是 [0, π/2] 区间内加密采样的角度数组
   xm = xm*pi/18-pi/2;
   xm = fliplr(-xm);


