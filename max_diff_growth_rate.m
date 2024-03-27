function [sig,lmax] = max_diff_growth_rate(ice_thickness, meltinput)

% semi-analytical growth rates for fig 6

[x,~,qq,~,bb,~,NN,~] = steady_background(0:1000,-0.02*(1000:-1:0) + 2000, -0.02*(1000:-1:0) + 2000+ice_thickness,10^-6*ones(size(0:1000)),meltinput/24/365/3600*ones(size(0:1000)));
xx = 10^3-x;

% Define constants
rhow=1000; % Density of water (kg/m3)
rhoi=917; % Bulk density of ice (kg/m3)
g=9.81; % Gravitational acceleration (m/s2)
ct = 0*7.5e-8; %Clapeyron slope, K/Pa
cw = 4.22e3; %specific heat capacity of water, J/kg/K
nu=1.787e-6; % Kinematic viscosity of water at 0 degC (m2/s)
n=3; % Rheological flow law exponent for ice
A=2.5e-25; %(Pa-3 s-1)
L=334e3; % Latent heat of fusion J/kg
om = 0.001; % 1/transition to turbulence Re
G = 0.05;
alpha = - 0.02;

Mbb = rhow*g/L*(1+om*abs(qq)/nu).^2./(1+2*om*abs(qq)/nu).*qq.^2*3./bb.^4*(12*nu)/g*(1 + ct*cw*rhow)+rhow*g/L*(1+om*abs(qq)/nu)./(1 + 2*om*abs(qq)/nu).*qq*3./bb*ct*cw*rhow*alpha;
Mhh = rhow*g/L*(2+3*om*abs(qq)/nu)./(1+2*om*abs(qq)/nu).*qq*(1 + ct*cw*rhow)-rhow*g/L*ct*cw*rhow*bb.^3*g/12/nu./(1+2*om*abs(qq)/nu)*alpha;
Qbb = (1 + om*abs(qq)/nu)./(1 +2*om*abs(qq)/nu)*3./bb.*qq;
KK = bb.^3*g./(12*nu*(1 +om*abs(qq)/nu));
mbar = rhow/L*(G + 12*nu*(1+ om*abs(qq)/nu).*qq.^2);
DD = bb.*mbar/rhoi;

sig0 = Mbb/rhoi - A*NN.^n;
sig00 = max(sig0);
sigx = max(gradient(sig0,xx));
sig = sig00 - 4*DD(1)^(1/4)*(1.0187/3)^(3/4)*sigx^(1/2)*(Qbb(1)*Mhh(1)/rhoi/KK(1)).^(1/4);

lmax = 2*pi*(sig00/1.0187)^(3/2)*(rhoi*KK(1)/(Qbb(1)*Mhh(1)))^(1/2)/sigx;

end