function [h_end,x,y] = backshoot_eigenfunction_meltdifflat(sig,k,xx,bb,NN,qq,hh,xmin)

% we are going towards the terminus, so xmax is the terminus and the sign of q must be positive
%ub = 10^-4*ones(size(xx));
ub = zeros(size(xx));

xmax = max(xx);
%xmin = 0.99*xmax;
options = odeset('RelTol',1e-11,'AbsTol',1e-12);

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
drag = 20;
lr = 2; % lateral spacing of bumps, m

Mbb = rhow*g/L*(1+om*abs(qq)/nu).^2./(1+2*om*abs(qq)/nu).*qq.^2*3./bb.^4*(12*nu)/g*(1 + ct*cw*rhow)+rhow*g/L*(1+om*abs(qq)/nu)./(1 + 2*om*abs(qq)/nu).*qq*3./bb*ct*cw*rhow*alpha;
Mhh = rhow*g/L*(2+3*om*abs(qq)/nu)./(1+2*om*abs(qq)/nu).*qq*(1 + ct*cw*rhow)-rhow*g/L*ct*cw*rhow*bb.^3*g/12/nu./(1+2*om*abs(qq)/nu)*alpha;
Qbb = (1 + om*abs(qq)/nu)./(1 +2*om*abs(qq)/nu)*3./bb.*qq;
Qhh = bb.^3*g./(12*nu*(1 +2*om*abs(qq)/nu));
KK = bb.^3*g./(12*nu*(1 +om*abs(qq)/nu));
UU = rhow*g/L*ub.^2*drag^2;
mbar = rhow/L*(G + 12*nu*(1+ om*abs(qq)/nu).*qq.^2);
DD = bb.*mbar/rhoi;

N = griddedInterpolant(xx,NN);
b = griddedInterpolant(xx,bb);
q = griddedInterpolant(xx,qq);
hbar = griddedInterpolant(xx,hh);
Mb = griddedInterpolant(xx,Mbb);
Mh = griddedInterpolant(xx,Mhh);
Qb = griddedInterpolant(xx,Qbb);
Qh = griddedInterpolant(xx,Qhh);
K = griddedInterpolant(xx,KK);
U = griddedInterpolant(xx,UU);
D = griddedInterpolant(xx,DD);
ubb = griddedInterpolant(xx,ub);

sig0 = @(x) Mb(x)/rhoi-A*N(x).^n-D(x)*k^2;%-ubb(x)/lr;
overlhs = @(x) Mh(x)./(Qb(x).*Mh(x) + rhoi*(Qh(x).*(sig-sig0(x))));
fun1 = @(x) Qb(x) + rhoi*Qh(x)./Mh(x).*(sig-sig0(x));
dfun1 = griddedInterpolant(xx,gradient(fun1(xx),xx));
fun2 = @(x) rhoi*Qh(x)./Mh(x)*A*n.*N(x).^(n-1)*rhow*g.*b(x)-U(x)/rhoi;
dfun2 = griddedInterpolant(0:10^4,gradient(fun2(0:10^4),0:10^4));

[x,y] = ode45(@bvpfcn,[xmin xmax],[-10^-12; 0],options);

h_end = y(end,1)/max(abs(y(:,1)));

subplot(2,1,1)
plot(x,y(:,1)/y(end,2))
subplot(2,1,2)
plot(x,y(:,2)/y(end,2))



function dydx = bvpfcn(x,y)
    h = y(1);
    bhat = y(2);
    dhdx = ((A*n*N(x)^(n-1)*rhow*g*b(x)-U(x)/rhoi)*h - (sig-sig0(x))*bhat)/(Mh(x)/rhoi);
    rhs = -sig*bhat-K(x)*k^2*h-D(x)/rhow*rhoi*k^2*bhat+Mb(x)/rhow*bhat - U(x)/rhow*h - Mh(x)/rhow*dhdx -dfun1(x)*bhat + fun2(x)*dhdx + dfun2(x)*h;
    dydx = [dhdx
        rhs.*overlhs(x)];
end

end