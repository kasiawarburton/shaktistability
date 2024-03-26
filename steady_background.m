function [x,hh,qq,hhx,bb,mm,NN,dzb] = steady_background(x_mesh,z_b,z_s,u_b,ieb)

% Define constants
rhow=1000; % Density of water (kg/m3)
rhoi=917; % Bulk density of ice (kg/m3)
g=9.81; % Gravitational acceleration (m/s2)
nu=1.787e-6; % Kinematic viscosity of water at 0 degC (m2/s)
n=3; % Rheological flow law exponent for ice
A=2.5e-25; %(Pa-3 s-1)
G=0.05; % Geothermal flux (W/m2)
om=0.001; % Parameter controlling transition to nonlinear resistance in basal system
L=334e3; % Latent heat of fusion J/kg
ct=7.5e-8; %Clapeyron slope (K/Pa)
cw=4.22e3; %specific heat capacity of water (J/kg/K)
drag = 20;
br=0.0; % Typical bedrock bump height (m)
lr=2.0; % Bedrock bump spacing (m)

ub = griddedInterpolant(x_mesh,u_b);
zb = griddedInterpolant(x_mesh,z_b);
dzbdx = griddedInterpolant(x_mesh,gradient(z_b,x_mesh));
pri = griddedInterpolant(x_mesh,rhoi*g.*(z_s-z_b));
input = griddedInterpolant(x_mesh,ieb);

xmax = max(x_mesh);
options = odeset('RelTol',1e-11,'AbsTol',1e-13);

    function q_end = shoot(q_start)
        [xx,y] = ode45(@bvpfcn,[0 xmax],[zb(0); q_start],options);
        q_end = y(end,2);
    end

% change bounds depending on ieb and size of domain
q_0 = fzero(@shoot,[10^-7, 10^-1],optimset('TolX',10^-11));
[x,yy] = ode45(@bvpfcn,[0 xmax],[zb(0); q_0],options);
hh = yy(:,1);
qq = yy(:,2);
hhx = hx(yy(:,1),x,yy(:,2));
bb = 1./invb(qq,hhx);
mm = 1./invm1(hh,x,qq,hhx);
NN = (pri(x)-pw(hh,x));
dzb = dzbdx(x);

if false
figure(3);
subplot(4,1,1);plot((x(end)-x)/10^3,hh,'linewidth',2);ylabel('h (m)','Interpreter','latex');hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'XTickLabel',[])
subplot(4,1,2);plot((x(end)-x)/10^3,qq,'linewidth',2);ylabel('q (m$^2$/s)','Interpreter','latex');hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'XTickLabel',[])
subplot(4,1,3);plot((x(end)-x)/10^3,bb,'linewidth',2);ylabel('b (m)','Interpreter','latex');hold on
set(gca,'TickLabelInterpreter','latex')
set(gca,'XTickLabel',[])
subplot(4,1,4);plot((x(end)-x)/10^3,NN/1e6,'linewidth',2);ylabel('N (MPa)','Interpreter','latex');hold on;
set(gca,'TickLabelInterpreter','latex')
xlabel('Distance (km)','Interpreter','latex')
end

function dydx = bvpfcn(x,y)
    h = y(1);
    q = y(2);
    if min(q) < 0
        dydx = [0
            -0.00000001];
    else
        dhdx = hx(h,x,q);
        m = 1./invm2(h,x,q,dhdx);
        dydx = [dhdx
             -m/rhoi - input(x)];
    end
end

    function PW = pw(h,x)
        PW = rhow*g*(h-zb(x));
    end

    function B = invb(q,dhdx)
        if q == 0
            B = 3;
        else
        B = abs(g*dhdx./(12*nu*(1+om/nu*abs(q)).*q)).^(1/3);
        end
    end

    function betaval = beta(q,dhdx)
        b = 1./invb(q,dhdx);
        betaval = zeros(size(b));
        betaval(b<br) = (br-b(b<br))./lr;
    end

    function M = invm1(h,x,q,dhdx)
        M = invb(q,dhdx)./(A*(pri(x)-pw(h,x)).^n - beta(q,dhdx).*ub(x).*invb(q,dhdx))/rhoi;
    end

    function M = invm2(h,x,q,dhdx)
        u = ub(x);
        taub=drag^2.*abs(pri(x)-pw(h,x)).*u;
        M = L./(G + u.*taub + abs(rhow*g*q.*dhdx)); % no PMP term
    end

    function HX = hx(h,x,q)
        HX = h;
        for i = 1: length(h)
            if q(i) == 0
                HX(i) = 0;
            else
                dHmax = 50000000;
                HX(i) = sign(q(i)).*abs(fzero(@(dhdx) invm1(h(i),x(i),q(i),dhdx) - invm2(h(i),x(i),q(i),dhdx),[0,dHmax*0.99]));
            end
        end
    end

end