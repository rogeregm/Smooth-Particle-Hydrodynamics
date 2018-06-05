%% Droplet Reference Solution
close all
clear all
clc
%Constants
g = -1*9.81;
gamma = 0.0728;
rho = 1000;
capil = 13.45%(g*rho)/gamma;
curv = 0.734%/0.1299;

%Domain
domain = linspace(0,1,10000);    
    
%Reset initial slope
dydx_0 = [curv;0;0;]% cos(curv*1*10^-3); sin(curv*1*10^-3)];
options = odeset('InitialSlope',dydx_0,'InitialStep',1*10^-3,'NonNegative',1,'NonNegative',3);
[s,Y] = ode15s(@(s,Y) odedrop(s,Y,curv,capil),domain,[curv;1;0],options);
    
%plotting
figure
plot(Y(:,2),Y(:,3),'.');
hold on
axis([0.995 1.005 0 0.01]);



%% With ODE45

g = 9.81;
gamma = 0.0728;
rho = 1000;

figure('Name','Drop Reference Solution')

dydx_0 = [curv; cos(curv*1*10^-3); sin(curv*1*10^-3)];
for nn=1:7
    n=1*10^nn;
    
    %Domain
%     domain = linspace(0,0.01,1000);
    domain = logspace(0,0.1,1000);

    %Physical Constants
    gravity =n*(1*10^-6)*g
    capil = (n*(1*10^-6)*g*rho)/gamma;
    curv = 1/0.1299;
    
    %Reset initial slope
    options = odeset('InitialSlope',dydx_0,'InitialStep',1*10^-3);
    [s,Y] = ode15s(@(s,Y) odedrop(s,Y,curv,capil),domain,[curv;1;0],options);
    
    %plotting
%     hold on
    figure
    plot(Y(:,2),Y(:,3));
end






%% Sparavigna
clear
curv = 0.734;
capil = 13.45;
%Initial Conditions @ s=0
theta_0 = 0;
r_0 = 0.734;
z_0 = 0;
x_0 = 0;
s_0 = 2*curv + capil*z_0 -r_0;

%Reference point
A_ref = 0.5*pi*0.1^2;
V_0 = pi*x_0^2*sin(theta_0);
V_ref=0.0089;

%Stepsize
ds = 0.001;

%
theta(1) = theta_0 + s_0*ds;
z(1) = sin(theta(1))*s_0;
x(1) = cos(theta(1))*s_0;
s(1) = 2*curv + capil*(z_0 + z(1)) -r_0;
V(1) = pi*x(1)^2*sin(theta(1));
n=2;
while V<V_ref
    theta(n) = theta_0 + s(n-1)*ds/2;
    x(n) = cos(theta(n))*ds/2;
    z(n) = sin(theta(n))*ds/2;
    s(n) = 2*curv + capil*(z(n-1)+z(n)) -r_0;%sin(theta(n))/(x_0 + x(n));
    V(n)= V(n-1)+pi*x(n)^2*sin(theta(n))
    
    n=n+1;
end

plot(x,z)



%%
clear
curv = -0.5;
capil = 13.45;
%Initial Conditions @ s=0
theta(1) = 0;
z(1) = 0;
x(1) = 0;
s(1) = curv;

domain = linspace(0,0.1,1000);
ds=1/length(domain);
V_ref=1;
% %
theta(2) = theta(1) + s(1)*ds
z(2) = sin(theta(2))*s(1)
x(2) = cos(theta(2))*s(1)
s(2) = 2*curv + capil*(z(1) + z(2)) -sin(theta(1) + theta(2))/(x(1) + x(2))

for n=3:length(domain)

    k1 = ds*[2*curv+capil* z(n-1)-sin(theta(n-1))/(x(n-1))                          ;cos(theta(n-1))         ;sin(theta(n-1))];
    k2 = ds*[2*curv+capil*(z(n-1)+k1(3)/2)-sin(theta(n-1)+k1(1)/2)/(x(n-1)+k1(2)/2) ;cos(theta(n-1)+k1(1)/2) ;sin(theta(n-1)+k1(1)/2)];
    k3 = ds*[2*curv+capil*(z(n-1)+k2(3)/2)-sin(theta(n-1)+k2(1)/2)/(x(n-1)+k2(2)/2) ;cos(theta(n-1)+k2(1)/2) ;sin(theta(n-1)+k2(1)/2)];
    k4 = ds*[2*curv+capil*(z(n-1)+k3(3)/2)-sin(theta(n-1)+k3(1))/(x(n-1)+k3(2))     ;cos(theta(n-1)+k3(1))   ;sin(theta(n-1)+k3(1))];
    
    theta(n) = theta(n-1) + (k1(1) + 2*k2(1) + 2*k3(1) + k4(1))/6;
    x(n) = x(n-1) + (k1(2) + 2* k2(2) + 2*k3(2) + k4(2))/6;
    z(n) = z(n-1) + (k1(3) + 2*-k2(3) + 2*k3(3) + k4(3))/6;
    V(n)= pi*x(n)^2*sin(theta(n))
    if V(n)>V_ref
        break
    end
end

% figure('Name','Sessile Drop')
hold on
plot(x,z,'.')


