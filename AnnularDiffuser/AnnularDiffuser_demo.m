%% AnnularDiffuser demo
% Author: Roberto Agromayor
% Date: 25/08/2018
clear all                               % Clear worskpace variables
close all                               % Close previous figures
clc                                     % Clear command window
addpath([pwd,'/FUNCTIONS'])             % Add path to the functions folder
plot_settings                           % Configure plot settings


%% Define input parameters
% Thermodynamic properties
fluid = 'air.mix';                      % Working fluid
T = 273.15+20;                          % Inlet static temperature
p = 101.235;                            % Inlet static pressure
d = refpropm('d','T',T,'p',p,fluid);    % Inlet density
a = refpropm('a','T',T,'p',p,fluid);    % Inlet speed of sound

% Friction factor
Cf = 0.010;                             % Mean friction coefficient

% Geometry (axial inlet)
R = 1.00;                               % Turbomachinery outlet radius
x = 0.70;                               % Turbomachinery outlet hub-to-tip ratio
H = 2*R*(1-x)/(1+x);                    % Turbomachinery outlet blade height

% Wall angles
phi = 30*pi/180;                        % Mean wall cant angle
div = 5*pi/180;                         % Divergence semi-angle
phi_1 = phi - div;                      % Inner cant angle
phi_2 = phi + div;                      % Outer cant angle

% Area ratio
AR = 5.00;                              % Area ratio

% Velocity vector
Ma_m = 0.30;                            % Inlet meridional Mach number
alpha = 30*pi/180;                      % Inlet flow angle
v_m = Ma_m*a;                           % Inlet meridional velocity
v_t = v_m*tan(alpha);                   % Inlet tangential velocity


%% Solve the flow in the diffuser
% Call the main function: computation_diffuser.m
[m,U,geometry,Cp] = AnnularDiffuser(AR,phi_1,phi_2,R,H,v_m,v_t,d,p,Cf,fluid);

% Solution along the diffuser
v_m = U(:,1);
v_t = U(:,2);
v = sqrt(v_m.^2+v_t.^2);
d = U(:,3);
p = U(:,4);

% Compute the enthalpy along the diffuser
h = zeros(length(m),1);
for i = 1:length(m)
    h(i) = refpropm('h','p',p(i),'d',d(i),fluid);
end
h_0 = h + v.^2/2;

% Retrieve the area ratio vector
AR_vec = geometry.AR;


%% Plot the pressure recovery coefficient distribution
figure1 = figure(1); ax_fig1 = gca;
hold on; axis square; box on
xlabel({' ';'$AR$ -- Area ratio '});
ylabel({'$C_{p}$ -- Pressure recovery coefficient';' '});
ax_fig1.XAxis.TickLabelFormat = '%.0f';
ax_fig1.YAxis.TickLabelFormat = '%.2f';
ax_fig1.XTick = 1:1:10;
ax_fig1.YTick = 0.00:0.20:1.00;
axis([1 AR 0 1])
plot(AR_vec,Cp,'k')


%% Plot the velocity evolution
figure2 = figure(2); ax_fig2 = gca;
hold on; axis square; box on
xlabel({' ';'$AR$ -- Area ratio '});
ylabel({'Velocity distribution (m/s)';' '});
ax_fig2.XAxis.TickLabelFormat = '%.0f';
ax_fig2.YAxis.TickLabelFormat = '%.0f';
ax_fig2.XTick = 1:1:10;
ax_fig2.YTick = 0.00:30:120;
axis([1 AR 0 120])

plot(AR_vec,v,'k')
plot(AR_vec,v_m,'b')
plot(AR_vec,v_t,'r')
legend('Velocity magnitude -- $v$','Meridional component -- $v_m$','Tangential component -- $v_{\theta}$','Location','Northeast')


%% Plot energy distribution
figure3 = figure(3); ax_fig3 = gca;
hold on; axis square; box on
xlabel({' ';'$AR$ -- Area ratio '});
ylabel({'Energy distribution (J/kg)';' '});
ax_fig3.XAxis.TickLabelFormat = '%.0f';
ax_fig3.YAxis.TickLabelFormat = '%.0f';
ax_fig3.XTick = 1:1:10;
ax_fig3.YTick = 0:1000:6000;
axis([1 AR 0 6000])

plot(AR_vec,h_0-h(1),'k')
plot(AR_vec,h-h(1),'b')
plot(AR_vec,v.^2/2,'r')
% plot(AR_vec,v_t,'r')
legend('Stagnation enthalpy -- $h_0$','Static enthalpy -- $h$','Kinetic energy -- $v^2/2$','Location','East')


%% Plot the diffuser geometry
figure4 = figure(4); ax_fig4 = gca;
hold on; axis image; box on
xlabel({' ';'$z$ -- Axial direction '});
ylabel({'$r$ -- Radial direction';' '});
ax_fig4.XAxis.TickLabelFormat = '%.0f';
ax_fig4.YAxis.TickLabelFormat = '%.0f';
ax_fig4.XTick = -10:1:10;
ax_fig4.YTick = -10:1:10;
axis([-0.5 3.5 0 4])
title({' '; '\textbf{Geometry of the diffuser}';' '})

z_mean = geometry.z;
r_mean = geometry.r;
z_inner = geometry.z_inner;
z_outer = geometry.z_outer;
r_inner = geometry.r_inner;
r_outer = geometry.r_outer;
z = [z_inner; z_outer(end:-1:1); z_inner(1)];
r = [r_inner; r_outer(end:-1:1); r_inner(1)];

plot(z,r,'k')
plot(z_mean,r_mean,'k--')