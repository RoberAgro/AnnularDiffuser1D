function [m,U,geometry,Cp,other] = computation_diffuser(AR,phi_1,phi_2,R,H,v_m,v_t,d,p,Cf,fluid)
%% Input parameters
phi = (phi_2 + phi_1)/2;            % Mean wall cant angle
div = (phi_2 - phi_1)/2;            % Wall divergence semi-angle
b_in = H/cos(phi);
r_in = R;
z_in = 0;
s = refpropm('s','p',p,'d',d,fluid);

%% Solution of the ODE system
% Define the initial conditions
U0 = [v_m v_t d p s];

% Integrate the ode system using ode45
% Use m and U and variables and the rest of inputs as extra parameters
% Integrate between 0 and infinity (the AR is the stopping criterion)
% Use the options RelTol and AbsTol to set the integration tolerance
options = odeset('Events',@(m,U)area_ratio(m,U,AR,phi,div,r_in,b_in),'RelTol',1e-6,'AbsTol',1e-6);
[m,U] = ode45(@(m,U)ode_diffuser(m,U,phi,div,r_in,b_in,z_in,Cf,fluid),[0,inf],U0,options);

% Rename the solution to the other calculations
v_m = U(:,1);
v_t = U(:,2);
d = U(:,3);
p = U(:,4);


%% Geometry
% Retrieve geometry
r = r_fun(r_in,phi,m);
z = z_fun(z_in,phi,m);
b = b_fun(b_in,div,m);
AR = (b.*r)/(b_in*r_in);

% Define the inner and outer surface vectors
z_outer = z - b/2.*sin(phi);
z_inner = z + b/2.*sin(phi);
r_outer = r + b/2.*cos(phi);
r_inner = r - b/2.*cos(phi);

% Store geometry
geometry.r = r;
geometry.z = z;
geometry.b = b;
geometry.AR = AR;
geometry.z_outer = z_outer;
geometry.z_inner = z_inner;
geometry.r_outer = r_outer;
geometry.r_inner = r_inner;


%% Compute the pressure recovery factor
a = v_t(1)/v_m(1);  % a = tan(alpha_in)
v_in = sqrt(v_m(1)^2+v_t(1)^2);
p_in = p(1);
d_in = d(1);
s_in = refpropm('s','p',p_in,'d',d_in,fluid);
h_in = refpropm('h','p',p_in,'d',d_in,fluid);
h0_in = h_in + v_in^2/2;
p0_in = refpropm('p','h',h0_in,'s',s_in,fluid);

% Different definitions
p = U(:,4);
Cp = (p-p_in)/(p0_in-p_in);
Cp_incompressible = (p-p_in)/(1/2*d_in*v_in^2)*1000;
Cp_ideal = 1 - (r_in./r).^2.*((b_in./b).^2 + a^2)/(1 + a^2);
other.Cp_incompressible =  Cp_incompressible;
other.Cp_ideal =  Cp_ideal;


end


function dUdm = ode_diffuser(m,U,phi,div,r_in,b_in,z_in,Cf,fluid)

% Rename variables
v_m = U(1);
v_t = U(2);
d   = U(3);
p   = U(4);
alpha = atan(v_t/v_m);
v = sqrt(v_m^2+v_t^2);

% Increment for finite differences
delta = 1e-5;

% Local geometry
r = r_fun(r_in,phi,m);              % Radius as a function of m
z = z_fun(z_in,phi,m);              % Axial distance as a function of m
b = b_fun(b_in,div,m);              % Channel width as a function of m

% Derivative of the area change (forward finite differences)
diff_br = (b_fun(b_in,div,m+delta)*r_fun(r_in,phi,m+delta) - b*r)/delta;

% Derivative of internal energy with respect to pressure (constant density)
e_1 = refpropm('u','p',p-delta,'d',d,fluid);
e_2 = refpropm('u','p',p+delta,'d',d,fluid);
diff_e = (e_2-e_1)/(2*delta)/1000;   % Convert from J/kPa to J/Pa

% % Ideal gas limit (check)
% k = refpropm('k','d',d,'p',p,fluid);
% Diff_e = 1/(d*(k-1))

% Speed of sound (avoid computations in the two phase region
a = refpropm('a','p',p,'d',d,fluid);

% Stress at the wall
tau_w = Cf*d*v^2/2;        % Skin friction coefficient

% Heat flux at the wall
q_w = 0;                   % Adiabatic wall

% Coefficient matrix A (pressure conversion from kPa to Pa)
A = [d           0          v_m       0;
     d*v_m       0            0       1*1000;
     0       d*v_m            0       0;
     0           0   -d*v_m*a^2   d*v_m*1000];
 
% Source term vector
S = zeros(4,1);
S(1) = -d*v_m/(b*r)*diff_br;
S(2) = +d*v_t*v_t/r*sin(phi) - 2*tau_w/b*cos(alpha);
S(3) = -d*v_t*v_m/r*sin(phi) - 2*tau_w/b*sin(alpha);
S(4) = 2*(tau_w*v + q_w)/b/diff_e;

% Obtain the slope of the solution by Gaussian elimination
dUdm = A\S;

% Check entropy generation
T = refpropm('T','d',d,'p',p,fluid);
sigma = 2/b*(tau_w*v);
dUdm(5) = sigma/(d*v_m)/T;      % ds/dm

end


function r = r_fun(r_in,phi,m)
r = r_in + sin(phi)*m;
end


function z = z_fun(z_in,phi,m)
z = z_in + cos(phi)*m;
end


function b = b_fun(b_in,div,m)
b = b_in + 2*tan(div)*m;
end


function [AR_check,isterminal,direction] = area_ratio(m,~,AR_prescribed,phi,div,r_in,b_in)

% Geometry
r = r_fun(r_in,phi,m);              % Radius as a function of m
b = b_fun(b_in,div,m);              % Channel width as a function of m
AR_current = (b*r)/(b_in*r_in);     % Current area ratio

% Stopping criterion
AR_check = AR_prescribed - AR_current;
isterminal = 1;   % stop the integration
direction = 0;    % negative direction

end



%% Extensions of the code to make it more general:

% Use general functions for the geometry as input (instead of linear funcs)
% If the geometry is provided as general functions the local angle phi has
% to be computed by differentiation (finite differences)

% If the parametrization variable is not the meridional coordinate (m) then
% it is necessary to integrate the arclength to compute m

% Provide an arbitrary variation of the skin friction coefficient or use an
% empirical correlation to compute it

% Implement the heat transfer model in the code
% Perhaps use the Chilton-Colburn analogy to compute the heat transfer
% coefficient. It is also necessary to compute the stagnation temperature
% and to prescribe a wall temperature distribution

