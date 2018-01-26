clc; clear;
addpath OptimTraj/
% Parameters
g = 9.81; % acceleration of gravity
m = .7146; % mass
J = diag([0.0019, 0.0019, 0.008]); % moment of inertia matrix in body frame
l = .15; % spar length
kF = 4.46e-06; % aerodynamic force coefficient
kM = 1.217e-7; % aerodynamic torque coefficient
sigmamax = 1e3; % maximum spin rate
% Initial time
t0 = 0;
% Initial state
o0 = [0; 0; 0];
theta0 = [0; 0; 0];
v0 = [0; 0; 0];
w0 = [0; 0; 0];
% Sample time
dt = (1/50);
% Final time
tf = 5;
o_des = [0;0;-1]; % hover position
physics_disturbance = [0;0;0]; % rudimentary wind

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up function handles                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Put symbolic EOMs in convenient form
syms o1 o2 o3 t1 t2 t3 v1 v2 v3 w1 w2 w3 ...
    u1 u2 u3 u4 real
theta = [t1; t2; t3]; %assigning variables for simplicity
v = [v1;v2;v3];
w = [w1;w2;w3];

% GetR
    c1 = cos(theta(1));
    s1 = sin(theta(1));
    c2 = cos(theta(2));
    s2 = sin(theta(2));
    c3 = cos(theta(3));
    s3 = sin(theta(3));
    R = [c1*c2 c1*s2*s3-s1*c3 c1*s2*c3+s1*s3;
         s1*c2 s1*s2*s3+c1*c3 s1*s2*c3-c1*s3;
         -s2 c2*s3 c2*c3];
% GetThetaDot
    thetadot = (1/c2)*[0 s3 c3; 0 c2*c3 -c2*s3; c2 s2*s3 s2*c3]*w;
f0 = R*[0;0;-kF*(u1^2 + u2^2 + u3^2 + u4^2)] + ...
    [0;0;m*g] + physics_disturbance; % forces in the inertial frame
T = [l*kF*(u1^2 - u2^2); ...
     l*kF*(u3^2 - u4^2); ...
     kM*(u1^2 + u2^2 - u3^2 - u4^2)]; %torques in the body frame

xdot = [v; % computes the motion derivatives
        thetadot;...
        f0/m;
        J\(T-cross(w,J*w))];

% Convert EOMs from symbolic to numeric
numf = matlabFunction(simplify(xdot),'vars',[o1 o2 o3 ...
                                                          t1 t2 t3 ...
                                                          v1 v2 v3 ...
                                                          w1 w2 w3...
                                                          u1 u2 u3 u4]);
                                                      
problem.func.dynamics = @(t,x,u)( numf(x(1,:), x(2,:), x(3,:),...
                                       x(4,:), x(5,:), x(6,:),...
                                       x(7,:), x(8,:), x(9,:),...
                                       x(10,:), x(11,:), x(12,:),...
                                       u(1,:), u(2,:), u(3,:), u(4,:) ));
                                   
                                   
                                   
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%          
%                           Optimize for:                                     %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


% Time
%   problem.func.pathObj = @(t,x,u)( sum((x(1,:)-x_f).^2,1) + sum((x(3,:)-z_f).^2,1))
% problem.func.bndObj = @(t0,x0,tF,xF)(tF-t0);
problem.func.pathObj = @(t,x,u)(t);
%problem.func.pathObj = @(t,x,u)( sum((x(1,:) - o_des(1)).^2,1) + ...
%    sum((x(2,:) - o_des(2)).^2,1) + sum((x(3,:) - o_des(3)).^2,1) );


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                     Set up problem bounds                               %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.bounds.initialTime.low = 0;
problem.bounds.initialTime.upp = 0;
problem.bounds.finalTime.low = 0;
problem.bounds.finalTime.upp = tf;

problem.bounds.initialState.low = [o0; v0; theta0; w0];
problem.bounds.initialState.upp = [o0; v0; theta0; w0];
problem.bounds.finalState.low = [o_des; 0; 0; 0;...
                                    0; 0; 0; 0; 0; 0];
problem.bounds.finalState.upp = [o_des; 0; 0; 0;...
                                    0; 0; 0; 0; 0; 0];

problem.bounds.control.low = [0; 0; 0; 0];
problem.bounds.control.upp = [sigmamax; sigmamax; sigmamax; sigmamax];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                    Initial guess at trajectory                          %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.guess.time = [0,tf];
problem.guess.state = [problem.bounds.initialState.low, problem.bounds.finalState.low];
problem.guess.control = [[m*g/4; m*g/4; m*g/4; m*g/4], [m*g/4; m*g/4; m*g/4; m*g/4]];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Solver options                                  %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',1e6);
problem.options.MaxFunctionEvaluations = 100000;
% problem.options.method = 'trapezoid';
% problem.options.method = 'hermiteSimpson';
% problem.options.method = 'rungeKutta';
% problem.options.method = 'chebyshev';

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                            Solve!                                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

soln = optimTraj(problem);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                        Display Solution                                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%% Unpack the simulation
t = linspace(soln.grid.time(1), soln.grid.time(end), soln.grid.time(end)*10*(1/dt));
x = soln.interp.state(t);
u = soln.interp.control(t);

clf;
figure(1)
subplot(2,1,1);
plot(t,x(1,:),t,x(2,:),t,x(3,:));
legend('x','y','z');
grid on
subplot(2,1,2)
plot(t,x(4,:),t,x(5,:),t,x(6,:));
legend('xdot','ydot','zdot');
grid on


figure(2)
subplot(2,1,1)
plot(t,x(7,:),t,x(8,:),t,x(9,:));
legend('\theta','\phi','\psi')
grid on
subplot(2,1,2);
plot(t,x(10,:),t,x(11,:),t,x(12,:));
legend('\thetadot','\phidot','\psidot')
grid on

figure(3)
plot(t,u);
legend('f1','f2','f3','f4');
grid on

save('traj.mat','t','x','u');