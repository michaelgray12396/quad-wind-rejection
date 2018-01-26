function [t, o, theta, odes] = simulate
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
    t1 = 10;
    % Tuning parameters
    Q = 10000*eye(12);
    R = eye(4);
    o_des = [0;0;-1]; % hover position
    wind = [1;0;0]; % rudimentary wind
    % Derive gain
    xe = [o_des;0;0;0;0;0;0;0;0;0];
    ue = [sqrt(m*g/(4*kF));sqrt(m*g/(4*kF));sqrt(m*g/(4*kF));sqrt(m*g/(4*kF))];
    [Ac,Bc] = control_design_init(J,m,g,kF,kM,l);
    
    % Create variables to keep track of time, state, input, and desired position
    ti = t0;
    xi = [o0; theta0; v0; w0];
    x_log = []; % does this have a purpose?

    % Simulate
    for i = 1:(t1/dt)
        tic
        % decide what traj point to use (find minimum to the time vector)
        % in traj, increase resolution of interpolation
        % make it solve faster, check with the one on the github, probably
        % turn up some tolerances
        % add logic such that after the trajectory, it goes back to normal
        % control
        % recreate an automation script
        % add wind
        % check what wind does to something not following the trajectory
        % check what wind does when perfectly considered within the planner
        % in a non-traj flight, estimate wind
        % check this estimate
        % find a way to pass that estimate into the traj planner
        % fix the visualize
        % put this all on the github
        % write up a thing for the meeting
        % hit up david again, ask about blade flapping and the drag model
        % also ask him what he thinks in general
        % consider that the drag model can be experimentally driven
        
        K = control_design_loop(Ac,Bc,Q,R,dt,xi(:,i)-xe,ue);
        u_desired = -K*(xi(:,i) - xe) + ue;
        ui = GetBoundedInputs(u_desired,sigmamax);
        u(:, i+1) = ui; % for the output log
        [t, x] = ode45(@(t, x) h(t, x, ui, g, m, J, kF, kM, l, wind), [ti(i) ti(i)+dt], xi(:,i));
        ti(i+1) = t(end);
        x = x';
        xi(:,i+1) = x(:,end);
        loop_time(i+1) = toc;
    end

    % package variables
    odes = zeros(size(xi));
    t = ti;
    o = xi(1:3,:);
    theta = xi(4:6,:);
    
    % plots
    figure(1) % plot positions
    plot(ti,o(1,:),ti,o(2,:),ti,o(3,:))
    xlabel('Time, s');
    ylabel('Position, m');
    legend('x','y','z')
    figure(2) % plot inputs
    plot(ti,u(1,:),ti,u(2,:),ti,u(3,:),ti,u(4,:))
    xlabel('Time, s')
    ylabel('Input, rad/s');
    legend('u1','u2','u3','u4')
    figure(3) % plot computation time
    plot(ti,loop_time)
    xlabel('Time, s')
    ylabel('Computation Time, s');
    
    % save all the variables in the workspace
    filename = 'data.mat';
    save(filename)
    
    % animate the flight
    % visualize(t, o, theta, odes, 'movie.mp4')
end

function xdot = h(~, x, u, g, m, J, kF, kM, l, wind)
    theta = x(4:6); %assigning variables for simplicity
    v = x(7:9);
    w = x(10:12);
    q = 0.1924; %rho * cd * SA * (1/2)
    
    drag = q*(wind - v).^2;
    f0 = GetR_ZYX(theta)*[0;0;-kF*(u(1)^2 + u(2)^2 + u(3)^2 + u(4)^2)] + ...
        [0;0;m*g] + drag; % forces in the inertial frame
    T = [l*kF*(u(1)^2 - u(2)^2); ...
         l*kF*(u(3)^2 - u(4)^2); ...
         kM*(u(1)^2 + u(2)^2 - u(3)^2 - u(4)^2)]; %torques in the body frame

    xdot = [v; % computes the motion derivatives
            GetThetaDot_ZYX(theta, w);...
            f0/m;
            J\(T-cross(w,J*w))];
end

function u = GetBoundedInputs(u_desired,sigmamax)
    u = zeros(length(u_desired),1); % Preallocate
    for i=1:length(u) % Filters each spinrate input to max and min values
        if u_desired(i) > sigmamax
            u(i) = sigmamax;
            disp('maxed')
        elseif u_desired(i) < 0
            u(i) = 0;
            disp('minned')
        else
            u(i) = u_desired(i);
        end
    end
end

function R = GetR_ZYX(theta)
    c1 = cos(theta(1));
    s1 = sin(theta(1));
    c2 = cos(theta(2));
    s2 = sin(theta(2));
    c3 = cos(theta(3));
    s3 = sin(theta(3));
    R = [c1*c2 c1*s2*s3-s1*c3 c1*s2*c3+s1*s3;
         s1*c2 s1*s2*s3+c1*c3 s1*s2*c3-c1*s3;
         -s2 c2*s3 c2*c3];
end

function thetadot = GetThetaDot_ZYX(theta, w)
    c2 = cos(theta(2));
    s2 = sin(theta(2));
    c3 = cos(theta(3));
    s3 = sin(theta(3));
    thetadot = (1/c2)*[0 s3 c3; 0 c2*c3 -c2*s3; c2 s2*s3 s2*c3]*w;
end

function [Ac, Bc] = control_design_init(J,m,g,kF,kM,l)
    syms o1 o2 o3 t1 t2 t3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4
    o = [o1;o2;o3];
    t = [t1;t2;t3];
    v = [v1;v2;v3];
    w = [w1;w2;w3];
    state = [o;t;v;w];
    input = [u1;u2;u3;u4];
    R = GetR_ZYX(t);
    what = [0 -w3 w2;w3 0 -w1;-w2 w1 0];
    f0 = R*[0;0;-kF*(u1^2 + u2^2 + u3^2 + u4^2)] + [0;0;m*g];
    T = [l*kF*(u1^2 - u2^2); ...
         l*kF*(u3^2 - u4^2); ...
         kM*(u1^2 + u2^2 - u3^2 - u4^2)];
    eoms = [v; GetThetaDot_ZYX(t,w); f0/m; J\(T-what*J*w)];
    Ac = matlabFunction(jacobian(eoms,state),'Vars',[t1;t2;t3;w1;w2;w3;u1;u2;u3;u4]);
    Bc = matlabFunction(jacobian(eoms,input),'Vars',[t1;t2;t3;w1;w2;w3;u1;u2;u3;u4]);
end

function K = control_design_loop(Ac,Bc,Qd,Rd,dt,x,u)
    [Ad,Bd] = ssdata(c2d(ss(Ac(x(4),x(5),x(6),x(10),x(11),x(12),u(1),u(2),u(3),u(4)), ...
        Bc(x(4),x(5),x(6),x(10),x(11),x(12),u(1),u(2),u(3),u(4)),[],[]),dt));
    K = dlqr(Ad,Bd,Qd,Rd);
end