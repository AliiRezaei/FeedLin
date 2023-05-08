%% Project Info

% Feedback Linearization for 2DOF robotic arm
% *** Computed Torque Control ***
% author : Ali Rezaei

clc
clear
close all

%% Manipulator Parameters

m1 = 1;    % mass of link 1
m2 = 1;    % mass of link 2
l1 = .6;    % length of link 1
l2 = .6;    % length of link 2
lc1 = l1/2; % length of center of mass link 1 
lc2 = l2/2; % length of center of mass link 2
I1 = 1/12 * m1 * l1^2;  % inertia moment link 1
I2 = 1/12 * m2 * l2^2;  % inertia moment link 2
gravity = 9.81;

%% Initial Values

q = [0, 0]';     % angels
w = [0, 0]';     % velocities
tau = [0, 0]';   % torques

%% Desired Values

q_des = [pi/3, pi/4]';  % desired angels
w_des = [0, 0]';    % desired velocities
alpha_des = [0, 0]';     % desired accelerations

%% PD Controller Parameters

kp = 16 * eye(2);     % Proportional  
kd = 8 * eye(2);       % Derivative

%% Dynamics and Control

dt = 0.01; % Time Step
SimTime = 5; % Simulation Time (s)
tspan = 0:dt:SimTime; % Time Span for solve eq.
MaxIter = numel(tspan); % Maximum iteration

Q = zeros(2,MaxIter);   % for save angels
Qdot = zeros(2,MaxIter);    % for save velocities
E = zeros(2,MaxIter);    % for save pos error
dE = zeros(2,MaxIter);    % for save vel error
Tau = zeros(2,MaxIter);   % for save torques

figure % for plot Robot animation      
xlabel('X')
ylabel('Y')
    axis(1.5 * [-1 1 -1 1])
    grid minor
qdes_deg = q_des * 180 / pi;
TitleInfo = sprintf('Robot Motion \n q_{1_{des}} = %3.2f , q_{2_{des}} = %3.2f',qdes_deg(1),qdes_deg(2));
title(TitleInfo)

% ----------------------Solve eq. of Motion-----------------------

for i = 1:MaxIter
    
    % ------------------------H(q) Matrix-------------------------
    
    H11 = m1 * lc1^2 + I1 + m2 * (l1^2 + ...
        lc2^2 + 2 * l1 * lc2 * cos(q(2))) + I2;
    H12 = m2 * l1 * lc2 *cos(q(2)) + m2 * lc2^2 + I2;
    H21 = H12;
    H22 = m2 * lc2^2 + I2;
    H = [H11 H12; H21 H22];
     
    % ------------------------g(q) Vector-------------------------
    
    g1 = m1 * lc1 * gravity * cos(q(1)) + m2 * gravity * ...
        (lc2 * cos(q(1) + q(2)) + l1 * cos(q(1)));
    g2 = m2 * lc2 * gravity * cos(q(1) + q(2));
    g = [g1; g2];
    
    % ----------------------C(q,qdot) Matrix----------------------
    
    h = m2 * l1 * lc2 * sin(q(2));
    C = [-h * w(2), -h * (w(1) + w(2)); h * w(2), 0];
    
    % --------------------Control Manipulator--------------------
    
    alpha = H^(-1) * (tau - C * w - g);
    w = alpha * dt + w;
    q = w * dt + q;
    
    e = (q_des - q); % pos error
    de = (w_des - w); % vel error

    V = alpha_des - kd * de - kp * e;
    tau = - H * V + C * w + g;
    
    % ------------------Save variables for plot------------------
    
    Q(:,i) = q;
    Qdot(:,i) = w;
    E(:,i) = e;
    dE(:,i) = de;
    Tau(1:2,i) = tau;
    
    % ---------------------Show Manipulator----------------------
    
    axis([-1.5 1.5 -1.5 1.5])
    link1 = line([0, l1 * cos(q(1))],[0, l1 * sin(q(1))]);
    set(link1,'LineWidth',3,'Marker','.','MarkerSize',30)
    link2 = line([l1 * cos(q(1)), l1 * cos(q(1)) + l2 * cos(q(1) + q(2))]...
        ,[l1 * sin(q(1)), l1 * sin(q(1)) + l2 * sin(q(1) + q(2))]);
    set(link2,'LineWidth',3,'Marker','.','MarkerSize',30)
    pause(dt)
    if i ~= MaxIter
        delete(link1);
        delete(link2);
    end
    hold on
    plot(l1 * cos(q(1)) + l2 * cos(q(1) + q(2)), l1 * sin(q(1)) + l2 * sin(q(1) + q(2))...
        ,'.r','MarkerSize',3)
    hold off
    
end

%% Plots

% -------------------Plot Positions---------------------

figure
plot(tspan,Q(1,:),'LineWidth',1.5)
hold on
grid on
plot(tspan,Q(2,:),'LineWidth',1.5)
xlabel('Time')
ylabel('q(t)')
title('Positions')
legend('q_1','q_2')

% --------------------Plot Velocities---------------------

figure
plot(tspan,Qdot(1,:),'LineWidth',1.5)
hold on
grid on
plot(tspan,Qdot(2,:),'LineWidth',1.5)
xlabel('Time')
ylabel('\omega(t)')
title('Angular Velocites')
legend('\omega_1','\omega_2')

% --------------------Plot Position Error---------------------

figure
plot(tspan,E(1,:),'LineWidth',1.5)
hold on
grid on
plot(tspan,E(2,:),'LineWidth',1.5)
xlabel('Time')
ylabel('e(t)')
title('Position Error')
legend('q_1 error','q_2 error')

% --------------------Plot Velocity Error---------------------

figure
plot(tspan,dE(1,:),'LineWidth',1.5)
hold on
grid on
plot(tspan,dE(2,:),'LineWidth',1.5)
xlabel('Time')
ylabel('de(t)/dt')
title('Velocity Error')
legend('w_1 error','w_2 error')

% --------------------Plot Velocity Error---------------------

figure
subplot(211)
plot(tspan,Tau(1,:),'LineWidth',1.5)
ylabel('\tau_1')
grid on
subplot(212)
plot(tspan,Tau(2,:),'LineWidth',1.5)
xlabel('Time')
ylabel('\tau_2')
grid on
suptitle('Motor Torques (N.m)')

%% Response Information

ResponseInfo.joint1 = stepinfo(Q(1,:),tspan);
ResponseInfo.joint2 = stepinfo(Q(2,:),tspan);

%% References

% Slotin, J.-J. E and Li W., applied non-linear control systems 
% Mark W. Spong, Seth Hutchinson, M. Vidyasagar, Robot Modeling and Control
