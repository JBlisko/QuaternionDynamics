clear all
close all
clc

%% Initialize quaternion
beta = pi/4; % Pitch
gamm = pi/4; % Yaw
dt = 0.1;
vec = [1 0 0]';

quat_tot = find_q_init(beta,gamm);
q_init = quat_tot;
quat_tot_exp = quat_tot;

rot_vec_0 = quat_rot_vec(vec,quat_tot); % Rotate initial orientation (assumed to parallel to x-axis)
rot_vec_bod_ode_0 = rot_vec_0;
rot_vec_bod_exp_0 = rot_vec_0;

%% Plotting axes
figure
hold on
%Initial vector
quiver3( 0,0,0, rot_vec_0(1), rot_vec_0(2), rot_vec_0(3), 'b' )
% x-axis
quiver3( 0,0,0, 1,0,0,'r' ,'LineWidth',2);
% y-axis
quiver3( 0,0,0, 0,1,0,'g' ,'LineWidth',2);
% z-axis
quiver3( 0,0,0, 0,0,1,'b' ,'LineWidth',2);

%% Solving quaternion dynamics

for n=1:11
    % Choose axis of rotation
    %omega = (1/dt)*[0 .5 0]'; % Use large rotations to test consistency between lab- and body-frame formulations
    omega = [0.1 0.1 0.1]'; % Use small rotation to test consistency between ode and exp methods
    
    % Lab frame formulation - system of ODEs
    quat_tot = quat_ode(omega,dt,quat_tot); % ODE method
    rot_vec_lab_ode = quat_rot_vec(vec,quat_tot);
    
    % Lab frame formulation - exponential approximation
    quat_tot_exp = quat_exp_approx(omega,dt,quat_tot_exp);
    rot_vec_lab_exp = quat_rot_vec(vec,quat_tot_exp);
    
    % Body frame formulation - system of ODEs
    quat_bod_ode = quat_ode_body(omega,dt);
    rot_vec_bod_ode = quat_rot_vec(rot_vec_bod_ode_0,quat_bod_ode);
    rot_vec_bod_ode_0 = rot_vec_bod_ode;
    
    % Body frame formulation - exponential approximation
    new_quat_bod_exp = quat_exp_approx_body(omega,dt);
    rot_vec_bod_exp = quat_rot_vec(rot_vec_bod_exp_0,new_quat_bod_exp);
    rot_vec_bod_exp_0 = rot_vec_bod_exp;
    
    % Plot using labe frame method
    quiver3( 0,0,0, rot_vec_lab_ode(1), rot_vec_lab_ode(2), rot_vec_lab_ode(3), 'r' );
    
    % Plot using body frame method
    quiver3( 0,0,0, rot_vec_bod_ode(1), rot_vec_bod_ode(2), rot_vec_bod_ode(3), 'k' );
end
grid on
axis equal
axis([-1  1.0 -1 1.0 -1 1.0])
xlabel( 'x' )
ylabel( 'y' )
zlabel( 'z' )
view( 135, 20 )
hold off

%% Functions
function [qinit] = find_q_init(beta,gamma)
    vec=[1 0 0]'; % WLOG we assume the initial orientation is along the x-axis
    rho = [cos(gamma)*cos(beta); sin(gamma)*cos(beta); -sin(beta)];
    
    qinit_0 = [1+dot(vec,rho);cross(vec,rho)];
    qinit = qinit_0/norm(qinit_0);
end

function [rot_vec] = quat_rot_vec(vec,quat)
    s = quat(1);
    r = quat(2:4);
    
    rot_vec = vec + 2*cross(r, s*vec + cross(r,vec));
end

% Note: This has to be in world frame formulation of quaternion dynamics.
% Otherwise, the order of quaternions is opposite and the dynamics becomes
% incorrect.

% Lab frame formulation - system of ODEs
function [new_quat] = quat_ode(omega,dt,quat)
    syms q0(t) q1(t) q2(t) q3(t)
    
    % Establish system of ODEs
    A = [0 -omega(1) -omega(2) -omega(3);...
        omega(1) 0 -omega(3) omega(2);...
        omega(2) omega(3) 0 -omega(1);...
        omega(3) -omega(2) omega(1) 0];
    Q = [q0; q1; q2; q3;];
    odes = diff(Q) == (1/2)*A*Q;
    
    % Solve system of equations
    C = Q(0) == quat;
    [q0Sol(t),q1Sol(t),q2Sol(t),q3Sol(t)] = dsolve(odes,C);
    
    % Evaluate new quaternion with timestep
    t=dt; % Timestep
    new_quat_0 = eval([q0Sol(t);q1Sol(t);q2Sol(t);q3Sol(t)]);
    new_quat = new_quat_0/norm(new_quat_0);
end

% Body frame formulation - system of ODEs
function [new_quat] = quat_ode_body(omega,dt)
    syms q0(t) q1(t) q2(t) q3(t)
    quat = [1; 0; 0; 0];
    
    % Establish system of ODEs
    A = [0 -omega(1) -omega(2) -omega(3);...
         omega(1) 0 omega(3) -omega(2);...
         omega(2) -omega(3) 0 omega(1);...
         omega(3) omega(2) -omega(1) 0];
    Q = [q0; q1; q2; q3;];
    odes = diff(Q) == (1/2)*A*Q;
    
    % Solve system of equations
    C = Q(0) == quat;
    [q0Sol(t),q1Sol(t),q2Sol(t),q3Sol(t)] = dsolve(odes,C);
    
    % Evaluate new quaternion with timestep
    t=dt; % Timestep
    new_quat = eval([q0Sol(t);q1Sol(t);q2Sol(t);q3Sol(t)]);
end

% Lab frame formulation - exponential approximation
function [new_quat] = quat_exp_approx(omega,dt,quat)
    A = [0 -omega(1) -omega(2) -omega(3);...
         omega(1) 0 -omega(3) omega(2);...
         omega(2) omega(3) 0 -omega(1);...
         omega(3) -omega(2) omega(1) 0];
    
    new_quat_0 = (eye(4)+(dt/2)*A)*quat;
    new_quat = new_quat_0/norm(new_quat_0);
end

% Body frame formulation - exponential approximation
function [new_quat] = quat_exp_approx_body(omega,dt)
    new_quat = [1;omega*dt/2];
end

% Multiply quaternions q1 and q2
function [q_mul] = mult_quat(q1,q2)
    q_mul_0 = [q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4);...
             q1(1)*q2(2) + q1(2)*q2(1) - q1(3)*q2(4) + q1(4)*q2(3);...
             q1(1)*q2(3) + q1(2)*q2(4) + q1(3)*q2(1) - q1(4)*q2(2);...
             q1(1)*q2(4) - q1(2)*q2(3) + q1(3)*q2(2) + q1(4)*q2(1);];
    
    q_mul = q_mul_0/norm(q_mul_0);
end
