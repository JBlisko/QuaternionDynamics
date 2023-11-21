clear all
close all
clc

%% Matrix Rotation
axlen=1.25; %Axes length
origin=[0 0 0]';

%Angles
a=0; %Yaw
b=pi/2; %Pitch
c=0.3; %Roll

%Rotation
R=[cos(a)*cos(b) cos(a)*sin(b)*sin(c)-sin(a)*cos(c) cos(a)*sin(b)*cos(c)+sin(a)*sin(c);...
    sin(a)*cos(b) sin(a)*sin(b)*sin(c)+cos(a)*cos(c) sin(a)*sin(b)*cos(c)-cos(a)*sin(c);...
    -sin(b) cos(b)*sin(c) cos(b)*cos(c)];
r_init=[1 0 0; 0 1 0; 0 0 1]'; %Initial basis
r_rot=R*r_init; %Rotated basis

figure
hold on
%quiver3(origin(1),origin(2),origin(3),r_rot(1),r_rot(2),r_rot(3),'LineWidth',2,'MaxHeadSize',1)
line([-axlen axlen], [0 0], [0 0],'Color','k')
line([0 0], [-axlen axlen], [0 0],'Color','k')
line([0 0], [0 0], [-axlen axlen],'Color','k')

rotlx=line([0 r_init(1,1)], [0 r_init(2,1)], [0 r_init(3,1)],'Color','r','LineWidth',2);
line([0 r_init(1,2)], [0 r_init(2,2)], [0 r_init(3,2)],'Color','g','LineWidth',2)
line([0 r_init(1,3)], [0 r_init(2,3)], [0 r_init(3,3)],'Color','b','LineWidth',2)

rotlxp=line([0 r_rot(1,1)], [0 r_rot(2,1)], [0 r_rot(3,1)],'Color','r','LineWidth',2,'LineStyle','--');
line([0 r_rot(1,2)], [0 r_rot(2,2)], [0 r_rot(3,2)],'Color','g','LineWidth',2,'LineStyle','--')
line([0 r_rot(1,3)], [0 r_rot(2,3)], [0 r_rot(3,3)],'Color','b','LineWidth',2,'LineStyle','--')
hold off
legend([rotlx rotlxp],'Initial basis','Rotated basis','Interpreter','Latex')
pbaspect([1 1 1])
axis([-axlen axlen -axlen axlen -axlen axlen])
xlabel('x-axis','Interpreter','Latex')
ylabel('y-axis','Interpreter','Latex')
zlabel('z-axis','Interpreter','Latex')
campos([13.6971  17.8505   12.9904]) %Rotate to be in the standard representation
grid on

%% Quaternion Rotation Matrix
%xax=[1 0 0]; %x-axis
%yax=[0 1 0]; %y-axis
%zax=[0 0 1]; %z-axis

init=[1 0 0; 0 1 0; 0 0 1];

alph=0.3; %Roll
beta=pi/2; %Pitch
gamm=0; %Yaw

%Quaternion representation of basis axes for rotation
q1=[cos(alph/2) sin(alph/2) 0 0];
q2=[cos(beta/2) 0 sin(beta/2) 0];
q3=[cos(gamm/2) 0 0 sin(gamm/2)];

Rx=[1-2*(q1(3)^2+q1(4)^2) 2*(q1(2)*q1(3)-q1(4)*q1(1)) 2*(q1(2)*q1(4)+q1(3)*q1(1));...
    2*(q1(2)*q1(3)+q1(4)*q1(1)) 1-2*(q1(2)^2+q1(4)^2) 2*(q1(3)*q1(4)-q1(2)*q1(1));...
    2*(q1(2)*q1(4)-q1(3)*q1(1)) 2*(q1(3)*q1(4)+q1(2)*q1(1)) 1-2*(q1(2)^2+q1(3)^2)];

Ry=[1-2*(q2(3)^2+q2(4)^2) 2*(q2(2)*q2(3)-q2(4)*q2(1)) 2*(q2(2)*q2(4)+q2(3)*q2(1));...
    2*(q2(2)*q2(3)+q2(4)*q2(1)) 1-2*(q2(2)^2+q2(4)^2) 2*(q2(3)*q2(4)-q2(2)*q2(1));...
    2*(q2(2)*q2(4)-q2(3)*q2(1)) 2*(q2(3)*q2(4)+q2(2)*q2(1)) 1-2*(q2(2)^2+q2(3)^2)];

Rz=[1-2*(q3(3)^2+q3(4)^2) 2*(q3(2)*q3(3)-q3(4)*q3(1)) 2*(q3(2)*q3(4)+q3(3)*q3(1));...
    2*(q3(2)*q3(3)+q3(4)*q3(1)) 1-2*(q3(2)^2+q3(4)^2) 2*(q3(3)*q3(4)-q3(2)*q3(1));...
    2*(q3(2)*q3(4)-q3(3)*q3(1)) 2*(q3(3)*q3(4)+q3(2)*q3(1)) 1-2*(q3(2)^2+q3(3)^2)];

rot=Rz*Ry*Rx*init;

figure
hold on
%quiver3(origin(1),origin(2),origin(3),r_rot(1),r_rot(2),r_rot(3),'LineWidth',2,'MaxHeadSize',1)
line([-axlen axlen], [0 0], [0 0],'Color','k')
line([0 0], [-axlen axlen], [0 0],'Color','k')
line([0 0], [0 0], [-axlen axlen],'Color','k')

lx=line([0 init(1,1)], [0 init(2,1)], [0 init(3,1)],'Color','r','LineWidth',2);
line([0 init(1,2)], [0 init(2,2)], [0 init(3,2)],'Color','g','LineWidth',2)
line([0 init(1,3)], [0 init(2,3)], [0 init(3,3)],'Color','b','LineWidth',2)

lxp=line([0 rot(1,1)], [0 rot(2,1)], [0 rot(3,1)],'Color','r','LineWidth',2,'LineStyle','--');
line([0 rot(1,2)], [0 rot(2,2)], [0 rot(3,2)],'Color','g','LineWidth',2,'LineStyle','--')
line([0 rot(1,3)], [0 rot(2,3)], [0 rot(3,3)],'Color','b','LineWidth',2,'LineStyle','--')
hold off
legend([lx lxp],'Initial basis','Rotated basis','Interpreter','Latex')
pbaspect([1 1 1])
axis([-axlen axlen -axlen axlen -axlen axlen])
xlabel('x-axis','Interpreter','Latex')
ylabel('y-axis','Interpreter','Latex')
zlabel('z-axis','Interpreter','Latex')
campos([13.6971  17.8505   12.9904]) %Rotate to be in the standard representation
grid on

%% Rotation Using Quaternion Multiplication
clear alph beta gamm init rot
alph=pi/2; %Roll
beta=pi/2; %Pitch
gamm=0; %Yaw

[initial_q] = find_q_init(beta,gamm);
[rot_x] = rotated_x(initial_q);

px=[1 0 0]';
py=[0 1 0]';
pz=[0 0 1]';
init=[1 0 0; 0 1 0; 0 0 1];
qx=[cos(alph/2) sin(alph/2) 0 0];
qy=[cos(beta/2) 0 sin(beta/2) 0];
qz=[cos(gamm/2) 0 0 sin(gamm/2)];
q=[cos(alph/2) sin(alph/2) 0 0; cos(beta/2) 0 sin(beta/2) 0; cos(gamm/2) 0 0 sin(gamm/2)];

%px=px+2*cross(qx(2:4),qx(1)*px+cross(qx(2:4),px));

for i=1:3
    px=px+2*cross(q(i,2:4)',q(i,1)'*px+cross(q(i,2:4)',px));
    py=py+2*cross(q(i,2:4)',q(i,1)'*py+cross(q(i,2:4)',py));
    pz=pz+2*cross(q(i,2:4)',q(i,1)'*pz+cross(q(i,2:4)',pz));
end

rot=[px(1) py(1) pz(1); px(2) py(2) pz(2); px(3) py(3) pz(3)];

figure
hold on
%quiver3(origin(1),origin(2),origin(3),r_rot(1),r_rot(2),r_rot(3),'LineWidth',2,'MaxHeadSize',1)
line([-axlen axlen], [0 0], [0 0],'Color','k')
line([0 0], [-axlen axlen], [0 0],'Color','k')
line([0 0], [0 0], [-axlen axlen],'Color','k')

lx=line([0 init(1,1)], [0 init(2,1)], [0 init(3,1)],'Color','r','LineWidth',2);
line([0 init(1,2)], [0 init(2,2)], [0 init(3,2)],'Color','g','LineWidth',2)
line([0 init(1,3)], [0 init(2,3)], [0 init(3,3)],'Color','b','LineWidth',2)

lxp=line([0 rot(1,1)], [0 rot(2,1)], [0 rot(3,1)],'Color','r','LineWidth',2,'LineStyle','--');
line([0 rot(1,2)], [0 rot(2,2)], [0 rot(3,2)],'Color','g','LineWidth',2,'LineStyle','--')
line([0 rot(1,3)], [0 rot(2,3)], [0 rot(3,3)],'Color','b','LineWidth',2,'LineStyle','--')
hold off
legend([lx lxp],'Initial basis','Rotated basis','Interpreter','Latex')
pbaspect([1 1 1])
axis([-axlen axlen -axlen axlen -axlen axlen])
xlabel('x-axis','Interpreter','Latex')
ylabel('y-axis','Interpreter','Latex')
zlabel('z-axis','Interpreter','Latex')
campos([13.6971  17.8505   12.9904]) %Rotate to be in the standard representation
grid on

%% Functions
function [qinit] = find_q_init(beta,gamma)
    vec=[1 0 0]';
    rho = [cos(gamma)*cos(beta); sin(gamma)*cos(beta); sin(beta)];
    
    qinit_0 = [1+dot(vec,rho);cross(vec,rho)];
    qinit = qinit_0/norm(qinit_0);
end

function [rot_x] = rotated_x(qinit)
    vecx = [1 0 0]';
    
    rot_x = vecx + 2*cross(qinit(2:4),qinit(1)*vecx+cross(qinit(2:4),vecx));
end
