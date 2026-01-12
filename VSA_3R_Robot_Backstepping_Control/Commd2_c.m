function [sys,x0,str,ts] = Commd2(t,x,u,flag)
switch flag,
case 0,
    [sys,x0,str,ts]=mdlInitializeSizes;
case 1,
    sys=mdlDerivatives(t,x,u);
case 3, 
    sys=mdlOutputs(t,x,u);
case {2,4,9}
    sys=[]; 
otherwise
    error(['Unhandled flag = ',num2str(flag)]);
end

function [sys,x0,str,ts]=mdlInitializeSizes
global c1 c2 eta epsilon 
sizes = simsizes;
sizes.NumContStates  = 1;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 13;
sizes.NumInputs      = 57;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0= 0;
str = [];
ts  = [];

eta=2 ;  epsilon=0.05 ;

c1=[-9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9 ;-9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9 ; 
    -9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9 ; -9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9 ] ;

c2=[-9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9 ; -9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9 ; -9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9 ; 
    -9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9 ; -9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9 ; -9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9 ; 
    -9,-36/5,-27/5,-18/5,-9/5,0,9/5,18/5,27/5,36/5,9  ] ;


function sys=mdlDerivatives(t,x,u)
global c1 c2 eta epsilon 

k01=[30 0 0 ; 0 30 0 ; 0 0 30]  ;

qdes=u(1:3) ; dqdes=u(4:6)  ; 
q=u(19:21) ;  q1=q(1) ;q2=q(2) ; q3=q(3) ;
dq=u(22:24) ; dq1=dq(1) ; dq2=dq(2) ;dq3=dq(3) ;
x1=u(25:27) ; the1=x1(1) ; the2=x1(2)  ; the3=x1(3) ;
x2= u(28:30) ; 

x2c=u(40:42) ; eps1=u(49:51) ; 
 
% % Etats
v=dqdes-k01*(q-qdes) ;
e=dq-v ;

kix2= x2-x2c ;
z1=e-eps1 ;   z3=kix2   ;   

% % ***RBF
% %  

rho_est= x(1) ;

gam=0.5 ; gam1=2 ;
z1t=z1' ;  z3t=z3' ;

Z1=[ q1; q2; q3; rho_est; ] ;
Z2=[ q1; q2; q3; the1; the2; the3; rho_est ] ;


for j=1:1:11
    S1(j)=exp(-norm(Z1-c1(:,j))^2/(eta^2));  %Hidden Layer
end 

for j=1:1:11
    S2(j)=exp(-norm(Z2-c2(:,j))^2/(eta^2));  %Hidden Layer
end 

U1=norm(S1)+1;

U2=norm(S2)+1;

U1_ch=[U1*tanh((z1(1)*U1)/epsilon) ,U1*tanh((z1(2)*U1)/epsilon) ,U1*tanh((z1(3)*U1)/epsilon)]' ;

U2_ch=[U2*tanh((z3(1)*U2)/epsilon) ,U2*tanh((z3(2)*U2)/epsilon) ,U2*tanh((z3(3)*U2)/epsilon)]' ;

% loi adaptative

sys=gam1*(z1t*U1_ch+z3t*U2_ch)-gam*rho_est ;


function sys=mdlOutputs(t,x,u)

global c1 c2 eta epsilon 

k01=[30 0 0 ; 0 30 0 ; 0 0 30]  ;
k02=[10 0 0 ; 0 10 0 ; 0 0 10]  ;

% % % Entrées

% % 
qdes=u(1:3) ; dqdes=u(4:6)  ; ddqdes=u(7:9) ;
x3des=u(10:12) ; dx3des=u(13:15); ddx3des=u(16:18) ;

q=u(19:21) ;  q1=q(1) ;q2=q(2) ; q3=q(3) ;
dq=u(22:24) ; dq1=dq(1) ; dq2=dq(2) ;dq3=dq(3) ;
x1=u(25:27) ; the1=x1(1) ; the2=x1(2)  ; the3=x1(3) ;
x2= u(28:30) ; 
x3=u(31:33) ; sig1=x3(1)  ; sig2=x3(2)  ; sig3=x3(3) ;
x4=u(34:36) ;

x1c=u(37:39) ; x2c=u(40:42) ; dx1c=u(43:45)  ; dx2c=u(46:48) ;
eps1=u(49:51) ; eps2=u(52:54) ;
Kdes=u(55:57) ; 
Rdes=[Kdes(1) 0 0 ; 0 Kdes(2) 0 ; 0 0 Kdes(3)] ;
    

% % Etats
v=dqdes-k01*(q-qdes) ;
e=dq-v ;
dv= ddqdes-k01*(dq-dqdes) ;

vsig=dx3des-k02*(x3-x3des) ;
esig=x4-vsig;
dvsig=ddx3des-k02*(x4-dx3des) ;

kix1= x1-x1c ;  kix2= x2-x2c ;

z1=e-eps1 ;   z2=kix1-eps2 ;
z3=kix2   ;   z4=esig+vsig  ;


 % Parametres robot

l2 =0.760; l3=0.930; d=0.28;
m1=5; m2=5; m3=3; lc1=0.22; lc2=0.51; lc3=0.67; 
Ix1=.34 ; Iy1=.36 ; Iz1=.31; Ix2=0.18; Iy2=1.32; Iz2=1.31; 
Ix3=0.07; Iy3=0.92; Iz3=0.93; g0=9.8;

% Inertie

M=zeros(3);

M(1,1)=4*lc3*m3*cos(q2)*cos(q3)*l3*sin(q2)*sin(q3)-2*Iy3*cos(q3)*cos(q2)*sin(q2)*sin(q3)-m3*lc3^2*cos(q3)^2-m3*lc3^2*cos(q2)^2-m3*l3^2*cos(q3)^2-m3*l3^2*cos(q2)^2-2*Ix3*cos(q3)^2*cos(q2)^2+2*Iy3*cos(q2)^2*cos(q3)^2+m3*cos(q2)^2*l2^2-4*lc3*m3*cos(q2)^2*cos(q3)^2*l3-2*m3*sin(q2)*sin(q3)*l3*cos(q2)*l2-2*m3*cos(q2)*cos(q3)*l3^2*sin(q2)*sin(q3)+2*lc3*m3*cos(q2)*l2*sin(q2)*sin(q3)-2*m3*lc3^2*sin(q2)*sin(q3)*cos(q2)*cos(q3)+cos(q2)^2*Iy2-Ix2*cos(q2)^2+m2*cos(q2)^2*lc2^2+m2*cos(q2)^2*l2^2-2*lc2*m2*cos(q2)^2*l2+Iy1+Ix3*cos(q3)^2+Ix3*cos(q2)^2-Iy3*cos(q3)^2-Iy3*cos(q2)^2+Ix2+Iy3-2*lc3*m3*cos(q2)^2*l2*cos(q3)+2*m3*cos(q2)^2*cos(q3)*l3*l2+m3*lc3^2+m3*l3^2+2*m3*lc3^2*cos(q2)^2*cos(q3)^2+2*lc3*m3*l3*cos(q3)^2+2*lc3*m3*l3*cos(q2)^2+2*m3*cos(q2)^2*cos(q3)^2*l3^2+2*Ix3*cos(q3)*cos(q2)*sin(q2)*sin(q3)-2*lc3*m3*l3 ;

M(1,2)= 0 ;

M(1,3)= 0 ;

M(2,1)= 0 ;

M(2,2)= -2*lc2*m2*l2+m2*lc2^2+m2*l2^2+Iz2+2*m3*cos(q3)*l3*l2-2*lc3*m3*l3+m3*l2^2+m3*l3^2-2*lc3*m3*l2*cos(q3)+m3*lc3^2+Iz3;

M(2,3)= Iz3+m3*lc3^2+m3*l3^2-2*lc3*m3*l3-m3*lc3*cos(q3)*l2+m3*cos(q3)*l3*l2 ;

M(3,1)= 0 ;

M(3,2)=Iz3-lc3*m3*l2*cos(q3)+m3*lc3^2-2*lc3*m3*l3+m3*l3^2+m3*l3*l2*cos(q3) ;

M(3,3)= m3*lc3^2+m3*l3^2-2*lc3*m3*l3+Iz3 ;

% Gravite

G=zeros(3,1);

G(1)= 0 ;

G(2)= g0*(m2*lc2*cos(q2)+m3*cos(q2)*l2+m3*lc3*cos(q2+q3)) ;

G(3)= g0*m3*lc3*cos(q2+q3) ;

% Coriolis

C=zeros(3);

C(1,1)= Iy3*dq2*sin(q3)*cos(q3)+Iy3*dq2*sin(q2)*cos(q2)+Iy3*dq3*sin(q3)*cos(q3)+Iy3*dq3*sin(q2)*cos(q2)-Ix3*dq2*sin(q3)*cos(q3)-Ix3*dq2*sin(q2)*cos(q2)-Ix3*dq3*sin(q3)*cos(q3)-Ix3*dq3*sin(q2)*cos(q2)-2*m3*dq2*sin(q2)*cos(q3)*l3*cos(q2)*l2-m3*dq3*l3*sin(q2)*cos(q3)*cos(q2)*l2+4*lc3*m3*dq2*sin(q2)*cos(q3)^2*l3*cos(q2)+2*lc3*m3*l2*dq2*sin(q2)*cos(q2)*cos(q3)+4*lc3*m3*dq3*l3*sin(q2)*cos(q3)^2*cos(q2)+4*lc3*m3*dq3*l3*cos(q3)*sin(q3)*cos(q2)^2+4*lc3*m3*dq2*cos(q3)*l3*sin(q3)*cos(q2)^2+m3*lc3*dq3*sin(q2)*cos(q3)*cos(q2)*l2+Ix2*cos(q2)*dq2*sin(q2)-cos(q1)^2*sin(q2)*Ix2*cos(q2)*dq2+cos(q1)^2*sin(q2)*Iy2*cos(q2)*dq2-Iy2*cos(q2)*dq2*sin(q2)+2*m3*lc3*dq2*cos(q2)^2*sin(q3)*l2+m3*lc3*dq3*cos(q2)^2*sin(q3)*l2-2*lc3*m3*dq2*cos(q3)*l3*sin(q3)-2*lc3*m3*dq2*cos(q2)*l3*sin(q2)-2*lc3*m3*dq3*l3*cos(q2)*sin(q2)-m3*dq3*l3*cos(q2)^2*sin(q3)*l2-2*lc3*m3*dq3*l3*cos(q3)*sin(q3)-2*m3*lc3^2*dq2*cos(q3)*sin(q3)*cos(q2)^2-2*m3*lc3^2*dq2*cos(q2)*sin(q2)*cos(q3)^2-2*m3*lc3^2*dq3*cos(q3)*sin(q3)*cos(q2)^2-2*m3*lc3^2*dq3*cos(q2)*sin(q2)*cos(q3)^2-2*m3*dq2*cos(q3)*l3^2*sin(q3)*cos(q2)^2-2*m3*dq2*cos(q2)*l3^2*sin(q2)*cos(q3)^2-2*m3*l2*dq2*sin(q3)*l3*cos(q2)^2-2*m3*dq3*l3^2*cos(q3)*sin(q3)*cos(q2)^2-2*m3*dq3*l3^2*cos(q2)*sin(q2)*cos(q3)^2+2*Iy3*dq2*sin(q3)*cos(q3)*cos(q1)^2*cos(q2)^2+2*Iy3*dq2*sin(q2)*cos(q2)*cos(q1)^2*cos(q3)^2+2*Iy3*dq3*sin(q3)*cos(q3)*cos(q1)^2*cos(q2)^2+2*Iy3*dq3*sin(q2)*cos(q2)*cos(q1)^2*cos(q3)^2-2*Ix3*dq2*sin(q3)*cos(q3)*cos(q1)^2*cos(q2)^2-2*Ix3*dq2*sin(q2)*cos(q2)*cos(q1)^2*cos(q3)^2-2*Ix3*dq3*sin(q3)*cos(q3)*cos(q1)^2*cos(q2)^2-2*Ix3*dq3*sin(q2)*cos(q2)*cos(q1)^2*cos(q3)^2+m3*lc3^2*dq2*cos(q3)*sin(q3)+m3*lc3^2*dq2*cos(q2)*sin(q2)+m3*lc3^2*dq3*cos(q3)*sin(q3)+m3*lc3^2*dq3*cos(q2)*sin(q2)+m3*dq2*cos(q3)*l3^2*sin(q3)+m3*dq2*cos(q2)*l3^2*sin(q2)+m3*l2*dq2*sin(q3)*l3+m3*dq3*l3^2*cos(q3)*sin(q3)+m3*dq3*l3^2*cos(q2)*sin(q2)-lc3*m3*l2*dq2*sin(q3)-m3*l2^2*dq2*sin(q2)*cos(q2)-2*Iy3*dq2*sin(q3)*cos(q3)*cos(q2)^2-Iy3*dq2*sin(q3)*cos(q3)*cos(q1)^2-2*Iy3*dq2*sin(q2)*cos(q2)*cos(q3)^2-Iy3*dq2*sin(q2)*cos(q2)*cos(q1)^2-2*Iy3*dq3*sin(q3)*cos(q3)*cos(q2)^2-Iy3*dq3*sin(q3)*cos(q3)*cos(q1)^2-2*Iy3*dq3*sin(q2)*cos(q2)*cos(q3)^2-Iy3*dq3*sin(q2)*cos(q2)*cos(q1)^2+2*Ix3*dq2*sin(q3)*cos(q3)*cos(q2)^2+Ix3*dq2*sin(q3)*cos(q3)*cos(q1)^2+2*Ix3*dq2*sin(q2)*cos(q2)*cos(q3)^2+Ix3*dq2*sin(q2)*cos(q2)*cos(q1)^2+2*Ix3*dq3*sin(q3)*cos(q3)*cos(q2)^2+Ix3*dq3*sin(q3)*cos(q3)*cos(q1)^2+2*Ix3*dq3*sin(q2)*cos(q2)*cos(q3)^2+Ix3*dq3*sin(q2)*cos(q2)*cos(q1)^2-m2*cos(q2)*lc2^2*dq2*sin(q2)-m2*cos(q2)*l2^2*dq2*sin(q2)+2*lc2*m2*cos(q2)*l2*dq2*sin(q2) ;

C(1,2)= -sin(q1)*Ix3*cos(q1)*dq3*cos(q3)^2-m2*cos(q2)*lc2^2*dq1*sin(q2)-m2*cos(q2)*l2^2*dq1*sin(q2)+Ix2*dq1*cos(q2)*sin(q2)-Ix2*dq2*sin(q1)*cos(q1)-Iy3*dq1*sin(q2)*cos(q2)*cos(q1)^2-2*Iy3*dq1*sin(q2)*cos(q2)*cos(q3)^2+2*sin(q1)*Ix3*cos(q1)*dq3*cos(q3)^2*cos(q2)^2-2*lc3*m3*dq1*cos(q3)*l3*sin(q3)-2*sin(q1)*Iy3*cos(q1)*dq3*cos(q2)^2*cos(q3)^2-2*lc3*m3*dq1*cos(q2)*l3*sin(q2)+2*m3*lc3*dq1*cos(q2)^2*sin(q3)*l2-2*m3*lc3^2*dq1*cos(q3)*sin(q3)*cos(q2)^2-2*m3*lc3^2*dq1*cos(q2)*sin(q2)*cos(q3)^2-2*m3*dq1*cos(q3)*l3^2*sin(q3)*cos(q2)^2-2*m3*dq1*cos(q2)*l3^2*sin(q2)*cos(q3)^2-2*m3*l2*dq1*sin(q3)*l3*cos(q2)^2+2*sin(q1)*Ix3*cos(q1)*dq2*cos(q3)^2*cos(q2)^2+2*Iy3*dq1*sin(q2)*cos(q2)*cos(q1)^2*cos(q3)^2-2*Ix3*dq1*sin(q2)*cos(q2)*cos(q1)^2*cos(q3)^2-2*sin(q1)*Iy3*cos(q1)*dq2*cos(q2)^2*cos(q3)^2+sin(q1)*Iy3*cos(q1)*dq2*cos(q3)^2+2*Iy3*dq1*sin(q3)*cos(q3)*cos(q1)^2*cos(q2)^2-2*Ix3*dq1*sin(q3)*cos(q3)*cos(q1)^2*cos(q2)^2+2*Ix3*dq1*sin(q2)*cos(q2)*cos(q3)^2+Ix3*dq1*sin(q2)*cos(q2)*cos(q1)^2+2*Ix3*dq1*sin(q3)*cos(q3)*cos(q2)^2+Ix3*dq1*sin(q3)*cos(q3)*cos(q1)^2+m3*lc3^2*dq1*cos(q2)*sin(q2)+m3*lc3^2*dq1*cos(q3)*sin(q3)-Iy3*dq1*sin(q3)*cos(q3)*cos(q1)^2-2*Iy3*dq1*sin(q3)*cos(q3)*cos(q2)^2-sin(q1)*Ix3*cos(q1)*dq2*cos(q2)^2-sin(q1)*Ix3*cos(q1)*dq2*cos(q3)^2-sin(q1)*Ix3*cos(q1)*dq3*cos(q2)^2-2*m3*dq1*sin(q2)*cos(q3)*l3*cos(q2)*l2+2*sin(q1)*Iy3*cos(q1)*dq2*sin(q2)*sin(q3)*cos(q3)*cos(q2)+2*sin(q1)*Iy3*cos(q1)*dq3*sin(q2)*sin(q3)*cos(q3)*cos(q2)+4*lc3*m3*dq1*cos(q3)*l3*sin(q3)*cos(q2)^2+2*lc3*m3*l2*dq1*sin(q2)*cos(q2)*cos(q3)+4*lc3*m3*dq1*sin(q2)*cos(q3)^2*l3*cos(q2)-2*sin(q1)*Ix3*cos(q1)*dq2*sin(q2)*sin(q3)*cos(q3)*cos(q2)-2*sin(q1)*Ix3*cos(q1)*dq3*sin(q2)*sin(q3)*cos(q3)*cos(q2)-lc3*m3*l2*dq1*sin(q3)-Ix3*dq1*sin(q3)*cos(q3)-Iy2*dq1*cos(q2)*sin(q2)+m3*dq1*cos(q3)*l3^2*sin(q3)+2*lc2*m2*cos(q2)*l2*dq1*sin(q2)-Ix2*dq1*cos(q2)*sin(q2)*cos(q1)^2+Ix2*dq2*sin(q1)*cos(q1)*cos(q2)^2+Iy2*dq1*cos(q2)*sin(q2)*cos(q1)^2-sin(q1)*Iy2*cos(q2)^2*cos(q1)*dq2+Iy3*dq1*sin(q3)*cos(q3)+Iy3*dq1*sin(q2)*cos(q2)-sin(q1)*Iy3*cos(q1)*dq3-sin(q1)*Iy3*cos(q1)*dq2-Ix3*dq1*sin(q2)*cos(q2)+m3*dq1*cos(q2)*l3^2*sin(q2)+m3*l2*dq1*sin(q3)*l3+sin(q1)*Iy3*cos(q1)*dq3*cos(q3)^2+sin(q1)*Iy3*cos(q1)*dq3*cos(q2)^2-m3*l2^2*dq1*sin(q2)*cos(q2)+sin(q1)*Iy3*cos(q1)*dq2*cos(q2)^2 ;

C(1,3)= 4*lc3*m3*l3*dq1*sin(q2)*cos(q3)^2*cos(q2)+4*lc3*m3*l3*dq1*cos(q2)^2*sin(q3)*cos(q3)-m3*l3*dq1*sin(q2)*cos(q3)*cos(q2)*l2+m3*lc3*dq1*sin(q2)*cos(q3)*cos(q2)*l2-Ix3*dq1*sin(q3)*cos(q3)-Ix3*dq1*sin(q2)*cos(q2)-sin(q1)*Iy3*cos(q1)*dq2-sin(q1)*Iy3*cos(q1)*dq3+Iy3*dq1*sin(q3)*cos(q3)+Iy3*dq1*sin(q2)*cos(q2)+m3*lc3*dq1*cos(q2)^2*sin(q3)*l2-m3*l3*dq1*cos(q2)^2*sin(q3)*l2-2*lc3*m3*l3*dq1*cos(q2)*sin(q2)-2*lc3*m3*l3*dq1*cos(q3)*sin(q3)-2*cos(q1)^2*sin(q2)*cos(q3)^2*Ix3*dq1*cos(q2)-2*m3*lc3^2*dq1*cos(q2)*sin(q2)*cos(q3)^2-2*m3*lc3^2*dq1*cos(q3)*sin(q3)*cos(q2)^2-2*m3*l3^2*dq1*cos(q3)*sin(q3)*cos(q2)^2-2*m3*l3^2*dq1*cos(q2)*sin(q2)*cos(q3)^2-2*cos(q1)^2*cos(q2)^2*sin(q3)*Ix3*dq1*cos(q3)+2*sin(q1)*Ix3*cos(q1)*dq2*cos(q3)^2*cos(q2)^2+2*sin(q1)*Ix3*cos(q1)*dq3*cos(q3)^2*cos(q2)^2+2*cos(q1)^2*sin(q3)*Iy3*dq1*cos(q3)*cos(q2)^2-2*sin(q1)*Iy3*cos(q1)*dq2*cos(q3)^2*cos(q2)^2+2*cos(q1)^2*sin(q2)*cos(q3)^2*Iy3*dq1*cos(q2)-2*sin(q1)*Iy3*cos(q1)*dq3*cos(q3)^2*cos(q2)^2+m3*lc3^2*dq1*cos(q2)*sin(q2)+m3*lc3^2*dq1*cos(q3)*sin(q3)+m3*l3^2*dq1*cos(q3)*sin(q3)+m3*l3^2*dq1*cos(q2)*sin(q2)+cos(q1)^2*sin(q3)*Ix3*dq1*cos(q3)+cos(q1)^2*sin(q2)*Ix3*dq1*cos(q2)+sin(q1)*Iy3*cos(q1)*dq2*cos(q3)^2+sin(q1)*Iy3*cos(q1)*dq2*cos(q2)^2+sin(q1)*Iy3*cos(q1)*dq3*cos(q3)^2+sin(q1)*Iy3*cos(q1)*dq3*cos(q2)^2-2*Iy3*dq1*sin(q3)*cos(q3)*cos(q2)^2-2*Iy3*dq1*sin(q2)*cos(q2)*cos(q3)^2+2*Ix3*dq1*sin(q3)*cos(q3)*cos(q2)^2+2*Ix3*dq1*sin(q2)*cos(q2)*cos(q3)^2-sin(q1)*Ix3*cos(q1)*dq2*cos(q3)^2-sin(q1)*Ix3*cos(q1)*dq2*cos(q2)^2-sin(q1)*Ix3*cos(q1)*dq3*cos(q3)^2-sin(q1)*Ix3*cos(q1)*dq3*cos(q2)^2-cos(q1)^2*sin(q3)*Iy3*dq1*cos(q3)-cos(q1)^2*sin(q2)*Iy3*dq1*cos(q2)-2*sin(q1)*Ix3*cos(q1)*dq3*sin(q2)*sin(q3)*cos(q3)*cos(q2)-2*sin(q1)*Ix3*cos(q1)*dq2*sin(q2)*sin(q3)*cos(q3)*cos(q2)+2*sin(q1)*Iy3*cos(q1)*dq2*sin(q2)*sin(q3)*cos(q3)*cos(q2)+2*sin(q1)*Iy3*cos(q1)*dq3*sin(q2)*sin(q3)*cos(q3)*cos(q2) ;

C(2,1)= -cos(q1)*Iz2*dq2*sin(q1)-m3*dq1*cos(q2)*l3^2*sin(q2)-m3*dq1*cos(q3)*l3^2*sin(q3)-m3*dq1*sin(q3)*l3*l2+m3*lc3*dq1*sin(q3)*l2-2*lc2*m2*sin(q2)*l2*dq1*cos(q2)+2*lc3*m3*dq1*sin(q3)*l3*cos(q3)+2*lc3*m3*dq1*sin(q2)*l3*cos(q2)+m2*sin(q2)*lc2^2*dq1*cos(q2)+m2*sin(q2)*l2^2*dq1*cos(q2)-4*lc3*m3*dq1*cos(q2)^2*cos(q3)*l3*sin(q3)-4*lc3*m3*dq1*cos(q2)*cos(q3)^2*l3*sin(q2)-2*lc3*m3*l2*dq1*cos(q2)*sin(q2)*cos(q3)-cos(q1)*Iz3*dq2*sin(q1)-cos(q1)*Iz3*dq3*sin(q1)+2*m3*dq1*cos(q2)*cos(q3)^2*l3^2*sin(q2)+2*m3*dq1*cos(q2)^2*cos(q3)*l3^2*sin(q3)+2*m3*l2*dq1*cos(q2)^2*sin(q3)*l3+2*m3*lc3^2*dq1*cos(q2)*cos(q3)^2*sin(q2)+2*m3*lc3^2*dq1*cos(q2)^2*cos(q3)*sin(q3)-2*m3*lc3*dq1*sin(q3)*l2*cos(q2)^2+2*m3*dq1*cos(q2)*cos(q3)*l3*sin(q2)*l2+m3*l2^2*dq1*cos(q2)*sin(q2)-m3*lc3^2*dq1*sin(q3)*cos(q3)-m3*lc3^2*dq1*sin(q2)*cos(q2) ;

C(2,2)= -cos(q1)*Iz2*dq1*sin(q1)+m3*lc3*sin(q3)*dq3*l2-m3*dq3*l3*sin(q3)*l2-Iz3*cos(q1)*dq1*sin(q1) ;

C(2,3)= -m3*l3*dq3*sin(q3)*l2-m3*l3*dq2*sin(q3)*l2-cos(q1)*Iz3*dq1*sin(q1)+m3*lc3*dq2*sin(q3)*l2+m3*lc3*dq3*sin(q3)*l2 ;

C(3,1)= -Iz3*dq2*sin(q1)*cos(q1)-Iz3*dq3*sin(q1)*cos(q1)-m3*lc3^2*dq1*sin(q3)*cos(q3)-m3*lc3^2*dq1*sin(q2)*cos(q2)-m3*l3^2*dq1*cos(q2)*sin(q2)-m3*l3^2*dq1*cos(q3)*sin(q3)+2*lc3*m3*dq1*cos(q2)*l3*sin(q2)+2*lc3*m3*dq1*cos(q3)*l3*sin(q3)-lc3*m3*dq1*cos(q2)*l2*sin(q2)*cos(q3)+m3*l3*dq1*cos(q2)*l2*sin(q2)*cos(q3)+2*m3*lc3^2*dq1*cos(q2)*cos(q3)^2*sin(q2)+2*m3*lc3^2*dq1*cos(q2)^2*cos(q3)*sin(q3)+m3*l3*dq1*cos(q2)^2*l2*sin(q3)+2*m3*l3^2*dq1*cos(q2)*cos(q3)^2*sin(q2)+2*m3*l3^2*dq1*cos(q2)^2*cos(q3)*sin(q3)-lc3*m3*dq1*cos(q2)^2*l2*sin(q3)-4*lc3*m3*dq1*cos(q2)*cos(q3)^2*l3*sin(q2)-4*lc3*m3*dq1*cos(q2)^2*cos(q3)*l3*sin(q3) ;

C(3,2)= -lc3*m3*dq2*l2*sin(q3)+m3*l3*dq2*l2*sin(q3)-cos(q1)*Iz3*dq1*sin(q1) ;

C(3,3)= -lc3*m3*dq2*l2*sin(q3)+m3*l3*dq2*l2*sin(q3)-cos(q1)*Iz3*dq1*sin(q1) ;


% Moteur 

 Jthe=[2,0,0 ; 0,2,0 ; 0,0,2] ;
 Jsig=[0.0002,0,0 ; 0,0.0002,0 ; 0,0,0.0002] ;
 Bthe=[0.001,0,0 ; 0,0.001,0 ; 0,0,0.001] ;
 Bsig=[0.001,0,0 ; 0,0.001,0 ; 0,0,0.001] ;
 Kthe=[0.26,0,0 ; 0,0.26,0 ; 0,0,0.26] ;
 Ksig=[0.26,0,0 ; 0,0.26,0 ; 0,0,0.26] ;

% % Couple_reaction

KR1= 2.43*10^-6 ;

n=0.006 ;
delta=0.015 ;

F= zeros(3,1);

F(1)= (KR1*sig1*(the1-q1)^2)/(delta-n*sig1)^3 ;

F(2)= (KR1*sig2*(the2-q2)^2)/(delta-n*sig2)^3 ;

F(3)= (KR1*sig3*(the3-q3)^2)/(delta-n*sig3)^3 ;

%  ***RBF
% % % 
rho_est=x(1) ;

Z1=[ q1; q2; q3; rho_est ] ;
Z2=[ q1; q2; q3; the1; the2; the3; rho_est ] ;


for j=1:1:11
    S1(j)=exp(-norm(Z1-c1(:,j))^2/(eta^2));  %Hidden Layer
end 

for j=1:1:11
    S2(j)=exp(-norm(Z2-c2(:,j))^2/(eta^2));  %Hidden Layer
end 

U1=norm(S1)+1;

U2=norm(S2)+1;

U1_ch=[U1*tanh((z1(1)*U1)/epsilon) ,U1*tanh((z1(2)*U1)/epsilon) ,U1*tanh((z1(3)*U1)/epsilon)]' ;

U2_ch=[U2*tanh((z3(1)*U2)/epsilon) ,U2*tanh((z3(2)*U2)/epsilon) ,U2*tanh((z3(3)*U2)/epsilon)]' ;

% Parametres Backstepping 
% % 

k1=[30 0 0 ; 0 30 0 ; 0 0 30] ;
k2=[30 0 0 ; 0 30 0 ; 0 0 30] ;
k3=[2 0 0 ; 0 1 0 ; 0 0 2] ;
k4=[2 0 0 ; 0 2 0 ; 0 0 2] ;


% fonctions package
% % 

f1=  C*v+ G + M*dv+ Rdes*q ;
 
f2= (Jthe\(Bthe*x2+Rdes*(x1-q))) ;
 
f3= (Jsig\(Bsig*x4 + Jsig*dvsig + F ))  ;

% commande virtuelles
% % 
alpha1=Rdes\(-k1*e+f1-rho_est*U1_ch) ;

alpha2=-k2*kix1-Rdes*z1+dx1c ;

% Commande
% % 
U_the= Kthe\(Jthe*(-k3*z3-z2+f2+ dx2c-rho_est*U2_ch)) ;

U_sig= Ksig\(Jsig*(-k4*z4+f3-dvsig)) ;


% % 

sys(1)=alpha1(1);
sys(2)=alpha1(2);
sys(3)=alpha1(3);

sys(4)=alpha2(1);
sys(5)=alpha2(2);
sys(6)=alpha2(3);

sys(7)=U_the(1) ;
sys(8)=U_the(2) ;
sys(9)=U_the(3) ;

sys(10)=U_sig(1);
sys(11)=U_sig(2);
sys(12)=U_sig(3);

sys(13)=rho_est ; 