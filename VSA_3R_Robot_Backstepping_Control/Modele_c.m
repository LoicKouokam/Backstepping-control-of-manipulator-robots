function [sys,x0,str,ts] = planVSA_2(t,x,u,flag)
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
sizes = simsizes;
sizes.NumContStates  = 18;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 18;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0= [0 0 0 0 0 0 0 0 0 0 0 0 2.83-0.55 2.83-0.55 2.83-0.55 0 0 0];
str = [];
ts  = [];

function sys=mdlDerivatives(t,x,u)
 
% Parametres robot

l2 =0.760; l3=0.930; d=0.28;
m1=5; m2=5; m3=3; lc1=0.22; lc2=0.51; lc3=0.67; 
Ix1=.34 ; Iy1=.36 ; Iz1=.31; Ix2=0.18; Iy2=1.32; Iz2=1.31; 
Ix3=0.07; Iy3=0.92; Iz3=0.93; g0=9.8;


% % etats  

q=x(1:3) ;  q1=q(1) ; q2=q(2) ; q3=q(3)  ;
dq=x(4:6) ; dq1=dq(1) ; dq2=dq(2) ;  dq3=dq(3);
the=x(7:9) ;   the1=the(1) ; the2=the(2) ; the3=the(3) ;
dthe=x(10:12) ;
sig=x(13:15) ; sig1=sig(1) ; sig2=sig(2) ; sig3=sig(3) ;
dsig=x(16:18) ;


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
 Jsig=[0.0002,0,0  ; 0,0.0002,0 ; 0,0,0.0002] ;
 Bthe=[0.001,0,0 ; 0,0.001,0 ; 0,0,0.001] ;
 Bsig=[0.001,0,0 ; 0,0.001,0 ; 0,0,0.001] ;
 Kthe=[0.26,0,0 ; 0,0.26,0 ; 0,0,0.26] ;
 Ksig=[0.26,0,0 ; 0,0.26,0 ; 0,0,0.26] ;

% % Rigidite

KE1= 1.62*10^-4 ;
% KE1 =0.001 ;
delta=0.015 ;
n=  0.006 ;

R= zeros(3);

R(1,1)=(KE1*sig1^2)/(delta-n*sig1)^2 ;

R(2,2)= (KE1*sig2^2)/(delta-n*sig2)^2 ;

R(3,3)= (KE1*sig3^2)/(delta-n*sig3)^2 ;

 
% % Couple_reaction

KR1= 2.43*10^-6 ;
 

F= zeros(3,1);

F(1)= (KR1*sig1*(the1-q1)^2)/(delta-n*sig1)^3 ;

F(2)= (KR1*sig2*(the2-q2)^2)/(delta-n*sig2)^3 ;

F(3)= (KR1*sig3*(the3-q3)^2)/(delta-n*sig3)^3 ;

% % 

diff_q=dq;

diff_the=dthe;

diff_sig=dsig;

uthe=u(1:3);
usig=u(4:6);
% % 

diff_diff_q= M\(R*(the-q)-C*diff_q-G) ;
diff_diff_the= Jthe\(Kthe*uthe-Bthe*diff_the-R*(the-q)) ;
diff_diff_sig= Jsig\(Ksig*usig-Bsig*diff_sig-F) ;


% %

sys(1)=diff_q(1);
sys(2)=diff_q(2);
sys(3)=diff_q(3);

sys(4)=diff_diff_q(1);
sys(5)=diff_diff_q(2);
sys(6)=diff_diff_q(3);

sys(7)=diff_the(1);
sys(8)=diff_the(2);
sys(9)=diff_the(3);

sys(10)=diff_diff_the(1);
sys(11)=diff_diff_the(2);
sys(12)=diff_diff_the(3);

sys(13)=diff_sig(1);
sys(14)=diff_sig(2);
sys(15)=diff_sig(3);

sys(16)=diff_diff_sig(1);
sys(17)=diff_diff_sig(2);
sys(18)=diff_diff_sig(3);


function sys=mdlOutputs(t,x,u)

sys(1)=x(1);
sys(2)=x(2);
sys(3)=x(3);
sys(4)=x(4);
sys(5)=x(5);
sys(6)=x(6);
sys(7)=x(7);
sys(8)=x(8);
sys(9)=x(9);
sys(10)=x(10);
sys(11)=x(11);
sys(12)=x(12);
sys(13)=x(13);
sys(14)=x(14);
sys(15)=x(15);
sys(16)=x(16);
sys(17)=x(17);
sys(18)=x(18);