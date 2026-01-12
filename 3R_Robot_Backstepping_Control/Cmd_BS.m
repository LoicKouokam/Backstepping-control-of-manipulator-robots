function [sys,x0,str,ts] = Modele(t,x,u,flag)

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
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 3;
sizes.NumInputs      = 15;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0= [];
str = [];
ts  = [];


function sys=mdlOutputs(t,x,u)

xr=u(1:3) ; dxr=u(4:6)  ; ddxr=u(7:9) ;
 
x1=u(10:12) ; x2=u(13:15) ;

q1=x1(1) ; q2=x1(2)  ;q3=x1(3) ;

dq1=x2(1) ; dq2=x2(2) ; dq3=x2(3) ;

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

% % 

c1=[10 0 0; 0 10 0 ;0 0 10] ; c2=[5 0 0;0 5 0 ; 0 0 5 ] ;

z1= x1-xr ;     alpha=-c1*z1+dxr ;

z2= x2-alpha ;

d_alpha=-c1*x2+c1*dxr+ddxr ;

tau=M*(-c2*z2-z1+d_alpha)+C*x2+G ;   


sys(1)=tau(1);
sys(2)=tau(2);
sys(3)=tau(3);