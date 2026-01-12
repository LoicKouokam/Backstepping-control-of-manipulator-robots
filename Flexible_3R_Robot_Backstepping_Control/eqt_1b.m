function [sys,x0,str,ts] = eqt_1b(t,x,u,flag)

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
sizes.NumContStates  = 6;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 6;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 0;
sys = simsizes(sizes);
x0= [0 0 0 0 0 0 ];
str = [];
ts  = [];

function sys=mdlDerivatives(t,x,u)
 
% Parametres robot

l2 =0.760; l3=0.930; d=0.28;
m1=5; m2=5; m3=3; lc1=0.22; lc2=0.51; lc3=0.67; 
Ix1=.34 ; Iy1=.36 ; Iz1=.31; Ix2=0.18; Iy2=1.32; Iz2=1.31; 
Ix3=0.07; Iy3=0.92; Iz3=0.93; g0=9.8;

% % etats  

q= x(1:3)  ;  q1=q(1)  ; q2=q(2) ;   q3=q(3)  ;

dq= x(4:6)  ;  dq1= dq(1) ; dq2=dq(2) ; dq3=dq(3) ;

the=u(1:3) ;

sig=u(4:6) ;  sig1=sig(1)  ;  sig2=sig(2)  ; sig3=sig(3) ;


% Inertie
M=zeros(3);

M(1,1)= -2*m3*lc3*l3-2*m3*sin(q2)*sin(q3)*l3^2*cos(q2)*cos(q3)-2*m3*lc3^2*sin(q2)*sin(q3)*cos(q2)*cos(q3)-4*m3*lc3*cos(q2)^2*cos(q3)^2*l3+cos(q2)^2*Iy2+2*m3*lc3^2*cos(q2)^2*cos(q3)^2+2*m3*cos(q2)^2*cos(q3)^2*l3^2+4*m3*lc3*sin(q2)*sin(q3)*cos(q2)*cos(q3)*l3+m3*l3^2+m3*lc3^2+2*m3*cos(q2)*l2*lc3*sin(q2)*sin(q3)-2*m3*cos(q2)*l2*sin(q2)*sin(q3)*l3-2*m3*cos(q2)^2*l2*lc3*cos(q3)+2*m3*cos(q2)^2*l2*cos(q3)*l3-2*m2*cos(q2)^2*lc2*l2+m2*cos(q2)^2*lc2^2+m2*cos(q2)^2*l2^2+m3*cos(q2)^2*l2^2+2*Iy3*cos(q2)^2*cos(q3)^2+2*Ix3*cos(q2)*sin(q3)*sin(q2)*cos(q3)-2*Iy3*cos(q2)*sin(q3)*sin(q2)*cos(q3)-m3*lc3^2*cos(q3)^2-m3*l3^2*cos(q2)^2-m3*lc3^2*cos(q2)^2-2*Ix3*cos(q2)^2*cos(q3)^2-m3*l3^2*cos(q3)^2+Ix3*cos(q3)^2-Ix2*cos(q2)^2-Iy3*cos(q3)^2-Iy3*cos(q2)^2+Ix3*cos(q2)^2+2*m3*lc3*cos(q3)^2*l3+2*m3*lc3*cos(q2)^2*l3+Iy3+Iy1+Ix2 ;
M(1,2)= 0 ;

M(1,3)= 0 ;

M(2,1)= 0 ;

M(2,2)= Iz3+Iz2+m3*lc3^2+m3*l3^2-2*m3*lc3*l3-2*m3*l2*lc3*cos(q3)+2*m3*l2*cos(q3)*l3-2*m2*lc2*l2+m2*lc2^2+m2*l2^2+m3*l2^2 ;

M(2,3)= -2*m3*lc3*l3+Iz3+m3*l3^2+m3*lc3^2-m3*l2*lc3*cos(q3)+m3*l2*cos(q3)*l3 ;

M(3,1)= 0 ;

M(3,2)=Iz3+m3*lc3^2+m3*l3^2-m3*l2*lc3*cos(q3)-2*m3*lc3*l3+m3*l2*cos(q3)*l3 ;

M(3,3)= Iz3-2*m3*lc3*l3+m3*lc3^2+m3*l3^2 ;
% Gravite

G=zeros(3,1);

G(1)= 0 ;

G(2)= g0*(m2*lc2*cos(q2)+m3*cos(q2)*l2+m3*lc3*cos(q2+q3)) ;

G(3)= g0*m3*lc3*cos(q2+q3) ;

% Coriolis

C=zeros(3);

C(1,1)= -2*dq2*m3*cos(q2)*l2*cos(q3)*l3*sin(q2)+2*dq2*m3*cos(q2)*l2*lc3*cos(q3)*sin(q2)+4*dq2*m3*lc3*cos(q2)*cos(q3)^2*l3*sin(q2)+4*dq2*m3*lc3*cos(q2)^2*sin(q3)*cos(q3)*l3-dq3*m3*cos(q2)*l2*cos(q3)*l3*sin(q2)+dq3*m3*cos(q2)*l2*lc3*cos(q3)*sin(q2)+4*dq3*m3*lc3*cos(q2)*cos(q3)^2*l3*sin(q2)+4*dq3*m3*lc3*cos(q2)^2*sin(q3)*cos(q3)*l3-2*dq2*m3*cos(q2)^2*l2*sin(q3)*l3-2*dq2*m3*lc3^2*cos(q2)^2*sin(q3)*cos(q3)+2*dq2*m3*cos(q2)^2*l2*lc3*sin(q3)-2*dq2*m3*lc3*cos(q2)*l3*sin(q2)-2*dq2*m3*cos(q2)^2*sin(q3)*l3^2*cos(q3)+2*dq2*m2*cos(q2)*lc2*l2*sin(q2)-2*dq2*m3*cos(q2)*cos(q3)^2*l3^2*sin(q2)-2*dq2*m3*lc3^2*cos(q2)*cos(q3)^2*sin(q2)-2*dq2*m3*lc3*sin(q3)*cos(q3)*l3-dq3*m3*cos(q2)^2*l2*sin(q3)*l3-2*dq3*m3*lc3^2*cos(q2)^2*sin(q3)*cos(q3)+dq3*m3*cos(q2)^2*l2*lc3*sin(q3)-2*dq3*m3*lc3*cos(q2)*l3*sin(q2)-2*dq3*m3*cos(q2)^2*sin(q3)*l3^2*cos(q3)-2*dq3*m3*cos(q2)*cos(q3)^2*l3^2*sin(q2)-2*dq3*m3*lc3^2*cos(q2)*cos(q3)^2*sin(q2)-2*dq3*m3*lc3*sin(q3)*cos(q3)*l3-dq2*cos(q2)*Iy2*sin(q2)+dq2*Iy3*cos(q2)*sin(q2)+dq2*Ix2*cos(q2)*sin(q2)-dq2*Ix3*cos(q2)*sin(q2)-dq2*Ix3*sin(q3)*cos(q3)+dq2*Iy3*sin(q3)*cos(q3)+dq3*Iy3*cos(q2)*sin(q2)-dq3*Ix3*cos(q2)*sin(q2)-dq3*Ix3*sin(q3)*cos(q3)+dq3*Iy3*sin(q3)*cos(q3)-dq2*m2*cos(q2)*l2^2*sin(q2)-2*dq2*Iy3*cos(q2)^2*sin(q3)*cos(q3)-dq2*m2*cos(q2)*lc2^2*sin(q2)+2*dq2*Ix3*cos(q2)^2*sin(q3)*cos(q3)+dq2*m3*l2*sin(q3)*l3+dq2*m3*l3^2*cos(q2)*sin(q2)+dq2*m3*lc3^2*cos(q2)*sin(q2)+2*dq2*Ix3*cos(q2)*cos(q3)^2*sin(q2)+dq2*m3*lc3^2*sin(q3)*cos(q3)-dq2*m3*l2*lc3*sin(q3)-dq2*m3*cos(q2)*l2^2*sin(q2)-2*dq2*Iy3*cos(q2)*cos(q3)^2*sin(q2)+dq2*m3*sin(q3)*l3^2*cos(q3)-2*dq3*Iy3*cos(q2)^2*sin(q3)*cos(q3)+2*dq3*Ix3*cos(q2)^2*sin(q3)*cos(q3)+dq3*m3*l3^2*cos(q2)*sin(q2)+dq3*m3*lc3^2*cos(q2)*sin(q2)+2*dq3*Ix3*cos(q2)*cos(q3)^2*sin(q2)+dq3*m3*lc3^2*sin(q3)*cos(q3)-2*dq3*Iy3*cos(q2)*cos(q3)^2*sin(q2)+dq3*m3*sin(q3)*l3^2*cos(q3) ;

C(1,2)= dq1*(2*m2*cos(q2)*lc2*l2*sin(q2)-2*m3*cos(q2)*cos(q3)^2*l3^2*sin(q2)-2*m3*lc3^2*cos(q2)*cos(q3)^2*sin(q2)-2*m3*lc3*cos(q2)*l3*sin(q2)-2*Iy3*cos(q2)^2*sin(q3)*cos(q3)-Ix3*sin(q3)*cos(q3)+2*Ix3*cos(q2)^2*sin(q3)*cos(q3)+Iy3*cos(q2)*sin(q2)-2*m3*cos(q2)*l2*cos(q3)*l3*sin(q2)+2*m3*cos(q2)*l2*lc3*cos(q3)*sin(q2)-2*m3*lc3*sin(q3)*cos(q3)*l3+m3*sin(q3)*l3^2*cos(q3)+m3*l2*sin(q3)*l3-m3*l2*lc3*sin(q3)+m3*lc3^2*sin(q3)*cos(q3)+Iy3*sin(q3)*cos(q3)+4*m3*lc3*cos(q2)*cos(q3)^2*l3*sin(q2)-2*m3*lc3^2*cos(q2)^2*sin(q3)*cos(q3)+2*m3*cos(q2)^2*l2*lc3*sin(q3)-2*m3*cos(q2)^2*l2*sin(q3)*l3-2*m3*cos(q2)^2*sin(q3)*l3^2*cos(q3)+2*Ix3*cos(q2)*cos(q3)^2*sin(q2)+m3*l3^2*cos(q2)*sin(q2)+m3*lc3^2*cos(q2)*sin(q2)-2*Iy3*cos(q2)*cos(q3)^2*sin(q2)-m3*cos(q2)*l2^2*sin(q2)-m2*cos(q2)*l2^2*sin(q2)-m2*cos(q2)*lc2^2*sin(q2)+4*m3*lc3*cos(q2)^2*sin(q3)*cos(q3)*l3-cos(q2)*Iy2*sin(q2)+Ix2*cos(q2)*sin(q2)-Ix3*cos(q2)*sin(q2));

C(1,3)= dq1*(m3*lc3^2*cos(q3)*sin(q3)+m3*l3^2*cos(q3)*sin(q3)+2*Ix3*cos(q2)^2*cos(q3)*sin(q3)+2*Ix3*sin(q2)*cos(q3)^2*cos(q2)-2*Iy3*sin(q2)*cos(q3)^2*cos(q2)-2*Iy3*cos(q2)^2*cos(q3)*sin(q3)-Ix3*cos(q3)*sin(q3)+Iy3*cos(q3)*sin(q3)-2*m3*lc3^2*cos(q2)^2*cos(q3)*sin(q3)-2*m3*cos(q2)^2*cos(q3)*l3^2*sin(q3)-2*m3*sin(q2)*cos(q3)^2*l3^2*cos(q2)-2*m3*lc3^2*sin(q2)*cos(q3)^2*cos(q2)+m3*cos(q2)^2*l2*lc3*sin(q3)-m3*cos(q2)^2*l2*sin(q3)*l3-2*m3*lc3*l3*cos(q3)*sin(q3)+m3*sin(q2)*l3^2*cos(q2)+m3*lc3^2*sin(q2)*cos(q2)+4*m3*lc3*sin(q2)*cos(q3)^2*cos(q2)*l3+4*m3*lc3*cos(q2)^2*cos(q3)*l3*sin(q3)+m3*cos(q2)*l2*lc3*sin(q2)*cos(q3)-m3*cos(q2)*l2*sin(q2)*cos(q3)*l3-2*m3*lc3*sin(q2)*cos(q2)*l3-Ix3*sin(q2)*cos(q2)+Iy3*sin(q2)*cos(q2)) ;

C(2,1)= -dq1*(-2*m3*lc3*sin(q3)*cos(q3)*l3+2*m3*cos(q2)^2*l2*lc3*sin(q3)-2*m3*cos(q2)^2*l2*sin(q3)*l3+4*m3*lc3*cos(q2)*cos(q3)^2*l3*sin(q2)+2*m3*cos(q2)*l2*lc3*cos(q3)*sin(q2)-2*m3*cos(q2)*l2*cos(q3)*l3*sin(q2)-2*m3*lc3^2*cos(q2)^2*sin(q3)*cos(q3)-2*m3*cos(q2)^2*sin(q3)*l3^2*cos(q3)-Ix3*cos(q2)*sin(q2)+Iy3*cos(q2)*sin(q2)-cos(q2)*Iy2*sin(q2)+Ix2*cos(q2)*sin(q2)-Ix3*cos(q3)*sin(q3)+Iy3*cos(q3)*sin(q3)+m3*lc3^2*cos(q2)*sin(q2)-m2*cos(q2)*lc2^2*sin(q2)-m2*cos(q2)*l2^2*sin(q2)-m3*cos(q2)*l2^2*sin(q2)+2*Ix3*cos(q2)*cos(q3)^2*sin(q2)-2*Iy3*cos(q2)*cos(q3)^2*sin(q2)+m3*l3^2*cos(q2)*sin(q2)-m3*l2*lc3*sin(q3)+m3*l2*sin(q3)*l3+4*m3*lc3*cos(q2)^2*sin(q3)*cos(q3)*l3+2*Ix3*cos(q2)^2*cos(q3)*sin(q3)-2*Iy3*cos(q2)^2*cos(q3)*sin(q3)+2*m2*cos(q2)*lc2*l2*sin(q2)-2*m3*lc3*cos(q2)*l3*sin(q2)-2*m3*lc3^2*cos(q2)*cos(q3)^2*sin(q2)-2*m3*cos(q2)*cos(q3)^2*l3^2*sin(q2)+m3*sin(q3)*l3^2*cos(q3)+m3*lc3^2*sin(q3)*cos(q3)) ;

C(2,2)= -dq3*m3*l2*sin(q3)*(-lc3+l3) ;

C(2,3)= -m3*l2*sin(q3)*(-lc3+l3)*(dq2+dq3) ;

C(3,1)= -dq1*(2*Ix3*cos(q3)*cos(q2)^2*sin(q3)+m3*cos(q3)*l3^2*sin(q3)+m3*lc3^2*cos(q3)*sin(q3)-2*Iy3*cos(q2)^2*cos(q3)*sin(q3)+2*Ix3*sin(q2)*cos(q3)^2*cos(q2)-2*Iy3*sin(q2)*cos(q3)^2*cos(q2)-Ix3*cos(q3)*sin(q3)+Iy3*cos(q3)*sin(q3)-2*m3*sin(q2)*cos(q3)^2*l3^2*cos(q2)-2*m3*lc3^2*sin(q2)*cos(q3)^2*cos(q2)-2*m3*lc3*cos(q3)*l3*sin(q3)-2*m3*lc3^2*cos(q2)^2*cos(q3)*sin(q3)-2*m3*cos(q2)^2*cos(q3)*l3^2*sin(q3)+m3*cos(q2)^2*l2*lc3*sin(q3)-m3*cos(q2)^2*l2*sin(q3)*l3+m3*sin(q2)*l3^2*cos(q2)+m3*lc3^2*sin(q2)*cos(q2)+m3*cos(q2)*l2*lc3*sin(q2)*cos(q3)-m3*cos(q2)*l2*sin(q2)*cos(q3)*l3+4*m3*lc3*cos(q2)^2*cos(q3)*l3*sin(q3)+4*m3*lc3*sin(q2)*cos(q3)^2*cos(q2)*l3-2*m3*lc3*sin(q2)*cos(q2)*l3-Ix3*sin(q2)*cos(q2)+Iy3*sin(q2)*cos(q2)) ;

C(3,2)= dq2*m3*l2*sin(q3)*(-lc3+l3) ;

C(3,3)= 0 ;

% % Rigidite

KE1= 1.62*10^-4 ;
delta=0.015 ; 
n= 0.006 ;

R= zeros(3);

R(1,1)=(KE1*sig1^2)/(delta-n*sig1)^2 ;

R(2,2)= (KE1*sig2^2)/(delta-n*sig2)^2 ;

R(3,3)= (KE1*sig3^2)/(delta-n*sig3)^2 ;


% % 


diff_diff_q= M\(R*(the-q)-C*dq-G) ;

% % 
sys(1)=dq(1);
sys(2)=dq(2);
sys(3)=dq(3);
sys(4)=diff_diff_q(1);
sys(5)=diff_diff_q(2);
sys(6)=diff_diff_q(3);


function sys=mdlOutputs(t,x,u)
sys(1)=x(1);
sys(2)=x(2);
sys(3)=x(3);

sys(4)=x(4);
sys(5)=x(5);
sys(6)=x(6);