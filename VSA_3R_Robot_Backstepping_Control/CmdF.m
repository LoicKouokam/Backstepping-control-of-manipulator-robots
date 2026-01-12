function [sys,x0,str,ts] = CmdF(t,x,u,flag)
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
sizes.NumContStates  = 12;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 12;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;
sys = simsizes(sizes);
x0  = [0,0,0,0,0,0,0,0,0,0,0,0];
str = [];
ts  = [0 0];

function sys=mdlDerivatives(t,x,u) 
%% CF parameters pulsation and atenuation
wn=[400 0 0; 0 400 0; 0 0 400] ;
zeta=[0.6 0 0; 0 0.6 0; 0 0 0.6] ;

%% command filter%%
%%% u=[alpha1,alpha2]; x=[x1c,dx1c,x2c,dx2c] (en realite dx1c=wn*dx1c);
alpha1=u(1:3)  ; alpha2=u(4:6) ;

x1c=x(1:3) ;  dx1c=x(4:6) ;

x2c=x(7:9) ;  dx2c=x(10:12) ;

diff_x1c=wn*dx1c ;
diff_diff_x1c=-2*zeta*wn*dx1c-wn*(x1c-alpha1) ;

diff_x2c=wn*dx2c ;
diff_diff_x2c=-2*zeta*wn*dx2c-wn*(x2c-alpha2) ;

sys(1)= diff_x1c(1) ;
sys(2)= diff_x1c(2) ;
sys(3)= diff_x1c(3) ;
sys(4)= diff_diff_x1c(1) ;
sys(5)= diff_diff_x1c(2) ;
sys(6)= diff_diff_x1c(3) ;

sys(7)= diff_x2c(1) ;
sys(8)= diff_x2c(2) ;
sys(9)= diff_x2c(3) ;
sys(10)= diff_diff_x2c(1) ;
sys(11)= diff_diff_x2c(2) ;
sys(12)= diff_diff_x2c(3) ;



%%


function sys=mdlOutputs(t,x,u)
% outputs y=[x1c,dx1c,x2c,dx2c];
wn=[400 0 0; 0 400 0; 0 0 400] ;

sys(1)=x(1);
sys(2)=x(2);
sys(3)=x(3);

sys(4:6)=(wn)*x(4:6);

sys(7)=x(7);
sys(8)=x(8);
sys(9)=x(9);

sys(10:12)=(wn)*x(10:12);


