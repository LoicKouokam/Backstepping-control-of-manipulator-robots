close all;


figure(1);

subplot(211);
plot(t2,error_q(:,1),'k',t2,error_q(:,2),'r:',t2,error_q(:,3),'g-','linewidth',2) ;
xlabel('time(s)');ylabel('q-qdes[rad]');
legend('joint 1','joint 2','joint 3');

subplot(212);
plot(t2,error_x3(:,1),'k',t2,error_x3(:,2),'g:',t2,error_x3(:,3),'b-','linewidth',2);
xlabel('time(s)');ylabel('x3-x3des[rad]');

figure(2);
subplot(211);
plot(t2,Uthe(:,1),'r',t2,Uthe(:,2),'g:',t2,Uthe(:,3),'b-','linewidth',2);
xlabel('time(s)');ylabel('Uthe');
legend('joint 1','joint 2','joint 3');

subplot(212);
plot(t2,Usig(:,1),'r',t2,Usig(:,2),'g:',t2,Usig(:,3),'b-','linewidth',2);
xlabel('time(s)');ylabel('Usig');
legend('joint 1','joint 2','joint 3');

figure(3) 
plot(t2,rho_est,'r','linewidth',2)
xlabel('time(s)');ylabel('rho');