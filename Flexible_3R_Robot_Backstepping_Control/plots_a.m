close all;


figure(1);

plot(t,error_q(:,1),'k',t,error_q(:,2),'r:',t,error_q(:,3),'g-','linewidth',2) ;
xlabel('time(s)');ylabel('q-qdes[rad]');
legend('joint 1','joint 2','joint 3');


figure(2)
plot(t,U(:,1),'r',t,U(:,2),'g:',t,U(:,3),'b-','linewidth',2);
xlabel('time(s)');ylabel('U');
legend('joint 1','joint 2','joint 3');