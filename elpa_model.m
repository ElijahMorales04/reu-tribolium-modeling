t0 = 0;
tfinal = 52;
p0 = [0; 0; 0; 50];

[t,p] = ode45(@ELPA_ODE,[t0 tfinal],p0);
plot(t,p)
title('E L P A Populations Over Time')
xlabel('t')
ylabel('Population')
legend('Eggs','Larvae','Pupae','Adutls')