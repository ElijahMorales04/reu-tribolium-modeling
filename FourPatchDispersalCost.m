%4-Cup Dispersal w/ LPA Model
%adding direct cost to dispersal

b = 0.00004;
mu_a = 0.0842;
mu_l = 0.6053;
c_ea = 0.0003;
c_el = 0.0179;
c_pa = 0.0;
gamma = 0.3; % chance of adult leaving any patch
Eps = 0.1; % chance of adult dying from dispersal

MaxT = 1000;

DF = zeros(MaxT+1,4);
DF(:,1) = 0:MaxT;
DF(1,:) = [0, 0, 0, 0.01];
DF2 = zeros(MaxT+1,4);
DF2(:,1) = 0:MaxT;
DF2(1,:) = [0, 0, 0, 0];
DF3 = zeros(MaxT+1,4);
DF3(:,1) = 0:MaxT;
DF3(1,:) = [0, 0, 0, 0];
DF4 = zeros(MaxT+1,4);
DF4(:,1) = 0:MaxT;
DF4(1,:) = [0, 0, 0, 0];

for i = 2:MaxT+1
    DF(i,2)= b*DF(i-1,4)*exp(-c_ea*DF(i-1,4)-c_el*DF(i-1,2));
    DF(i,3)= (1-mu_l)*DF(i-1,2);
    DF(i,4) = DF(i-1,3)*exp(-c_pa*DF(i-1,4))+((1-mu_a)*(1-gamma))*DF(i-1,4)+(DF2(i-1,4)+DF3(i-1,4))*gamma*(1-mu_a)*(1-Eps);
    
    DF2(i,2)= b*DF2(i-1,4)*exp(-c_ea*DF2(i-1,4)-c_el*DF2(i-1,2));
    DF2(i,3)= (1-mu_l)*DF2(i-1,2);
    DF2(i,4) = DF2(i-1,3)*exp(-c_pa*DF2(i-1,4))+((1-mu_a)*(1-gamma))*DF2(i-1,4)+(DF4(i-1,4)+DF(i-1,4))*gamma*(1-mu_a)*(1-Eps);
    
    DF3(i,2)= b*DF3(i-1,4)*exp(-c_ea*DF3(i-1,4)-c_el*DF3(i-1,2));
    DF3(i,3)= (1-mu_l)*DF3(i-1,2);
    DF3(i,4) = DF3(i-1,3)*exp(-c_pa*DF3(i-1,4))+((1-mu_a)*(1-gamma))*DF3(i-1,4)+(DF(i-1,4)+DF4(i-1,4))*gamma*(1-mu_a)*(1-Eps);
    
    DF4(i,2)= b*DF4(i-1,4)*exp(-c_ea*DF4(i-1,4)-c_el*DF4(i-1,2));
    DF4(i,3)= (1-mu_l)*DF4(i-1,2);
    DF4(i,4) = DF4(i-1,3)*exp(-c_pa*DF4(i-1,4))+((1-mu_a)*(1-gamma))*DF4(i-1,4)+(DF3(i-1,4)+DF2(i-1,4))*gamma*(1-mu_a)*(1-Eps);
end

figure(2);
subplot(2,2,1);
x1 = DF(:,1);
y1 = DF(:,2:4);
plot(x1,y1)
title('Patch 1')
xlabel('Time(Weeks)')
ylabel('Population')
ax = gca;
comm.ax = ax.YLim;

subplot(2,2,2); 
x2 = DF2(:,1);
y2 = DF2(:,2:4);
plot(x2,y2)
title('Patch 2')
ax = gca;
ax.YLim = comm.ax;

subplot(2,2,3);
x3 = DF3(:,1);
y3 = DF3(:,2:4);
plot(x3,y3)
title('Patch 3')
ax = gca;
ax.YLim = comm.ax;

subplot(2,2,4); 
x4 = DF4(:,1);
y4 = DF4(:,2:4);
plot(x4,y4)
title('Patch 4')
ax = gca;
ax.YLim = comm.ax;
