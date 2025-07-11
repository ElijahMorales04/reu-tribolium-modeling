%5-Cup Dispersal w/ LPA Model
%5th cup connects all 4 cups, making an X shape

b = 11.68;
mu_a = 0.8;
mu_l = 0.5129;
c_ea = 0.0110;
c_el = 0.0110;
c_pa = 0.0178;
gamma1 = 0.001; % chance of adult leaving cup 1
gamma2 = 0.001; % chance of adult leaving cup 2
gamma3 = 0.001; % chance of adult leaving cup 3
gamma4 = 0.001; % chance of adult leaving cup 4
gamma5 = 0.001; % chance of adult leaving cup 5

MaxT = 52;

DF = zeros(MaxT+1,4)
DF(:,1) = 0:MaxT
DF(1,:) = [0, 0, 0, 0]
DF2 = zeros(MaxT+1,4)
DF2(:,1) = 0:MaxT
DF2(1,:) = [0, 0, 0, 0]
DF3 = zeros(MaxT+1,4)
DF3(:,1) = 0:MaxT
DF3(1,:) = [0, 0, 0, 0]
DF4 = zeros(MaxT+1,4)
DF4(:,1) = 0:MaxT
DF4(1,:) = [0, 0, 0, 0]
DF5 = zeros(MaxT+1,4)
DF5(:, 1) = 0:MaxT
DF5(1,:) = [0, 0, 0, 50]


for i = 2:MaxT+1
    
    DF(i,2)= b*DF(i-1,4)*exp(-c_ea*DF(i-1,4)-c_el*DF(i-1,2));
    DF(i,3)= (1-mu_l)*DF(i-1,2);
    DF(i,4) = DF(i-1,3)*exp(-c_pa*DF(i-1,4))+(1-mu_a)*(1-gamma1)*DF(i-1,4)+ gamma5*DF5(i-1,4);
    
    DF2(i,2)= b*DF2(i-1,4)*exp(-c_ea*DF2(i-1,4)-c_el*DF2(i-1,2));
    DF2(i,3)= (1-mu_l)*DF2(i-1,2);
    DF2(i,4) = DF2(i-1,3)*exp(-c_pa*DF2(i-1,4))+(1-mu_a)*(1-gamma2)*DF2(i-1,4)+ gamma5*DF5(i-1,4);
    
    DF3(i,2)= b*DF3(i-1,4)*exp(-c_ea*DF3(i-1,4)-c_el*DF3(i-1,2));
    DF3(i,3)= (1-mu_l)*DF3(i-1,2);
    DF3(i,4) = DF3(i-1,3)*exp(-c_pa*DF3(i-1,4))+(1-mu_a)*(1-gamma3)*DF3(i-1,4)+ gamma5*DF5(i-1,4);
    
    DF4(i,2)= b*DF4(i-1,4)*exp(-c_ea*DF4(i-1,4)-c_el*DF4(i-1,2));
    DF4(i,3)= (1-mu_l)*DF4(i-1,2);
    DF4(i,4) = DF4(i-1,3)*exp(-c_pa*DF4(i-1,4))+(1-mu_a)*(1-gamma4)*DF4(i-1,4)+ gamma5*DF5(i-1,4);

    DF5(i,2)= b*DF5(i-1,4)*exp(-c_ea*DF5(i-1,4)-c_el*DF5(i-1,2));
    DF5(i,3)= (1-mu_l)*DF5(i-1,2);
    DF5(i,4) = DF5(i-1,3)*exp(-c_pa*DF5(i-1,4))+(1-mu_a)*(1-gamma5)*DF5(i-1,4)+ gamma1*DF(i-1,4)+ gamma2*DF2(i-1,4)+ gamma3*DF3(i-1,4)+ gamma4*DF4(i-1,4);

end

subplot(3,2,1);
x1 = DF(:,1);
y1 = DF(:,2:4);
plot(x1,y1)
title('Patch 1')
xlabel('Time(Weeks)')
ylabel('Population')
ax = gca;
comm.ax = ax.YLim;

subplot(3,2,2); 
x2 = DF2(:,1);
y2 = DF2(:,2:4);
plot(x2,y2)
title('Patch 2')
ax = gca;
ax.YLim = comm.ax;

subplot(3,2,3);
x3 = DF3(:,1);
y3 = DF3(:,2:4);
plot(x3,y3)
title('Patch 3')
ax = gca;
ax.YLim = comm.ax;

subplot(3,2,4); 
x4 = DF4(:,1);
y4 = DF4(:,2:4);
plot(x4,y4)
title('Patch 4')
ax = gca;
ax.YLim = comm.ax;

subplot(3,2,5); 
x5 = DF5(:,1);
y5 = DF5(:,2:4);
plot(x5,y5)
title('Patch 5')
ax = gca;
ax.YLim = comm.ax;

fig = gcf;
lgnd = legend('Larvae','Pupae','Adults')
lgnd.Position(1) = 0.65;
lgnd.Position(2) = 0.2;