%4-Cup Dispersal w/ LPA Model
%5th Patch connects all 4 patches, making an X shape
%adding direct cost to dispersal

b = 20;
mu_a = 0.0842;
mu_l = 0.6053;
c_ea = 0.0179;
c_el = 0.0003;
c_pa = 0.0;
gammaR = 0.3; % chance of adult leaving any patch
gammaM = 0.1; % chance of mutant adult leaving any patch
Eps = 0.1; % chance of adult dying from dispersal

MaxT = 200;

DF1 = zeros(MaxT+1,4);
DF1(:,1) = 0:MaxT;
DF1(1,:) = [0, 0, 0, 0];
DF2 = zeros(MaxT+1,4);
DF2(:,1) = 0:MaxT;
DF2(1,:) = [0, 0, 0, 0];
DF3 = zeros(MaxT+1,4);
DF3(:,1) = 0:MaxT;
DF3(1,:) = [0, 0, 0, 0];
DF4 = zeros(MaxT+1,4);
DF4(:,1) = 0:MaxT;
DF4(1,:) = [0, 0, 0, 0];
DF5 = zeros(MaxT+1,4)
DF5(:,1) = 0:MaxT
DF5(1,:) = [0, 0, 0, 50];

DF1M = zeros(MaxT+1,4);
DF1M(:,1) = 0:MaxT;
DF1M(1,:) = [0, 0, 0, 0];
DF2M = zeros(MaxT+1,4);
DF2M(:,1) = 0:MaxT;
DF2M(1,:) = [0, 0, 0, 0];
DF3M = zeros(MaxT+1,4);
DF3M(:,1) = 0:MaxT;
DF3M(1,:) = [0, 0, 0, 0];
DF4M = zeros(MaxT+1,4);
DF4M(:,1) = 0:MaxT;
DF4M(1,:) = [0, 0, 0, 0];
DF5M = zeros(MaxT+1,4)
DF5M(:,1) = 0:MaxT
DF5M(1,:) = [0, 0, 0, 1];

for i = 2:MaxT+1
    DF1(i,2)= b*DF1(i-1,4)*exp(-c_ea*(DF1(i-1,4)+DF1M(i-1,4))-c_el*(DF1(i-1,2)+DF1M(i-1,2)));
    DF1(i,3)= (1-mu_l)*DF1(i-1,2);
    DF1(i,4) = DF1(i-1,3)*exp(-c_pa*(DF1(i-1,4)+DF1M(i-1,4)))+((1-mu_a)*(1-gammaR))*DF1(i-1,4)+(DF5(i-1,4))*gammaR*(1-mu_a)*(1-Eps);
    
    DF2(i,2)= b*DF2(i-1,4)*exp(-c_ea*(DF2(i-1,4)+DF2M(i-1,4))-c_el*(DF2(i-1,2)+DF2M(i-1,2)));
    DF2(i,3)= (1-mu_l)*DF2(i-1,2);
    DF2(i,4) = DF2(i-1,3)*exp(-c_pa*(DF2(i-1,4)+DF2M(i-1,4)))+((1-mu_a)*(1-gammaR))*DF2(i-1,4)+(DF5(i-1,4))*gammaR*(1-mu_a)*(1-Eps);
    
    DF3(i,2)= b*DF3(i-1,4)*exp(-c_ea*(DF3(i-1,4)+DF3M(i-1,4))-c_el*(DF3(i-1,2)+DF3M(i-1,2)));
    DF3(i,3)= (1-mu_l)*DF3(i-1,2);
    DF3(i,4) = DF3(i-1,3)*exp(-c_pa*(DF3(i-1,4)+DF3M(i-1,4)))+((1-mu_a)*(1-gammaR))*DF3(i-1,4)+(DF5(i-1,4))*gammaR*(1-mu_a)*(1-Eps);
    
    DF4(i,2)= b*DF4(i-1,4)*exp(-c_ea*(DF4(i-1,4)+DF4M(i-1,4))-c_el*(DF4(i-1,2)+DF4M(i-1,2)));
    DF4(i,3)= (1-mu_l)*DF4(i-1,2);
    DF4(i,4) = DF4(i-1,3)*exp(-c_pa*(DF4(i-1,4)+DF4M(i-1,4)))+((1-mu_a)*(1-gammaR))*DF4(i-1,4)+(DF5(i-1,4))*gammaR*(1-mu_a)*(1-Eps);
    
    DF5(i,2)= b*DF5(i-1,4)*exp(-c_ea*(DF5(i-1,4)+DF5M(i-1,4))-c_el*(DF5(i-1,2)+DF5M(i-1,2)));
    DF5(i,3)= (1-mu_l)*DF5(i-1,2);
    DF5(i,4) = DF5(i-1,3)*exp(-c_pa*(DF5(i-1,4)+DF5M(i-1,4)))+((1-mu_a)*(1-gammaR))*DF5(i-1,4)+(1/4*DF1(i-1,4)+1/4*DF2(i-1,4)+1/4*DF3(i-1,4)+1/4*DF4(i-1,4))*gammaR*(1-mu_a)*(1-Eps);
    
    DF1M(i,2)= b*DF1M(i-1,4)*exp(-c_ea*(DF1(i-1,4)+DF1M(i-1,4))-c_el*(DF1(i-1,2)+DF1M(i-1,2)));
    DF1M(i,3)= (1-mu_l)*DF1M(i-1,2);
    DF1M(i,4) = DF1M(i-1,3)*exp(-c_pa*(DF1(i-1,4)+DF1M(i-1,4)))+((1-mu_a)*(1-gammaM))*DF1M(i-1,4)+(DF5M(i-1,4))*gammaM*(1-mu_a)*(1-Eps);
    
    DF2M(i,2)= b*DF2M(i-1,4)*exp(-c_ea*(DF2(i-1,4)+DF2M(i-1,4))-c_el*(DF2(i-1,2)+DF2M(i-1,2)));
    DF2M(i,3)= (1-mu_l)*DF2M(i-1,2);
    DF2M(i,4) = DF2M(i-1,3)*exp(-c_pa*(DF2(i-1,4)+DF2M(i-1,4)))+((1-mu_a)*(1-gammaM))*DF2M(i-1,4)+(DF5M(i-1,4))*gammaM*(1-mu_a)*(1-Eps);
    
    DF3M(i,2)= b*DF3M(i-1,4)*exp(-c_ea*(DF3(i-1,4)+DF3M(i-1,4))-c_el*(DF3(i-1,2)+DF3M(i-1,2)));
    DF3M(i,3)= (1-mu_l)*DF3M(i-1,2);
    DF3M(i,4) = DF3M(i-1,3)*exp(-c_pa*(DF3(i-1,4)+DF3M(i-1,4)))+((1-mu_a)*(1-gammaM))*DF3M(i-1,4)+(DF5M(i-1,4))*gammaM*(1-mu_a)*(1-Eps);
    
    DF4M(i,2)= b*DF4M(i-1,4)*exp(-c_ea*(DF4(i-1,4)+DF4M(i-1,4))-c_el*(DF4(i-1,2)+DF4M(i-1,2)));
    DF4M(i,3)= (1-mu_l)*DF4M(i-1,2);
    DF4M(i,4) = DF4M(i-1,3)*exp(-c_pa*(DF4(i-1,4)+DF4M(i-1,4)))+((1-mu_a)*(1-gammaM))*DF4M(i-1,4)+(DF5M(i-1,4))*gammaM*(1-mu_a)*(1-Eps);
    
    DF5M(i,2)= b*DF5M(i-1,4)*exp(-c_ea*(DF5(i-1,4)+DF5M(i-1,4))-c_el*(DF5(i-1,2)+DF5M(i-1,2)));
    DF5M(i,3)= (1-mu_l)*DF5M(i-1,2);
    DF5M(i,4) = DF5M(i-1,3)*exp(-c_pa*(DF5(i-1,4)+DF5M(i-1,4)))+((1-mu_a)*(1-gammaM))*DF5M(i-1,4)+(1/4*DF1M(i-1,4)+1/4*DF2M(i-1,4)+1/4*DF3M(i-1,4)+1/4*DF4M(i-1,4))*gammaM*(1-mu_a)*(1-Eps);
end

% figure(2);
% subplot(3,2,1);
% x1 = DF1(:,1);
% y1 = DF1(:,2:4);
% plot(x1,y1)
% title('Patch 1')
% xlabel('Time(Weeks)')
% ylabel('Population')
% ax = gca;
% comm.ax = ax.YLim;
% 
% subplot(3,2,2); 
% x2 = DF2(:,1);
% y2 = DF2(:,2:4);
% plot(x2,y2)
% title('Patch 2')
% ax = gca;
% ax.YLim = comm.ax;
% 
% subplot(3,2,3);
% x3 = DF3(:,1);
% y3 = DF3(:,2:4);
% plot(x3,y3)
% title('Patch 3')
% ax = gca;
% ax.YLim = comm.ax;
% 
% subplot(3,2,4); 
% x4 = DF4(:,1);
% y4 = DF4(:,2:4);
% plot(x4,y4)
% title('Patch 4')
% ax = gca;
% ax.YLim = comm.ax;
% 
% subplot(3,2,5); 
% x5 = DF5(:,1);
% y5 = DF5(:,2:4);
% plot(x5,y5)
% title('Patch 5')
% ax = gca;
% ax.YLim = comm.ax;
% 
% fig = gcf;
% lgnd = legend('Larvae','Pupae','Adults')
% lgnd.Position(1) = 0.65;
% lgnd.Position(2) = 0.2;

figure(1);
y1p = zeros(size(DF1));
subplot(3,2,1);
x1p = DF1M(:,1); %Proportion of population are mutants.
y1p(:,2:4) = DF1M(:,2:4)./(DF1M(:,2:4) + DF1(:,2:4));
y1p(isnan(y1p)) = 0;
plot(x1p,y1p(:,2:4), linewidth = 2)
title('Patch 1')

y2p = zeros(size(DF1));
subplot(3,2,2);
x2p = DF2M(:,1); %Proportion of population are mutants.
y2p(:,2:4) = DF2M(:,2:4)./(DF2M(:,2:4) + DF2(:,2:4));
y2p(isnan(y2p)) = 0;
plot(x2p,y2p(:,2:4), linewidth = 2)
title('Patch 2')

y3p = zeros(size(DF1));
subplot(3,2,3);
x3p = DF3M(:,1); %Proportion of population are mutants.
y3p(:,2:4) = DF3M(:,2:4)./(DF3M(:,2:4) + DF3(:,2:4));
y3p(isnan(y3p)) = 0;
plot(x3p,y3p(:,2:4), linewidth = 2)
title('Patch 3')
ylabel('Proportion of Mutants')

y4p = zeros(size(DF1));
subplot(3,2,4);
x4p = DF4M(:,1); %Proportion of population are mutants.
y4p(:,2:4) = DF4M(:,2:4)./(DF4M(:,2:4) + DF4(:,2:4));
y4p(isnan(y4p)) = 0;
plot(x4p,y4p(:,2:4), linewidth = 2)
title('Patch 4')

y5p = zeros(size(DF1));
subplot(3,2,5);
x5p = DF5M(:,1); %Proportion of population are mutants.
y5p(:,2:4) = DF5M(:,2:4)./(DF5M(:,2:4) + DF5(:,2:4));
y5p(isnan(y5p)) = 0;
plot(x5p,y5p(:,2:4), linewidth = 2)
title('Patch 5')
xlabel('Time(bi-weekly)')
lgnd = legend('Larvae','Pupae','Adults')

lgnd = legend('Larvae','Pupae','Adults')
lgnd.Position(1) = 0.65;
lgnd.Position(2) = 0.2;
lgnd.FontSize = 14;