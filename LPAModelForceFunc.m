%LPA Model. Changing medium, killing percentage of eggs.

b = 20;
mu_a = 0.0842;
mu_l = 0.6053;
mu_g = 0.359;
c_ea = 0.0179;
c_el = 0.0003;
c_pa = 0.0;

IC = [0, 0, 50];
MaxT = 52;

DF = zeros(MaxT+1,4);
DF(:,1) = 0:MaxT;
DF(1,:) = [0, 0, 0, 50];

for i = 2:MaxT+1
    if mod(i,4) == 0;
    
        DF(i,2)= (1-mu_g)*b*DF(i-1,4)*exp(-c_ea*DF(i-1,4)-c_el*DF(i-1,2));
        DF(i,3)= (1-mu_l)*DF(i-1,2);
        DF(i,4) = DF(i-1,3)*exp(-c_pa*DF(i-1,4))+DF(i-1,4)*(1-mu_a);
    else
        DF(i,2)= b*DF(i-1,4)*exp(-c_ea*DF(i-1,4)-c_el*DF(i-1,2));
        DF(i,3)= (1-mu_l)*DF(i-1,2);
        DF(i,4) = DF(i-1,3)*exp(-c_pa*DF(i-1,4))+DF(i-1,4)*(1-mu_a);
    end
end

plot(DF(:,1),DF(:,2), 'LineWidth',2,'Color','b')
hold on
plot(DF(:,1),DF(:,3), 'LineWidth',2,'Color','r')
plot(DF(:,1),DF(:,4), 'LineWidth',2,'Color','k')
title('LPA Model w/ Forcing Function')
xlabel('Time(Weeks)')
ylabel('Population')
lgnd = legend('Larvae','Pupae','Adults')
lgnd.FontSize = 15