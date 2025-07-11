%LPA Model.

b = 0.215;
mu_a = 0.0842;
mu_l = 0.6053;
c_ea = 0.0110;
c_el = 0.0093;
c_pa = 0.0178;

IC = [0, 0, 0.01];
MaxT = 40000;

DF = zeros(MaxT+1,4);
DF(:,1) = 0:MaxT;
DF(1,:) = [0, 0, 0, 0.01];

for i = 2:MaxT+1
    DF(i,2)= b*DF(i-1,4)*exp(-c_ea*DF(i-1,4)-c_el*DF(i-1,2));
    DF(i,3)= (1-mu_l)*DF(i-1,2);
    DF(i,4) = DF(i-1,3)*exp(-c_pa*DF(i-1,4))+DF(i-1,4)*(1-mu_a);
end

plot(DF(:,1),DF(:,2:4), linewidth=2)
title('b = 0.215, R_0 > 1')
xlabel('Time(Weeks)')
ylabel('Population')
lgnd = legend('Larvae','Pupae','Adults')
lgnd.FontSize = 15