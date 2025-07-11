%LPAA Model with noise.

b = 11.6772;
mu_a = 0.1108;
mu_l = 0.5129;
mu_p = rand;
c_1 = 0.0099;
c_2 = 0.0028;

IC = [0, 0, 0, 50];
n = 0;
MaxT = 52;

DF = zeros(MaxT+1,5)
DF(:,1) = 0:MaxT
DF(1,:) = [0, 0, 0, 0, 50]

for i = 2:MaxT
    DF(i,2) = b*DF(i-1,5)*exp(-c_1*DF(i-1,5));
    DF(i,3) = (1-mu_l)*DF(i-1,2);
    DF(i,4) = (1-mu_p)*DF(i-1,3);
    DF(i,5) = DF(i-1,4)*exp(-c_2*DF(i-1,5))+DF(i-1,5)*(1-mu_a);
end

plot(DF(:,1),DF(:,2:5))