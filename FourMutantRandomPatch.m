%4-Cup Dispersal w/ LPA Model
%adding direct cost to dispersal

bakari = 1000;
b = 20;
mu_a = 0.0842;
mu_l = 0.6053;
c_ea = 0.0179;
c_el = 0.0003;
c_pa = 0.0;
gammaR = 0.3; % chance of adult leaving any patch
gammaA = rand; % chance of mutant adult leaving any patch
gammaB = rand; % chance of mutant adult leaving any patch
gammaC = rand; % chance of mutant adult leaving any patch

Eps = 0.1; % chance of adult dying from dispersal
k = 4;

index = 1;

MaxT = 20;

DF1i = zeros(MaxT,12);
DF1i(1,1) = 67;
DF1i(1,5) = 26;
DF1i(1,9) = 236;

DF2i = zeros(MaxT,12);
DF2i(1,1) = 67;
DF2i(1,5) = 26;
DF2i(1,9) = 236;

DF3i = zeros(MaxT,12);
DF3i(1,1) = 67;
DF3i(1,5) = 26;
DF3i(1,9) = 236;

DF4i = zeros(MaxT,12);
DF4i(1,1) = 67;
DF4i(1,5) = 26;
DF4i(1,9) = 236;

for i = 0:MaxT
    done = 'a';
    % disp([num2str(i/MaxT),'%'])
    for j = i*20+1:20*(i+1)
        %PATCH 1
        DF1i(j+1,1)= b*DF1i(j,9)*exp(-c_ea*(sum(DF1i(j,2*k+1:3*k))-c_el*(sum(DF1i(j,1:1*k)))));
        DF1i(j+1,2)= b*DF1i(j,10)*exp(-c_ea*(sum(DF1i(j,2*k+1:3*k))-c_el*(sum(DF1i(j,1:1*k)))));
        DF1i(j+1,3)= b*DF1i(j,11)*exp(-c_ea*(sum(DF1i(j,2*k+1:3*k))-c_el*(sum(DF1i(j,1:1*k)))));
        DF1i(j+1,4)= b*DF1i(j,12)*exp(-c_ea*(sum(DF1i(j,2*k+1:3*k))-c_el*(sum(DF1i(j,1:1*k)))));

        DF1i(j+1,5)= (1-mu_l)*(DF1i(j,1));
        DF1i(j+1,6)= (1-mu_l)*(DF1i(j,2));
        DF1i(j+1,7)= (1-mu_l)*(DF1i(j,3));
        DF1i(j+1,8)= (1-mu_l)*(DF1i(j,4));

        DF1i(j+1,9) = DF1i(j,5)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF1i(j,9)+1/2*(DF2i(j,2*k+1)+DF3i(j,2*k+1))*gammaR*(1-mu_a)*(1-Eps);
        DF1i(j+1,10) = DF1i(j,6)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF1i(j,10)+1/2*(DF2i(j,2*k+2)+DF3i(j,2*k+2))*gammaA*(1-mu_a)*(1-Eps);
        DF1i(j+1,11) = DF1i(j,7)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF1i(j,11)+1/2*(DF2i(j,2*k+3)+DF3i(j,2*k+3))*gammaB*(1-mu_a)*(1-Eps);
        DF1i(j+1,12) = DF1i(j,8)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF1i(j,12)+1/2*(DF2i(j,2*k+4)+DF3i(j,2*k+4))*gammaC*(1-mu_a)*(1-Eps);

        %PATCH 2
        DF2i(j+1,1)= b*DF2i(j,9)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));
        DF2i(j+1,2)= b*DF2i(j,10)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));
        DF2i(j+1,3)= b*DF2i(j,11)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));
        DF2i(j+1,4)= b*DF2i(j,12)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));

        DF2i(j+1,5)= (1-mu_l)*(DF2i(j,1));
        DF2i(j+1,6)= (1-mu_l)*(DF2i(j,2));
        DF2i(j+1,7)= (1-mu_l)*(DF2i(j,3));
        DF2i(j+1,8)= (1-mu_l)*(DF2i(j,4));

        DF2i(j+1,9) = DF2i(j,5)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF2i(j,9)+1/2*(DF4i(j,2*k+1)+DF1i(j,2*k+1))*gammaR*(1-mu_a)*(1-Eps);
        DF2i(j+1,10) = DF2i(j,6)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF2i(j,10)+1/2*(DF4i(j,2*k+2)+DF1i(j,2*k+2))*gammaA*(1-mu_a)*(1-Eps);
        DF2i(j+1,11) = DF2i(j,7)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF2i(j,11)+1/2*(DF4i(j,2*k+3)+DF1i(j,2*k+3))*gammaB*(1-mu_a)*(1-Eps);
        DF2i(j+1,12) = DF2i(j,8)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF2i(j,12)+1/2*(DF4i(j,2*k+4)+DF1i(j,2*k+4))*gammaC*(1-mu_a)*(1-Eps);

        %PATCH 3
        DF3i(j+1,1)= b*DF3i(j,9)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));
        DF3i(j+1,2)= b*DF3i(j,10)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));
        DF3i(j+1,3)= b*DF3i(j,11)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));
        DF3i(j+1,4)= b*DF3i(j,12)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));

        DF3i(j+1,5)= (1-mu_l)*(DF3i(j,1));
        DF3i(j+1,6)= (1-mu_l)*(DF3i(j,2));
        DF3i(j+1,7)= (1-mu_l)*(DF3i(j,3));
        DF3i(j+1,8)= (1-mu_l)*(DF3i(j,4));

        DF3i(j+1,9) = DF3i(j,5)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF3i(j,9)+1/2*(DF4i(j,2*k+1)+DF1i(j,2*k+1))*gammaR*(1-mu_a)*(1-Eps);
        DF3i(j+1,10) = DF3i(j,6)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF3i(j,10)+1/2*(DF4i(j,2*k+2)+DF1i(j,2*k+2))*gammaA*(1-mu_a)*(1-Eps);
        DF3i(j+1,11) = DF3i(j,7)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF3i(j,11)+1/2*(DF4i(j,2*k+3)+DF1i(j,2*k+3))*gammaB*(1-mu_a)*(1-Eps);
        DF3i(j+1,12) = DF3i(j,8)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF3i(j,12)+1/2*(DF4i(j,2*k+4)+DF1i(j,2*k+4))*gammaC*(1-mu_a)*(1-Eps);

        %PATCH 4
        DF4i(j+1,1)= b*DF4i(j,9)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));
        DF4i(j+1,2)= b*DF4i(j,10)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));
        DF4i(j+1,3)= b*DF4i(j,11)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));
        DF4i(j+1,4)= b*DF4i(j,12)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));

        DF4i(j+1,5)= (1-mu_l)*(DF4i(j,1));
        DF4i(j+1,6)= (1-mu_l)*(DF4i(j,2));
        DF4i(j+1,7)= (1-mu_l)*(DF4i(j,3));
        DF4i(j+1,8)= (1-mu_l)*(DF4i(j,4));

        DF4i(j+1,9) = DF4i(j,5)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF4i(j,9)+1/2*(DF2i(j,2*k+1)+DF3i(j,2*k+1))*gammaR*(1-mu_a)*(1-Eps);
        DF4i(j+1,10) = DF4i(j,6)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF4i(j,10)+1/2*(DF2i(j,2*k+2)+DF3i(j,2*k+2))*gammaA*(1-mu_a)*(1-Eps);
        DF4i(j+1,11) = DF4i(j,7)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF4i(j,11)+1/2*(DF2i(j,2*k+3)+DF3i(j,2*k+3))*gammaB*(1-mu_a)*(1-Eps);
        DF4i(j+1,12) = DF4i(j,8)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF4i(j,12)+1/2*(DF2i(j,2*k+4)+DF3i(j,2*k+4))*gammaC*(1-mu_a)*(1-Eps);
    end
    randPatchExtinct = randi(4);
    if randPatchExtinct == 1;
        DF1i(20*(i+1)+1,:) = 0;
    elseif randPatchExtinct == 2;
        DF2i(20*(i+1)+1,:) = 0;
    elseif randPatchExtinct == 3;
        DF3i(20*(i+1)+1,:) = 0;
    elseif randPatchExtinct == 4;
        DF4i(20*(i+1)+1,:) = 0;
    end
    
    while done == 'a';
        randPatch = randi(4);
        while randPatchExtinct == randPatch;
            randPatch = randi(4);
        end
            
        if randPatch == 1 && index < 4;
            if index == 1;
                DF1i(20*(i+1)+1,10) = DF1i(20*(i+1),10)+bakari;
                index = index+1;
                done = 'b';
            elseif index == 2;
                DF1i(20*(i+1)+1,11) = DF1i(20*(i+1),11)+bakari;
                index = index+1;
                done = 'b';
            elseif index == 3;
                DF1i(20*(i+1)+1,12) = DF1i(20*(i+1),12)+bakari;
                index = index+1;
                done = 'b';
            end
            
        elseif randPatch == 2 && index < 4;
            if index == 1;
                DF2i(20*(i+1)+1,10) = DF2i(20*(i+1),10)+bakari;
                index = index+1;
                done = 'b';
            elseif index == 2;
                DF2i(20*(i+1)+1,11) = DF2i(20*(i+1),11)+bakari;
                index = index+1;
                done = 'b';
            elseif index == 3;
                DF2i(20*(i+1)+1,12) = DF2i(20*(i+1),12)+bakari;
                index = index+1;
                done = 'b';
            end
            
        elseif randPatch == 3 && index < 4;
            if index == 1;
                DF3i(20*(i+1)+1,10) = DF3i(20*(i+1),10)+bakari;
                index = index+1;
                done = 'b';
            elseif index == 2;
                DF3i(20*(i+1)+1,11) = DF3i(20*(i+1),11)+bakari;
                index = index+1;
                done = 'b';
            elseif index == 3;
                DF3i(20*(i+1)+1,12) = DF3i(20*(i+1),12)+bakari;
                index = index+1;
                done = 'b';
            end
            
        elseif randPatch == 4 && index < 4;
            if index == 1;
                DF4i(20*(i+1)+1,10) = DF4i(20*(i+1),10)+bakari;
                index = index+1;
                done = 'b';
            elseif index == 2;
                DF4i(20*(i+1)+1,11) = DF4i(20*(i+1),11)+bakari;
                index = index+1;
                done = 'b';
            elseif index == 3;
                DF4i(20*(i+1)+1,12) = DF4i(20*(i+1),12)+bakari;
                index = index+1;
                done = 'b';
            end
        else
            done = 'b';
        end
    end
    
end


%Proportion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = linspace(1,421,421);
LPAr1sum = DF1i(:,1) + DF1i(:,5) + DF1i(:,9);
LPAm1asum = DF1i(:,2) + DF1i(:,6) + DF1i(:,10);
LPAm1bsum = DF1i(:,3) + DF1i(:,7) + DF1i(:,11);
LPAm1csum = DF1i(:,4) + DF1i(:,8) + DF1i(:,12);
subplot(2,2,1);
LPAr1prop = LPAr1sum./(LPAr1sum+LPAm1asum+LPAm1bsum+LPAm1csum);
LPAm1aprop = LPAm1asum./(LPAr1sum+LPAm1asum+LPAm1bsum+LPAm1csum);
LPAm1bprop = LPAm1bsum./(LPAr1sum+LPAm1asum+LPAm1bsum+LPAm1csum);
LPAm1cprop = LPAm1csum./(LPAr1sum+LPAm1asum+LPAm1bsum+LPAm1csum);
patch(x,LPAr1prop,'red')

LPAr2sum = DF2i(:,1) + DF2i(:,5) + DF2i(:,9);
LPAm2asum = DF2i(:,2) + DF2i(:,6) + DF2i(:,10);
LPAm2bsum = DF2i(:,3) + DF2i(:,7) + DF2i(:,11);
LPAm2csum = DF2i(:,4) + DF2i(:,8) + DF2i(:,12);
subplot(2,2,2);
LPAr2prop = LPAr2sum./(LPAr2sum+LPAm2asum+LPAm2bsum+LPAm2csum);
LPAm2aprop = LPAm2asum./(LPAr2sum+LPAm2asum+LPAm2bsum+LPAm2csum);
LPAm2bprop = LPAm2bsum./(LPAr2sum+LPAm2asum+LPAm2bsum+LPAm2csum);
LPAm2cprop = LPAm2csum./(LPAr2sum+LPAm2asum+LPAm2bsum+LPAm2csum);

LPAr3sum = DF3i(:,1) + DF3i(:,5) + DF3i(:,9);
LPAm3asum = DF3i(:,2) + DF3i(:,6) + DF3i(:,10);
LPAm3bsum = DF3i(:,3) + DF3i(:,7) + DF3i(:,11);
LPAm3csum = DF3i(:,4) + DF3i(:,8) + DF3i(:,12);
subplot(2,2,3);
LPAr3prop = LPAr3sum./(LPAr3sum+LPAm3asum+LPAm3bsum+LPAm3csum);
LPAm3aprop = LPAm3asum./(LPAr3sum+LPAm3asum+LPAm3bsum+LPAm3csum);
LPAm3bprop = LPAm3bsum./(LPAr3sum+LPAm3asum+LPAm3bsum+LPAm3csum);
LPAm3cprop = LPAm3csum./(LPAr3sum+LPAm3asum+LPAm3bsum+LPAm3csum);

LPAr4sum = DF4i(:,1) + DF4i(:,5) + DF4i(:,9);
LPAm4asum = DF4i(:,2) + DF4i(:,6) + DF4i(:,10);
LPAm4bsum = DF4i(:,3) + DF4i(:,7) + DF4i(:,11);
LPAm4csum = DF4i(:,4) + DF4i(:,8) + DF4i(:,12);
subplot(2,2,4);
LPAr4prop = LPAr4sum./(LPAr4sum+LPAm4asum+LPAm4bsum+LPAm4csum);
LPAm4aprop = LPAm4asum./(LPAr4sum+LPAm4asum+LPAm4bsum+LPAm4csum);
LPAm4bprop = LPAm4bsum./(LPAr4sum+LPAm4asum+LPAm4bsum+LPAm4csum);
LPAm4cprop = LPAm4csum./(LPAr4sum+LPAm4asum+LPAm4bsum+LPAm4csum);

lgnd = legend(strcat('GammaR = ',num2str(gammaR)),strcat('GammaA = ',num2str(gammaA)),...
    strcat('GammaB = ',num2str(gammaB)),strcat('GammaC = ',num2str(gammaC)),'NumColumns',4)
lgnd.Position(1) = 0.002;
lgnd.Position(2) = 0.475;