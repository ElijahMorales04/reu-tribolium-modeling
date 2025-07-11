%4-Cup Dispersal w/ LPA Model
%adding direct cost to dispersal

bakari = 1000;
b = 20;
mu_a = 0.0842;
mu_l = 0.6053;
c_ea = 0.0179;
c_el = 0.0003;
c_pa = 0.0;
gammaR = 0.8; % chance of adult leaving any patch
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

DF5i = zeros(MaxT,12);
DF5i(1,1) = 67;
DF5i(1,5) = 26;
DF5i(1,9) = 236;

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

        DF1i(j+1,9) = DF1i(j,5)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF1i(j,9)+1/4*DF5i(j,2*k+1)*gammaR*(1-mu_a)*(1-Eps);
        DF1i(j+1,10) = DF1i(j,6)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF1i(j,10)+1/4*DF5i(j,2*k+2)*gammaA*(1-mu_a)*(1-Eps);
        DF1i(j+1,11) = DF1i(j,7)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF1i(j,11)+1/4*DF5i(j,2*k+3)*gammaB*(1-mu_a)*(1-Eps);
        DF1i(j+1,12) = DF1i(j,8)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF1i(j,12)+1/4*DF5i(j,2*k+4)*gammaC*(1-mu_a)*(1-Eps);

        %PATCH 2
        DF2i(j+1,1)= b*DF2i(j,9)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));
        DF2i(j+1,2)= b*DF2i(j,10)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));
        DF2i(j+1,3)= b*DF2i(j,11)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));
        DF2i(j+1,4)= b*DF2i(j,12)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));

        DF2i(j+1,5)= (1-mu_l)*(DF2i(j,1));
        DF2i(j+1,6)= (1-mu_l)*(DF2i(j,2));
        DF2i(j+1,7)= (1-mu_l)*(DF2i(j,3));
        DF2i(j+1,8)= (1-mu_l)*(DF2i(j,4));

        DF2i(j+1,9) = DF2i(j,5)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF2i(j,9)+1/4*DF5i(j,2*k+1)*gammaR*(1-mu_a)*(1-Eps);
        DF2i(j+1,10) = DF2i(j,6)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF2i(j,10)+1/4*DF5i(j,2*k+2)*gammaA*(1-mu_a)*(1-Eps);
        DF2i(j+1,11) = DF2i(j,7)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF2i(j,11)+1/4*DF5i(j,2*k+3)*gammaB*(1-mu_a)*(1-Eps);
        DF2i(j+1,12) = DF2i(j,8)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF2i(j,12)+1/4*DF5i(j,2*k+4)*gammaC*(1-mu_a)*(1-Eps);

        %PATCH 3
        DF3i(j+1,1)= b*DF3i(j,9)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));
        DF3i(j+1,2)= b*DF3i(j,10)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));
        DF3i(j+1,3)= b*DF3i(j,11)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));
        DF3i(j+1,4)= b*DF3i(j,12)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));

        DF3i(j+1,5)= (1-mu_l)*(DF3i(j,1));
        DF3i(j+1,6)= (1-mu_l)*(DF3i(j,2));
        DF3i(j+1,7)= (1-mu_l)*(DF3i(j,3));
        DF3i(j+1,8)= (1-mu_l)*(DF3i(j,4));

        DF3i(j+1,9) = DF3i(j,5)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF3i(j,9)+1/4*DF5i(j,2*k+1)*gammaR*(1-mu_a)*(1-Eps);
        DF3i(j+1,10) = DF3i(j,6)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF3i(j,10)+1/4*DF5i(j,2*k+2)*gammaA*(1-mu_a)*(1-Eps);
        DF3i(j+1,11) = DF3i(j,7)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF3i(j,11)+1/4*DF5i(j,2*k+3)*gammaB*(1-mu_a)*(1-Eps);
        DF3i(j+1,12) = DF3i(j,8)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF3i(j,12)+1/4*DF5i(j,2*k+4)*gammaC*(1-mu_a)*(1-Eps);

        %PATCH 4
        DF4i(j+1,1)= b*DF4i(j,9)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));
        DF4i(j+1,2)= b*DF4i(j,10)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));
        DF4i(j+1,3)= b*DF4i(j,11)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));
        DF4i(j+1,4)= b*DF4i(j,12)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));

        DF4i(j+1,5)= (1-mu_l)*(DF4i(j,1));
        DF4i(j+1,6)= (1-mu_l)*(DF4i(j,2));
        DF4i(j+1,7)= (1-mu_l)*(DF4i(j,3));
        DF4i(j+1,8)= (1-mu_l)*(DF4i(j,4));

        DF4i(j+1,9) = DF4i(j,5)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF4i(j,9)+1/4*DF5i(j,2*k+1)*gammaR*(1-mu_a)*(1-Eps);
        DF4i(j+1,10) = DF4i(j,6)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF4i(j,10)+1/4*DF5i(j,2*k+2)*gammaA*(1-mu_a)*(1-Eps);
        DF4i(j+1,11) = DF4i(j,7)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF4i(j,11)+1/4*DF5i(j,2*k+3)*gammaB*(1-mu_a)*(1-Eps);
        DF4i(j+1,12) = DF4i(j,8)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF4i(j,12)+1/4*DF5i(j,2*k+4)*gammaC*(1-mu_a)*(1-Eps);
    
        %PATCH 5
        DF5i(j+1,1)= b*DF5i(j,9)*exp(-c_ea*(sum(DF5i(j,2*k+1:3*k))-c_el*(sum(DF5i(j,1:1*k)))));
        DF5i(j+1,2)= b*DF5i(j,10)*exp(-c_ea*(sum(DF5i(j,2*k+1:3*k))-c_el*(sum(DF5i(j,1:1*k)))));
        DF5i(j+1,3)= b*DF5i(j,11)*exp(-c_ea*(sum(DF5i(j,2*k+1:3*k))-c_el*(sum(DF5i(j,1:1*k)))));
        DF5i(j+1,4)= b*DF5i(j,12)*exp(-c_ea*(sum(DF5i(j,2*k+1:3*k))-c_el*(sum(DF5i(j,1:1*k)))));

        DF5i(j+1,5)= (1-mu_l)*(DF5i(j,1));
        DF5i(j+1,6)= (1-mu_l)*(DF5i(j,2));
        DF5i(j+1,7)= (1-mu_l)*(DF5i(j,3));
        DF5i(j+1,8)= (1-mu_l)*(DF5i(j,4));

        DF5i(j+1,9) = DF5i(j,5)*exp(-c_pa*(sum(DF5i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF5i(j,9)+(DF1i(j,2*k+1)+DF2i(j,2*k+1)+DF3i(j,2*k+4)+DF4i(j,2*k+4))*gammaR*(1-mu_a)*(1-Eps);
        DF5i(j+1,10) = DF5i(j,6)*exp(-c_pa*(sum(DF5i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF5i(j,10)+(DF1i(j,2*k+2)+DF2i(j,2*k+2)+DF3i(j,2*k+4)+DF4i(j,2*k+4))*gammaA*(1-mu_a)*(1-Eps);
        DF5i(j+1,11) = DF5i(j,7)*exp(-c_pa*(sum(DF5i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF5i(j,11)+(DF1i(j,2*k+3)+DF2i(j,2*k+3)+DF3i(j,2*k+4)+DF4i(j,2*k+4))*gammaB*(1-mu_a)*(1-Eps);
        DF5i(j+1,12) = DF5i(j,8)*exp(-c_pa*(sum(DF5i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF5i(j,12)+(DF1i(j,2*k+4)+DF2i(j,2*k+4)+DF3i(j,2*k+4)+DF4i(j,2*k+4))*gammaC*(1-mu_a)*(1-Eps);
    end
    randPatchExtinct = randi(5);
    if randPatchExtinct == 1;
        DF1i(20*(i+1)+1,:) = 0;
    elseif randPatchExtinct == 2;
        DF2i(20*(i+1)+1,:) = 0;
    elseif randPatchExtinct == 3;
        DF3i(20*(i+1)+1,:) = 0;
    elseif randPatchExtinct == 4;
        DF4i(20*(i+1)+1,:) = 0;
    elseif randPatchExtinct == 5;
        DF5i(20*(i+1)+1,:) = 0;
    end
    
    while done == 'a';
        randPatch = randi(5);
        while randPatchExtinct == randPatch;
            randPatch = randi(5);
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
            
        elseif randPatch == 5 && index < 4;
            if index == 1;
                DF5i(20*(i+1)+1,10) = DF5i(20*(i+1),10)+bakari;
                index = index+1;
                done = 'b';
            elseif index == 2;
                DF5i(20*(i+1)+1,11) = DF5i(20*(i+1),11)+bakari;
                index = index+1;
                done = 'b';
            elseif index == 3;
                DF5i(20*(i+1)+1,12) = DF5i(20*(i+1),12)+bakari;
                index = index+1;
                done = 'b';
            end
        else
            done = 'b';
        end
    end
    
end

%Proportion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = [1:421 421:-1:1 0]';
bottom = zeros([1 421])';
LPAr1sum = DF1i(:,1) + DF1i(:,5) + DF1i(:,9);
LPAm1asum = DF1i(:,2) + DF1i(:,6) + DF1i(:,10);
LPAm1bsum = DF1i(:,3) + DF1i(:,7) + DF1i(:,11);
LPAm1csum = DF1i(:,4) + DF1i(:,8) + DF1i(:,12);
subplot(3,2,1);
LPAr1prop = LPAr1sum./(LPAr1sum+LPAm1asum+LPAm1bsum+LPAm1csum);
LPAr1prop(isnan(LPAr1prop)) = 0;
LPAm1aprop = LPAm1asum./(LPAr1sum+LPAm1asum+LPAm1bsum+LPAm1csum);
LPAm1aprop(isnan(LPAm1aprop)) = 0;
LPAm1bprop = LPAm1bsum./(LPAr1sum+LPAm1asum+LPAm1bsum+LPAm1csum);
LPAm1bprop(isnan(LPAm1bprop)) = 0;
LPAm1cprop = LPAm1csum./(LPAr1sum+LPAm1asum+LPAm1bsum+LPAm1csum);
LPAm1cprop(isnan(LPAm1cprop)) = 0;
yr1 = [LPAr1prop+LPAm1aprop+LPAm1cprop+LPAm1bprop; LPAm1aprop(421:-1:1,1); 1];
patch(x,yr1,'blue')
hold on
ym1a = [LPAm1aprop+LPAm1cprop+LPAm1bprop; LPAm1bprop(421:-1:1,1); LPAm1aprop(1,1)];
patch(x,ym1a,'red')
ym1b = [LPAm1cprop+LPAm1bprop; LPAm1cprop(421:-1:1,1); LPAm1bprop(1,1)];
patch(x,ym1b,'yellow')
ym1c = [LPAm1cprop; bottom(421:-1:1,1); LPAm1cprop(1,1)];
patch(x,ym1c,'magenta')
title('Patch 1')
ax= gca;
ax.XLim = [0, 421];
currentax = ax.XLim;
hold off


LPAr2sum = DF2i(:,1) + DF2i(:,5) + DF2i(:,9);
LPAm2asum = DF2i(:,2) + DF2i(:,6) + DF2i(:,10);
LPAm2bsum = DF2i(:,3) + DF2i(:,7) + DF2i(:,11);
LPAm2csum = DF2i(:,4) + DF2i(:,8) + DF2i(:,12);
subplot(3,2,2);
LPAr2prop = LPAr2sum./(LPAr2sum+LPAm2asum+LPAm2bsum+LPAm2csum);
LPAr2prop(isnan(LPAr2prop)) = 0;
LPAm2aprop = LPAm2asum./(LPAr2sum+LPAm2asum+LPAm2bsum+LPAm2csum);
LPAm2aprop(isnan(LPAm2aprop)) = 0;
LPAm2bprop = LPAm2bsum./(LPAr2sum+LPAm2asum+LPAm2bsum+LPAm2csum);
LPAm2bprop(isnan(LPAm2bprop)) = 0;
LPAm2cprop = LPAm2csum./(LPAr2sum+LPAm2asum+LPAm2bsum+LPAm2csum);
LPAm2cprop(isnan(LPAm2cprop)) = 0;
yr2 = [LPAr2prop+LPAm2aprop+LPAm2cprop+LPAm2bprop; LPAm2aprop(421:-1:1,1); 1];
patch(x,yr2,'blue')
hold on
ym2a = [LPAm2aprop+LPAm2cprop+LPAm2bprop; LPAm2bprop(421:-1:1,1); LPAm2aprop(1,1)];
patch(x,ym2a,'red')
ym2b = [LPAm2cprop+LPAm2bprop; LPAm2cprop(421:-1:1,1); LPAm2bprop(1,1)];
patch(x,ym2b,'yellow')
ym2c = [LPAm2cprop; bottom(421:-1:1,1); LPAm2cprop(1,1)];
patch(x,ym2c,'magenta')
title('Patch 2')
ax= gca;
ax.XLim = currentax;
hold off


LPAr3sum = DF3i(:,1) + DF3i(:,5) + DF3i(:,9);
LPAm3asum = DF3i(:,2) + DF3i(:,6) + DF3i(:,10);
LPAm3bsum = DF3i(:,3) + DF3i(:,7) + DF3i(:,11);
LPAm3csum = DF3i(:,4) + DF3i(:,8) + DF3i(:,12);
ax3 = subplot(3,2,3);
LPAr3prop = LPAr3sum./(LPAr3sum+LPAm3asum+LPAm3bsum+LPAm3csum);
LPAr3prop(isnan(LPAr3prop)) = 0;
LPAm3aprop = LPAm3asum./(LPAr3sum+LPAm3asum+LPAm3bsum+LPAm3csum);
LPAm3aprop(isnan(LPAm3aprop)) = 0;
LPAm3bprop = LPAm3bsum./(LPAr3sum+LPAm3asum+LPAm3bsum+LPAm3csum);
LPAm3bprop(isnan(LPAm3bprop)) = 0;
LPAm3cprop = LPAm3csum./(LPAr3sum+LPAm3asum+LPAm3bsum+LPAm3csum);
LPAm3cprop(isnan(LPAm3cprop)) = 0;
yr3 = [LPAr3prop+LPAm3aprop+LPAm3cprop+LPAm3bprop; LPAm3aprop(421:-1:1,1); 1];
patch(x,yr3,'blue')
hold on
ym3a = [LPAm3aprop+LPAm3cprop+LPAm3bprop; LPAm3bprop(421:-1:1,1); LPAm3aprop(1,1)];
patch(x,ym3a,'red')
ym3b = [LPAm3cprop+LPAm3bprop; LPAm3cprop(421:-1:1,1); LPAm3bprop(1,1)];
patch(x,ym3b,'yellow')
ym3c = [LPAm3cprop; bottom(421:-1:1,1); LPAm3cprop(1,1)];
patch(x,ym3c,'magenta')
title('Patch 3')
ax= gca;
ax.XLim = currentax;
ax3.YLabel.String = 'Proportion';
hold off

LPAr4sum = DF4i(:,1) + DF4i(:,5) + DF4i(:,9);
LPAm4asum = DF4i(:,2) + DF4i(:,6) + DF4i(:,10);
LPAm4bsum = DF4i(:,3) + DF4i(:,7) + DF4i(:,11);
LPAm4csum = DF4i(:,4) + DF4i(:,8) + DF4i(:,12);
subplot(3,2,4);
LPAr4prop = LPAr4sum./(LPAr4sum+LPAm4asum+LPAm4bsum+LPAm4csum);
LPAr4prop(isnan(LPAr4prop)) = 0;
LPAm4aprop = LPAm4asum./(LPAr4sum+LPAm4asum+LPAm4bsum+LPAm4csum);
LPAm4aprop(isnan(LPAm4aprop)) = 0;
LPAm4bprop = LPAm4bsum./(LPAr4sum+LPAm4asum+LPAm4bsum+LPAm4csum);
LPAm4bprop(isnan(LPAm4bprop)) = 0;
LPAm4cprop = LPAm4csum./(LPAr4sum+LPAm4asum+LPAm4bsum+LPAm4csum);
LPAm4cprop(isnan(LPAm4cprop)) = 0;
yr4 = [LPAr4prop+LPAm4aprop+LPAm4cprop+LPAm4bprop; LPAm4aprop(421:-1:1,1); 1];
patch(x,yr4,'blue')
hold on
ym4a = [LPAm4aprop+LPAm4cprop+LPAm4bprop; LPAm4bprop(421:-1:1,1); LPAm4aprop(1,1)];
patch(x,ym4a,'red')
ym4b = [LPAm4cprop+LPAm4bprop; LPAm4cprop(421:-1:1,1); LPAm4bprop(1,1)];
patch(x,ym4b,'yellow')
ym4c = [LPAm4cprop; bottom(421:-1:1,1); LPAm4cprop(1,1)];
patch(x,ym4c,'magenta')
title('Patch 4')
ax= gca;
ax.XLim = currentax;
hold off

LPAr5sum = DF5i(:,1) + DF5i(:,5) + DF5i(:,9);
LPAm5asum = DF5i(:,2) + DF5i(:,6) + DF5i(:,10);
LPAm5bsum = DF5i(:,3) + DF5i(:,7) + DF5i(:,11);
LPAm5csum = DF5i(:,4) + DF5i(:,8) + DF5i(:,12);
ax5 = subplot(3,2,5);
LPAr5prop = LPAr5sum./(LPAr5sum+LPAm5asum+LPAm5bsum+LPAm5csum);
LPAr5prop(isnan(LPAr5prop)) = 0;
LPAm5aprop = LPAm5asum./(LPAr5sum+LPAm5asum+LPAm5bsum+LPAm5csum);
LPAm5aprop(isnan(LPAm5aprop)) = 0;
LPAm5bprop = LPAm5bsum./(LPAr5sum+LPAm5asum+LPAm5bsum+LPAm5csum);
LPAm5bprop(isnan(LPAm5bprop)) = 0;
LPAm5cprop = LPAm5csum./(LPAr5sum+LPAm5asum+LPAm5bsum+LPAm5csum);
LPAm5cprop(isnan(LPAm5cprop)) = 0;
yr5 = [LPAr5prop+LPAm5aprop+LPAm5cprop+LPAm5bprop; LPAm5aprop(421:-1:1,1); 1];
patch(x,yr5,'blue')
hold on
ym5a = [LPAm5aprop+LPAm5cprop+LPAm5bprop; LPAm5bprop(421:-1:1,1); LPAm5aprop(1,1)];
patch(x,ym5a,'red')
ym5b = [LPAm5cprop+LPAm5bprop; LPAm5cprop(421:-1:1,1); LPAm5bprop(1,1)];
patch(x,ym5b,'yellow')
ym5c = [LPAm5cprop; bottom(421:-1:1,1); LPAm5cprop(1,1)];
patch(x,ym5c,'magenta')
title('Patch 5')
ax= gca;
ax.XLim = currentax;
ax5.XLabel.String = 'Time(Bi-weekly)';
hold off

lgnd = legend(strcat('\gamma_{R} = ',num2str(gammaR)),strcat('\gamma_{A} = ',num2str(gammaA)),...
    strcat('\gamma_{B} = ',num2str(gammaB)),strcat('\gamma_{C} = ',num2str(gammaC)),'NumColumns',1)
lgnd.Position(1) = 0.6;
lgnd.Position(2) = 0.15;
