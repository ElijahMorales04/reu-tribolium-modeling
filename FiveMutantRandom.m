%4-Cup Dispersal w/ LPA Model
%adding direct cost to dispersal

bakari = 0.01;
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

        DF1i(j+1,9) = DF1i(j,5)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF1i(j,9)+1/4*(DF5i(j,2*k+1))*gammaR*(1-mu_a)*(1-Eps);
        DF1i(j+1,10) = DF1i(j,6)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF1i(j,10)+1/4*(DF5i(j,2*k+2))*gammaA*(1-mu_a)*(1-Eps);
        DF1i(j+1,11) = DF1i(j,7)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF1i(j,11)+1/4*(DF5i(j,2*k+3))*gammaB*(1-mu_a)*(1-Eps);
        DF1i(j+1,12) = DF1i(j,8)*exp(-c_pa*(sum(DF1i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF1i(j,12)+1/4*(DF5i(j,2*k+4))*gammaC*(1-mu_a)*(1-Eps);

        %PATCH 2
        DF2i(j+1,1)= b*DF2i(j,9)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));
        DF2i(j+1,2)= b*DF2i(j,10)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));
        DF2i(j+1,3)= b*DF2i(j,11)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));
        DF2i(j+1,4)= b*DF2i(j,12)*exp(-c_ea*(sum(DF2i(j,2*k+1:3*k))-c_el*(sum(DF2i(j,1:1*k)))));

        DF2i(j+1,5)= (1-mu_l)*(DF2i(j,1));
        DF2i(j+1,6)= (1-mu_l)*(DF2i(j,2));
        DF2i(j+1,7)= (1-mu_l)*(DF2i(j,3));
        DF2i(j+1,8)= (1-mu_l)*(DF2i(j,4));

        DF2i(j+1,9) = DF2i(j,5)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF2i(j,9)+1/4*(DF5i(j,2*k+1))*gammaR*(1-mu_a)*(1-Eps);
        DF2i(j+1,10) = DF2i(j,6)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF2i(j,10)+1/4*(DF5i(j,2*k+2))*gammaA*(1-mu_a)*(1-Eps);
        DF2i(j+1,11) = DF2i(j,7)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF2i(j,11)+1/4*(DF5i(j,2*k+3))*gammaB*(1-mu_a)*(1-Eps);
        DF2i(j+1,12) = DF2i(j,8)*exp(-c_pa*(sum(DF2i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF2i(j,12)+1/4*(DF5i(j,2*k+4))*gammaC*(1-mu_a)*(1-Eps);

        %PATCH 3
        DF3i(j+1,1)= b*DF3i(j,9)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));
        DF3i(j+1,2)= b*DF3i(j,10)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));
        DF3i(j+1,3)= b*DF3i(j,11)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));
        DF3i(j+1,4)= b*DF3i(j,12)*exp(-c_ea*(sum(DF3i(j,2*k+1:3*k))-c_el*(sum(DF3i(j,1:1*k)))));

        DF3i(j+1,5)= (1-mu_l)*(DF3i(j,1));
        DF3i(j+1,6)= (1-mu_l)*(DF3i(j,2));
        DF3i(j+1,7)= (1-mu_l)*(DF3i(j,3));
        DF3i(j+1,8)= (1-mu_l)*(DF3i(j,4));

        DF3i(j+1,9) = DF3i(j,5)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF3i(j,9)+1/4*(DF5i(j,2*k+1))*gammaR*(1-mu_a)*(1-Eps);
        DF3i(j+1,10) = DF3i(j,6)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF3i(j,10)+1/4*(DF5i(j,2*k+2))*gammaA*(1-mu_a)*(1-Eps);
        DF3i(j+1,11) = DF3i(j,7)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF3i(j,11)+1/4*(DF5i(j,2*k+3))*gammaB*(1-mu_a)*(1-Eps);
        DF3i(j+1,12) = DF3i(j,8)*exp(-c_pa*(sum(DF3i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF3i(j,12)+1/4*(DF5i(j,2*k+4))*gammaC*(1-mu_a)*(1-Eps);

        %PATCH 4
        DF4i(j+1,1)= b*DF4i(j,9)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));
        DF4i(j+1,2)= b*DF4i(j,10)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));
        DF4i(j+1,3)= b*DF4i(j,11)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));
        DF4i(j+1,4)= b*DF4i(j,12)*exp(-c_ea*(sum(DF4i(j,2*k+1:3*k))-c_el*(sum(DF4i(j,1:1*k)))));

        DF4i(j+1,5)= (1-mu_l)*(DF4i(j,1));
        DF4i(j+1,6)= (1-mu_l)*(DF4i(j,2));
        DF4i(j+1,7)= (1-mu_l)*(DF4i(j,3));
        DF4i(j+1,8)= (1-mu_l)*(DF4i(j,4));

        DF4i(j+1,9) = DF4i(j,5)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaR))*DF4i(j,9)+1/4*(DF5i(j,2*k+1))*gammaR*(1-mu_a)*(1-Eps);
        DF4i(j+1,10) = DF4i(j,6)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaA))*DF4i(j,10)+1/4*(DF5i(j,2*k+2))*gammaA*(1-mu_a)*(1-Eps);
        DF4i(j+1,11) = DF4i(j,7)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaB))*DF4i(j,11)+1/4*(DF5i(j,2*k+3))*gammaB*(1-mu_a)*(1-Eps);
        DF4i(j+1,12) = DF4i(j,8)*exp(-c_pa*(sum(DF4i(j,2*k+1:3*k))))+((1-mu_a)*(1-gammaC))*DF4i(j,12)+1/4*(DF5i(j,2*k+4))*gammaC*(1-mu_a)*(1-Eps);
    
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


figure(2);
%PATCH 1
subplot(3,2,1);
x = linspace(1,421,421);
Lr1 = DF1i(:,1);
plot(x,Lr1,'r-', linewidth=2)
title('Patch 1')

hold on
Lm1a = DF1i(:,2);
plot(x,Lm1a,'r--',linewidth=2)
Lm1b = DF1i(:,3);
plot(x,Lm1b,'r-.',linewidth=2)
Lm1c = DF1i(:,4);
plot(x,Lm1c,'r:',linewidth=2)

Pr1 = DF1i(:,5);
plot(x,Pr1,'b-', linewidth=2)

Pm1a = DF1i(:,6);
plot(x,Pm1a,'b--',linewidth=2)
Pm1b = DF1i(:,7);
plot(x,Pm1b,'b-.',linewidth=2)
Pm1c = DF1i(:,8);
plot(x,Pm1c,'b:',linewidth=2)

Ar1 = DF1i(:,9);
plot(x,Ar1,'k-', linewidth=2)

Am1a = DF1i(:,10);
plot(x,Am1a,'k--',linewidth=2)
Am1b = DF1i(:,11);
plot(x,Am1b,'k-.',linewidth=2)
Am1c = DF1i(:,12);
plot(x,Am1c,'k:',linewidth=2)
hold off


%PATCH 2
subplot(3,2,2);
x = linspace(1,421,421);
Lr2 = DF1i(:,1);
plot(x,Lr2,'r-', linewidth=2)
title('Patch 2')

hold on
Lm2a = DF1i(:,2);
plot(x,Lm2a,'r--',linewidth=2)
Lm2b = DF1i(:,3);
plot(x,Lm2b,'r-.',linewidth=2)
Lm2c = DF1i(:,4);
plot(x,Lm2c,'r:',linewidth=2)

Pr2 = DF1i(:,5);
plot(x,Pr2,'b-', linewidth=2)

Pm2a = DF1i(:,6);
plot(x,Pm2a,'b--',linewidth=2)
Pm2b = DF1i(:,7);
plot(x,Pm2b,'b-.',linewidth=2)
Pm2c = DF1i(:,8);
plot(x,Pm2c,'b:',linewidth=2)

Ar2 = DF1i(:,9);
plot(x,Ar2,'k-', linewidth=2)

Am2a = DF1i(:,10);
plot(x,Am2a,'k--',linewidth=2)
Am2b = DF1i(:,11);
plot(x,Am2b,'k-.',linewidth=2)
Am3c = DF1i(:,12);
plot(x,Am3c,'k:',linewidth=2)
hold off

%PATCH 3
ax3 = subplot(3,2,3);
x = linspace(1,421,421);
Lr3 = DF1i(:,1);
plot(x,Lr3,'r-', linewidth=2)
title('Patch 3')
ax3.YLabel.String = 'Population';

hold on
Lm3a = DF1i(:,2);
plot(x,Lm3a,'r--',linewidth=2)
Lm3b = DF1i(:,3);
plot(x,Lm3b,'r-.',linewidth=2)
Lm3c = DF1i(:,4);
plot(x,Lm3c,'r:',linewidth=2)

Pr3 = DF1i(:,5);
plot(x,Pr3,'b-', linewidth=2)

Pm3a = DF1i(:,6);
plot(x,Pm3a,'b--',linewidth=2)
Pm3b = DF1i(:,7);
plot(x,Pm3b,'b-.',linewidth=2)
Pm3c = DF1i(:,8);
plot(x,Pm3c,'b:',linewidth=2)

Ar3 = DF1i(:,9);
plot(x,Ar3,'k-', linewidth=2)

Am3a = DF1i(:,10);
plot(x,Am3a,'k--',linewidth=2)
Am3b = DF1i(:,11);
plot(x,Am3b,'k-.',linewidth=2)
Am3c = DF1i(:,12);
plot(x,Am3c,'k:',linewidth=2)
hold off

%PATCH 4
subplot(3,2,4);
x = linspace(1,421,421);
Lr4 = DF1i(:,1);
plot(x,Lr4,'r-', linewidth=2)
title('Patch 4')

hold on
Lm4a = DF1i(:,2);
plot(x,Lm4a,'r--',linewidth=2)
Lm4b = DF1i(:,3);
plot(x,Lm4b,'r-.',linewidth=2)
Lm4c = DF1i(:,4);
plot(x,Lm4c,'r:',linewidth=2)

Pr4 = DF1i(:,5);
plot(x,Pr4,'b-', linewidth=2)

Pm4a = DF1i(:,6);
plot(x,Pm4a,'b--',linewidth=2)
Pm4b = DF1i(:,7);
plot(x,Pm4b,'b-.',linewidth=2)
Pm4c = DF1i(:,8);
plot(x,Pm4c,'b:',linewidth=2)

Ar4 = DF1i(:,9);
plot(x,Ar4,'k-', linewidth=2)

Am4a = DF1i(:,10);
plot(x,Am4a,'k--',linewidth=2)
Am4b = DF1i(:,11);
plot(x,Am4b,'k-.',linewidth=2)
Am4c = DF1i(:,12);
plot(x,Am4c,'k:',linewidth=2)
hold off

%PATCH 5
ax5 = subplot(3,2,5);
x = linspace(1,421,421);
Lr5 = DF5i(:,1);
plot(x,Lr5,'r-', linewidth=2)
title('Patch 5')
ax5.XLabel.String = 'Time(Bi-weekly)';

hold on
Lm5a = DF5i(:,2);
plot(x,Lm5a,'r--',linewidth=2)
Lm5b = DF5i(:,3);
plot(x,Lm5b,'r-.',linewidth=2)
Lm5c = DF5i(:,4);
plot(x,Lm5c,'r:',linewidth=2)

Pr5 = DF5i(:,5);
plot(x,Pr5,'b-', linewidth=2)

Pm5a = DF5i(:,6);
plot(x,Pm5a,'b--',linewidth=2)
Pm5b = DF5i(:,7);
plot(x,Pm5b,'b-.',linewidth=2)
Pm5c = DF5i(:,8);
plot(x,Pm5c,'b:',linewidth=2)

Ar5 = DF5i(:,9);
plot(x,Ar5,'k-', linewidth=2)

Am5a = DF5i(:,10);
plot(x,Am5a,'k--',linewidth=2)
Am5b = DF5i(:,11);
plot(x,Am5b,'k-.',linewidth=2)
Am5c = DF5i(:,12);
plot(x,Am5c,'k:',linewidth=2)
hold off


% lgnd = legend('Larvae','LarveM1','LarvaeM2','LarvaeM3','Pupae','PupaeM1','PupaeM2','PupaeM3','Adults','AdultsM1','AdultsM2','AdultsM3','NumColumns',4);
% lgnd.Position(1) = 0.48;
% lgnd.Position(2) = 0.15;

