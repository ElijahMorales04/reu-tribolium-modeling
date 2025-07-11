%4-Cup Dispersal w/ LPA Model
%adding direct cost to dispersal

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

index1 = 1;
index2 = 1;
index3 = 1;
index4 = 1;

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
    disp([num2str(i/MaxT),'%'])
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
    randPatch = randi(4);
    if randPatch == 1;
        DF1i(20*(i+1)+1,:) = 0;
    elseif randPatch == 2;
        DF2i(20*(i+1)+1,:) = 0;
    elseif randPatch == 3;
        DF3i(20*(i+1)+1,:) = 0;
    elseif randPatch == 4;
        DF4i(20*(i+1)+1,:) = 0;
    end
    
    while done == 'a';
        randPatch = randi(4);
        if randPatch == 1 && DF1i(20*(i+1)+1,1)~=0 && index1~=4;
            if index1 == 1;
                DF1i(20*(i+1)+1,10) = DF1i(20*(i+1),10)+0.01;
                index1 = index1+1;
                done = 'b';
            elseif index1 == 2;
                DF1i(20*(i+1)+1,11) = DF1i(20*(i+1),11)+0.01;
                index1 = index1+1;
                done = 'b';
            elseif index1 == 3;
                DF1i(20*(i+1)+1,12) = DF1i(20*(i+1),12)+0.01;
                index1 = index1+1;
                done = 'b';
            end
            
        elseif randPatch == 2 && DF2i(20*(i+1)+1,1)~=0 && index2~=4;
            if index2 == 1;
                DF2i(20*(i+1)+1,10) = DF2i(20*(i+1),10)+0.01;
                index2 = index2+1;
                done = 'b';
            elseif index2 == 2;
                DF2i(20*(i+1)+1,11) = DF2i(20*(i+1),11)+0.01;
                index2 = index2+1;
                done = 'b';
            elseif index2 == 3;
                DF2i(20*(i+1)+1,12) = DF2i(20*(i+1),12)+0.01;
                index2 = index2+1;
                done = 'b';
            end
            
        elseif randPatch == 3 && DF3i(20*(i+1)+1,1)~=0 && index3~=4;
            if index3 == 1;
                DF3i(20*(i+1)+1,10) = DF3i(20*(i+1),10)+0.01;
                index3 = index3+1;
                done = 'b';
            elseif index3 == 2;
                DF3i(20*(i+1)+1,11) = DF3i(20*(i+1),11)+0.01;
                index3 = index3+1;
                done = 'b';
            elseif index3 == 3;
                DF3i(20*(i+1)+1,12) = DF3i(20*(i+1),12)+0.01;
                index3 = index3+1;
                done = 'b';
            end
            
        elseif randPatch == 4 && DF4i(20*(i+1)+1,1)~=0 && index4~=4;
            if index4 == 1;
                DF4i(20*(i+1)+1,10) = DF4i(20*(i+1),10)+0.01;
                index4 = index4+1;
                done = 'b';
            elseif index4 == 2;
                DF4i(20*(i+1)+1,11) = DF4i(20*(i+1),11)+0.01;
                index4 = index4+1;
                done = 'b';
            elseif index4 == 3;
                DF4i(20*(i+1)+1,12) = DF4i(20*(i+1),12)+0.01;
                index4 = index4+1;
                done = 'b';
            end
        else
            done = 'b';
        end
    end
    
end


figure(1);
%PATCH 1
subplot(2,2,1);
x = linspace(1,421,421);
Lr1 = DF1i(:,1);
plot(x,Lr1, linewidth=2)

hold on
Lm1a = DF1i(:,2);
plot(x,Lm1a,linewidth=2)
Lm1b = DF1i(:,3);
plot(x,Lm1b,linewidth=2)
Lm1c = DF1i(:,4);
plot(x,Lm1c,linewidth=2)

Pr1 = DF1i(:,5);
plot(x,Pr1, linewidth=2)

Pm1a = DF1i(:,6);
plot(x,Pm1a,linewidth=2)
Pm1b = DF1i(:,7);
plot(x,Pm1b,linewidth=2)
Pm1c = DF1i(:,8);
plot(x,Pm1c,linewidth=2)

Ar1 = DF1i(:,9);
plot(x,Ar1, linewidth=2)

Am1a = DF1i(:,10);
plot(x,Am1a,linewidth=2)
Am1b = DF1i(:,11);
plot(x,Am1b,linewidth=2)
Am1c = DF1i(:,12);
plot(x,Am1c,linewidth=2)
hold off


%PATCH 2
subplot(2,2,2);
x = linspace(1,421,421);
Lr2 = DF1i(:,1);
plot(x,Lr2, linewidth=2)

hold on
Lm2a = DF1i(:,2);
plot(x,Lm2a,linewidth=2)
Lm2b = DF1i(:,3);
plot(x,Lm2b,linewidth=2)
Lm2c = DF1i(:,4);
plot(x,Lm2c,linewidth=2)

Pr2 = DF1i(:,5);
plot(x,Pr2, linewidth=2)

Pm2a = DF1i(:,6);
plot(x,Pm2a,linewidth=2)
Pm2b = DF1i(:,7);
plot(x,Pm2b,linewidth=2)
Pm2c = DF1i(:,8);
plot(x,Pm2c,linewidth=2)

Ar2 = DF1i(:,9);
plot(x,Ar2, linewidth=2)

Am2a = DF1i(:,10);
plot(x,Am2a,linewidth=2)
Am2b = DF1i(:,11);
plot(x,Am2b,linewidth=2)
Am3c = DF1i(:,12);
plot(x,Am3c,linewidth=2)
hold off

%PATCH 3
subplot(2,2,3);
x = linspace(1,421,421);
Lr3 = DF1i(:,1);
plot(x,Lr3, linewidth=2)

hold on
Lm3a = DF1i(:,2);
plot(x,Lm3a,linewidth=2)
Lm3b = DF1i(:,3);
plot(x,Lm3b,linewidth=2)
Lm3c = DF1i(:,4);
plot(x,Lm3c,linewidth=2)

Pr3 = DF1i(:,5);
plot(x,Pr3, linewidth=2)

Pm3a = DF1i(:,6);
plot(x,Pm3a,linewidth=2)
Pm3b = DF1i(:,7);
plot(x,Pm3b,linewidth=2)
Pm3c = DF1i(:,8);
plot(x,Pm3c,linewidth=2)

Ar3 = DF1i(:,9);
plot(x,Ar3, linewidth=2)

Am3a = DF1i(:,10);
plot(x,Am3a,linewidth=2)
Am3b = DF1i(:,11);
plot(x,Am3b,linewidth=2)
Am3c = DF1i(:,12);
plot(x,Am3c,linewidth=2)
hold off

%PATCH 4
subplot(2,2,4);
x = linspace(1,421,421);
Lr4 = DF1i(:,1);
plot(x,Lr4, linewidth=2)

hold on
Lm4a = DF1i(:,2);
plot(x,Lm4a,linewidth=2)
Lm4b = DF1i(:,3);
plot(x,Lm4b,linewidth=2)
Lm4c = DF1i(:,4);
plot(x,Lm4c,linewidth=2)

Pr4 = DF1i(:,5);
plot(x,Pr4, linewidth=2)

Pm4a = DF1i(:,6);
plot(x,Pm4a,linewidth=2)
Pm4b = DF1i(:,7);
plot(x,Pm4b,linewidth=2)
Pm4c = DF1i(:,8);
plot(x,Pm4c,linewidth=2)

Ar4 = DF1i(:,9);
plot(x,Ar4, linewidth=2)

Am4a = DF1i(:,10);
plot(x,Am4a,linewidth=2)
Am4b = DF1i(:,11);
plot(x,Am4b,linewidth=2)
Am4c = DF1i(:,12);
plot(x,Am4c,linewidth=2)
hold off

lgnd = legend('Larvae','Pupae','Adults','Larvae','Pupae','Adults','Larvae','Pupae','Adults','Larvae','Pupae','Adults','NumColumns',3);
lgnd.Position(1) = 0.55;
lgnd.Position(2) = 0.85;
