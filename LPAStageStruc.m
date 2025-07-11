%LPA Model with Projected Matrix.
%Parameters were changed from original LPA.

b = 1.4286;
mu_a = 0.006;
mu_l = 0.0432;
c_ea = 0.0013;
c_el = 0.0002;
c_pa = 0.0008;

MaxT = 29;

ProjMatrix = zeros(29);
ProjMatrix(1,29) = b;

P = zeros(29,100);
P(29,1) = 50;

MaxN = 100;
 
for n = 1:MaxN
    ProjMatrix(2,1) = exp(-c_ea*P(29,n)-c_el*sum(P(4:11,n)));
    ProjMatrix(29,29) = (1-mu_a);
    for i = 3:29
        if i < 5; 
            ProjMatrix(i,i-1) =exp(-c_ea*P(29,n)-c_el*sum(P(4:11,n)));
        elseif i < 13;
            ProjMatrix(i,i-1) =(1-mu_l);
        elseif i < 29;
            ProjMatrix(i,i-1) = exp(-c_pa*P(29,n));
        elseif i == 29;
            ProjMatrix(i,i-1) = exp(-c_pa*P(29,n));
        end
        P(:,n+1) = ProjMatrix*P(:,n);
    end
    
end
figure(2)
[T,S] = meshgrid(1:MaxN+1,1:29);
Z = T.*S;
surf(T,S,P,Z)
colormap
xlabel('Time(weeks')
ylabel('Stage')
zlabel('Population')
hold off
LPAMat = zeros(MaxN,4);
LPAMat(:,1) = 1:MaxN;
LPAMat(1,4) = 50;
for j = 2:MaxN
%     [T,S] = meshgrid(1:MaxN+1,1:29);
%     surf(T,S,P)
    LPAMat(j,2) = sum(P(1:14,j));
    LPAMat(j,3) = sum(P(15:28,j));
    LPAMat(j,4) = sum(P(29,j));
end
figure(1)
plot(LPAMat(:,1),LPAMat(:,2:4), linewidth = 2)
title('Stage Structured Model')
xlabel('Time(Days)')
ylabel('Population')
lgnd = legend('Larvae','Pupae','Adutls')
lgnd.Location = 'northeast';