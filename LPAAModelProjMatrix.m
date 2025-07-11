%LPAA Model with Projected Matrix.
%Parameters were changed from original LPAA.

b = 20.0;
mu_a = 0.0842;
mu_l = 0.6053;
c_ea = .0099;
c_el = .003;
c_aa = 0.0028;

MaxT = 20;

DF = zeros(MaxT+1,4);
DF(1,:) = [ 0, 0, 0, 50];


ProjMatrix = [0, 0, 0, b*exp(-c_ea*DF(1,4));
              (1-mu_l), 0, 0, 0;
              0, (1-mu_p), 0, 0;
              0, 0, exp(-c_aa*DF(1,4)),(1-mu_a)];
          
for i = 2:MaxT+1
    DF(i,:) = ProjMatrix*DF(i-1,:)';
    
    ProjMatrix = [0, 0, 0, b*exp(-c_ea*DF(i,4));
                 (1-mu_l), 0, 0, 0;
                  0, (1-mu_p), 0, 0;
                  0, 0, exp(-c_aa*DF(i,4)),(1-mu_a)];
end

plot(1:MaxT+1,DF(:,1:4), linewidth = 2)
xlabel('Time (2-week period)', ...
       'Interpreter', 'latex', ...
       'FontSize', 22)

ylabel('Number of beetles', ...
       'Interpreter', 'latex', ...
       'FontSize', 22)

lgnd = legend('Larvae', 'Pupae', 'Adults(Callow)', 'Adults(Mature)');
lgnd.Interpreter = 'latex';
lgnd.Location = 'northeast';