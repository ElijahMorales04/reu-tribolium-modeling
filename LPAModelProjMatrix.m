%LPA Model with Projected Matrix.
%Parameters were changed from original LPA.

b = 20.0;
mu_a = 0.0842;
mu_l = 0.6053;
c_ea = .0179;
c_el = .003;
c_pa = 0;

MaxT = 52;

DF = zeros(MaxT+1,3);
DF(1,:) = [ 0, 0, 50];


ProjMatrix = [0, 0, b*exp(-c_ea*DF(1,3)-c_el*DF(1,1));
              (1-mu_l), 0, 0;
              0,exp(-c_pa*DF(1,3)),(1-mu_a)];
          
for i = 2:MaxT+1
    DF(i,:) = ProjMatrix*DF(i-1,:)';
    
    ProjMatrix = [0, 0, b*exp(-c_ea*DF(i,3)-c_el*DF(i,1));
                 (1-mu_l), 0, 0;
                  0, exp(-c_pa*DF(i,3)), (1-mu_a)];
end

plot(1:MaxT+1,DF(:,1:3), linewidth = 2)
xlabel('Time (2-week period)', ...
       'Interpreter', 'latex', ...
       'FontSize', 22)

ylabel('Number of beetles', ...
       'Interpreter', 'latex', ...
       'FontSize', 22)

lgnd = legend('Larvae', 'Pupae', 'Adults');
lgnd.Interpreter = 'latex';
lgnd.Location = 'northeast';