clear 
t(1)    =0  ;  
L(1)    =1  ; % larvae
P(1)    =1  ; %pupae 
A(1)    =50 ;% adults 

b       =11.6772   ; %number of eggs
c_ea    =0.0110   ; %cannibalism rate eggs by adults
c_el    =0.0093   ; %cannibalisim rate eggs by larvae
c_pa    =0.0162   ; %cannibalism rate pupae by adults
mu_l    =0.5129   ; %survival rate for larvae
mu_a    =0.1108   ; %survial rate for adults
mu_p    =0.7;           

i=2;
while (L(i-1) > 0 && t(i-1)<100000)         
     r1 =   b*A(i-1).*exp(-c_ea * A(i-1)-c_el * L(i-1)); % larvae birth rate    
     r2 =   mu_l * L(i-1); % larvae death rate
     r3 =   (1-mu_l)*L(i-1); %pupae birth rate
     r4 =   mu_p*P(i-1); % pupae death rate
     r5 =   P(i-1).*exp(-c_pa*A(i-1)); %Adults birth rate
     r6 =   mu_a * A(i-1);  %Adults death rate
     
     total_rate=r1+r2+r3+r4+r5+r6;    
     
     t(i)=t(i-1)-log(rand(1)/total_rate);  % generates  random variable    
     u=total_rate*rand(1);
        
     if 0<u & u<r1       % larvae birth rate         
         L(i)   =   L(i-1)+1;          
         P(i)   =   P(i-1);
         A(i)   =   A(i-1);
         
       
     elseif r1<u & u<r1+r2 % larvae death rate
         L(i)   =   L(i-1)-1;          
         P(i)   =   P(i-1);
         A(i)   =   A(i-1);
         
     elseif r1+r2<u & u<r1+r2+r3 %pupae birth rate
         L(i)   =   L(i-1);          
         P(i)   =   P(i-1)+1;
         A(i)   =   A(i-1);
      
     elseif  r1+r2+r3<u & u<r1+r2+r3+r4  %pupae death rate
         L(i)   =   L(i-1);          
         P(i)   =   P(i-1)-1;
         A(i)   =   A(i-1);
     
     elseif r1+r2+r3+r4<u & u<r1+r2+r3+r4+r5 %adults birth rate
         L(i)   =   L(i-1);          
         P(i)   =   P(i-1);
         A(i)   =   A(i-1)+1;
     else                               %adults death rate
         L(i)   =   L(i-1);          
         P(i)   =   P(i-1);
         A(i)   =   A(i-1)-1;
     end  
     i =i+1;
    
end
plot(t,L,'-k','DisplayName','Larvae')
hold on
plot(t,P,'-r','DisplayName','Pupae')
plot(t,A,'-b','DisplayName','Adults')
xlabel('Time')
ylabel('Population')
legend