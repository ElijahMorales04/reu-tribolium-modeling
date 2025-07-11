function dpdt = ELPA_ODE(t,p)
lambda      =   13;
mu_ea       =   0.1;
mu_el       =   0.05;
alpha_e     =   0.9;
mu_la       =   0.05;
alpha_l     =   0.85;
mu_pa       =   0.04;
alpha_p     =   0.5;
delta       =   0.4;

dpdt = [lambda * p(4)- mu_ea * p(4) .* p(1) - mu_el * p(2) .* p(1) - alpha_e * p(1);
        alpha_e * p(1) - 0.02 * p(2) .* (mu_el * p(1) - mu_la * p(4) ) - alpha_l * p(2);
        alpha_l * p(2) - mu_pa * p(4) .* p(3) - alpha_p * p(3);
        alpha_p * p(3) + 0.02 *  p(4) .* (mu_ea * p(1) + mu_el * p(2) + mu_pa * p(3)) - delta * p(4)        ];
end


