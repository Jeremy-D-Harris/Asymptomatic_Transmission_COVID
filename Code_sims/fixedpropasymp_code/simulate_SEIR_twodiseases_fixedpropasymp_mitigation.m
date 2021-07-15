function dydt = simulate_SEIR_twodiseases_fixedpropasymp_mitigation(t,y,params)

% simulate the SEIR model with p the proportion asymptomatic incidence

% get contact rates over time
[beta_a, beta_s]= mitigation_function(t,params);

% parameters to local variables
gamma_a = params.gamma_a; gamma_s = params.gamma_s;
gamma_e = params.gamma_e;

p = params.p;

S = y(1); E_a = y(2); E_s = y(3); I_a = y(4); I_s = y(5); R_a = y(6); R_s = y(7);

dSdt = -beta_a*S*I_a - beta_s*S*I_s;

dEadt = p*(beta_a*S*I_a + beta_s*S*I_s) - gamma_e*E_a;

dEsdt = (1-p)*(beta_a*S*I_a + beta_s*S*I_s) - gamma_e*E_s;

dIadt = gamma_e*E_a - gamma_a*I_a;

dIsdt = gamma_e*E_s - gamma_s*I_s;

dRadt = gamma_a*I_a;

dRsdt = gamma_s*I_s;

dydt = [dSdt; dEadt; dEsdt; dIadt; dIsdt; dRadt; dRsdt];
