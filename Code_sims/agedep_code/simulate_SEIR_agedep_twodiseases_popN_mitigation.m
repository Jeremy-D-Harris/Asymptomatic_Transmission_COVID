function dydt = simulate_SEIR_agedep_twodiseases_popN_mitigation(t,y,params)

[beta_a,beta_s]=mitigation_function(t,params);
M_tilde = params.M_tilde;
% M = get_contact_matrix(t,params);
% parameters to local variables
% [beta_a, beta_s]= mitigation_function(t,params);
gamma_a = params.gamma_a; gamma_s = params.gamma_s;
gamma_e = params.gamma_e;
% M = params.M_tilde;



sigma = params.sigma;
prob_symp = params.prob_symp;
prob_asymp = 1-prob_symp;
N = params.N;


% p_aa = params.p_aa; p_as = params.p_as;

% I_a_init = eigen_prop_asymp*1e-6*ones(N,1);
% I_s_init = (1-eigen_prop_asymp)*1e-6*ones(N,1);
% R_a_init = 0*ones(N,1);
% R_s_init = 0*ones(N,1);
% pop_tot = I_a_init+I_s_init+R_a_init+R_s_init;
% S_init = ones(N,1)-pop_tot;
y_matrix = reshape(y,[N, 7]);
S = y_matrix(:,1);
E_a = y_matrix(:,2); E_s = y_matrix(:,3);
I_a = y_matrix(:,4); I_s = y_matrix(:,5);
R_a = y_matrix(:,6); R_s = y_matrix(:,7);


q_a_vector = prob_asymp.*sigma.*S;
q_s_vector = prob_symp.*sigma.*S;

q_susceptible = q_a_vector + q_s_vector;

dSdt = -beta_a*q_susceptible.*(M_tilde*I_a) -beta_s*q_susceptible.*(M_tilde*I_s);

dEadt = beta_a*q_a_vector.*(M_tilde*I_a) + beta_s*q_a_vector.*(M_tilde*I_s) - gamma_e*E_a;

dEsdt = beta_a*q_s_vector.*(M_tilde*I_a) + beta_s*q_s_vector.*(M_tilde*I_s) - gamma_e*E_s;

dIadt = gamma_e*E_a - gamma_a*I_a;

dIsdt = gamma_e*E_s - gamma_s*I_s;

dRadt = gamma_a*I_a;

dRsdt = gamma_s*I_s;

dydt = [dSdt; dEadt; dEsdt; dIadt; dIsdt; dRadt; dRsdt];

