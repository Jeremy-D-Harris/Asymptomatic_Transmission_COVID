function r_calc = get_r_SEIR_agedep(params)

% parameters to local variables
beta_a=params.beta_a; beta_s = params.beta_s;
gamma_a = params.gamma_a; gamma_s = params.gamma_s;
gamma_e = params.gamma_e;

M = params.M;
sigma = params.sigma;
prob_symp = params.prob_symp;
prob_asymp = 1-prob_symp;
N = params.N;

diag_gamma_a = diag(gamma_a*ones(N,1));
diag_gamma_s = diag(gamma_s*ones(N,1));
diag_gamma_e = diag(gamma_e*ones(N,1));

% y_matrix = reshape(y,[N, 5]);
% S = y_matrix(:,1); 
% I_a = y_matrix(:,2); I_s = y_matrix(:,3); 
% R_a = y_matrix(:,4); R_s = y_matrix(:,5);


q_a_vector = diag(prob_asymp.*sigma);
q_s_vector = diag(prob_symp.*sigma);

q_susceptible = q_a_vector + q_s_vector;

dSdt = [zeros(N,N) zeros(N,N) zeros(N,N) -beta_a*q_susceptible*M -beta_s*q_susceptible*M zeros(N,2*N)];

dEadt = [zeros(N,N)  -diag_gamma_e zeros(N,N) beta_a*q_a_vector*M  beta_s*q_a_vector*M zeros(N,2*N)];

dEsdt = [zeros(N,N)  zeros(N,N) -diag_gamma_e  beta_a*q_s_vector*M beta_s*q_s_vector*M zeros(N,2*N)];

dIadt = [zeros(N,N)  diag_gamma_e zeros(N,N) -diag_gamma_a zeros(N,N) zeros(N,2*N)];

dIsdt = [zeros(N,N)  zeros(N,N) diag_gamma_e  zeros(N,N) -diag_gamma_s zeros(N,2*N)];

dRadt = [zeros(N,N) zeros(N,N) zeros(N,N) diag_gamma_a zeros(N,N) zeros(N,2*N)];

dRsdt = [zeros(N,N) zeros(N,N) zeros(N,N) zeros(N,N) diag_gamma_s zeros(N,2*N)];

A = [dSdt; dEadt; dEsdt; dIadt; dIsdt; dRadt; dRsdt];


% S = y(1); I_a = y(2); I_s = y(3); R_a = y(4); R_s = y(5);
% A = [0 -beta_a -beta_s 0 0; 0 p*beta_a-gamma_a p*beta_s 0 0; 0 (1-p)*beta_a (1-p)*beta_s-gamma_s 0 0; 0 gamma_a 0 0 0; 0 0 gamma_s 0 0];

[eigen_directions, eigen_values] = eig(A); % get eigenvalues/eigenvectors
[val ind] = max(real(diag(eigen_values)));
r_calc = val;

