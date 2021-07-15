function R0_calc = get_R0_SEIR_agedep(params)


% parameters to local variables
beta_a=params.beta_a; beta_s=params.beta_s;
gamma_a = params.gamma_a; gamma_s = params.gamma_s;
gamma_e = params.gamma_e;

M = params.M;
sigma = params.sigma;
prob_symp = params.prob_symp;
prob_asymp = 1-prob_symp;
N = params.N;


% for n=1:length(params.t_span)
%
%     t = params.t_span(n);

%     M = get_contact_matrix(t,params);
%     [beta_a,beta_s]=mitigation_agestratified_simple(t,params);

diag_gamma_a = diag(gamma_a*ones(N,1));
diag_gamma_s = diag(gamma_s*ones(N,1));
diag_gamma_e = diag(gamma_e*ones(N,1));

%     this_y = y_traj(n,:);
%     y_matrix = reshape(this_y,[N, 5]);
%     S = y_matrix(:,1);
% I_a = y_matrix(:,2); I_s = y_matrix(:,3);

q_a_vector = diag(prob_asymp.*sigma);
q_s_vector = diag(prob_symp.*sigma);

% q_susceptible = q_a_vector + q_s_vector;

% dSdt = [- beta_a*q_susceptible*M - beta_s*q_susceptible*M zeros(N,4*N)];
% Ea(1), Es(2), Ia(3), Is(4), Ra(5), Rs(6)
T1 = [zeros(N,N) zeros(N,N) beta_a*q_a_vector*M beta_s*q_a_vector*M];
T2 = [zeros(N,N) zeros(N,N) beta_a*q_s_vector*M beta_s*q_s_vector*M];
T3 = [zeros(N,N) zeros(N,N) zeros(N,N) zeros(N,N)];
T4 = [zeros(N,N) zeros(N,N) zeros(N,N) zeros(N,N)];

S1 = [-diag_gamma_e zeros(N,N) zeros(N,N) zeros(N,N)];
S2 = [zeros(N,N) -diag_gamma_e zeros(N,N) zeros(N,N)];
S3 = [diag_gamma_e zeros(N,N) -diag_gamma_a zeros(N,N)];
S4 = [zeros(N,N)  diag_gamma_e zeros(N,N) -diag_gamma_s];


T = [T1; T2; T3; T4];

Sigma = [S1; S2; S3; S4];

NGM = -T*(inv(Sigma));

eigen_values = eig(NGM);

R0_calc = max(eigen_values);

% end

