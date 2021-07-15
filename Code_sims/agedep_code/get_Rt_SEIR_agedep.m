function Rt_calc = get_Rt_SEIR_agedep(params,y_traj)


% parameters to local variables
gamma_a = params.gamma_a; gamma_s = params.gamma_s;
gamma_e = params.gamma_e;

M_tilde = params.M_tilde;
sigma = params.sigma;
prob_symp = params.prob_symp;
prob_asymp = 1-prob_symp;
N = params.N;

% p_a = params.p_aa; p_s = params.p_as;

for count=1:length(params.t_span)
    
    t = params.t_span(count);
    
%     M = get_contact_matrix(t,params);
    %     [beta_a,beta_s]=mitigation_agestratified_simple(t,params);
    [beta_a,beta_s]=mitigation_function(t,params);
    
    diag_gamma_a = diag(gamma_a*ones(N,1));
    diag_gamma_s = diag(gamma_s*ones(N,1));
    diag_gamma_e = diag(gamma_e*ones(N,1));
    
    this_y = y_traj(count,:);
    y_matrix = reshape(this_y,[N, 7]);
    S = y_matrix(:,1);
    % I_a = y_matrix(:,2); I_s = y_matrix(:,3);
    
    q_a_vector = diag(prob_asymp.*sigma.*S);
    q_s_vector = diag(prob_symp.*sigma.*S);
    
    % q_susceptible = q_a_vector + q_s_vector;
    
    % dSdt = [- beta_a*q_susceptible*M - beta_s*q_susceptible*M zeros(N,4*N)];
    
    % Ea(1), Es(2), Ia(3), Is(4), Ra(5), Rs(6)
    T1 = [zeros(N,N) zeros(N,N) beta_a*q_a_vector*M_tilde beta_s*q_a_vector*M_tilde];
    T2 = [zeros(N,N) zeros(N,N) beta_a*q_s_vector*M_tilde beta_s*q_s_vector*M_tilde];
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
    
    Rt_calc(count) = max(eigen_values);
    
end

