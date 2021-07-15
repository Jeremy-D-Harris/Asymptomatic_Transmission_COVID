function [asymp_transmission,symp_transmission] = get_transmission_SEIR_agedep(params,y_traj)

% params.t_span;

% parameters to local variables
% [beta_a, beta_s]= mitigation_function(t,params);
% delta_a = params.delta_a; delta_s = params.delta_s;
% delta_e = params.delta_e;

M_tilde = params.M_tilde;
sigma = params.sigma;
prob_symp = params.prob_symp;
prob_asymp = 1-prob_symp;
N = params.N;

asymp_transmission = [];
symp_transmission = [];
for count=1:length(params.t_span)
    
    t = params.t_span(count);
    
%     [beta_a,beta_s]=mitigation_agestratified_simple(t,params);
    [beta_a,beta_s]=mitigation_function(t,params);
    
    this_y = y_traj(count,:);
    y_matrix = reshape(this_y,[N, 7]);
    S = y_matrix(:,1);
    I_a = y_matrix(:,4); I_s = y_matrix(:,5);
    
    q_a_vector = prob_asymp.*sigma.*S;
    q_s_vector = prob_symp.*sigma.*S;
    
    q_susceptible = q_a_vector + q_s_vector;
    
    asymp_transmission(count,:) = beta_a*q_susceptible.*(M_tilde*I_a);
    symp_transmission(count,:) = beta_s*q_susceptible.*(M_tilde*I_s);
    
    
end

