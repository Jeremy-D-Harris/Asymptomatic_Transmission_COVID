function [i_a,i_s] = get_incidence_SEIR_agedep(params,y_traj)


sigma = params.sigma;
prob_symp = params.prob_symp;
prob_asymp = 1-prob_symp;
N = params.N;

for count=1:length(params.t_span)
    
    t = params.t_span(count);
    
    [beta_a,beta_s]=mitigation_function(t,params);
    M_tilde = params.M_tilde;
    
    this_y = y_traj(count,:);
    y_matrix = reshape(this_y,[N, 7]);
    S = y_matrix(:,1); 
    I_a = y_matrix(:,4); I_s = y_matrix(:,5); 
    
    q_a_vector = prob_asymp.*sigma.*S;
    q_s_vector = prob_symp.*sigma.*S;

    q_susceptible = q_a_vector + q_s_vector;
    
    i_a(count,:) = beta_a*q_a_vector.*(M_tilde*I_a)+beta_s*q_a_vector.*(M_tilde*I_s);
    i_s(count,:) = beta_a*q_s_vector.*(M_tilde*I_a)+beta_s*q_s_vector.*(M_tilde*I_s);
    
end