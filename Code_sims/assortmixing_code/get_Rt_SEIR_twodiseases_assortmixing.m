function Rt_calc = get_Rt_SEIR_twodiseases_assortmixing(params,y)

% calculate effective reproduction number, R_t

% parameters to local variables
gamma_a = params.gamma_a; gamma_s = params.gamma_s;
gamma_e = params.gamma_e;

p_a = params.p_aa; p_s = params.p_as;

for n=1:length(params.t_span)

t = params.t_span(n);
S = y(n,1); 

% get contact rates
[beta_a, beta_s]= mitigation_function(t,params);

% transmissions
T = [0 0 p_a*beta_a*S p_s*beta_s*S; 0 0 (1-p_a)*beta_a*S  (1-p_s)*beta_s*S; 0 0 0 0; 0 0 0 0];

% transitions
Sigma = [-gamma_e 0 0 0; 0 -gamma_e 0 0; gamma_e 0 -gamma_a 0; 0 gamma_e 0 -gamma_s];

NGM = -T*(inv(Sigma)); % next generation matrix

eigen_values = eig(NGM);

Rt_calc(n) = max(eigen_values);

end

