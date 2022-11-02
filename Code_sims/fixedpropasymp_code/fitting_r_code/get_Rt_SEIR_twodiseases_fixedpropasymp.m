function Rt_calc = get_Rt_SEIR_twodiseases_fixedpropasymp(params,y)

% calculate effective reproduction number, R_t

% parameters to local variables
gamma_a = params.gamma_a; gamma_s = params.gamma_s;
gamma_e = params.gamma_e;

p = params.p;

for count=1:length(params.t_span)

t = params.t_span(count);
S = y(count,1);

% get contact rates
[beta_a, beta_s]= mitigation_function(t,params);

% transmissions
T = [0 0 p*beta_a*S p*beta_s*S; 0 0 (1-p)*beta_a*S  (1-p)*beta_s*S; 0 0 0 0; 0 0 0 0];

% transitions
Sigma = [-gamma_e 0 0 0; 0 -gamma_e 0 0; gamma_e 0 -gamma_a 0; 0 gamma_e 0 -gamma_s];

NGM = -T*(inv(Sigma)); % next generation matrix

eigen_values = eig(NGM);

Rt_calc(count) = max(eigen_values);

end
