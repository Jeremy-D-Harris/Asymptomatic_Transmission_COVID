function R0_calc = get_R0_SEIR_twodiseases_fixedpropasymp(params,y)

% calculate basic reproduction number, R0

% parameters to local variables
beta_a=params.beta_a; beta_s = params.beta_s;
gamma_a = params.gamma_a; gamma_s = params.gamma_s;
gamma_e = params.gamma_e;

p = params.p;

% transmissions
T = [0 0 p*beta_a p*beta_s; 0 0 (1-p)*beta_a  (1-p)*beta_s; 0 0 0 0; 0 0 0 0]; % transmission

% transitions
Sigma = [-gamma_e 0 0 0; 0 -gamma_e 0 0; gamma_e 0 -gamma_a 0; 0 gamma_e 0 -gamma_s];

NGM = -T*(inv(Sigma)); % next generation matrix

eigen_values = eig(NGM);

R0_calc = max(eigen_values);

