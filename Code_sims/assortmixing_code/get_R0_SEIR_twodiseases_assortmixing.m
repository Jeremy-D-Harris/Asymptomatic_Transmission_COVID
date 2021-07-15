function R0_calc = get_R0_SEIR_twodiseases_assortmixing(params)

% calculate basic reproduction number, R0

% parameters to local variables
beta_a=params.beta_a; beta_s = params.beta_s;
gamma_a = params.gamma_a; gamma_s = params.gamma_s;
gamma_e = params.gamma_e;

p_aa = params.p_aa; p_as = params.p_as;

% transmissions
T = [0 0 p_aa*beta_a p_as*beta_s; 0 0 (1-p_aa)*beta_a  (1-p_as)*beta_s; 0 0 0 0; 0 0 0 0];

% transitions
Sigma = [-gamma_e 0 0 0; 0 -gamma_e 0 0; gamma_e 0 -gamma_a 0; 0 gamma_e 0 -gamma_s];

NGM = -T*(inv(Sigma)); % next generation matrix

eigen_values = eig(NGM);

eigen_values(imag(eigen_values) ~= 0) = NaN;

R0_calc = max(eigen_values);


