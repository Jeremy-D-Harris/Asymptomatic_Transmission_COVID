function eigen_vector = get_eigendirection_SEIR_twodiseases_fixedpropasymp(params)

% parameters to local variables
beta_a=params.beta_a; beta_s = params.beta_s;
gamma_a = params.gamma_a; gamma_s = params.gamma_s;
gamma_e = params.gamma_e;

p = params.p;

% S = y(1); E_a = y(2); E_s = y(3); I_a = y(4); I_s = y(5); R_a = y(6); R_s = y(7);
A = [0 0 0 -beta_a -beta_s 0 0; 0 -gamma_e 0 p*beta_a p*beta_s 0 0; 0 0 -gamma_e (1-p)*beta_a (1-p)*beta_s 0 0;...
    0 gamma_e 0 -gamma_a 0 0 0; 0 0 gamma_e 0 -gamma_s 0 0; 0 0 0 0 gamma_a 0 0; 0 0 0 0 0 gamma_s 0];

[eigen_directions, eigen_values] = eig(A); % get eigenvalues/eigenvectors

[val, ind] = max(diag(eigen_values));
eigen_vector = eigen_directions(:,ind); % corresponds to max eigenvalue
