function [beta_a, beta_s]= mitigation_function(t,params)

% mitigation function to reduce contact rates

% slope_down = (1-params.mitigation_level)/(params.t_m2 - params.t_m1);

mu = params.mitigation_level;
% lam = 1/2;
% fifth order polynomial with:
% g(t1) = 1 and g(t2) = mitigation level,
% additionally, g'(t1) = 0 and g'(t2) = 0

% % Interpolation step
% y_interp = @(s) interp1(nominales(:,1),nominales(:,3),s);
% i_interp = @(t,s) feval(y_interp,t-s);
%
% conv_i_g = integral(@(s) g(s).*i_interp(t,s),lb(t),ub(t),'RelTol',0,'AbsTol',Tolintegral);


t1=params.t_m1; t2=params.t_m2;


% t1=params.t_m1; t2=params.t_m2;
% b = 2*(t1+t2);
% c = t2^2+4*t1*t2+t1^2;
% d = 2*t1*t2*(t1+t2);
% e = t1^2*t2^2;
% k1 = t1^5/5-b*t1^4/4+c*t1^3/3-d*t1^2/2+e*t1;
% k2 = t2^5/5-b*t2^4/4+c*t2^3/3-d*t2^2/2+e*t2;
% a = (1-mu)/(k1-k2);
% f = (mu*k1-k2)/(k1-k2);
% mit_g = @(t) a*(t^5/5-b*t^4/4+c*t^3/3-d*t^2/2+e*t)+f;

% t1=params.t_m1; t2=params.t_m2;
% b = t1+t2;
% c = t1*t2;
% k1 = t1^3/3-b*t1^2/2+c*t1;
% k2 = t2^3/3-b*t2^2/2+c*t2;
% a = (1-mu)/(k1-k2);
% d = (mu*k1-k2)/(k1-k2);
% mit_g = @(t) a*t^3/3-a*b*t^2/2+a*c*t+d;

if t >= params.t_m1 && t < params.t_m2
    
    h_prime = @(s) (s-t1).^2.*(s-t2).^3;
    lb=@(t) t1;
    ub=@(t) t;
    H = @(t) integral(@(s) h_prime(s),lb(t),ub(t),'RelTol',10^-6,'AbsTol',10^-6);
    
    % feval(H,t1);
    % feval(H,t2);
    a = (mu-1)/H(t2);
    c = 1;
    
    g_mit = @(t) a*H(t)+c;
    
    beta_a = params.beta_a*g_mit(t);
    beta_s = params.beta_s*g_mit(t);
    
elseif t >= params.t_m2
    
    beta_a=params.beta_a*params.mitigation_level;
    beta_s=params.beta_s*params.mitigation_level;
    
else
    
    beta_a=params.beta_a;
    beta_s=params.beta_s;
    
end