function [beta_a, beta_s]= mitigation_function(t,params)

% mitigation function to reduce contact rates

% slope_down = (1-params.mitigation_level)/(params.t_m2 - params.t_m1);

mu = params.mitigation_level;

% g(t1) = 1 and g(t2) = mitigation level,
% additionally, g'(t1) = 0 and g'(t2) = 0

t1=params.t_m1; t2=params.t_m2;
b = 2*(t1+t2); 
c = t2^2+4*t1*t2+t1^2;
d = 2*t1*t2*(t1+t2);
e = t1^2*t2^2;
k1 = t1^5/5-b*t1^4/4+c*t1^3/3-d*t1^2/2+e*t1;
k2 = t2^5/5-b*t2^4/4+c*t2^3/3-d*t2^2/2+e*t2;
a = (1-mu)/(k1-k2);
f = (mu*k1-k2)/(k1-k2);
mit_g = @(t) a*(t^5/5-b*t^4/4+c*t^3/3-d*t^2/2+e*t)+f;

if t >= params.t_m1 && t < params.t_m2
    
    beta_a = params.beta_a*mit_g(t);
    beta_s = params.beta_s*mit_g(t);
    
elseif t >= params.t_m2
    
    beta_a=params.beta_a*params.mitigation_level;
    beta_s=params.beta_s*params.mitigation_level;
    
else
    
    beta_a=params.beta_a;
    beta_s=params.beta_s;
    
end