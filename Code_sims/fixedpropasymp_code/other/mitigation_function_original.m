function [beta_a, beta_s]= mitigation_function_original(t,params)

% mitigation function to reduce contact rates

slope_down = (1-params.mitigation_level)/(params.t_m2 - params.t_m1);

if t >= params.t_m1 && t < params.t_m2
    
    beta_a = params.beta_a*(1-slope_down*(t-params.t_m1));
    beta_s = params.beta_s*(1-slope_down*(t-params.t_m1));
    
elseif t >= params.t_m2
    
    beta_a=params.beta_a*params.mitigation_level;
    beta_s=params.beta_s*params.mitigation_level;
    
else
    
    beta_a=params.beta_a;
    beta_s=params.beta_s;
    
end