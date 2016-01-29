function params = fit_psych(psi_obs);
    
    init = [0.5,1,1];
    
    params = fminsearch(@(x)(psi_cost(x,psi_obs)),init);        
end

function cost = psi_cost(params,psi_obs);
  
    x = linspace(0,1,length(psi_obs))';
    
    psi_fit = psi_func(x,params);
        
    cost = sum((psi_fit-psi_obs).^2);
end

