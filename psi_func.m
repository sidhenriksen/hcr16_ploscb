function y = psi_func(x,params)
    g = params(1);
    a = params(2);
    b = params(3);
    y = (1-(1-g)*exp(-(x/a).^b));
end