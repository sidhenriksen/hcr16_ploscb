function A = fractional_area(psi)
        % Computes the fractional area as given in Doi et al. (2013)

        [~,i] = min(abs(psi-0.5));
        
        A_num = 2*sum(0.5-psi(1:i));
        A_denom = sum(abs(0.5-psi));
        
        A = A_num/A_denom;
end