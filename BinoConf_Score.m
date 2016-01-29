function [Lo Up]=BinoConf_Score(m,n,varargin)
%[Lo Up]=BinoConf_Score(m,n,alpha)
%m=hits, n=total of trials, alpha=0.05 for CI 95%, alpha = 0.32 for 68% CI corresponding to +/-1SD of normal
%Binomial confidence interval
%Score confidence interval Edwin B. Wilson (1927)
%According to Agresti and Coull (1998) this is the best confidence interval,
%it yields coverage probabilities close to nominal confidence levels, even
%for very small sample sizes.
%
%Agresti, A. and Coull, B. A. (1998). Approximate is better than "exact"
%for interval estimation of binomial proportions", The American Statistician, 52(2), 119-12
%Wilson, E. B. (1927). Probable inference, the law of succession,
%and statistical inference. J. Amer. Statist. Assoc. 22 209–212.
%

%28-07-2008
%Dr. Ignacio Serrano-Pedraza

if nargin==2
    alpha = 0.05;
else
    alpha = varargin{1};
end

p=m./n;
z=norminv(alpha/2);
A=p+((z^2)./(2*n));
B=z*sqrt((p.*(1-p)+z^2./(4*n))./n);
Lo=(A+B)./(1+z^2./n);
Up=(A-B)./(1+z^2./n);