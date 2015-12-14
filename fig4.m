% This will recreate Figure 4 from Henriksen, Cumming, & Read (2016).

clear all; close all; clc;
silent = 0;
bem = BEMunit('silent',silent);
bem.Nx = 292; bem.Ny = 292;
bem.deg_per_pixel = 0.03;
bem.dx = 0.09;
bem = bem.update();

antibem = bem;
antibem.dx = bem.dx*-1;
antibem = antibem.update();

correlated_rds = pairedRDS(); correlated_rds.dx = bem.dx/bem.deg_per_pixel;
correlated_rds.Nx = 292; correlated_rds.Ny = 292;
correlated_rds.dotsize=3; % pixels

% Create halfmatched stereogram
halfmatched_rds = correlated_rds;
halfmatched_rds.correlation = 0;

sxs = linspace(0.025,0.30,21);
densities = logspace(log10(0.04),log10(4),11);

correlated_bem = zeros(length(sxs),length(densities));
halfmatched_bem = zeros(length(sxs),length(densities));

correlated_antibem = zeros(length(sxs),length(densities));
halfmatched_antibem = zeros(length(sxs),length(densities));

N = 2000; % number of images to average across

bootstrap_mode = 1;

load_bootstrap = bootstrap_mode && N ==0;

for j = 1:length(sxs);
    fprintf('Running size %i of %i\n',j,length(sxs))    

    sx=bem.subunits(1).rf_params.right.sx;
    % Rescale the RF of the model units
    bem = bem.rescale(sxs(j)/sx);
    antibem = antibem.rescale(sxs(j)/sx);
    

    for k = 1:length(densities);
        if j < 5 && k < 4
            N = 4000;
        else
            N = 2000;
        end
        correlated_rds.density = densities(k);
        halfmatched_rds.density = densities(k);                

        % Compute correlated and half-matched response for bem unit
        if load_bootstrap
            bem = bem.load_bootstrap(correlated_rds);
        end
        Cc_bem = bem.simulate_spatial(correlated_rds,N,bootstrap_mode);
        
        if load_bootstrap
            bem = bem.load_bootstrap(halfmatched_rds);
        end
        Chm_bem = bem.simulate_spatial(halfmatched_rds,N,bootstrap_mode);
        
                
        % Compute the same for its antineuron
        if load_bootstrap
            antibem = antibem.load_bootstrap(halfmatched_rds);
        end
        
        Chm_antibem = antibem.simulate_spatial(halfmatched_rds,N,bootstrap_mode);

        correlated_bem(j,k) = mean(Cc_bem.^2);
        halfmatched_bem(j,k) = mean(Chm_bem.^2);      
        halfmatched_antibem(j,k) = mean(Chm_antibem.^2);
    end
    
    
end

R_bem = halfmatched_bem./correlated_bem;
R_antibem = halfmatched_antibem./correlated_bem;

figure();
imagesc(1:11,sxs./0.075,R_bem-R_antibem);
set(gca,'ydir','norm','xtick',1:2:11,'xticklabels',round(densities(1:2:11),2))
ylabel('\sigma / r')
xlabel('Density')