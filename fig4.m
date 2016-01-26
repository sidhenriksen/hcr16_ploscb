% This will recreate Figure 4 from Henriksen, Cumming, & Read (2016).
function fig4()
    %%  Set run-time parameters
    run_parallel = 1;
    silent = 1;
    
    %% Create BEM unit and its antineuron
    bem = BEMunit('silent',silent);
    bem.memory_threshold=bem.memory_threshold*0.5;
    bem.Nx = 292; bem.Ny = 292;
    bem.deg_per_pixel = 0.03;
    bem.dx = 0.10;
    bem.outputNL = @(x)(x.^2);
    for j = 1:length(bem.subunits);
        bem.subunits(j).rf_params.left.f=3;
        bem.subunits(j).rf_params.right.f=3;
    end
    bem = bem.update();
    

    %% Create stimuli
    % First correlated random dot stereogram
    correlated_rds = pairedRDS(); correlated_rds.dx = bem.dx/bem.deg_per_pixel;
    correlated_rds.Nx = 292; correlated_rds.Ny = 292;
    correlated_rds.dotsize=3; % pixels

    % Create halfmatched stereogram
    halfmatched_rds = correlated_rds;
    halfmatched_rds.correlation = 0;
    
    uncorrelated_rds = halfmatched_rds;
    uncorrelated_rds.toggle_match = 0;

    all_rds = {correlated_rds,halfmatched_rds,uncorrelated_rds};
    
    %% Stimulus parameters
    sxs = linspace(0.025,0.30,21);
    densities = logspace(log10(0.04),log10(4),11);

    resps = zeros(length(sxs),length(densities),3);

    %% 
    N = 2e4; % number of images to average across
    bootstrap_mode = 1; % only set this to 1 if you want to save responses to disk

        
    %% Run simulations; iterate over all RF sizes and all densities
    for j = 1:21
        fprintf('Running size %i of %i\n',j,length(sxs))    
        
                        
        sx=bem.subunits(1).rf_params.right.sx;
        % Rescale the RF of the model units
        bem = bem.rescale(sxs(j)/sx);
        
        
        
        for k = 1:length(densities);                        
            
            for i = 1:length(all_rds);
                rds = all_rds{i};
                rds.density = densities(k);
                                
                if bootstrap_mode
                    bem = bem.load_bootstrap(rds);
                end

                if ~bem.check_bootstrap(rds) && bootstrap_mode
                    Nb = (N==0)*3e3 + N; % this will just be N if N >0, and 3k if N == 0. 
                    C = bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel);
                else                
                    C = bem.simulate_spatial(rds,N,bootstrap_mode,run_parallel);
                end
                
                resps(j,k,i) = mean(C);
            end
            
        end
    end

    %% Compute normalised responses and plot data

    c=resps(:,:,1);
    hm=resps(:,:,2);
    u=resps(:,:,3);
    
    R = (hm-u)./(c-u);

    a=figure();
    imagesc(1:11,sxs./0.075,R);
    set(gca,'ydir','norm','xtick',1:2:11,'xticklabel',round(densities(1:2:11),2),'fontsize',20)        
    ylabel('\sigma / r','fontsize',22)
    xlabel('Density','fontsize',22)
   
    set(gcf,'color','white')
    
    colormap('hot');
    c=colorbar;
    ylabel(c,'Normalized response');
    set(c,'ticks',0:0.1:0.4);
    savefig(a,'fig4.fig');
    
end