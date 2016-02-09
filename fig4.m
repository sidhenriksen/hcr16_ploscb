% This will recreate Figure 4 from Henriksen, Cumming, & Read (2016).
function fig4()
    %%  Set run-time parameters
    run_parallel = 1;
    silent = 1;
    
    %% Create BEM unit and set parameters
    bem = BEMunit('silent',silent);
    bem.memory_threshold=bem.memory_threshold*0.5;
    bem.Nx = 292; bem.Ny = 292;
    bem.deg_per_pixel = 0.03;
    bem.dx = 0.09;
    bem.outputNL = @(x)(x.^2);
    for j = 1:length(bem.subunits);
        bem.subunits(j).rf_params.left.f=3;
        bem.subunits(j).rf_params.right.f=3;
    end
    
    %% Set temporal kernel properties
    bem.temporal_kernel = 'gamma-cosine';
    bem.tk.tau = 0.035;
    bem.tk.omega = 4;
    bem.dt = 1/1000;
    
    bem = bem.update();
    

    %% Create correlated, half-matched and uncorrelated stimulus objects
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
    sxs = linspace(0.025,0.30,21); % RF sizes
    densities = logspace(log10(0.04),log10(4),11);
    freqs = [1,10:10:100]; % dot pattern refresh rates
    duration = 10; % stimulus duration

    %% Pre-allocate space for the figure A and B matrices
    respsA = zeros(length(sxs),length(densities),3);
    respsB = zeros(length(sxs),length(freqs),3);
    
    %% 
    bootstrap_mode = 1; % only set this to 1 if you want to save responses to disk
    
    % Figure A
    N = 2e4; % number of images to average across
    
    % Figure B
    N_trials= 5e3;      
    
    %% Run simulations; iterate over all RF sizes and all densities
    for j = 1:21
        sx=bem.subunits(1).rf_params.right.sx;
        % Rescale the RF of the model units
        bem = bem.rescale(sxs(j)/sx);        
        
        fprintf('%i. Running - 4a ',j);        
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
                
                respsA(j,k,i) = mean(C);
            end
            
        end
        
        
        k = 5;
        fprintf('- 4b');
        for i = 1:length(all_rds);
            rds = all_rds{i};
            rds.density = densities(k);
            
            bem = bem.load_bootstrap(rds);
            
            for f = 1:length(freqs);
                freq = freqs(f);
                n_frames = duration*freq;
                
                
                Cs = zeros(1,N_trials);
                parfor rep = 1:N_trials;
                    C = bem.simulate_spatiotemporal(rds,n_frames,duration,bootstrap_mode);
                    Cs(rep) = mean(C);
                end
                
                respsB(j,f,i) = mean(Cs);
            end
        end
        fprintf('. Done.\n');
        
    end

    %% Compute normalised responses and plot data

    % Normalized responses for 4a
    c=respsA(:,:,1);
    hm=respsA(:,:,2);
    u=respsA(:,:,3);   
    Ra = (hm-u)./(c-u);

    % Normalized responses for 4b
    c=respsB(:,:,1);
    hm=respsB(:,:,2);
    u=respsB(:,:,3);
    Rb = (hm-u)./(c-u);

    xa = 1:length(densities);
    figure();
    subplot(1,2,1);
    imagesc(xa,sxs./0.075,Ra');
    xlabel('Density','fontsize',20)
    ylabel('\sigma/r','fontsize',20);    
    cbar1=colorbar;
    
    set(cbar1,'ticks',0:0.1:0.4,'fontsize',16);
    
    subplot(1,2,2);
    imagesc(freqs,sxs./0.075,Rb')
    xlabel('Dot pattern refresh rate (Hz)','fontsize',20)
    ylabel('\sigma/r','fontsize',20);   
    
    cbar2=colorbar;
    set(cbar2,'ticks',0:0.1:0.4,'fontsize',16);
    ylabel(cbar2,'Normalized response','fontsize',20);
    set_plot_params(gcf)
    
    
    subplot(1,2,1);
    set(gca,'clim',[0,0.4],'ydir','norm','xtick',1:3:length(densities),'xticklabels',round(densities(1:3:length(densities)),2));
    subplot(1,2,2);
    set(gca,'clim',[0,0.4],'ydir','norm','xtick',10:30:100);
    
    colormap('hot');
    
    
    savefig(gcf,'fig4.fig');
end