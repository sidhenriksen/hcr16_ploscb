% This will recreate Figure 5 from Henriksen, Cumming, & Read (2016).

function fig5()

    %% Set run-time parameters
    bootstrap_mode = 1;
    silent = 1;
        
    %% Create the BEM neuron
    bem = BEMunit('silent',silent);
    bem.Nx = 292; bem.Ny = 292;
    bem.deg_per_pixel = 0.03;
    bem.dx = 0.1;
    bem.temporal_kernel = 'gamma-cosine';
    bem.outputNL = @(x)(x.^2);
    for j = 1:length(bem.subunits);
        bem.subunits(j).rf_params.left.f=3;
        bem.subunits(j).rf_params.right.f=3;
    end    
    bem = bem.update();    
    bem.dt = 1/200;
    
    sxs = linspace(0.025,0.30,21);    

    %% Set stimulus parameters
    
    densities = logspace(log10(0.04),log10(4),11);
    density = densities(5);
    duration = 10; % 10 second trials
    freqs = [1,10:10:100]; % rds pattern refresh rates
    N = 5e3; % number of trials
    
    %% Create the stimulus; half-matched and correlated RDS objects
    halfmatched_rds = pairedRDS(); halfmatched_rds.dx = bem.dx/bem.deg_per_pixel;
    halfmatched_rds.Nx = 292; halfmatched_rds.Ny = 292;
    halfmatched_rds.dotsize=3; % pixels
    halfmatched_rds.correlation = 0;
    halfmatched_rds.density = density;
    correlated_rds = halfmatched_rds; correlated_rds.correlation = 1;
    
    uncorrelated_rds = halfmatched_rds;
    uncorrelated_rds.toggle_match = 0;
    
    all_rds = {correlated_rds,halfmatched_rds,uncorrelated_rds};

    mother_bem = bem; % this one doesn't change; makes it easier to rescale

    cells_to_run = [2,5,10];

    
    %% Figure 5a
    Ma = zeros(length(cells_to_run),length(freqs)); %init matrix

    fprintf('Running simulations for Figure 5a.\n');
    for j = 1:length(cells_to_run);
        fprintf('Cell %i of %i\n',j,length(cells_to_run));
        k = cells_to_run(j);
        sx = sxs(k);

        bem = mother_bem.rescale(sx/0.1);

        bem = bem.compute_normalization_constant(halfmatched_rds,0,bootstrap_mode);                                                                        
        cbem = bem.load_bootstrap(correlated_rds);
        bem = bem.load_bootstrap(halfmatched_rds);
        ubem = bem.load_bootstrap(uncorrelated_rds);
        


        fprintf('Frequencies: ')
        for f = 1:length(freqs);
            fprintf('%i ',f);
            freq = freqs(f);
            n_frames = ceil(freq*duration);
            HM_bem = zeros(1,N);
            C_bem = zeros(1,N);
            U_bem = zeros(1,N);
            
            parfor i = 1:N;
                HM_bem(i) = mean(bem.simulate_spatiotemporal(halfmatched_rds,n_frames,duration,bootstrap_mode));                
                C_bem(i) = mean(cbem.simulate_spatiotemporal(correlated_rds,n_frames,duration,bootstrap_mode));
                U_bem(i) = mean(ubem.simulate_spatiotemporal(uncorrelated_rds,n_frames,duration,bootstrap_mode));
            end
            

            Ma(j,f) = (mean(HM_bem)-mean(U_bem))./(mean(C_bem)-mean(U_bem));
        end
        fprintf('. Done.\n');

    end
    

    %% Figure 5b
    fs_to_run = [1,10,100];
    Mb = zeros(length(sxs),length(fs_to_run));

    fprintf('Running simulations for Figure 5b.\n')
    for j = 1:length(sxs);
        fprintf('Running cell %i of %i.\n',j,length(sxs));
        sx = bem.subunits(1).rf_params.left.sx;
        bem = bem.rescale(sxs(j)/sx);
        
        if bootstrap_mode            
            if ~bem.check_bootstrap(correlated_rds);
                Nb=5e3;
                bem.simulate_spatial(correlated_rds,Nb,bootstrap_mode);
                bem.simulate_spatial(halfmatched_rds,Nb,bootstrap_mode);
            end
                
            cbem = bem.load_bootstrap(correlated_rds);
            bem = bem.load_bootstrap(halfmatched_rds);
            ubem = bem.load_bootstrap(uncorrelated_rds);
        end

        for f = 1:length(fs_to_run);
            HM_bem = zeros(1,N);            
            C_bem = zeros(1,N);
            U_bem = zeros(1,N);
            
            freq = fs_to_run(f);
            parfor k = 1:N;
                HM_bem(k) = mean(bem.simulate_spatiotemporal(halfmatched_rds,freq*duration,duration,bootstrap_mode));                
                C_bem(k) = mean(cbem.simulate_spatiotemporal(correlated_rds,freq*duration,duration,bootstrap_mode));                
                U_bem(k) = mean(ubem.simulate_spatiotemporal(uncorrelated_rds,freq*duration,duration,bootstrap_mode));                
            end


            Mb(j,f) = (mean(HM_bem)-mean(U_bem))/(mean(C_bem)-mean(U_bem));

        end
    end
    fprintf('Done.\n\n');



    %% Plot data
    sxs_pix = sxs/bem.deg_per_pixel;
    dotsize_pix = round(correlated_rds.dotsize);

    fig=figure();
    subplot(1,2,1); hold on;
    leg_a = cell(1,length(cells_to_run));
    for j = 1:length(cells_to_run);
        rrfsize = round(sxs_pix(cells_to_run(j))/dotsize_pix,1); %relative rf size
        plot(freqs,Ma(j,:),'-','linewidth',3);
        leg_a{j} = sprintf('\\sigma/r=%.1f',rrfsize);
    end
    legend_a=legend(leg_a);
    xlabel('Pattern refresh rate (Hz)','fontsize',16)
    ylabel('Normalized response','fontsize',16);

    subplot(1,2,2); hold on;
    leg_b = cell(1,length(fs_to_run));
    for j = 1:length(fs_to_run);    
        plot(sxs_pix(3:end)/dotsize_pix,Mb(3:end,j),'-','linewidth',3);
        leg_b{j} = sprintf('%i Hz',fs_to_run(j));
    end
    legend_b = legend(leg_b);
    xlabel('RF size/dot size','fontsize',16);
    ylabel('Normalized response','fontsize',16);
    savefig(fig,'fig5.fig');

end