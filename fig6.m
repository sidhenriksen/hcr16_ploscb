function fig6();
    % This function will recreate Figure 6 from Henriksen, Cumming, & Read
    % (2016). These are time-consuming simulations, and will first
    % pre-compute the spatial responses of cells to a large number of
    % random dot stereograms. A large number of trials are then simulated
    % by computing the spatiotemporal response with appropriate noise
    % levels. Finally, a decision model is created.
    %
    % The pre-computed spatial responses will take up some space (roughly
    % 400 MB). Contact Sid @ sid.henriksen@gmail.com for any questions.
    %
    % All code is released under GPL v2, copies of which can be found online.

    %% Define run-time parameters
    run_bootstrap = 1;
    bootstrap_mode = 1;
    build_bootstrap_arrays = 1;

    save_bootstrap_arrays = 0;


    %% Define the parameters of the mother BEMunit; all the other energy model
    %% units are derived from this  mother unit.
    silent = 1;
    bem = BEMunit('silent',silent,'x0',0,'y0',0);
    bem.Nx = 292; bem.Ny=292;
    bem.deg_per_pixel = 0.03;
    bem.temporal_kernel = 'gamma-cosine';
    bem.outputNL = @(x)(x.^2); % squaring output nonlinearity

    % Set temporal properties of bem unit
    bem.tk.tau = 0.035;
    bem.tk.omega = 4;
    bem.dt=1/1000;


    f = 0.3125; % cycles per SD
    %f=0.3;
    % Set frequency of Gabor
    for j = 1:length(bem.subunits);
        bem.subunits(j).rf_params.left.f=f/bem.subunits(j).rf_params.left.sx;
        bem.subunits(j).rf_params.right.f=f/bem.subunits(j).rf_params.right.sx;
    end
    bem = bem.update();


    %% Create the daughter energy model units.
    % There are four types of cells: 
    % near_fine, far_fine, coarse_fine, coarse_far
    % near_fine has negative disparity, small dx/RF.

    % function to map dx to rf size
    dx_to_sx = @(x)(0.023 + abs(x)*0.41);
    sx = bem.subunits(1).rf_params.left.sx;

    fine_near_bem = bem; fine_near_bem.dx = -0.03;
    fine_sx = dx_to_sx(fine_near_bem.dx);
    fine_near_bem = fine_near_bem.rescale(fine_sx/sx);

    coarse_near_bem = bem; coarse_near_bem.dx=-0.48;
    coarse_sx = dx_to_sx(coarse_near_bem.dx);
    coarse_near_bem = coarse_near_bem.rescale(coarse_sx/sx);

    % These are exactly the same except that their disparities are the opposite
    % sign.
    fine_far_bem = fine_near_bem; 
    fine_far_bem.dx = fine_near_bem.dx*-1;
    fine_far_bem = fine_far_bem.update();

    coarse_far_bem = coarse_near_bem; 
    coarse_far_bem.dx = coarse_near_bem.dx*-1;
    coarse_far_bem = coarse_far_bem.update();



    %% Create stimulus generator object
    rds = pairedRDS(); rds.dotsize=3;
    rds.Nx = bem.Nx; rds.Ny = bem.Ny;


    %% Stimulus parameters
    dxs = [-0.48,-0.03,0.03,0.48]./bem.deg_per_pixel;
    correlation_levels = -1:0.2:1; 
    freqs = [85/16,85/4,85/2]; % rds pattern refresh rates
    duration = 1.5; % stimulus duration in seconds

    nt = round(duration/bem.dt); % number of time points
    n_bootstrap = 4e4; % this many samples to save
    run_parallel = 1;
    %% Precompute the monocular responses; this lets us compute a large
    %% number of dynamic RDS responses very quickly
    if run_bootstrap          
        for j = 1:length(correlation_levels);
            rds.correlation = correlation_levels(j);
            fprintf('Correlation %i of %i\n',j,length(correlation_levels));
            for k = 1:length(dxs);
                rds.dx = dxs(k);

                temp_bem = coarse_near_bem.load_bootstrap(rds);
                if isfield(temp_bem.subunits(1),'V_L');
                    n_temp = length(temp_bem.subunits(1).V_L);
                else
                    n_temp = 0;
                end
                if n_temp < n_bootstrap
                    Nb = n_bootstrap-n_temp;

                    % the seed guarantees that they see the same image sequences
                    seed = randi(1e9,1);
                    fine_far_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel,seed);
                    fine_near_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel,seed);
                    coarse_far_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel,seed);
                    coarse_near_bem.simulate_spatial(rds,Nb,bootstrap_mode,run_parallel,seed);
                end
            end
        end
    end


    %% Build bootstrap arrays
    % First we build the bootstrap arrays. We make the decision model by
    % resampling responses of cells to a sequence of RDSs. The rationale for
    % doing this is that generating the stimuli, computing the cells'
    % spatiotemporal response and generating the decision would take forever
    % with 200 cells. Instead, we get the distribution of responses for each
    % cell per trial (given some noise model), assume that their responses 
    % are independent (as they see
    % independent parts of the image), and take random picks from this
    % distribution. Once we've built the bootstrap arrays, we can quickly find 
    % the optimal number of cells + noise combination that fits the observed
    % data.
    n_bootstrap_trials = 4e4; % number of samples in bootstrap distribution
    conditions = CombVec(freqs,correlation_levels,dxs);
    noise_levels = [0,linspace(30,100,40)];

    n_bootstrap=4e4;
    if build_bootstrap_arrays
        ff_resp = zeros(size(conditions,2),n_bootstrap_trials,length(noise_levels));
        fn_resp = zeros(size(conditions,2),n_bootstrap_trials,length(noise_levels));
        cn_resp = zeros(size(conditions,2),n_bootstrap_trials,length(noise_levels));
        cf_resp = zeros(size(conditions,2),n_bootstrap_trials,length(noise_levels));

        norm_constants = zeros(size(conditions,2),4);
        completed_freqs = zeros(3,4);
        for j = 1:size(conditions,2);
            fprintf('Condition %i of %i\n',j,size(conditions,2));
            freq = conditions(1,j);
            c = conditions(2,j);
            dx = conditions(3,j);
            n_frames = round(duration*freq);
            rds.correlation = c;
            rds.dx = dx;


            % these indices ensure that the only difference between near and
            % far cells (neuron-antineuron pairs) is the difference in
            % disparity
            fine_indices = randi(n_bootstrap,1,n_frames*n_bootstrap_trials);
            fine_indices = reshape(fine_indices,[n_frames,n_bootstrap_trials]);

            coarse_indices = randi(n_bootstrap,1,n_frames*n_bootstrap_trials);
            coarse_indices = reshape(coarse_indices,[n_frames,n_bootstrap_trials]);

            % Load the data for the subunits
            bems = {fine_far_bem,fine_near_bem,coarse_far_bem,coarse_near_bem};                
            for i = 1:length(bems);
                bem = bems{i};

                if c == -1 && j <= 3;                                
                    fprintf('freq: %.2f\n', freq);
                    crds = rds;
                    crds.correlation = 1;

                    crds.dx = bem.dx/bem.deg_per_pixel;

                    cbem = bem.load_bootstrap(crds);

                    this_f = conditions(1,:) == freq;

                    resps = zeros(1,n_bootstrap_trials);
                    for k = 1:n_bootstrap_trials;
                        resps(k) = mean(cbem.simulate_spatiotemporal(crds,n_frames,duration,bootstrap_mode));
                    end
                    norm_constants(this_f,i) = mean(resps);
                end


                bem = bem.load_bootstrap(rds);

                current_norm = norm_constants(j,i);
                parfor k = 1:n_bootstrap_trials;

                    if i < 3
                        idx=fine_indices(:,k);
                    else
                        idx=coarse_indices(:,k);
                    end

                    C = bem.simulate_spatiotemporal(rds,n_frames,duration,bootstrap_mode,idx);
                    C = C/current_norm;

                    Cn = zeros(1,length(noise_levels));

                    for n_i = 1:length(noise_levels);                    
                        Cn(n_i) = sum(C + randn(size(C)).*sqrt(C).*noise_levels(n_i));                    
                    end

                    switch i 
                        case 1
                            ff_resp(j,k,:) = Cn;
                        case 2
                            fn_resp(j,k,:) = Cn;
                        case 3
                            cf_resp(j,k,:) = Cn;
                        case 4
                            cn_resp(j,k,:) = Cn;
                    end
                end

            end

        end

        if save_bootstrap_arrays
            save('fig6_bootstrap_data.mat','ff_resp','fn_resp','cf_resp','cn_resp')
        end
    else
        load('fig6_bootstrap_data.mat');
    end




    all_n_cells = [40]; all_ps = [1.0];

    N = 2e4;


    %% Run the decision model
    fprintf('Running conditions...')
    sub_dims = [length(freqs),length(correlation_levels),length(dxs)];
    Psi = zeros([length(all_n_cells),length(noise_levels),sub_dims]);
    for j = 1:size(conditions,2);   
        [f,g,k] = ind2sub(sub_dims,j);


        for n_i = 1:length(noise_levels);                
            % Get the bootstrap distribution for the current condition
            current_ff = squeeze(ff_resp(j,:,n_i));
            current_fn = squeeze(fn_resp(j,:,n_i));
            current_cf = squeeze(cf_resp(j,:,n_i));
            current_cn = squeeze(cn_resp(j,:,n_i));

            for nk = 1:length(all_n_cells);

                n_cells = all_n_cells(nk);



                idx = randi(length(current_ff),[n_cells,N]);
                ff = current_ff(idx) - current_fn(idx);                        
                fn = -ff;

                idx = randi(length(current_ff),[n_cells,N]);            
                cf = current_cf(idx) - current_cn(idx);            
                cn = -cf;

                near = [fn;cn];

                psi = sum(near) > 0;

                Psi(nk,n_i,f,g,k) = mean(psi);

            end
        end
    end
    fprintf('Done.\n');








    %i =1; k = 6;


    clear PsiM
    PsiM(:,:,:,2) = squeeze(Psi(i,:,:,:,1)+(1-Psi(i,:,:,:,4)))/2;
    PsiM(:,:,:,1) = squeeze(Psi(i,:,:,:,2)+(1-Psi(i,:,:,:,3)))/2;
    brown = [0.6,0.3,0.1];
    %

    psi_fine = squeeze(PsiM(j,2,:,1)); psi_coarse = squeeze(PsiM(j,2,:,2));
    psi_slow = squeeze(PsiM(j,1,:,1)); psi_fast = squeeze(PsiM(j,3,:,1));

    psi_func  = @(x,params)(1-(1-params(1))*exp(-(x/params(2)).^params(3)));
    x = linspace(0,1,501); x2 = (x-0.5)*2;
    params_fine = fit_psych(psi_fine); params_coarse = fit_psych(psi_coarse);
    params_slow = fit_psych(psi_slow); params_fast = fit_psych(psi_fast);
    fit_fine = psi_func(x,params_fine);
    fit_coarse = psi_func(x,params_coarse);
    fit_slow = psi_func(x,params_slow);
    fit_fast = psi_func(x,params_fast);

    fig=figure(); 
    sub1=subplot(1,2,1); hold on;
    plot(correlation_levels,psi_fine,'g o','markerfacecolor','g','linewidth',3,'markersize',8);
    plot(correlation_levels,psi_coarse,'o','color',brown,'markerfacecolor',brown,'linewidth',3,'markersize',8);

    L1=legend('\Deltax=\pm 0.03^\circ','\Deltax=\pm 0.48^\circ');
    set(L1,'location','southeast');
    plot(x2,fit_fine,'g -','linewidth',2);
    plot(x2,fit_coarse,'-','linewidth',2,'color',brown);

    xlabel('Correlation');
    ylabel('Proportion correct');
    xlim([-1,1]); ylim([0,1]);
    plot([-1,1],[0.5,0.5],'k --')
    set(gca,'ytick',0:0.25:1);

    red = [0.8,0.1,0.1];
    sub2=subplot(1,2,2); hold on;
    plot(correlation_levels,psi_slow,'o' ,'color',red,'markerfacecolor',red,'linewidth',3,'markersize',8);
    plot(correlation_levels,psi_fast,'b o','markerfacecolor','b','linewidth',3,'markersize',8);
    L2=legend('5.3 Hz', '42.5 Hz');
    set(L2,'location','southeast')
    plot(x2,fit_slow,'-','linewidth',2,'color',red);
    plot(x2,fit_fast,'b -','linewidth',2)


    xlabel('Correlation');
    ylabel('Proportion correct');
    xlim([-1,1]); ylim([0,1]);
    set(gca,'ytick',0:0.25:1);
    plot([-1,1],[0.5,0.5],'k --')


    set_plot_params(gcf);



    A_dx = [fractional_area(fit_fine),fractional_area(fit_coarse)];
    a = figure(); hold on;    
    plot([0.03,0.48],A_dx,'k o -','markerfacecolor','k','linewidth',2,...
        'markersize',8);

    xlabel('\Deltax');
    ylabel(sprintf('Fractional\n area'));
    ylim([0.0,1.0]); xlim([0.03-0.05,0.48+0.05]);
    set(gca,'ytick',0:0.5:1,'xtick',[0.03,0.48]);

    %set_plot_params(a);

    inset_size = 0.15;
    inset_handle = gca;
    new_handle = copyobj(inset_handle,fig);
    close(a);
    ax=get(sub1,'Position');
    set(new_handle,'Position',[1.8*ax(1)+ax(3)-inset_size 0.1*ax(2)+ax(4)-inset_size inset_size*0.5 inset_size])



    subplot(1,2,2);
    % Make inset for subplot 2
    a = figure(); hold on;    
    A_hz = [fractional_area(fit_slow),fractional_area(fit_fast)];
    plot([85/16,85/2],A_hz,'k o -','markerfacecolor','k','linewidth',2,...
        'markersize',8);

    xlabel('Hz');
    ylabel(sprintf('Fractional\n area'));
    ylim([0.2,0.5]); 
    set(gca,'ytick',[0.2,0.35,0.5],'xtick',round([85/16,85/2],2),'xscale','linear');
    xlim([0,85/2+85/16])
    %set_plot_params(a);

    inset_size = 0.15;
    inset_handle = gca;
    new_handle = copyobj(inset_handle,fig);
    close(a);
    ax=get(sub2,'Position');
    set(new_handle,'Position',[1.2*ax(1)+ax(3)-inset_size 0.1*ax(2)+ax(4)-inset_size inset_size*0.5 inset_size])
end