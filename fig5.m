
% Generating Figure 5 from Henriksen, Cumming, & Read (2016).
% Tuning curves for the different cells.
% Note: This takes a while to run.
function fig5(N)

    if nargin < 1
        N = 500; 
    end

    bootstrap_mode = 0;
    run_parallel = 1;
    
    dx_to_sx = @(x)(0.023 + abs(x)*0.41);
    dx = 0.03;    
    bem = BEMunit('dx',dx);
    bem.deg_per_pixel = 0.03/2;
    bem.Nx = 292*2; bem.Ny=292*2;
    bem.outputNL = @(x)(x.^2);
    bem = bem.update();
    

    mother_bem = bem;
    m_sx = mother_bem.subunits(1).rf_params.left.sx;

    rds = RDS(); 
    rds.Nx = bem.Nx; rds.Ny = bem.Ny;
    dxs = [-0.48,-0.03,0.03,0.48];
    stim_dxs = [-68:4:-12,-10:1:10,12:4:68];
    tcs = zeros(length(dxs),length(stim_dxs));
    
    
    for k = 1:length(dxs);
        bem = mother_bem;
        bem.dx = dxs(k);
        bem = bem.rescale(dx_to_sx(dxs(k))/m_sx);
        

        for j = 1:length(stim_dxs);
            rds.dx = stim_dxs(j);
            if bootstrap_mode && ~N
                bem = bem.load_bootstrap(rds);
            end
            
            tcs(k,j) = mean(bem.simulate_spatial(rds,N,bootstrap_mode,run_parallel));

        end
    end

    
    fig=figure(); hold on;
    red=[0.8,0.1,0.1];
    green=[0.1,0.8,0.1];
    cols = {red,green,green,red};
    ls={'-','-','--','--'}';
    for j = 1:length(dxs);
     plot(stim_dxs*bem.deg_per_pixel,tcs(j,:)./max(tcs(j,:)),'linewidth',3,...
         'linestyle',ls{j},'color',cols{j});
    end

    xlabel('Disparity (deg)','fontsize',16);
    ylabel('Normalized response','fontsize',16);
    xlim([-1,1]);
    ylim([0,1]);
    set_plot_params(fig)
    savefig(fig,'fig5.fig')

end