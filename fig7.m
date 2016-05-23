% This recreates Figure 7 from Henriksen, Cumming, & Read (2016)

run_bootstrap = 0;
build_arrays = 1;
bootstrap_mode = 1;

%% Define the parameters of the mother BEMunit; all the other energy model
%% units are derived from this  mother unit.
silent = 1;
bem = BEMunit('silent',silent,'x0',0,'y0',0,'bootstrap_dir','/sid/Modelling/hcr16_ploscb/fig7_data/');
bem.Nx = 292; bem.Ny=292;
bem.deg_per_pixel = 0.03;
bem.temporal_kernel = 'gamma-cosine';
bem.outputNL = @(x)(x.^2); % squaring output nonlinearity

% Set temporal properties of bem unit
bem.tk.tau = 0.035;
bem.tk.omega = 4;
bem.dt = 1/1000;


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

% The reason for this peculiar scaling is that the experiment was done with
% dotsize multiples of 0.075. In our simulation, a pixel width is 0.03 deg,
% so in order to get this to work, we rescale the RF slightly, such
% that the model is effectively seeing dots with a raidus of 0.025 deg
% rather than the 0.03 deg.
fine_near_bem = bem; fine_near_bem.dx = -0.03;
fine_sx = dx_to_sx(fine_near_bem.dx);
fine_near_bem = fine_near_bem.rescale(fine_sx/sx * 0.03/0.025);

coarse_near_bem = bem; coarse_near_bem.dx=-0.48;
coarse_sx = dx_to_sx(coarse_near_bem.dx);
coarse_near_bem = coarse_near_bem.rescale(coarse_sx/sx * 0.03/0.025);

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
rds.correlation = 0;
crds = rds;
crds.correlation = 1;

%% Stimulus parameters
dotsizes = [1,2,3];
dxs = [-0.03,0.03]./bem.deg_per_pixel;
freq = 85/4;
duration = 0.5; % stimulus duration in seconds
n_frames = ceil(freq*duration);

nt = round(duration/bem.dt); % number of time points
n_bootstrap = 2e4; % this many samples to save
run_parallel = 1;

n_bootstrap_trials = 2e4;

bems = {fine_far_bem,fine_near_bem,coarse_far_bem,coarse_near_bem};
%% 

if run_bootstrap
    fprintf('Running bootstrap.\n');
    for ds = 1:length(dotsizes);        
        norms = zeros(1,length(bems));

        rds.dotsize = dotsizes(ds);
        crds.dotsize = dotsizes(ds);

        idx = randi(1e9);


        %% Compute the spatial responses
        for i = 1:length(bems);
            for k = 1:length(dxs)
                rds.dx = dxs(k);
                bem = bems{i};        
                bem.simulate_spatial(rds,n_bootstrap,bootstrap_mode,run_parallel,idx);
            end

            crds.dx = bem.dx/bem.deg_per_pixel;
            bem.simulate_spatial(crds,n_bootstrap,bootstrap_mode,run_parallel);

        end
    end
end



%%

noise_levels = [0,linspace(10,50,30)];
if build_arrays
    ff_resp = zeros(length(dotsizes),length(dxs),n_bootstrap_trials,length(noise_levels));
    fn_resp = zeros(length(dotsizes),length(dxs),n_bootstrap_trials,length(noise_levels));
    cf_resp = zeros(length(dotsizes),length(dxs),n_bootstrap_trials,length(noise_levels));
    cn_resp = zeros(length(dotsizes),length(dxs),n_bootstrap_trials,length(noise_levels));
    %%
    norms = zeros(1,length(bems));
    for ds = 1:length(dotsizes)
        fprintf('Running dotsize %i of %i\n',ds,length(dotsizes));
        rds.dotsize = dotsizes(ds);
        crds.dotsize = dotsizes(ds);
        for i = 1:length(bems);
            bem = bems{i};

            
            
            crds.dx = bem.dx/bem.deg_per_pixel;                
            bem = bem.load_bootstrap(crds);
            Cs = zeros(1,n_bootstrap_trials);            
            parfor j = 1:n_bootstrap_trials;
                C = bem.simulate_spatiotemporal(crds,n_frames,duration,bootstrap_mode);
                Cs(j) = mean(C);
            end
            norms(i) = mean(Cs);

            current_norm = norms(i);

            for k = 1:length(dxs);
                rds.dx = dxs(k);

                bem = bem.load_bootstrap(rds);

                fine_indices = randi(n_bootstrap,1,n_frames*n_bootstrap_trials);
                fine_indices = reshape(fine_indices,[n_frames,n_bootstrap_trials]);

                coarse_indices = randi(n_bootstrap,1,n_frames*n_bootstrap_trials);
                coarse_indices = reshape(coarse_indices,[n_frames,n_bootstrap_trials]);
                parfor trial = 1:n_bootstrap_trials;

                    if i < 3;
                        idx = fine_indices(:,trial);
                    else
                        idx = coarse_indices(:,trial);
                    end

                    C = bem.simulate_spatiotemporal(rds,n_frames,duration,bootstrap_mode,idx);
                    C = C/current_norm;
                    Cn = zeros(1,length(noise_levels));

                    for n_i = 1:length(noise_levels);
                        Cn(n_i) = sum(C + randn(size(C)).*sqrt(C).*noise_levels(n_i));
                    end

                    switch i 
                        case 1
                            ff_resp(ds,k,trial,:) = Cn;
                        case 2
                            fn_resp(ds,k,trial,:) = Cn;
                        case 3
                            cf_resp(ds,k,trial,:) = Cn;
                        case 4
                            cn_resp(ds,k,trial,:) = Cn;
                    end
                end

            end
        end
    end
end

%%
n_cells = 40;
    
Psi = zeros(length(dotsizes),length(dxs),length(noise_levels));
N = 5e3;
for ds = 1:length(dotsizes);
    for k = 1:length(dxs);
        for n_i = 1:length(noise_levels);
            current_ff = squeeze(ff_resp(ds,k,:,n_i));
            current_fn = squeeze(fn_resp(ds,k,:,n_i));
            current_cf = squeeze(cf_resp(ds,k,:,n_i));
            current_cn = squeeze(cn_resp(ds,k,:,n_i));


            idx = randi(length(current_ff),[n_cells,N]);
            ff = current_ff(idx) - current_fn(idx);                        
            fn = -ff;

            idx = randi(length(current_ff),[n_cells,N]);            
            cf = current_cf(idx) - current_cn(idx);            
            cn = -cf;

            near = [fn;cn];
            
            psi = sum(near)>0;

            Psi(ds,k,n_i) = mean(psi);
        end
    end
end


%% Generate plots
k=8;

x = dotsizes*0.025;
P = dotsize_psych_analysis();
P2 = squeeze(P(1,1,:,:));
cols = rand([3,size(P2,2)]);
P2m = mean(P2,2);
figure(); 
subplot(1,2,1); hold on;
for j = 1:size(P2,2);
    plot(x,P2(:,j),'o -','markerfacecolor',cols(:,j),'color',cols(:,j),'linewidth',2,...
        'markersize',6);
end
clear xlim ylim;
xlim([0.015,0.085]); ylim([0,1]);
xlabel('Dot size (deg)');
L1= legend('AD','DS','KM','SH');
set(L1,'location','southeast');
set(gca,'xtick',x,'ytick',0:0.25:1,'linewidth',2)
plot([0,0.1],[0.5,0.5],'k --','linewidth',2);
ylabel('Proportion correct');

subplot(1,2,2);
PsiM = (Psi(:,1,:) + (1-Psi(:,2,:)))/2;
ExpN = 160;
[Lexp,Uexp] = BinoConf_Score(P2m*ExpN,ExpN);
Lexp = P2m-Lexp; Uexp = Uexp-P2m;
hold on;

E1 = errorbar(x,P2m,Lexp,Uexp);
set(E1,'color','k','markersize',6,'marker','o','linestyle','None',...
'markerfacecolor','k','linewidth',3);
a=plot(x,P2m,'k - o','markersize',6,'markerfacecolor','k','linewidth',2);
b=plot(x,PsiM(:,k),'m o -','markersize',6,'markerfacecolor','m','linewidth',2);
L2=legend([a,b],'Average of 4 subjects','Model');


xlim([0.015,0.085]); ylim([0,1]);
set(gca,'xtick',x,'ytick',0:0.25:1,'linewidth',2);
xlabel('Dot size (deg)');


set(L2,'location','southeast')
plot([0,0.1],[0.5,0.5],'k --','linewidth',2);
set_plot_params(gcf,'fs_scalar',0.4)

savefig(gcf,'fig7.fig')
%end