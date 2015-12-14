clear all;

% If you want to save the responses to disk (note: this takes up some
% space)
run_bootstrap = 0;
bootstrap_mode = 1;
compute_responses = 1;

% Set parameters
silent = 1;
bem = BEMunit('silent',silent);
bem.Nx = 292; bem.Ny=292;
bem.deg_per_pixel = 0.03;
bem.temporal_kernel = 'gamma-cosine';
bem.outputNL = @(x)(x.^2); % squaring output nonlinearity

% Set temporal properties of bem unit
bem.tk.tau = 0.035;
bem.tk.omega = 0.004;

bem = bem.update();

% function to map dx to rf size
dx_to_sx = @(x)(0.023 + abs(x)*0.4);


sx = bem.subunits(1).rf_params.left.sx;

% There are four types of cells: 
% near_fine, far_fine, coarse_fine, coarse_far
% near_fine has negative disparity, small dx/RF.


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



% Create stimulus generator object
rds = pairedRDS();
rds.Nx = bem.Nx; rds.Ny = bem.Ny;



dxs = [-0.48,-0.03,0.03,0.48]./bem.deg_per_pixel;
correlation_levels = -1:0.2:1;

N = 10; % number of trials

duration = 1.5;

n_cells = 50;

freqs = [85/16,85/4,85/2];


if run_bootstrap
    n_bootstrap = 2500; % this many samples to save
    for j = 1:length(correlation_levels);
        rds.correlation = correlation_levels(j);
        for k = 1:length(dxs);
            rds.dx = dxs(k);
            fine_far_bem.simulate_spatial(rds,n_bootstrap,bootstrap_mode);
            fine_near_bem.simulate_spatial(rds,n_bootstrap,bootstrap_mode);
            coarse_far_bem.simulate_spatial(rds,n_bootstrap,bootstrap_mode);
            coarse_near_bem.simulate_spatial(rds,n_bootstrap,bootstrap_mode);

        end
    end
end
% 
run_norm = ~fine_far_bem.check_normalization_constant(rds);
if run_norm;
    fprintf('Calculating normalisation constants... ');
    fine_far_bem = fine_far_bem.compute_normalization_constant(rds);
    fine_near_bem = fine_near_bem.compute_normalization_constant(rds);
    coarse_near_bem = coarse_near_bem.compute_normalization_constant(rds);
    coarse_far_bem = coarse_far_bem.compute_normalization_constant(rds);
    fprintf('Done.\n');
end


%ff = fine far; fn = fine near; cn = coarse near; cf = coarse far
ff_resp = zeros(length(correlation_levels),length(dxs),length(freqs),N);
fn_resp = zeros(length(correlation_levels),length(dxs),length(freqs),N);
cn_resp = zeros(length(correlation_levels),length(dxs),length(freqs),N);
cf_resp = zeros(length(correlation_levels),length(dxs),length(freqs),N);
Psi = zeros(length(correlation_levels),length(dxs),length(freqs));

noise = 0;%7.333;
duration = 1.5;

n_bootstrap = 1000;

n_cells = 50;
K = 30; % top 30




rnorm = ones(4,length(freqs));


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
for j = 1:length(correlation_levels);
    fprintf('Correlation %i of %i\n',j,length(correlation_levels));
    rds.correlation = correlation_levels(j);
    for k = 1:length(dxs);
         % load bootstrap
        rds.dx = dxs(k);
        
        fine_far_bem = fine_far_bem.load_bootstrap(rds);
        fine_near_bem = fine_near_bem.load_bootstrap(rds);
        coarse_far_bem = coarse_far_bem.load_bootstrap(rds);
        coarse_near_bem = coarse_near_bem.load_bootstrap(rds);
         
        for f = 1:length(freqs);
            n_frames = round(duration*freqs(f))+1;
            
            
            parfor i = 1:n_bootstrap;
                C_ff = fine_far_bem.simulate_spatiotemporal(rds,n_frames,...
                    duration,bootstrap_mode)./rnorm(1,f);
                C_fn = fine_near_bem.simulate_spatiotemporal(rds,n_frames,...
                    duration,bootstrap_mode)./rnorm(2,f);
                
                C_cf = coarse_far_bem.simulate_spatiotemporal(rds,n_frames,...
                    duration,bootstrap_mode)./rnorm(3,f);
                C_cn = coarse_near_bem.simulate_spatiotemporal(rds,n_frames,...
                    duration,bootstrap_mode)./rnorm(4,f);
                
                % Noise variance proportional to its current value
                ff_resp(j,k,f,i) = sum(C_ff + noise .* randn(size(C_ff)).*sqrt(C_ff));
                fn_resp(j,k,f,i) = sum(C_fn + noise .* randn(size(C_fn)).*sqrt(C_fn));
                cf_resp(j,k,f,i) = sum(C_cf + noise .* randn(size(C_cf)).*sqrt(C_cf));
                cn_resp(j,k,f,i) = sum(C_cn + noise .* randn(size(C_cn)).*sqrt(C_cn));

            end
            
                       
        end
        
    end   
end
fprintf('Finished building bootstrap arrays.\n');

N = 500;
dims = size(ff_resp);
Psi = zeros(dims(1:3));
for j = 1:length(correlation_levels);
    for k = 1:length(dxs);
        for f = 1:length(freqs);
            
            current_ff = squeeze(ff_resp(j,k,f,:));
            current_fn = squeeze(fn_resp(j,k,f,:));
            current_cf = squeeze(cf_resp(j,k,f,:));
            current_cn = squeeze(cn_resp(j,k,f,:));
            psi = zeros(1,N);
            for i = 1:N
                ff = randsample(current_ff,n_cells,1)-randsample(current_fn,n_cells,1);
                fn = -ff;
                cf = randsample(current_cf,n_cells,1)-randsample(current_cn,n_cells,1);
                cn = -cf;
                
                far = [ff;cf];
                near = [fn;cn];
                far_sorted = sort(far,'descend');
                near_sorted = sort(near,'descend');
                psi(i) = sum(near_sorted(1:K)) > sum(far_sorted(1:K));                                
            end
            Psi(j,k,f) = mean(psi);
            
        end
    end
end
PsiM = zeros(length(correlation_levels),2,length(freqs));
PsiM(:,1,:) = (squeeze(Psi(:,2,:) + 1-Psi(:,3,:) ))/2;
PsiM(:,2,:) = (squeeze(Psi(:,1,:) + 1-Psi(:,4,:) ))/2;


brown = [0.6,0.3,0.1];
figure(); subplot(1,2,1); hold on;
plot(correlation_levels,PsiM(:,1,2),'g o -','markerfacecolor','g','linewidth',3,'markersize',8);
plot(correlation_levels,PsiM(:,2,2),'o -','color',brown,'markerfacecolor',brown,'linewidth',3,'markersize',8);
xlabel('Correlation');
ylabel('Proportion correct');
xlim([-1,1]); ylim([0,1]);

red = [0.8,0.1,0.1];
subplot(1,2,2); hold on;
plot(correlation_levels,PsiM(:,1,1),'b o -','markerfacecolor','b','linewidth',3,'markersize',8);
plot(correlation_levels,PsiM(:,1,3),'o -' ,'color',red,'markerfacecolor',red,'linewidth',3,'markersize',8);
xlabel('Correlation');
ylabel('Proportion correct');
xlim([-1,1]); ylim([0,1]);