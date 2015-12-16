clear all;

% If you want to save the responses to disk (note: this takes up some
% space)
run_bootstrap = 1;
bootstrap_mode = 1;
compute_responses = 1;
build_bootstrap_arrays = 1;

% Set parameters
silent = 1;
bem = BEMunit('silent',silent);
bem.Nx = 292; bem.Ny=292;
bem.deg_per_pixel = 0.03;
bem.temporal_kernel = 'gamma-cosine';
bem.outputNL = @(x)(x.^2); % squaring output nonlinearity

% Set temporal properties of bem unit
bem.tk.tau = 0.035;
bem.tk.omega = 4;

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
    n_bootstrap = 5000; % this many samples to save
    for j = 1:length(correlation_levels);
        rds.correlation = correlation_levels(j);
        for k = 1:length(dxs);
            rds.dx = dxs(k);
            
            temp_bem = fine_far_bem.load_bootstrap(rds);
            n_temp = length(temp_bem.subunits(1).V_L);
            if n_temp < n_bootstrap
                N = n_bootstrap-n_temp;
                fine_far_bem.simulate_spatial(rds,N,bootstrap_mode);
                fine_near_bem.simulate_spatial(rds,N,bootstrap_mode);
                coarse_far_bem.simulate_spatial(rds,N,bootstrap_mode);
                coarse_near_bem.simulate_spatial(rds,N,bootstrap_mode);
            end

        end
    end
end

run_norm = ~fine_far_bem.check_normalization_constant(rds);
if run_norm;
    fprintf('Calculating normalisation constants... ');
    fine_far_bem = fine_far_bem.compute_normalization_constant(rds,0,bootstrap_mode);
    fine_near_bem = fine_near_bem.compute_normalization_constant(rds,0,bootstrap_mode);
    coarse_near_bem = coarse_near_bem.compute_normalization_constant(rds,0,bootstrap_mode);
    coarse_far_bem = coarse_far_bem.compute_normalization_constant(rds,0,bootstrap_mode);
    fprintf('Done.\n');
end

noise = 2.73;
duration = 1.5;

n_bootstrap = 1000;


nt = round(duration/bem.dt);
%ff = fine far; fn = fine near; cn = coarse near; cf = coarse far
ff_resp = zeros(length(correlation_levels),length(dxs),length(freqs),N,nt);
fn_resp = zeros(length(correlation_levels),length(dxs),length(freqs),N,nt);
cn_resp = zeros(length(correlation_levels),length(dxs),length(freqs),N,nt);
cf_resp = zeros(length(correlation_levels),length(dxs),length(freqs),N,nt);
Psi = zeros(length(correlation_levels),length(dxs),length(freqs));


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
if build_bootstrap_arrays;
    for j = 1:length(correlation_levels);
        fprintf('Correlation %i of %i\n',j,length(correlation_levels));
        rds.correlation = correlation_levels(j);
        for k = 1:length(dxs);
             % load bootstrap
            rds.dx = dxs(k);

            bems = {fine_far_bem,fine_near_bem,coarse_far_bem,coarse_near_bem};

            for i = 1:length(bems);
                bem = bems{i};
                try
                    bem = bem.load_bootstrap(rds);
                catch
                    fprintf('File corrupt; remaking bootstrap file.\n');
                    fname = num2str(bem.get_identifier());
                    in_dir = bem.bootstrap_dir;
                    delete([in_dir,fname,'.mat']);
                    delete([in_dir,fname,'.idx'])
                    bem.simulate_spatial(rds,n_bootstrap,bootstrap_mode);
                    bem = bem.load_bootstrap(rds);
                end

                switch i
                    case 1
                        fine_far_bem = bem;
                    case 2
                        fine_near_bem = bem;
                    case 3
                        coarse_far_bem = bem;
                    case 4
                        coarse_near_bem = bem;
                end

            end


            for f = 1:length(freqs);
                n_frames = round(duration*freqs(f));


                parfor i = 1:n_bootstrap;
                    C_ff = fine_far_bem.simulate_spatiotemporal(rds,n_frames,...
                        duration,bootstrap_mode);
                    C_fn = fine_near_bem.simulate_spatiotemporal(rds,n_frames,...
                        duration,bootstrap_mode);                
                    C_cf = coarse_far_bem.simulate_spatiotemporal(rds,n_frames,...
                        duration,bootstrap_mode);
                    C_cn = coarse_near_bem.simulate_spatiotemporal(rds,n_frames,...
                        duration,bootstrap_mode);


                    % Noise variance proportional to its current value
                    ff_resp(j,k,f,i,:) = C_ff;
                    fn_resp(j,k,f,i,:) = C_fn;
                    cf_resp(j,k,f,i,:) = C_cf;
                    cn_resp(j,k,f,i,:) = C_cn;

                end


            end

        end   
    end
    save('fig7_bootstrap_arrays.mat','ff_resp','fn_resp','cf_resp','cn_resp');
    fprintf('Finished building bootstrap arrays.\n');
else
    fprintf('Loading bootstrap arrays... ')
    load('fig7_bootstrap_arrays.mat');
    fprintf('Done.\n');
end


all_n_cells = 40:10:60; all_p = [0.25,0.5]; all_noise_levels = linspace(0,5,6);
conds = CombVec(all_n_cells,all_p,all_noise_levels);

N = 500;
dims = size(ff_resp);
Psi = zeros([size(conds,2),dims(1:3)]);



old_noise = -1;
fprintf('Running condition: ')
for c_i = 1:size(conds,2)
    fprintf('%i ',c_i);
    n_cells = conds(1,c_i);
    p = conds(2,c_i);
    noise = conds(3,c_i);
    
    K = round(n_cells*p);
    
    % Only update when noise level changes...
    if noise ~= old_noise
        all_ff = squeeze(mean(ff_resp + sqrt(ff_resp).*randn(size(ff_resp)) .* noise,5));
        all_fn = squeeze(mean(fn_resp + sqrt(fn_resp).*randn(size(fn_resp)) .* noise,5));
        all_cf = squeeze(mean(cf_resp + sqrt(cf_resp).*randn(size(cf_resp)) .* noise,5));
        all_cn = squeeze(mean(cn_resp + sqrt(cn_resp).*randn(size(cn_resp)) .* noise,5));
        old_noise = noise;
    end


    for j = 1:length(correlation_levels);
        for k = 1:length(dxs);
            for f = 1:length(freqs);
                current_ff = squeeze(all_ff(j,k,f,:));
                current_fn = squeeze(all_fn(j,k,f,:));
                current_cf = squeeze(all_cf(j,k,f,:));
                current_cn = squeeze(all_cn(j,k,f,:));
                
                psi = zeros(1,N);
                parfor i = 1:N
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
                Psi(c_i,j,k,f) = mean(psi);
            
            end
        end
    end
end
fprintf('Done.\n');

data = load('doi_data.mat');
x = linspace(0,1,9);
x2 = linspace(0,1,11);
doi_dx(1,:) = interp1(x,data.doi_dx(1,:),x2);
doi_dx(2,:) = interp1(x,data.doi_dx(2,:),x2);
doi_hz(1,:) = interp1(x,data.doi_hz(1,:),x2);
doi_hz(2,:) = interp1(x,data.doi_hz(2,:),x2);


n_conds = size(conds,2);
allSS = zeros(1,n_conds);
for j = 1:n_conds;
    temp_Psi = squeeze(Psi(j,:,:,:));
    current_Psi(:,1,:) = (temp_Psi(:,2,:)+(1-temp_Psi(:,3,:)))/2;
    current_Psi(:,2,:) = (temp_Psi(:,1,:)+(1-temp_Psi(:,4,:)))/2;
    obs_dx(2,:) = current_Psi(:,1,2);
    obs_dx(1,:) = current_Psi(:,2,2);
    obs_hz(1,:) = current_Psi(:,1,1);
    obs_hz(2,:) = current_Psi(:,1,3);
    SS = sum((obs_dx(:)-doi_dx(:)).^2) + ...
         sum((obs_hz(:)-doi_hz(:)).^2);
     
    allSS(j) = SS;
end
[~,b] = sort(allSS);
idx=b(5); idx = 4;
PsiM = zeros(length(correlation_levels),2,length(freqs));
PsiM(:,1,:) = (squeeze(Psi(idx,:,2,:) + 1-Psi(idx,:,3,:) ))/2;
PsiM(:,2,:) = (squeeze(Psi(idx,:,1,:) + 1-Psi(idx,:,4,:) ))/2;


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