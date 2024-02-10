%% sample random quadratic models and compute rho for Gaussian distributions according to
% 1. Logit equilibrium
% 2. Replicator with Gaussian mutation
% 3. OU process

% assume centroid has converged (gradient vanished)
% note: the intransitivity is invariant to the scale of the covariance


%% clear
clear

%% set experimental parameters
num_dimensions = 38;
min_dimension = 2;
max_dimension = 50; % Default: 100
dimension_range = unique(round(10.^(linspace(log(min_dimension)/log(10),log(max_dimension)/log(10),num_dimensions))));
num_dimensions = length(dimension_range);

num_trials = 1.6*10^3;

quantiles = [0.05,0.1,0.2,0.4,0.6,0.8,0.9,0.95];

%% set distributional parameters for Hessian blocks
% this matches the distribution of H_xx to H_xy as closely as possible. May
% be better to always use a Gaussian model for H_xy since it only enters
% via its expected Frobenius norm, and we can scale/conditionin on a fixed
% norm without skewing the ensuing distribution

norm_ratio_xx_to_xy = 1; % does not apply if we use the GP model, then the ratio is sqrt(1 - 1/T) automatically
                         % note: this is not epsilon_o, this is
                         % sqrt(1/epsilon_o). Use 1 for 1, 2 for 1/4, etc

singular_value_parameters.fudge = 'off';

% delta
% singular_value_parameters.mode = 'delta'; % delta

% power law decay -- THEORY DOES NOT MATCH HERE
% singular_value_parameters.mode = 'power'; % power law
% singular_value_parameters.power_b = 0.4; % absolute value of the power used for power law decay, should be less than 0.5, else, rho goes to zero

% exponential/geometric decay
% singular_value_parameters.mode = 'geometric'; % power law decay (deterministic)
% singular_value_parameters.geo_base = 1/2; % base of power law decay
% singular_value_parameters.geo_floor = 0; % lower floor for singular values

% uniform
% singular_value_parameters.mode = 'uni'; % uniform

% pareto --  surprisingly marginal match to theory... theory is ok
% singular_value_parameters.mode = 'par'; % pareto
% singular_value_parameters.par_k = 2.5;
% singular_value_parameters.par_floor = 1; % lower floor for singular values, set to one since this is a scale parameter that does not influence expected epsilon or rho
% singular_value_parameters.fudge = 'on'; % fudges the estimate for the iso case by a correction factor (chosen empirically)

% exponential
% singular_value_parameters.mode = 'exp'; % exponential
% singular_value_parameters.exp_floor = 0.8; % lower floor for singular values

% gamma
% singular_value_parameters.mode = 'gamma';
% singular_value_parameters.gamma_d = 4; % shape parameter for the gamma distribution

% semi-circle
% singular_value_parameters.mode = 'semicircle'; % semi-circle, truncated at zero

% GP distribution
singular_value_parameters.mode = 'GP'; % based on GP model for F, uses a GOE for H_xx conditioned on a maximum eigenvalue, uses a skew-symmetric Gaussian for H_xy and D from Marchenko-Pastur
singular_value_parameters.lambda_max = -10^(-1);
singular_value_parameters.lambda_scaling = 'proportional'; % 'fixed' = lambda floor is fixed with T, 'proportional' = lambda floor scaled with sqrt(T)

%% set H_xy mode
singular_value_parameters.xy_mode = 'Gaussian';
%singular_value_parameters.xy_mode = 'matched';


%% Set distributional parameters for D
% modes for diffusion matrix. Use I for comparison to analysis, Marchenko
% for general case (GP), anything else defaults to copy the model used for
% H_xx
%singular_value_parameters.D_mode = 'I'; % sets the diffusion matriox to an identity
singular_value_parameters.D_mode = 'Marchenko'; % sets the diffusion matriox to an identity



%% preallocate
% rhos.logit(trial,dimension_index) = nan([num_trials,num_dimensions]);
% rhos.replicator(trial,dimension_index) = nan([num_trials,num_dimensions]);
% rhos.OU(trial,dimension_index) = nan([num_trials,num_dimensions]);

rhos.mean.logit = nan([num_dimensions,1]);
rhos.mean.replicator = nan([num_dimensions,1]);
rhos.mean.OU = nan([num_dimensions,1]);

rhos.median.logit = nan([num_dimensions,1]);
rhos.median.replicator = nan([num_dimensions,1]);
rhos.median.OU = nan([num_dimensions,1]);

rhos.std.logit = nan([num_dimensions,1]);
rhos.std.replicator = nan([num_dimensions,1]);
rhos.std.OU = nan([num_dimensions,1]);

rhos.quantiles.logit = nan([8,num_dimensions]); % 5%, 10%, 20%, 40%, 60%, 80%, 90%, 95%
rhos.quantiles.replicator = nan([8,num_dimensions]);
rhos.quantiles.OU = nan([8,num_dimensions]);

epsilons.predict.iso = nan([num_dimensions,1]);
epsilons.predict.logit = nan([num_dimensions,1]);
epsilons.predict.rep = nan([num_dimensions,1]);

rhos.predict.iso = nan([num_dimensions,1]);
rhos.predict.logit = nan([num_dimensions,1]);
rhos.predict.rep = nan([num_dimensions,1]);

%% loop over dimension
for dimension_index = 1:num_dimensions
    num_traits = dimension_range(dimension_index);
    fprintf('\n Dimension: %d', num_traits)

    %% loop over trials
    for trial = 1:num_trials
        %% sample H_xx and H_xy
        [H_xx,H_xy,D] = sample_random_matrices(norm_ratio_xx_to_xy,singular_value_parameters,num_traits);

        %% isovariant
        ss_cov.iso = eye(num_traits,num_traits);


        %% Logit: solve for steady state covariance (note: we only need this up to scaling)
        ss_cov.logit = inv(-H_xx);

        %% Replicator: solve for steady state covariance (note: we only need this up to scaling)
        [U,S,~] = svd(-H_xx);
        s_sqrt = sqrt(diag(S));
        H_xx_sqrt = U*diag(s_sqrt)*U';
        H_xx_sqrt_inv = U*diag(1./s_sqrt)*U';
        M = H_xx_sqrt*D*H_xx_sqrt;
        [U,S,~] = svd(M);
        M_sqrt = U*sqrt(S)*U';
        ss_cov.replicator = H_xx_sqrt_inv*M_sqrt*H_xx_sqrt_inv;


        %% OU: solve for steady state covariance (note: we only need this up to scaling)
        J = H_xx + H_xy;
        J_ksum_J = kron(H_xx,eye(num_traits,num_traits)) + kron(eye(num_traits,num_traits),H_xx);
        d = reshape(D,[num_traits^2,1]);
        ss_cov_vec = -J_ksum_J\d;
        ss_cov.OU = reshape(ss_cov_vec,[num_traits,num_traits]);

        %% compute rho
        epsilon.iso = trace(-H_xy*ss_cov.iso*H_xy*ss_cov.iso)/trace(H_xx*ss_cov.iso*H_xx*ss_cov.iso);
        epsilon.logit = trace(-H_xy*ss_cov.logit*H_xy*ss_cov.logit)/trace(H_xx*ss_cov.logit*H_xx*ss_cov.logit);
        epsilon.replicator = trace(-H_xy*ss_cov.replicator*H_xy*ss_cov.replicator)/trace(H_xx*ss_cov.replicator*H_xx*ss_cov.replicator);
        epsilon.OU = trace(-H_xy*ss_cov.OU*H_xy*ss_cov.OU)/trace(H_xx*ss_cov.OU*H_xx*ss_cov.OU);

        rho_samples.iso(trial) = (1/2)*(1/(1 + epsilon.iso));
        rho_samples.logit(trial) = (1/2)*(1/(1 + epsilon.logit));
        rho_samples.replicator(trial) = (1/2)*(1/(1 + epsilon.replicator));
        rho_samples.OU(trial) = (1/2)*(1/(1 + epsilon.OU));

        %% print
        if mod(trial,10^2) == 0
            if trial == 100
                fprintf(' Trials Completed: ')
            end
            fprintf(' %d,',trial)
        end

    end

    %% compute statistics
    rhos.mean.iso(dimension_index) = mean(rho_samples.iso);
    rhos.mean.logit(dimension_index) = mean(rho_samples.logit);
    rhos.mean.replicator(dimension_index) = mean(rho_samples.replicator);
    rhos.mean.OU(dimension_index) = mean(rho_samples.OU);

    rhos.median.iso(dimension_index) = median(rho_samples.iso);
    rhos.median.logit(dimension_index) = median(rho_samples.logit);
    rhos.median.replicator(dimension_index) = median(rho_samples.replicator);
    rhos.median.OU(dimension_index) = median(rho_samples.OU);

    rhos.std.iso(dimension_index) = std(rho_samples.iso);
    rhos.std.logit(dimension_index) = std(rho_samples.logit);
    rhos.std.replicator(dimension_index) = std(rho_samples.replicator);
    rhos.std.OU(dimension_index) = std(rho_samples.OU);

    rhos.quantiles.iso(:,dimension_index) = quantile(rho_samples.iso,quantiles);
    rhos.quantiles.logit(:,dimension_index) = quantile(rho_samples.logit,quantiles);
    rhos.quantiles.replicator(:,dimension_index) = quantile(rho_samples.replicator,quantiles);
    rhos.quantiles.OU(:,dimension_index) = quantile(rho_samples.OU,quantiles);

    %% print statistics
    fprintf('\n Rho Isovariant: %f +- %f, Median: %f, Quantiles: [ ',...
        rhos.mean.iso(dimension_index), rhos.std.iso(dimension_index), rhos.median.iso(dimension_index))
    fprintf('%g, ', rhos.quantiles.iso(1:end-1, dimension_index))
    fprintf('%g ]', rhos.quantiles.iso(end, dimension_index))

    fprintf('\n Rho Logit: %f +- %f, Median: %f, Quantiles: [ ',...
        rhos.mean.logit(dimension_index), rhos.std.logit(dimension_index), rhos.median.logit(dimension_index))
    fprintf('%g, ', rhos.quantiles.logit(1:end-1, dimension_index))
    fprintf('%g ]', rhos.quantiles.logit(end, dimension_index))

    fprintf('\n Rho Replicator: %f +- %f, Median: %f, Quantiles: [ ',...
        rhos.mean.replicator(dimension_index), rhos.std.replicator(dimension_index), rhos.median.replicator(dimension_index))
    fprintf('%g, ', rhos.quantiles.replicator(1:end-1, dimension_index))
    fprintf('%g ]', rhos.quantiles.replicator(end, dimension_index))

    fprintf('\n Rho OU: %f +- %f, Median: %f, Quantiles: [ ',...
        rhos.mean.OU(dimension_index), rhos.std.OU(dimension_index), rhos.median.OU(dimension_index))
    fprintf('%g, ', rhos.quantiles.OU(1:end-1, dimension_index))
    fprintf('%g ]', rhos.quantiles.OU(end, dimension_index))

    fprintf('\n\n')

    %% compute analytical predictions
    epsilon_predictions = predict_epsilon(num_traits,norm_ratio_xx_to_xy,singular_value_parameters);
    epsilons.predict.iso(dimension_index) = epsilon_predictions.iso;
    epsilons.predict.logit(dimension_index) = epsilon_predictions.logit;
    epsilons.predict.rep(dimension_index) = epsilon_predictions.rep;

    rhos.predict.iso(dimension_index) = 0.5/(1 + epsilon_predictions.iso);
    rhos.predict.logit(dimension_index) = 0.5/(1 + epsilon_predictions.logit);
    rhos.predict.rep(dimension_index) = 0.5/(1 + epsilon_predictions.rep);

end


%% display
% figure(1)
% clf
% 
% % epsilons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,3,1)
% 
% hold on
% plot(dimension_range,(0.5./rhos.median.logit)-1,'b-','Linewidth',1)
% 
% transparency = 0.8;
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.logit(1,:),fliplr(rhos.quantiles.logit(8,:))]) - 1,'b','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.logit(2,:),fliplr(rhos.quantiles.logit(7,:))]) - 1,'b','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.logit(3,:),fliplr(rhos.quantiles.logit(6,:))]) - 1,'b','facealpha',transparency*0.1,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.logit(4,:),fliplr(rhos.quantiles.logit(5,:))]) -1,'b','facealpha',transparency*0.2,'Linestyle','none')
% 
% axis tight
% %set(gca,'xscale','log')
% set(gca,'xscale','linear')
% set(gca,'yscale','log')
% grid on
% 
% %xlabel('Num Traits','FontSize',16,'interpreter','latex')
% ylabel('$\epsilon$','FontSize',16,'interpreter','latex')
% title('Logit','Fontsize',16,'interpreter','latex')
% 
% 
% 
% subplot(2,3,2)
% 
% hold on
% plot(dimension_range,(0.5./rhos.median.OU)-1,'m-','Linewidth',1)
% 
% transparency = 0.8;
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.OU(1,:),fliplr(rhos.quantiles.OU(8,:))]) - 1,'m','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.OU(2,:),fliplr(rhos.quantiles.OU(7,:))]) - 1,'m','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.OU(3,:),fliplr(rhos.quantiles.OU(6,:))]) - 1,'m','facealpha',transparency*0.1,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.OU(4,:),fliplr(rhos.quantiles.OU(5,:))]) -1,'m','facealpha',transparency*0.2,'Linestyle','none')
% 
% axis tight
% %set(gca,'xscale','log')
% set(gca,'xscale','linear')
% set(gca,'yscale','log')
% grid on
% 
% %xlabel('Num Traits','FontSize',16,'interpreter','latex')
% %ylabel('$\epsilon$','FontSize',16,'interpreter','latex')
% title('OU','Fontsize',16,'interpreter','latex')
% 
% 
% 
% subplot(2,3,3)
% 
% hold on
% plot(dimension_range,(0.5./rhos.median.replicator)-1,'r-','Linewidth',1)
% 
% transparency = 0.8;
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.replicator(1,:),fliplr(rhos.quantiles.replicator(8,:))]) - 1,'r','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.replicator(2,:),fliplr(rhos.quantiles.replicator(7,:))]) - 1,'r','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.replicator(3,:),fliplr(rhos.quantiles.replicator(6,:))]) - 1,'r','facealpha',transparency*0.1,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.replicator(4,:),fliplr(rhos.quantiles.replicator(5,:))]) -1,'r','facealpha',transparency*0.2,'Linestyle','none')
% 
% axis tight
% %set(gca,'xscale','log')
% set(gca,'xscale','linear')
% set(gca,'yscale','log')
% grid on
% 
% %xlabel('Num Traits','FontSize',16,'interpreter','latex')
% %ylabel('$\epsilon$','FontSize',16,'interpreter','latex')
% title('replicator','Fontsize',16,'interpreter','latex')
% 
% 
% 
% % rhos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2,3,4)
% 
% hold on
% plot(dimension_range,rhos.median.logit,'b-','Linewidth',1)
% 
% transparency = 0.8;
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.logit(1,:),fliplr(rhos.quantiles.logit(8,:))],'b','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.logit(2,:),fliplr(rhos.quantiles.logit(7,:))],'b','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.logit(3,:),fliplr(rhos.quantiles.logit(6,:))],'b','facealpha',transparency*0.1,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.logit(4,:),fliplr(rhos.quantiles.logit(5,:))],'b','facealpha',transparency*0.2,'Linestyle','none')
% 
% axis tight
% %set(gca,'xscale','log')
% set(gca,'xscale','linear')
% ylim([0,0.5])
% grid on
% 
% xlabel('Num Traits','FontSize',16,'interpreter','latex')
% ylabel('Correlation $\rho = 1/(2(1 + \epsilon))$','FontSize',16,'interpreter','latex')
% title('Logit','Fontsize',16,'interpreter','latex')
% 
% subplot(2,3,5)
% 
% hold on
% plot(dimension_range,rhos.median.OU,'m-','Linewidth',1)
% 
% transparency = 0.8;
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.OU(1,:),fliplr(rhos.quantiles.OU(8,:))],'m','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.OU(2,:),fliplr(rhos.quantiles.OU(7,:))],'m','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.OU(3,:),fliplr(rhos.quantiles.OU(6,:))],'m','facealpha',transparency*0.1,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.OU(4,:),fliplr(rhos.quantiles.OU(5,:))],'m','facealpha',transparency*0.2,'Linestyle','none')
% 
% axis tight
% % set(gca,'xscale','log')
% set(gca,'xscale','linear')
% ylim([0,0.5])
% grid on
% 
% xlabel('Num Traits','FontSize',16,'interpreter','latex')
% %ylabel('Correlation $\rho$','FontSize',16,'interpreter','latex')
% title('OU','Fontsize',16,'interpreter','latex')
% 
% subplot(2,3,6)
% 
% hold on
% plot(dimension_range,rhos.median.replicator,'r-','Linewidth',1)
% 
% transparency = 0.8;
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.replicator(1,:),fliplr(rhos.quantiles.replicator(8,:))],'r','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.replicator(2,:),fliplr(rhos.quantiles.replicator(7,:))],'r','facealpha',transparency*0.05,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.replicator(3,:),fliplr(rhos.quantiles.replicator(6,:))],'r','facealpha',transparency*0.1,'Linestyle','none')
% fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.replicator(4,:),fliplr(rhos.quantiles.replicator(5,:))],'r','facealpha',transparency*0.2,'Linestyle','none')
% 
% axis tight
% %set(gca,'xscale','log')
% set(gca,'xscale','linear')
% ylim([0,0.5])
% grid on
% 
% xlabel('Num Traits','FontSize',16,'interpreter','latex')
% %ylabel('Correlation $\rho$','FontSize',16,'interpreter','latex')
% title('replicator','Fontsize',16,'interpreter','latex')



%% paper figure
figure(2)
clf

% epsilons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1)

transparency = 0.8;

hold on
plot(dimension_range,(0.5./rhos.median.logit)-1,'b-','Linewidth',1)
plot(dimension_range,(0.5./rhos.median.OU)-1,'r-','Linewidth',1)
plot(dimension_range,(0.5./rhos.median.replicator)-1,'m-','Linewidth',1)
plot(dimension_range,(0.5./rhos.median.iso)-1,'k-','Linewidth',1)

plot(dimension_range,epsilons.predict.logit,'b--','Linewidth',1)
plot(dimension_range,epsilons.predict.rep,'m--','Linewidth',1)
plot(dimension_range,epsilons.predict.iso,'k--','Linewidth',1)

fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.iso(1,:),fliplr(rhos.quantiles.iso(8,:))]) - 1,'k','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.iso(2,:),fliplr(rhos.quantiles.iso(7,:))]) - 1,'k','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.iso(3,:),fliplr(rhos.quantiles.iso(6,:))]) - 1,'k','facealpha',transparency*0.1,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.iso(4,:),fliplr(rhos.quantiles.iso(5,:))]) -1,'k','facealpha',transparency*0.2,'Linestyle','none')

fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.logit(1,:),fliplr(rhos.quantiles.logit(8,:))]) - 1,'b','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.logit(2,:),fliplr(rhos.quantiles.logit(7,:))]) - 1,'b','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.logit(3,:),fliplr(rhos.quantiles.logit(6,:))]) - 1,'b','facealpha',transparency*0.1,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.logit(4,:),fliplr(rhos.quantiles.logit(5,:))]) -1,'b','facealpha',transparency*0.2,'Linestyle','none')

fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.OU(1,:),fliplr(rhos.quantiles.OU(8,:))]) - 1,'r','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.OU(2,:),fliplr(rhos.quantiles.OU(7,:))]) - 1,'r','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.OU(3,:),fliplr(rhos.quantiles.OU(6,:))]) - 1,'r','facealpha',transparency*0.1,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.OU(4,:),fliplr(rhos.quantiles.OU(5,:))]) -1,'r','facealpha',transparency*0.2,'Linestyle','none')

fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.replicator(1,:),fliplr(rhos.quantiles.replicator(8,:))]) - 1,'m','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.replicator(2,:),fliplr(rhos.quantiles.replicator(7,:))]) - 1,'m','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.replicator(3,:),fliplr(rhos.quantiles.replicator(6,:))]) - 1,'m','facealpha',transparency*0.1,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],(0.5./[rhos.quantiles.replicator(4,:),fliplr(rhos.quantiles.replicator(5,:))]) -1,'m','facealpha',transparency*0.2,'Linestyle','none')


axis tight
set(gca,'xscale','linear')
set(gca,'yscale','log')
grid on

set(gca,'Fontsize',16)

xlabel('Num Traits','FontSize',20,'interpreter','latex')
ylabel('$\epsilon$','FontSize',20,'interpreter','latex')
%title('Logit','Fontsize',16,'interpreter','latex')

l = legend('Logit','OU','Replicator + D','Iso');
% l = legend('Logit','Replicator + D','Iso');
set(l,'location','best','FontSize',18,'interpreter','latex')



% rhos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)

hold on
plot(dimension_range,rhos.median.iso,'k-','Linewidth',1)
plot(dimension_range,rhos.median.replicator,'m-','Linewidth',1)
plot(dimension_range,rhos.median.OU,'r-','Linewidth',1)
plot(dimension_range,rhos.median.logit,'b-','Linewidth',1)


plot(dimension_range,rhos.predict.iso,'k--','Linewidth',1)
plot(dimension_range,rhos.predict.rep,'m--','Linewidth',1)
plot(dimension_range,rhos.predict.logit,'b--','Linewidth',1)



fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.iso(1,:),fliplr(rhos.quantiles.iso(8,:))],'k','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.iso(2,:),fliplr(rhos.quantiles.iso(7,:))],'k','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.iso(3,:),fliplr(rhos.quantiles.iso(6,:))],'k','facealpha',transparency*0.1,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.iso(4,:),fliplr(rhos.quantiles.iso(5,:))],'k','facealpha',transparency*0.2,'Linestyle','none')

fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.logit(1,:),fliplr(rhos.quantiles.logit(8,:))],'b','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.logit(2,:),fliplr(rhos.quantiles.logit(7,:))],'b','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.logit(3,:),fliplr(rhos.quantiles.logit(6,:))],'b','facealpha',transparency*0.1,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.logit(4,:),fliplr(rhos.quantiles.logit(5,:))],'b','facealpha',transparency*0.2,'Linestyle','none')

fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.OU(1,:),fliplr(rhos.quantiles.OU(8,:))],'r','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.OU(2,:),fliplr(rhos.quantiles.OU(7,:))],'r','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.OU(3,:),fliplr(rhos.quantiles.OU(6,:))],'r','facealpha',transparency*0.1,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.OU(4,:),fliplr(rhos.quantiles.OU(5,:))],'r','facealpha',transparency*0.2,'Linestyle','none')

fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.replicator(1,:),fliplr(rhos.quantiles.replicator(8,:))],'m','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.replicator(2,:),fliplr(rhos.quantiles.replicator(7,:))],'m','facealpha',transparency*0.05,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.replicator(3,:),fliplr(rhos.quantiles.replicator(6,:))],'m','facealpha',transparency*0.1,'Linestyle','none')
fill([dimension_range,fliplr(dimension_range)],[rhos.quantiles.replicator(4,:),fliplr(rhos.quantiles.replicator(5,:))],'m','facealpha',transparency*0.2,'Linestyle','none')

axis tight
set(gca,'xscale','linear')
%set(gca,'xscale','linear')
ylim([0,0.5])
grid on

set(gca,'Fontsize',16)

xlabel('Num Traits','FontSize',20,'interpreter','latex')
ylabel('Correlation $\rho = 1/(2(1 + \epsilon))$','FontSize',20,'interpreter','latex')
%title('Logit','Fontsize',16,'interpreter','latex')

l = legend('Iso','Replicator + D','OU','Logit');
% l = legend('Iso','Replicator + D','Logit');
set(l,'location','best','FontSize',18,'interpreter','latex')



%% save
results.parameters.num_dimensions = num_dimensions;
results.parameters.dimension_range = [min(dimension_range),max(dimension_range)];
results.parameters.num_trials = num_trials;
results.parameters.norm_ratio_xx_to_xy = norm_ratio_xx_to_xy;
results.parameters.singular_value_parameters = singular_value_parameters;
results.parameters.quantiles = quantiles;

results.rhos = rhos;


save('SS_Correlation_Quad_f_Gauss_SS_result','results')






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% sample matrices
function [A,B,D] = sample_random_matrices(ratio,singular_value_parameters,num_traits)

%% sample random H_xx
U = orth(randn(num_traits));
if strcmp(singular_value_parameters.mode,'delta')
    v = ones([num_traits,1]);
elseif strcmp(singular_value_parameters.mode,'power')
    v = (1:num_traits).^(-singular_value_parameters.power_b);
elseif strcmp(singular_value_parameters.mode,'geometric')
    v = singular_value_parameters.geo_base.^(1:num_traits) + singular_value_parameters.geo_floor;
elseif strcmp(singular_value_parameters.mode,'uni')
    v = rand([num_traits,1]);
elseif strcmp(singular_value_parameters.mode,'par')
    v = singular_value_parameters.par_floor./((1 - rand([num_traits,1])).^(1/singular_value_parameters.par_k));
elseif strcmp(singular_value_parameters.mode,'exp')
    v = exprnd(1,1,num_traits) + singular_value_parameters.exp_floor;
elseif strcmp(singular_value_parameters.mode,'gamma')
    v = gamrnd(singular_value_parameters.gamma_d,1,num_traits,1);
elseif strcmp(singular_value_parameters.mode,'semicircle')
    M = randn([num_traits,num_traits]);
    M = 0.5*(M + M');
    v = abs(eig(M));
elseif strcmp(singular_value_parameters.mode,'GP')
    if strcmp(singular_value_parameters.lambda_scaling,'fixed')
        lambda_bound = -singular_value_parameters.lambda_max/sqrt(num_traits);
    elseif strcmp(singular_value_parameters.lambda_scaling,'proportional')
        lambda_bound = -singular_value_parameters.lambda_max;
    end
    v = Bounded_Eigenvalue_Sampler_GOE(num_traits,lambda_bound,1,1);
    v = sqrt(2)*squeeze(v); % outputs of ^^ are normalized to have variance 1 off the diagonal, want variance 2
end
S = diag(v);
A = -U*S*U'; % Make matrix symmetric & negative-definite


%% sample H_xy
if strcmp(singular_value_parameters.xy_mode,'Gaussian')
    % sample
    B = randn([num_traits,num_traits]);
    B = (B - B')/sqrt(2); % currently scales to Frobenius norm T(T-1)

    % scale to the correct expected norm (up to (1 - 1/T) factor)
    expected_xx_eig_squared = compute_xx_eig_sqr_expectation(num_traits,singular_value_parameters);
    % expected xx norm is T^2*E[lambda^2]

    norm_scaling = sqrt(expected_xx_eig_squared/((num_traits-1)*ratio^2));
    B = norm_scaling*B;

else
    if ~strcmp(singular_value_parameters.mode,'GP')
        U = orth(randn(num_traits));
        if strcmp(singular_value_parameters.mode,'delta')
            v = ones([num_traits,1]);
        elseif strcmp(singular_value_parameters.mode,'power')
            v = (1:num_traits).^(-singular_value_parameters.power_b);
        elseif strcmp(singular_value_parameters.mode,'geometric')
            v = singular_value_parameters.geo_base.^(1:num_traits) + singular_value_parameters.geo_floor;
        elseif strcmp(singular_value_parameters.mode,'uni')
            v = rand([num_traits,1]);
        elseif strcmp(singular_value_parameters.mode,'par')
            v = singular_value_parameters.par_floor./((1 - rand([num_traits,1])).^(1/singular_value_parameters.par_k));
        elseif strcmp(singular_value_parameters.mode,'exp')
            v = exprnd(1,1,num_traits) + singular_value_parameters.exp_floor;
        elseif strcmp(singular_value_parameters.mode,'gamma')
            v = gamrnd(singular_value_parameters.gamma_d,1,num_traits,1);
        elseif strcmp(singular_value_parameters.mode,'semicircle')
            M = randn([num_traits,num_traits]);
            M = 0.5*(M + M');
            v = eig(M);
        end
        for i = 1:floorDiv(num_traits,2)
            v(2*i) = v(2*i-1); % Set even-indexed singular values to be equal to preceding odd-indexed s.v.
        end
        if mod(num_traits,2) == 1
            v(num_traits) = 0; % If odd, set last singular value to 0
        end
        S = diag(v);
        w = zeros(num_traits);
        for i = 1:floorDiv(num_traits,2) % Create rotation matrix w
            w(2*i-1,2*i) = 1;
            w(2*i,2*i-1) = -1;
        end
        S = w*S;
        B = U*S*U'; %Generate random skew symmetric matrix B with decaying singular values
        B = (1/ratio)*norm(A,'fro')/norm(B, 'fro')*B; % fix norm ratio NOTE: THIS CONDITIONS ON THE RATIO, rather than enforcing the ratio in expectation. May violate theory, does not properly account for dimensional effect (1 - 1/T) scaling
    else
        B = randn([num_traits,num_traits]);
        B = (B - B');
    end
end

%% sample D
if strcmp(singular_value_parameters.D_mode,'I')
    D = eye([num_traits,num_traits]);
elseif strcmp(singular_value_parameters.D_mode,'Marchenko')
    D = randn([num_traits,2*num_traits]);
    D = D*D'/(2*num_traits);
else
    if ~strcmp(singular_value_parameters.mode,'GP') % GP defaults to Marchenko
        U = orth(randn(num_traits));
        if strcmp(singular_value_parameters.mode,'delta')
            v = ones([num_traits,1]);
        elseif strcmp(singular_value_parameters.mode,'power')
            v = (1:num_traits).^(-singular_value_parameters.power_b);
        elseif strcmp(singular_value_parameters.mode,'geometric')
            v = singular_value_parameters.geo_base.^(1:num_traits) + singular_value_parameters.geo_floor;
        elseif strcmp(singular_value_parameters.mode,'uni')
            v = rand([num_traits,1]);
        elseif strcmp(singular_value_parameters.mode,'par')
            v = singular_value_parameters.par_floor./((1 - rand([num_traits,1])).^(1/singular_value_parameters.par_k));
        elseif strcmp(singular_value_parameters.mode,'exp')
            v = exprnd(1,1,num_traits) + singular_value_parameters.exp_floor;
        elseif strcmp(singular_value_parameters.mode,'gamma')
            v = gamrnd(singular_value_parameters.gamma_d,1,num_traits,1);
        elseif strcmp(singular_value_parameters.mode,'semicircle')
            M = randn([num_traits,num_traits]);
            M = 0.5*(M + M');
            v = abs(eig(M));
        end
        S = diag(v);
        D = U*S*U'; % Make matrix symmetric & positive definite
        D = D/norm(D, 'fro'); % Normalize
    else
        D = randn([num_traits,2*num_traits]);
        D = D*D'/(2*num_traits);
    end
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function epsilon = predict_epsilon(T,norm_ratio,parameters)

% find epsilon naught (isovariant)
epsilon_o = norm_ratio^(-2);

% eigenvalue stats
if strcmp(parameters.mode,'delta')
    eig_stat.logit = 1;
    eig_stat.rep = 1;
elseif strcmp(parameters.mode,'power')
    b = parameters.power_b;
    eig_stat.logit = 1/((1 - 2*b)*(1+b)^2);
    eig_stat.rep = (1-b)/((1 - 2*b)*(1 + b/2)^2);
elseif strcmp(parameters.mode,'geometric')
    r = parameters.geo_base;
    eig_stat.logit = (1/(T^2*(T-1)))*2/((1-r^2)*(1 - r^(-1))*(1 - r^(-2)))*r^(-2*T);
    eig_stat.rep =(1/(T*(T-1)))*2/((1+r)*(1 - r^(-1/2))*(1 - r^(-1)))*r^(-T) ;

    eig_stat.logit = eig_stat.logit/(1 - 1/T); % correct for convention used below
    eig_stat.rep = eig_stat.rep/(1 - 1/T); % correct for convention used below
elseif strcmp(parameters.mode,'uni')
    %%%%%
elseif strcmp(parameters.mode,'par')
    k = parameters.par_k;
    eig_stat.logit = k^3/((k-2)*(k+1)^2);
    eig_stat.rep = k^2*(k-1)/((k-2)*(k+0.5)^2);
elseif strcmp(parameters.mode,'exp')
    floor = parameters.exp_floor;
    eig_stat.logit = (floor^2 + 2*(floor + 1))*exp(2*floor)*expint(floor)^2;
    eig_stat.rep = pi*(floor^2 + 2*(floor + 1))/(floor+1)*exp(2*floor)*(1 - erf(sqrt(floor)))^2;
    % eig_stat.logit = 2*(log(floor) + 0.577)^2; % asymptotic for small floor
    % eig_stat.rep = pi; % asymptotic for small floor
elseif strcmp(parameters.mode,'gamma')
    d = parameters.gamma_d;
    eig_stat.logit = d*(d+1)/(d - 1)^2;
    eig_stat.rep = (d+1)*(gamma(d - 0.5)/gamma(d))^2;
end

%% predict
if strcmp(parameters.fudge,'on')
    if strcmp(parameter.mode,'iso')
        epsilon.iso = (1 - 1/T)*epsilon_o*(k/(k-0.5)); % WARNING: this is a fudge factor
    end
else
    epsilon.iso = (1 - 1/T)*epsilon_o;
end
epsilon.logit = (1 - 1/T)*epsilon_o*eig_stat.logit;
epsilon.rep = (1 - 1/T)*epsilon_o*eig_stat.rep;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function expected_xx_eig_squared = compute_xx_eig_sqr_expectation(T,parameters)
    % eigenvalue stats
if strcmp(parameters.mode,'delta')
    expected_xx_eig_squared = 1;
elseif strcmp(parameters.mode,'power')
    b = parameters.power_b;
    v = (1:T).^(-b);
    expected_xx_eig_squared = sum(v.^2)/T;
elseif strcmp(parameters.mode,'geometric')
    v = parameters.geo_base.^(1:T) + parameters.geo_floor;
    expected_xx_eig_squared = sum(v.^2)/T;
elseif strcmp(parameters.mode,'uni')
    %%%%%
elseif strcmp(parameters.mode,'par')
    k = parameters.par_k;
    floor = parameters.par_floor;
    expected_xx_eig_squared = (k/(k-2))*floor^2;
elseif strcmp(parameters.mode,'exp')
    floor = parameters.exp_floor;
    expected_xx_eig_squared = floor^2 + 2*floor + 2;
elseif strcmp(parameters.mode,'gamma')
    d = parameters.gamma_d;
    expected_xx_eig_squared = d*(d+1);
end
end