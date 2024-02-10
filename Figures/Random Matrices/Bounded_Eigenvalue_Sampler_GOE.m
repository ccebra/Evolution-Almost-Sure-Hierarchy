function lambdas = Bounded_Eigenvalue_Sampler_GOE(T,lambda_bound,n_real,n_time_samples)
% samples eigenvalues for a GOE matrix conditioned on a lower bound
% (normalized so variance off diagonal is 1)

%% pick run time parameters
dt = 10^(-3);
burn_in = 3;

%% pick stride
stride_time = 0.1;
stride = round(stride_time/dt);

%% set stop time
max_time = burn_in + stride_time*n_time_samples;
steps = round(max_time/dt);

%% pick Gibbs sampling (boundary correction) parameters
Gibbs_refinement = 40;

%% preallocate
lambdas = nan([n_real,round((max_time - burn_in)/stride_time),T]);

%% loop over realizations
for j = 1:n_real
    %% initialize chain
    A = randn([T,T])/sqrt(T);
    A = (A + A')/sqrt(2);
    lambda = abs(eig(A)) + lambda_bound;
    lambda = sort(lambda,'descend');
    stride_count = 0;

    %% run
    for k = 1:steps

        %% compute distance and interaction (can we replace with a fast multipole expansion?)
        dlambdas = lambda - lambda';
        dlambdas(dlambdas == 0) = inf; % eliminates self interactions

        interaction = sum(1./dlambdas,2)/T;
        interaction_cieling = 10;
        interaction(abs(interaction) > interaction_cieling) =...
            sign(interaction(abs(interaction) > interaction_cieling))*interaction_cieling; 

        %% compute forces
        force = (-lambda + interaction);

        %% draw noise
        z = randn(T,1)/sqrt(T);

        %% update
        lambda = lambda + force*dt + z*sqrt(dt);
        lambda(lambda <= lambda_bound) = lambda_bound;
        lambda = sort(lambda,'descend');
        
        %% use Gibbs sampling to enforce boundary condition
        bound_set = find(lambda == lambda_bound);
        if ~isempty(bound_set)
            if bound_set(1) == 1
                lambda_min = 10;
            else
                lambda_min = lambda(bound_set(1) - 1);
            end

            for h = 1:length(bound_set)
                % set up quadrature nodes for cdf
                mus = linspace(lambda_bound,lambda_min,Gibbs_refinement+1)';

                % approximate cdf
                dist = abs(mus - lambda');
                dist(:,bound_set(1) - 1 + h) = 1; % drop self comparison
                us = -T/2*mus.^2 + sum(log(dist),2);
                us = us - mean(us(2:Gibbs_refinement));
                ps = exp(us);
                ps(1) = 0;
                ps(end) = 0;
                ps = ps/sum(ps);
                cdf = cumsum(ps);

                % sample vai inverse cdf (linear interp)
                z = rand;
                intercept = find(cdf < z,1,'last');
                mu = mus(intercept) + ((z - cdf(intercept))/(cdf(intercept+1) - cdf(intercept)))*(mus(intercept+1) - mus(intercept));

                % update
                lambda(bound_set(1) - 1 + h) = mu;
                lambda_min = mu;
            end
        end

        %% store
        if k > burn_in/dt
            if mod(k,stride) == 0
                %% store
                stride_count = stride_count + 1;
                lambdas(j,stride_count,:) = sort(lambda,'descend');
            end
        end
    end
   

   %% scale
   lambdas = sqrt(2*T)*lambdas;
end
