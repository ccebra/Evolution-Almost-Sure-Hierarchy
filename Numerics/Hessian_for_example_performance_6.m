function [g,H] = Hessian_for_example_performance_6(z,trait_pairs,alpha,linear_amplitude,phase)

%% Inputs
% 1. z: a row vector storing the traits
% 3. trait_pairs: an array whose columns correspond to a particular frequency,
% whose rows correspond to a particular pair to include in the sum
% 4. alpha: a matrix of the same size as pairs with the coefficient for the
% associated term in the sum
% 5. linear_amplitude: size of the linear term
% 6. phase: a matrix with the same size as alpha, corresponding to a phase
% shift in the sin term

%% Outputs
% 1. f: value of performance function
% 2. g: gradient of f(x,y) at x = y = z with respect to traits of
% 3. H: a struct with components H.xx and H.xy which are the diagonal and
% off diagonal blocks of the Hessian H at x = y = z

%% get dimensions
[n_trait_pairs,n_freq] = size(alpha);
[n_traits] = length(z);

%% get trait pairs
is = floor((trait_pairs - 1)/n_traits) + 1;
js = trait_pairs - n_traits*(is - 1);

%% set x = y = z (Note: this step is not necessary, but helps us keep track of the partial derivatives)
x = z; 
y = z; 

%% initialize
g = linear_amplitude*ones([n_traits,1]);
H.xx = sparse([],[],[],n_traits,n_traits);
H.xy = sparse([],[],[],n_traits,n_traits);

%% loop over trait pairs and frequency
for frequency = 1:n_freq
    for trait_pair = 1:n_trait_pairs
        phi = phase(trait_pair,frequency);
        i = is(trait_pair,frequency);
        j = js(trait_pair,frequency);
        
        %% gradient
        dfdx = (alpha(trait_pair,frequency)*2*pi/frequency)*...
            [cos(2*pi*frequency*x(i) - phi)*cos(2*pi*frequency*y(j) - phi);...
            sin(2*pi*frequency*y(i) - phi)*sin(2*pi*frequency*x(j) - phi)];
        g(i) = g(i) + dfdx(1);
        g(j) = g(j) + dfdx(2);
        
        %% diagonal block of the Hessian
        dfdxdx = -alpha(trait_pair,frequency)*(2*pi)^2*sin(2*pi*frequency*x(i) - phi)*cos(2*pi*frequency*y(j) - phi);
        if i ~= j %(terms cancel if i = j)
            H.xx(i,i) = H.xx(i,i) +  dfdxdx;
            H.xx(j,j) = H.xx(j,j) - dfdxdx;
        end
        
        %% off-diagonal block of the Hessian
        dfdxdy = -alpha(trait_pair,frequency)*(2*pi)^2*cos(2*pi*frequency*x(i) - phi)*sin(2*pi*frequency*y(j) - phi);
        if i ~= j %(terms cancel if i = j, and zero on diagonal, want to preserve structure)
            H.xy(i,j) = H.xy(i,j) + dfdxdy;
            H.xy(j,i) = H.xy(j,i) - dfdxdy;
        end
    end
end


end