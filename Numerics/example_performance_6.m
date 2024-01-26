function f = example_performance_6(x,y,trait_pairs,alpha,linear_amplitude,phase)

%% Inputs
% 1. x: a matrix whose rows correspond to the traits of a particular
% competitor
% 2. y: a matrix whose rows correspond to the traits of a particular
% competitor
% 3. trait_pairs: an array whose columns correspond to a particular frequency,
% whose rows correspond to a particular pair to include in the sum
% 4. alpha: a matrix of the same size as pairs with the coefficient for the
% associated term in the sum
% 5. linear_amplitude: size of the linear term
% 6. phase: a matrix with the same size as alpha, corresponding to a phase
% shift in the sin term

%% get dimensions
[n_trait_pairs,n_freq] = size(alpha);
[~,n_traits] = size(x);

%% get endpoints
is = floor((trait_pairs - 1)/n_traits) + 1;
js = trait_pairs - n_traits*(is - 1);

%% loop over traits
advantages = (x - y);
f = linear_amplitude*sum(advantages,2);
for frequency = 1:n_freq
    for trait_pair = 1:n_trait_pairs
        phi = phase(trait_pair,frequency);
        i = is(trait_pair,frequency);
        j = js(trait_pair,frequency);
        f = f + alpha(trait_pair,frequency)/frequency^2*(sin(2*pi*frequency*x(:,i) - phi).*cos(2*pi*frequency*y(:,j) - phi)...
            - sin(2*pi*frequency*y(:,i) - phi).*cos(2*pi*frequency*x(:,j) - phi));
    end
end


end