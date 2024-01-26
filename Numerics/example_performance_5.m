function f = example_performance_5(x,y,alpha,linear_amplitude,phase)

%% get dimensions
[~,~,n] = size(alpha);
[pairs,T] = size(x);

%% loop over traits
advantages = (x - y);
f = linear_amplitude*sum(advantages,2);
for frequency = 1:n
    for i = 1:T
        for j = i:T
            phi = phase(i,j,frequency);
            f = f + alpha(i,j,frequency)/frequency^2*(sin(2*pi*frequency*x(:,i) - phi).*cos(2*pi*frequency*y(:,j) - phi) - sin(2*pi*frequency*y(:,i) - phi).*cos(2*pi*frequency*x(:,j) - phi));
        end
    end
end


end