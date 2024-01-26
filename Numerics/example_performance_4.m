function f = example_performance_4(x,y,alpha,linear_amplitude,phase)

%% get dimensions
[~,n] = size(alpha);
[pairs,~] = size(x);

%% loop over traits
advantages = (x - y)';
p = linear_amplitude*advantages;
for frequency = 1:n
    phi = repmat(phase(:,frequency),1,pairs);
    p = p + diag(alpha(:,frequency)/frequency^2)*(sin(2*pi*frequency*x' - phi).*cos(2*pi*frequency*y' - phi) - sin(2*pi*frequency*y' - phi).*cos(2*pi*frequency*x' - phi));
end

%% get performance
f = sum(p); %automatically fair since sin is an odd function, is a sum of sines in each advantage

end