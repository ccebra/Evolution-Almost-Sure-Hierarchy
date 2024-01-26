function [g,H] = Hessian_for_performance_chicken_moran_24_individuals(f,x,tol)
%% Inputs
% 1. f is the function we are taking the numerical gradient and Hessian of
% 2. x is the trait (probability of defecting)
% 3. tol is the tolerance (amount of difference between numerical
% approximation of gradient and Hessian at a scale and half that scale that
% we tolerate)

%% outputs
% 1. g: gradient
% 2. H: Hessian

%% gradient
d_0 = 10^-2;
d_1 = 0.5*10^-2;
stop = 0;

while stop == 0
    g_0 = (f(x+d_0,x)-f(x-d_0,x))/(2*d_0);
    g_1 = (f(x+d_1,x)-f(x-d_1,x))/(2*d_1);
    if abs(g_0-g_1) < tol
        stop = 1;
    end
    d_0 = d_0/2;
    d_1 = d_1/2;
end
g = g_1;

%% Hessian
d_0 = 10^-2;
d_1 = 0.5*10^-2;
stop = 0;
while stop == 0
    h_0 = (f(x+d_0,x)-2*f(x,x)+f(x-d_0,x))/(d_0^2);
    h_1 = (f(x+d_1,x)-2*f(x,x)+f(x-d_1,x))/(d_1^2);
    if abs(h_0-h_1) < tol
        stop = 1;
    end
    d_0 = d_0/2;
    d_1 = d_1/2;
end
H.xx = h_1;
H.xy=0;
