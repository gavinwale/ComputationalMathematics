clear;
clc;

a = -15;
b  = 15;
f = @(x) 2*x + 1;
e = .1;
n = 0;
x = (a+b) / 2;

stop = log2((b-a)/e);

% Check signs are opposite
if f(a) * f(b) >= 0
    error('The function will not change sign in the given interval')
end

% Until incremented n is greater than our stopping point
while n < stop
    
    if (f(a)*f(x)) < 0
        b = x;
    else
        a = x;
    end

    x = (a+b) / 2;

    n = n + e; %If incremented by e, much more accurate than n++

end

fprintf('the solution is %f\n', x);