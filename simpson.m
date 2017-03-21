function [res] = simpson(a, b, x, f)
% Composite Simpson's Rule applied to a given array of x and f(x) on
% interval [a, b]. x values must be equally spaced as required for
% Simpson's Rule.
    i_start = find(x == a) + 1;
    i_end = find(x == b) - 1;
    h = x(3) - x(1);
    
    int = [];
    for i = i_start:2:i_end
        int = [int, (f(i-1) + 4*f(i) + f(i+1))*(h/6)];
    end
    res = sum(int);
end