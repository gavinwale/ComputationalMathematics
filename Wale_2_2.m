%Clear the workspace
clear;
clc;

%Hard coded function
L = @(x, y) [y^2 - x^2 - 3, x^2 + y - 3];

%Jacobian matrix J(x,y) of above function
J = @(x, y) [-2*x, 2*y; 2*x, 1];

%Hard coded initial guess
x0 = 5;
y0 = 4;

err = 0.0001;

%Maximum number of iterations
N = 37;

xx = x0;
yy = y0;
XXYY = [xx, yy];

%For each loop
for i = 1:N
    [XXYY] = XXYY - L(XXYY(1,1), XXYY(1,2)) / J(XXYY(1,1), XXYY(1,2));
    %fprintf('%d %.4f %.4f\n', i, XXYY(1,1), XXYY(1,2))

    % Check the err between the previous iteration and the current for both
    % values. Ellipsis used for long logical statement...
    if (i > 1)
        if (abs(XXYY(1,1) - prev1)/abs(XXYY(1,1)) < err ...
            && abs(XXYY(1,2) - prev2)/abs(XXYY(1,2)) < err)
            fprintf('Solution is %f %f at %d iterations\n', ...
                XXYY(1,1), XXYY(1,2), i);
            return
        end
    end

    prev1 = XXYY(1,1);
    prev2 = XXYY(1,2);
end

fprintf('Convergence not achieved for %d iterations', i);

