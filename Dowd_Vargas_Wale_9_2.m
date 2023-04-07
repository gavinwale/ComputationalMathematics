%{
    Solves a differential equation with 1 Neumann boundary
    condition and 1 Dirichlet boundary condition. The Neumann
    condition affects the shape of the solution curve and changes
    behavior at the endpoint. A linear system is solved using
    the tridiagonal function- plots visual solution to BVP.

    @authors Gavin Wale, John Dowd, Chris Vargas
%}
clear;
clc;
clf;

% Initialize functions
p = @(x) 2;
q = @(x) 1;
r = @(x) x.^2;

% Initialize script variables
n = 100;
a = 0;
b = 1;
A = 0;
G = 2;
% Step size
h = (b - a) / (n + 1);

% Initialize the tridiagonal matrix and solution vector
main_diag = zeros(n, 1);
upper_diag = zeros(n - 1, 1);
lower_diag = zeros(n - 1, 1);
rhs = zeros(n, 1);

% Fill in the tridiagonal matrix and the solution vector
for i = 1:n
    x = a + i * h;
    main_diag(i) = 2 / h^2 - q(x);
    rhs(i) = r(x);
    if i < n
        upper_diag(i) = -1 / h^2 - p(x) / (2 * h);
        lower_diag(i) = -1 / h^2 + p(x) / (2 * h);
    end
end
% Add and handle boundary conditions
rhs(1) = rhs(1) - A * (-1 / h^2 + p(a + h) / (2 * h));
rhs(n) = rhs(n) - G * (1 - p(b - h) / (2 * h));

% Solve tridiagonal matrix equation
yIn = tridiag_solver(main_diag, upper_diag, lower_diag, rhs);
% Concatenate boundary conditions
y = [A; yIn];
y = [y; yIn(end) + 2 * h * G];

x = linspace(a, b, length(y))';
% Plot solution
plot(x, y);
xlabel('x');
ylabel('y');
title('Approximate BVP Solution');
grid on;

%{
   Implements Thomas algorithm for solving tridiagonal systems of
    linear equations.

   @param - main_diag, upper_diag, lower_diag, rhs
%}
function x = tridiag_solver(main_diag, upper_diag, lower_diag, rhs)
    % Size of system
    n = length(main_diag);
    for i = 1:(n - 1)
        % Divides lower diag by main diag
        factor = lower_diag(i) / main_diag(i);
        % Updates main diag of next row
        main_diag(i + 1) = main_diag(i + 1) - factor * upper_diag(i);
        % Update rhs
        rhs(i + 1) = rhs(i + 1) - factor * rhs(i);
    end

    x = zeros(n, 1);
    % Compute last element by dividing rhs by main diag of last row
    x(n) = rhs(n) / main_diag(n);

    % Reverse iteration
    for i = (n - 1):-1:1
        % Compute solutions for x
        x(i) = (rhs(i) - upper_diag(i) * x(i + 1)) / main_diag(i);
    end
end
