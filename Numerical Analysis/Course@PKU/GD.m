% This script implement several algorithms, such as conjugate gradient,
% steepest descent

% initialization 

n = 50;
A = zeros(2 * n, 2 * n);
for i = 1:2 * n
    for j = 1:2 * n
        if i == j
            A(i, j) = 2 * i;
            D(i, j) = 2 * i;
        elseif abs(i - j) == 2
            A(i, j) = 0.5 * i;
            if i > j
                L(i, j) = 0.5 * i;
            end
        elseif abs(i - j) == 4
            A(i, j) = 0.25 * i;
            if i > j
                L(i, j) = 0.25 * i;
            end
        end
    end
end
b = ones(2 * n, 1) * pi;
dim = 2 * n;
epsilon = 10e-10;
max_iter = 10e5;
r_SD = zeros(dim, max_iter);
r_CG = zeros(dim, max_iter);

% steepest descent
iter_SD = 2;
x = zeros(2 * n, 1);
r_SD(:, 1) = A * x - b;
while abs(max(r_SD(:, iter_SD - 1))) > epsilon
    % step size t
    t = r_SD(:, iter_SD - 1) .* r_SD(:, iter_SD - 1) / r_SD(:, iter_SD - 1) .* (A * r_SD(iter_SD - 1));
    x = x + t * r_SD(:, iter_SD - 1);
    r_SD(:,iter_SD) = A * x - b;
    iter_SD = iter_SD + 1;
end

% conjugate gradient
iter_CG = 2;
x = zeros(2 * n, 1);
r_SD(:, 1) = - b;
d = -b;
while abs(max(r_SD(:, iter_CG - 1))) > epsilon
    % find the next conjugate direction
    beta = ((A * x - b)' * d) / (d' * A * d);
    d_new = - A * x + b + beta * d;
    d = d_new;
    
    alpha = ((A * x - b)' * d) / (d' * A * d);
    x = x + alpha * d;
    iter_CG = iter_CG + 1;
end