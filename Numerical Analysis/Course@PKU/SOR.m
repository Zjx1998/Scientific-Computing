% In this script, we consider the convergence rate of Jacobi, Gauss-Seidel,
% SOR iteration. Furthermore, we try to find the best parameter w which
% guarantees the best convergence rate for SOR iteration. And check the
% result with the one predicted by theory.

% initialization

tic;
n = 200;
A = zeros(2 * n, 2 * n);
D = zeros(2 * n, 2 * n);
L = zeros(2 * n, 2 * n);
U = zeros(2 * n, 2 * n);
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
U = L';
w = 1.1;
epsilon = 10e-10;
iter = 0;
% we calculate the inverse matrix for one time and use it in every
% iteration
Inv_SOR = (D + w * L)^(-1);
Inv_GS = (D + L)^(-1);
Inv_J = D^(-1);
T_SOR = Inv_SOR * ((1 - w) * D - w * U);
T_GS = Inv_GS * U;
T_J = Inv_J * (L + U);

% iteration of SOR
x_old = zeros(2 * n, 1);
while true
    x_new = Inv_SOR * (((1 - w) * D - w * U) * x_old + w * b);
    if abs(max(x_new - x_old)) < epsilon
        break;
    end
    iter = iter + 1;
    x_old = x_new;
end
iter_SOR = iter;

% iteration of Gauss-Seidel 
iter = 0;
x_old = zeros(2 * n, 1);
while true
    x_new = Inv_GS * (-U * x_old + b);
    if abs(max(x_new - x_old)) < epsilon
        break;
    end
    iter = iter + 1;
    x_old = x_new;
end
iter_GS = iter;


% iteration of Jacobi
iter = 0;
x_old = zeros(2 * n, 1);
while true
    x_new = Inv_J * ((L + U) * x_old + b);
    if abs(max(x_new - x_old)) < epsilon
        break;
    end
    iter = iter + 1;
    x_old = x_new;
end
iter_J = iter;


disp([iter_SOR iter_GS iter_J]);
disp([norm(max(eig(T_SOR))) norm(max(eig(T_GS))) norm(max(eig(T_J)))]);