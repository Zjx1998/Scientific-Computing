% This script implement the quasi-Newton method to find the zero for a
% function

epoch = 10000;
J = 0;
J_ = 0;
epsilon = 0.01;
x = zeros(dim, epoch);

% quasi-Newton
iter = 3;
while abs(max(F(x(:, iter - 1)))) > epsilon
    incre = (F(x(:, iter - 1)) - F(x(:, iter - 2)) ...
            - J_* (x(:, iter - 1) - x(:, iter - 2))) * (x(:, iter - 1) - x(:, iter - 2))'  ...
             / norm(x(:, iter - 1) - x(:, iter - 2));
    J_new = J + incre;
    J_new_ = J_ - J_ * incre * J_ / (1 + (x(:, iter - 1) - x(:, iter - 2))' ...
    * J_ * (F(x(:, iter - 1)) - F(x(:, iter - 2)) - J_* (x(:, iter - 1) - x(:, iter - 2))) ...
     / norm(x(:, iter - 1) - x(:, iter - 2)));
    disp([max(eig(J_new * J_new_)) min(eig(J_new * J_new_))]);
    J = J_new;
    J_ = J_new_;
    x(:, iter) = x(:, iter - 1) - J_ * F(x(:, iter - 1));
end

% Newton method