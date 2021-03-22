% In this script, we analyze different discretizations of the ODE to find
% there local/global truncate error, we solve the ODE \frac{dx}{dt} = f(t, x)

epoch = 10000;
h = 0.001;
w_Euler_for = zeros(1, epoch);
w_Euler_back = zeros(1, epoch);
w_RK = zeros(1, epoch);

% Forward Euler
t = 0;
for i = 1:epoch
    w_Euler_for(i + 1) = w_Euler_for(i) + h * f(t, w_Euler_for(i));
    t = t + h;
end

% Backward Euler
t = 0;
for i = 1:epoch
    w_new = Proximal(w_Euler_back(i), t);
    w_Euler_back(i + 1) = w_new;
end


% Runge-Kutta scheme
t = 0;
for i = 1:epoch
    k_1 = h * f(t, w_RK(i));
    k_2 = h * f(t + h/2, w_RK(i) + k_1/2);
    k_3 = h * f(t + h/2, w_RK(i) + k_2/2);
    k_4 = h * f(t + h, w_RK(i) + k_3);
    w_RK(i + 1) = w_RK(i) + (k_1 + 2 * k_2 + k_3 * 2 + k_4)/6;
    t = t + h;
end