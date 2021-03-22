% Für eine Einführung zu die Methode Adams-Moulton, 
% sieht https://web.mit.edu/10.001/Web/Course_Notes/Differential_Equations_Notes/node6.html
% y_{n + 1} = y_n + h/2 * (f(y_n, t_n) + f(y_{n + 1}, t_{n + 1}))
w = zeros(1, 21);
h = 0.01;
epsilon = 10e-2;
w(1) = 1;
w(2) = 1 - log(1 - exp(1) * h);
w(3) = 1 - log(1 - exp(1) * 2 * h);

for i = 1:18
    temp_1 = w(i + 2);
    temp_2 = w(i + 1);
    temp_3 = w(i);
    temp = temp_1;
    
    % using fixed point iteration to solve the implicit scheme
    while true
        temp_old = temp;
        temp = temp_1 + h * ([9 19 -5 1] * exp([temp temp_1 temp_2 temp_3]')) / 24;
        if abs(temp - temp_old) < epsilon
            break;
        end
    end
    w(i + 3) = temp;
end