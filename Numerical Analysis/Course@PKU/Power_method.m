% Initialization of the matrix
A = diag([4 3.9 2 1]);
O = rand(4,4);
A = O * A * O^(-1);
disp(eig(A))


% Initialization of the vector
x = rand(4, 1);
x = x/max(x);
epsilon = 10e-7;
epoch = 10e8;
err = zeros(epoch, 1);
result = zeros(epoch, 1);


% Iteration
for i = 1:epoch
    y = A * x;
    result(i) = max(y);
    y = y/result(i);
    error = max(abs(x - y));
    if  error > epsilon
        err(i) = error;
        x = y;
    else
        break;
    end
end

disp(i)
plot(1:i,log(err(1:i)));
xlabel('iteration');
ylabel('err');
title('\lambda_2/\lambda_1 = 3.9/4');
print(figure(1), "Power_method_1.jpg" , '-djpeg', '-r500');
%plot(1:50,log(err(1:50)));
    