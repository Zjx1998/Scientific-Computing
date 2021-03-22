A = [5 -1 0 0 0;
     -1 4.5 0.2 0 0;
     0 0.2 1 -0.4 0;
     0 0 -0.4 3 1;
     0 0 0 1 3;];

for i = 2:5
    for j = 1:i - 1
        P = eye(5);
        if A(i, i) == A(j, j)
            P(i, i) = sqrt(2)/2;
            P(j, j) = sqrt(2)/2;
            P(i, j) = sqrt(2)/2;
            P(j, i) = -sqrt(2)/2;
        else
            c = 2 * A(i, j) * sign(A(i, i) - A(j, j));
            b = abs(A(i, i) - A(j, j));
            P(i, i) = sqrt(0.5 * (1 + b / sqrt(b^2 + c^2)));
            P(j, j) = P(i, i);
            P(i, j) = c / (2 * P(i, i) * sqrt(b^2 + c^2));
            P(j, i) = -P(i, j);
        end
        A = P * A * P';
        disp(P * P');
    end
end