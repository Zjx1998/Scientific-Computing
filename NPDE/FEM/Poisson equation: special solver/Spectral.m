% We test the convergence rate of the spectral method, which should be
% exponential on N. However, we do not observe this.
tic;
%profile viewer
N   = pow2(6);
% This part of code may be simplfied using meshgrid, but I do not know.
x   = repmat(reshape(linspace(0, 1, N + 1), [N + 1, 1, 1]), [1, N + 1, N + 1]);
y   = repmat(reshape(linspace(0, 1, N + 1), [1, N + 1, 1]), [N + 1, 1, N + 1]);
z   = repmat(reshape(linspace(0, 1, N + 1), [1, 1, N + 1]), [N + 1, N + 1, 1]);

% Generate the bench-mark solution by calculating truncation error up to 20
u   = poiunit(30, x, y, z);

start   = 8;
k       = 10;
err     = zeros(k);
for i = start:start + k - 1
    ui  = poiunit(i, x, y, z);
    err(i - start + 1) = max(max(max(abs(ui - u))));
end
plot(start:start + k - 1, err);
toc;