function u=poiunit(N,x,y,z)
%POIUNIT Solution to Poisson's equation -div*(grad(u)) = 1 on the unit square.
%
%  Syntax: U = POIUNIT(N,X,Y,Z) where N is number of terms in Fourier expansion,
%          and X,Y,Z are the coordinates to evaluate u. The equation is solved
%          in 1-D, 2-D, or 3-D, depending on the number of input parameters.
% The Fourier transform on unit square is written as 
% f(x) = a_0 + \sum_{k = 1}^{\infty} a_k\cos(k \pi x) + b_k\sin(k \pi x)

%   Per-Olof Persson <persson@mit.edu>, November 2006.

u   = zeros(size(x));
if nargin == 2
% 1D problem, the Fourier transform of 1 on unit square is given by 
% a_i = b_{2i} = 0, b_{2i + 1} = 1 / (2i + 1) / \pi, similar for the high
% dimension case
    for i = 1:2:N
        u   = u + 2^2 / pi^3 / i^3 * sin(i*pi*x);
    end
    
elseif nargin==3
% 2D problem
    for i = 1:2:N
        for j = 1:2:N
            u = u + 2^4 / pi^4 / (i*j) / (i^2+j^2) * sin(i*pi*x) .* sin(j*pi*y);
        end
    end
elseif nargin==4
% 3D problem
    for i = 1:2:N
        for j = 1:2:N
            for k = 1:2:N
                u = u + 2^6/pi^5/(i*j*k)/(i^2+j^2+k^2)*sin(i*pi*x).*sin(j*pi*y).*sin(k*pi*z);
            end
        end
    end
end

