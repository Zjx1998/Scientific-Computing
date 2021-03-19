function u=laplacefft(u0,bndcond)
%LAPLACEFFT Solve Laplace's Equation using Fourier's Method.
%   U = LAPLACEFFT(U0) solves Laplace's equation on a rectangle with Dirichlet
%   boundary conditions on all four boundaries ("given values"). The size of
%   the domain and the boundary conditions are given by U0. Note: Only the
%   boundary values of U0 are actually used.
%
%   U = LAPLACEFFT(U0,'DIRICHLET') is equivalent to the above.
%
%   U = LAPLACEFFT(U0,'NEUMANN') uses Neumann conditions ("insulation") on the
%   top and bottom boundaries. Note: Only U0(:,1) and U0(:,end) are used in
%   this case.
%
%   The solver might run faster if the dimensions of U0 are 2^p+1 for
%   integer p (as in the example below).
%
%   Example:
%      h=1/2^5;
%      [xx,yy]=meshgrid(0:h:2,0:h:1);
%      U0=xx.^2-yy.^2;     % Bnd conds AND exact solution
%      U=laplacefft(U0);
%      contourf(xx,yy,U),view(2),axis equal,colorbar
%      error=max(abs(U(:)-U0(:)))

%   References:
%      G. Strang, "Introduction to Applied Mathematics", Wellesley-
%      Cambridge Press, 1986. (Section 5.5)
%
%      W. Press, S. Teukolsky, W. Vetterling, B. Flannery, "Numerical
%      Recipes in C", 2nd Edition, Cambridge University Press, 1992.
%      (Section 19.4).

%   Per-Olof Persson <persson@mit.edu>, March 2004.

dirichlet=1;
if nargin>=2
  switch upper(bndcond)
   case 'DIRICHLET'
   case 'NEUMANN'
    dirichlet=0;
   otherwise
    error('The second input argument must be ''DIRICHLET or ''NEUMANN''.');
  end
end

[M,N]=size(u0);
if dirichlet % Dirichlet on all four boundaries
  f=zeros(M-2,N-2);
  f(:,1)=f(:,1)-u0(2:M-1,1);
  f(:,N-2)=f(:,N-2)-u0(2:M-1,N);
  f(1,:)=f(1,:)-u0(1,2:N-1);
  f(M-2,:)=f(M-2,:)-u0(M,2:N-1);

  K=2*(cos(pi*(1:M-2)'*ones(1,N-2)/(M-1))+cos(pi*ones(M-2,1)*(1:N-2)/(N-1))-2);
  fhat=dst(dst(f)')';
  uhat=fhat./K;
  u1=idst(idst(uhat')');
  
  u=u0;
  u(2:M-1,2:N-1)=u1;
else % Dirichlet on left and right, Neumann on top and bottom
  f=zeros(M,N-2);
  f(:,1)=f(:,1)-u0(:,1);
  f(:,N-2)=f(:,N-2)-u0(:,N);

  K=2*(cos(pi*(0:M-1)'*ones(1,N-2)/(M-1))+cos(pi*ones(M,1)*(1:N-2)/(N-1))-2);
  fhat=dst(dct(f)')';
  uhat=fhat./K;
  u1=idst(idct(uhat)')';

  u=u0;
  u(:,2:N-1)=u1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real Fourier Transforms - These operate column-wise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y=dst(x)
% Discrete Sine Transform DST-I
[M,N]=size(x);
y=[zeros(1,N);x;zeros(1,N);-flipud(x)];
yy=fft(y);
y=real(yy(2:M+1,:)/(-2*i));

function y=idst(x)
% Inverse Discrete Sine Transform IDST-I
M=size(x,1);
y=2/(M+1)*dst(x);

function y=dct(x)
% Discrete Cosine Transform DCT-I
[M,N]=size(x);
y=[x;flipud(x(2:M-1,:))];
yy=fft(y);
y=real(yy(1:M,:)/2);

function y=idct(x)
% Inverse Discrete Cosine Transform IDCT-I
M=size(x,1);
y=2/(M-1)*dct(x);