% FEMPOISSON - Solve Poisson's equation with linear finite elements
% Per-Olof Persson <persson@mit.edu>, November 2006.

% Generate mesh for unit square
m=11; n=11;
[x,y]=ndgrid((0:m-1)/(m-1),(0:n-1)/(n-1));
p=[x(:),y(:)];
t=[1,2,m+2; 1,m+2,m+1];
t=kron(t,ones(m-1,1))+kron(ones(size(t)),(0:m-2)');
t=kron(t,ones(n-1,1))+kron(ones(size(t)),(0:n-2)'*m);
e=[1:m,m+1:m:m*n,2*m:m:m*n,m*n-m+1:m*n-1];

% Assemble K and F
N=size(p,1);
K=sparse(N,N);
F=zeros(N,1);
for ielem=1:size(t,1)
  el=t(ielem,:);
  
  Q=[ones(3,1),p(el,:)];
  Area=abs(det(Q))/2;
  c=inv(Q);
    
  Kh=Area*(c(2,:)'*c(2,:)+c(3,:)'*c(3,:));
  Fh=Area/3;
  
  K(el,el)=K(el,el)+Kh;
  F(el)=F(el)+Fh;
end

% Implement Dirichlet boundary conditions
K(e,:)=0; K(:,e)=0; F(e)=0;
K(e,e)=speye(length(e),length(e));

% Solve
U=K\F;

% Plot
trisurf(t,p(:,1),p(:,2),0*p(:,1),U,'edgecolor','k','facecolor','interp');
view(2),axis equal,colorbar