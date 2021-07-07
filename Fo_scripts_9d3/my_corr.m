function [R]=my_corr(x);
x=x(:)';
N=length(x);
zrs(N-1)=0;
x1=[x zrs];
x2=[zrs x];

for k=0:N-1
	R(k+1)=x1(1:N-k)*x2((N+k):(2*N-1))';
end
R=R/R(1);