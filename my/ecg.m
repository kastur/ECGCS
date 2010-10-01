N=1000;
res=1e-3;
f = 0.9/res/2;
t=0:res:res*(N-1);
x=sin(2*pi*f*t)';
K=200;
F=fft(eye(N))/sqrt(N);
Fi=F';
rand('state', 0);
randn('state', 0);
A = randn(K,N);
A = orth(A')';
B=A*Fi;
z=A*x;
y0 = F*A'*z;
size(y0)
size(B)
size(z)
yp = l1eq_pd(y0, B, [], z);
xp = Fi*yp;
plot(1:N, x-real(xp));
