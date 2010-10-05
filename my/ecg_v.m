N = 64;
K = 40;
res = 1/128.0;
f = 0.9/res/2;
t = 0:res:res*(N-1);
%x = sin(2*pi*f*t)';
x = w(1,1:64)';

wav = daubcqf(10);
idwt_fun = @(x) midwt(x,wav,6);
dwt_fun = @(x) mdwt(x,wav,6);
F = zeros(N);
Fi = zeros(N);
I = eye(N);
for i = 1:N, F(:,i) = dwt_fun(I(:,i)); end
for i = 1:N, Fi(:,i) = idwt_fun(I(:,i)); end

rand('state', 0);
randn('state', 0);
A = randn(K, N);
A = orth(A')';


B = A*Fi;
z = A*x;
y0 = F*A'*z;
size(y0)
size(B)
size(z)
%yp = GPSR_BB(z, B, 1.0);
%yp = l1_ls(B, z, 5.0);
yp = reconstructAmp(B, z);
xp = Fi*yp;

figure;
plot(1:N, x-real(xp));

figure; hold on;
plot(real(xp), 'r');
plot(x);
hold off;

figure; hold on;
stem(mdwt(x, daubcqf(10), 6));
stem(mdwt(xp, daubcqf(10), 6), 'r');
hold off;
