% Testing our concepts of CS

%
% Our data lies in vector x, which is Nx1
% We take random projections of x, using the matrix A
% A is a KxN matrix where each row represents a random 
% linear combination of the data values in x. The resulting
% measurement vector z is Kx1
%         z = A x
% Note that each row of A could be binary with only a single
% one. This represents random time domain sampling

% If x was sparse by itself, then we could use this 
% optimization formulation to reconstruct it:
%         x* = arg min ||x||_1   s.t.  A x = z
%
% x0 is initialized using: x0 = A^-1 z

% If instead x is sparse in some other domain, we imagine
% a vector y that is the sparse one. Imagine that y is 
% related to the data x using the orthogonal transformation:
%         y = F x
% To reconstruct x, we first find the sparse vector y using:
%         y* = arg min ||y||_1   s.t.  A F^-1 y = z
% and     x* = F^-1 y*
%
% y0 is initialized using: y0 = F A^-1 z

N=1000;
res=1e-3;

% N points over 1 sec
% Sampling rate is N sps
% Frequency is f Hz
% Nyquist rate is 2f sps
f = 0.9/res/2;
t=0:res:res*(N-1);
x=sin(2*pi*f*t)';

% Average CS sampling rate is K sps
K=200;

% Hopefully K << 2f

% Fourier domain conversion
F=fft(eye(N))/sqrt(N);
%Fi=inv(F); % For the fourier basis
Fi=F';

rand('state', 0);
randn('state', 0);
% Random projection measurement matrix
A = randn(K,N);
A = orth(A')';

%{
% Uniform random sampling
q = randperm(N);
OMEGA = q(1:K);
P = eye(N);
A=P(:,OMEGA)';

% L3 patent based sampling
% Extend the sampling instants K-1 time units into the past
Pt= P(:, OMEGA);
for i=1:K
    ind=OMEGA(i);
    lf=max(1, ind-1+1);
    Pt(lf:ind, i)=1;
end
A=Pt';

% A2I measurement matrix
% +1/-1 with probability 0.5
rand('state', 0);
randn('state', 0);
A=(1-2*(rand(K,N)>0.5))/sqrt(N);

% A2I matrix with multiply and integrate
D=floor(N/K);
D=20;
At=(1-2*(rand(K,D)>0.5))/sqrt(D);
A=zeros(K,N);
for i=1:K
    lf=(i-1)*D+1;
    rg=min(N, lf+D-1);
    A(i, lf:rg)=At(i, :);
end
%}
%%


% Combined matrix
B=A*Fi;

% Random projections of time domain samples
z=A*x;

% Apply L1 magic
% initial guess = min energy
% pinv(A) = A'
%y0 = F*pinv(A)*z;
y0 = F*A'*z;
yp = l1eq_pd(y0, B, [], z);
%yp = l1_ls(B,y0,1);

% Reconstruct the signal back in time domain
xp = Fi*yp;

plot(1:N, x-real(xp));

disp(sprintf('Nyquist sampling rate = %f sps', f*2));
disp(sprintf('CS sampling rate = %f sps', K));
disp(sprintf('CS gain = %f', 2*f/K));
enorm=norm(x-xp);
disp(sprintf('Error norm = %f', enorm));

%% How exact is the signal?
%Convert to int16 and see if CRC would match
x16=int16(x*2^15);
xp16=int16(real(xp)*2^15);
disp(sprintf('Sum of absolute difference = %d', sum(abs(xp16-x16))));

%% Plot the original signal

Fs=1e3;
f=-(N/2)*Fs/N:Fs/N:(N/2-1)*Fs/N;

%cfigure(17.36, 12);
subplot(2,1,1);
plot(t, x);
xlabel('Time (s)');
title('Original signal - Sine wave 450Hz sampled at 1Ksps', 'FontSize', 16);

subplot(2,1,2);
plot(f, abs(fftshift(fft(x))));
xlabel('Frequency (Hz)');

set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-r0', 'figs/csconcept_orig.png', '-dpng');
print(gcf, '-r0', 'figs/csconcept_orig.eps', '-depsc2');

%% Plot the reconstructed signal

Fs=1e3;
f=-(N/2)*Fs/N:Fs/N:(N/2-1)*Fs/N;

%cfigure(17.36, 12);
subplot(2,1,1);
plot(t, xp);
xlabel('Time (s)');
title('Reconstructed Signal from 180 measurements', 'FontSize', 16);

subplot(2,1,2);
plot(f, abs(fftshift(fft(xp))));
xlabel('Frequency (Hz)');

set(gcf, 'PaperPositionMode', 'auto');
print(gcf, '-r0', 'figs/csconcept_reco.png', '-dpng');
print(gcf, '-r0', 'figs/csconcept_reco.eps', '-depsc2');
