function w = ecg_get()
%%
load ecg128.mat
[rr,rs] = rrextract(ecg128', 128);


%%
w = zeros(length(rr), 64);
for i = 1:length(rr),
  rri = rr(i);
  x = ecg128((rri-32):(rri+31));
  y = w(i);
  w(i,:) = x;
end

%%
plot(w');

%%
end
