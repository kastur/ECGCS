function w = ecg_get()
%%
load ecg128.mat
[rr,rs] = rrextract(ecg128', 128);


%%
w = zeros(length(rr), 40);
for i = 1:length(rr),
  rri = rr(i);
  x = ecg128((rri-20):(rri+19));
  y = w(i);
  w(i,:) = x;
end

%%
plot(w');

%%
end