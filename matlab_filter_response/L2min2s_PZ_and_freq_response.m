% Title: PZ plot and frequency response of the L2min2s filter
% Author Affiliation: Hamad Bin Khalifa University, 34110, Qatar
% Code function: plot the pz, impulse response, and frequency response for
% the L2min2s filter with a given OSR
% Contact: Bo Wang, bwang@hbku.edu.qa

format long

%%%%filter weight function
M = 13;
l2min2s=zeros(1,2*M-1);
for i=1:M
    l2min2s(1,i) = (M^2+M)/2-(M-i+1)^2/2+(M-i+1)/2;
end
for i=M+1:2*M-1
    l2min2s(1,i) = (M^2+M)/2-(i-M+1)^2/2+(i-M+1)/2;
end

%%%%impulse response
figure(1)
subplot(3,1,1), stem(1:2*M-1,l2min2s), axis([0,2*M-2,0,M*(M+1)/2]+1),
xlabel('Bit index'),ylabel('Weight'),
title('L2min2s Filter: Impulse Response (w/o normalization)');

%%%%frequency response
subplot(3,1,2:3),
[Hz,f] = freqz(l2min2s/sum(l2min2s),1,'whole',1e8); %1e8 is adjustable for better simulation accuracy.
plot(f/pi/2,20*log10(abs(Hz))), xlabel('\omega/(2\pi)'),ylabel('Gain, dB'),
title('Magnitude RespoMse');
axis([0,1,-100,10]), grid,

%%%%pole-zero plot
figure(2),
zplane(l2min2s);
title('L2min2s Filter: pole-zero plot');