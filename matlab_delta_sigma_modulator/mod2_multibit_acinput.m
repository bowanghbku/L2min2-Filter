% Title: simplified 2nd-order single-bit incremental modulator
% Author Affiliation: Hamad Bin Khalifa University, 34110, Qatar
% Code function: basic matlab simulation of the output quantization error using
% different reconstruction filters for a single-bit 2nd-order modulator, including CoI1,
% CoI2, sinc2, and the proposed L2min2, L2min2s filters
% Pure differential NTF is used for the modulator model
% Contact: Bo Wang, bwang@hbku.edu.qa

%% simulation parameter setup
vin             = 2/3;                           %input sin wave amplitude
input_scaling   = 1;
bw              = 1e3;                           %data converter bandwidth;
OSR             = 1024;                          %oversampling ratio
N               = 2^14*OSR;                      %number of data points, 64*OSR results 0.5dB std in snr prediction
fs              = 2*bw*OSR;                      %sampling frequency
ts              = 1/fs;                          %sampling period
Mcycle          = 4879;                          %relative input frequency
fin             = Mcycle/N*fs;                   %coherent sampling, place fin in one of the fft bin
fbin            = Mcycle + 1;                    %input bin location, the first bin is dc
fbin_bw         = (bw/(fs/N))+1;                 %in-band bins
t               = (0:N-1)*ts;                    %sampling time vector
vinsp           = vin*sin(2*pi*fin*t);           %modulator input voltage vector
q_level         = 17;                            %number of internal quantizer levels (>2 and preferred to be 2^n+1 for easier implementation, with n being the DAC bit)
LSB             = 2/(q_level-1);                 %LSB of the midtread quantizer
noise_seed      = normrnd(0,64e-6,[40000,1]);    %example thermal noise seed

%% modulator parameters
%signal is scaled by 8 to limit the integrator outputs, the reader can change or modify the below code if other scaling or topology will be used
signal_scaling = 8;
a1 = 1/signal_scaling;
a2 = 0;
a3 = 0;
b1 = 1/signal_scaling;
b2 = 2/signal_scaling;
c1 = 0;
%check the meaning of a1,a2,a3,b1,b2,c1 in Fig. 11 of the below paper
%B. Wang, M. -K. Law, S. B. Belhaouari and A. Bermak, "Near-Optimal Decoding of Incremental Delta-Sigma ADC Output," in IEEE Transactions on Circuits and Systems I: Regular Papers, vol. 67, no. 11, pp. 3670-3680, Nov. 2020, doi: 10.1109/TCSI.2020.3010991.

%% decimated output storage
deci_out_coi2           = zeros(1,N/OSR);      %decimated output for incremental operation
deci_out_coi3           = zeros(1,N/OSR);      %decimated output for incremental operation
deci_out_sinc2          = zeros(1,N/OSR);      %decimated output for incremental operation
deci_out_sinc3          = zeros(1,N/OSR);      %decimated output for incremental operation
deci_out_l2min2         = zeros(1,N/OSR);      %decimated output for incremental operation
deci_out_l2min2s        = zeros(1,N/OSR);      %decimated output for incremental operation

%% filter weighting functions
coefficient_l2min2 = OSR*(OSR+1)/2 - ((1:1:OSR)-1).*(1:1:OSR)/2;
M = floor((OSR+1)/2);
wn = zeros(1,M);
for i=1:M
    wn(i) = M*(M+1)/2 - sum(0:1:i-1);
end
wn = [flip(wn),wn(2:end)];
coefficient_l2min2s = wn;

coefficient_tmp     = ones(1,floor((OSR+1)/2));
coefficient_sinc2   = conv(coefficient_tmp,coefficient_tmp);

coefficient_tmp     = ones(1,floor((OSR+2)/3));
coefficient_sinc3   = conv(conv(coefficient_tmp,coefficient_tmp),coefficient_tmp);

coefficient_tmp     = conv(1:1:OSR,ones(1,OSR));
coefficient_coi3    = flip(coefficient_tmp(1:OSR));
    
%% simulate the modulator with different input
for itotal  =  1:N/OSR               %total number of incremental conversion points
    q  = zeros(1,OSR);               %stores the quantizer output
    U1 = zeros(1,OSR);
    U2 = zeros(1,OSR);
    Y  = zeros(1,OSR);
    for i = 1:OSR
        noise_in = 0;%noise_seed(unidrnd(40000),1); %uncomment if thermal noise simulation is needed
        if (i==1)
            U1(1,i) = vinsp(i+(itotal-1)*OSR)*a1+ a1*noise_in;
            U2(1,i) = 0;
            Y(1,i) = signal_scaling*(vinsp(i+(itotal-1)*OSR)*a3+c1*U1(1,i)+U2(1,i));
            q(1,i) = 0;
        else
            U1(1,i) = U1(1,i-1)+a1*vinsp(i+(itotal-1)*OSR)-b1*q(1,i-1)*LSB + a1*noise_in;
            U2(1,i) = U2(1,i-1)+U1(1,i-1)+a2*vinsp(i)-b2*q(1,i-1)*LSB+a2*noise_in;
            Y(1,i) = signal_scaling*(a3*vinsp(i+(itotal-1)*OSR)+c1*U1(1,i)+U2(1,i)+a3*noise_in);
            q(1,i) = min(max(round(Y(1,i)/LSB),-(q_level-1)/2),(q_level-1)/2);  %multibit dac
        end
    end  
    deci_out_coi2(1,itotal)     = sum(q.*[OSR:-1:1])/sum([OSR:-1:1])*LSB;                %output reconstruction
    deci_out_coi3(1,itotal)     = sum(q(1:length(coefficient_coi3)).*coefficient_coi3)/sum(coefficient_coi3)*LSB;                 
    deci_out_sinc2(1,itotal)    = sum(q(1:length(coefficient_sinc2)).*coefficient_sinc2)/sum(coefficient_sinc2)*LSB;            
    deci_out_sinc3(1,itotal)    = sum(q(1:length(coefficient_sinc3)).*coefficient_sinc3)/sum(coefficient_sinc3)*LSB;     
    deci_out_l2min2(1,itotal)   = sum(q(1:length(coefficient_l2min2)).*coefficient_l2min2)/sum(coefficient_l2min2)*LSB;
    deci_out_l2min2s(1,itotal)  = sum(q(1:length(coefficient_l2min2s)).*coefficient_l2min2s)/sum(coefficient_l2min2s)*LSB;
end

%% freq. domain analysis of the reconstructed output
N_incre = N/OSR;                            %number of output samples
fnyq = fs/OSR;                              %nyquist frequency

%% min4 term window function
w=0.35875-0.48829*cos(2*pi*(0:N_incre-1)/N_incre)+0.14128*cos(4*pi*(0:N_incre-1)/N_incre)-0.01168*cos(6*pi*(0:N_incre-1)/N_incre);

%% for coi2
Y = deci_out_coi2.*w;
Y = fft(Y)/(norm(w,1)/N_incre);
Y = Y(1,1:N_incre/2+1);
Y(1) = Y(1)*1/N_incre;
Y(end) = Y(end)*1/N_incre;
Y(2:end-1) = Y(2:end-1)*2/N_incre;
Y = abs(Y);
f = fnyq*(0:(N_incre/2))/N_incre;
plot(f/2/bw,func_dbv(Y),'k') % plot the psd
hold on
%calculate sqnr
NBW = (norm(w,2)/norm(w,1))^2;
psig = sum(Y(1,fbin-3:fbin+3).^2)/2/(NBW*N_incre);
ptotal = sum(Y(1,1:fbin_bw).^2)/2/(NBW*N_incre);
pnoise = ptotal - psig;
psig = func_dbp(psig);
pnoise = func_dbp(pnoise);
sqnr_coi2 = psig - pnoise

%% for coi3
Y = deci_out_coi3.*w;
Y = fft(Y)/(norm(w,1)/N_incre);
Y = Y(1,1:N_incre/2+1);
Y(1) = Y(1)*1/N_incre;
Y(end) = Y(end)*1/N_incre;
Y(2:end-1) = Y(2:end-1)*2/N_incre;
Y = abs(Y);
f = fnyq*(0:(N_incre/2))/N_incre;
plot(f/2/bw,func_dbv(Y),'g')
%calculate sqnr
NBW = (norm(w,2)/norm(w,1))^2;
psig = sum(Y(1,fbin-3:fbin+3).^2)/2/(NBW*N_incre);
ptotal = sum(Y(1,1:fbin_bw).^2)/2/(NBW*N_incre);
pnoise = ptotal - psig;
psig = func_dbp(psig);
pnoise = func_dbp(pnoise);
sqnr_coi3 = psig - pnoise

%% for sinc2
Y = deci_out_sinc2.*w;
Y = fft(Y)/(norm(w,1)/N_incre);
Y = Y(1,1:N_incre/2+1);
Y(1) = Y(1)*1/N_incre;
Y(end) = Y(end)*1/N_incre;
Y(2:end-1) = Y(2:end-1)*2/N_incre;
Y = abs(Y);
f = fnyq*(0:(N_incre/2))/N_incre;
plot(f/2/bw,func_dbv(Y),'b')
%calculate sqnr
NBW = (norm(w,2)/norm(w,1))^2;
psig = sum(Y(1,fbin-3:fbin+3).^2)/2/(NBW*N_incre);
ptotal = sum(Y(1,1:fbin_bw).^2)/2/(NBW*N_incre);
pnoise = ptotal - psig;
psig = func_dbp(psig);
pnoise = func_dbp(pnoise);
sqnr_sinc2 = psig - pnoise

%% for sinc3
Y = deci_out_sinc3.*w;
Y = fft(Y)/(norm(w,1)/N_incre);
Y = Y(1,1:N_incre/2+1);
Y(1) = Y(1)*1/N_incre;
Y(end) = Y(end)*1/N_incre;
Y(2:end-1) = Y(2:end-1)*2/N_incre;
Y = abs(Y);
f = fnyq*(0:(N_incre/2))/N_incre;
plot(f/2/bw,func_dbv(Y),'y')
%calculate sqnr
NBW = (norm(w,2)/norm(w,1))^2;
psig = sum(Y(1,fbin-3:fbin+3).^2)/2/(NBW*N_incre);
ptotal = sum(Y(1,1:fbin_bw).^2)/2/(NBW*N_incre);
pnoise = ptotal - psig;
psig = func_dbp(psig);
pnoise = func_dbp(pnoise);
sqnr_sinc3 = psig - pnoise

%% for l2min2s
Y = deci_out_l2min2s.*w;
Y = fft(Y)/(norm(w,1)/N_incre);
Y = Y(1,1:N_incre/2+1);
Y(1) = Y(1)*1/N_incre;
Y(end) = Y(end)*1/N_incre;
Y(2:end-1) = Y(2:end-1)*2/N_incre;
Y = abs(Y);
f = fnyq*(0:(N_incre/2))/N_incre;
plot(f/2/bw,func_dbv(Y),'c')
%calculate sqnr
NBW = (norm(w,2)/norm(w,1))^2;
psig = sum(Y(1,fbin-3:fbin+3).^2)/2/(NBW*N_incre);
ptotal = sum(Y(1,1:fbin_bw).^2)/2/(NBW*N_incre);
pnoise = ptotal - psig;
psig = func_dbp(psig);
pnoise = func_dbp(pnoise);
sqnr_l2min2s = psig - pnoise

%% for l2min2
Y = deci_out_l2min2.*w;
Y = fft(Y)/(norm(w,1)/N_incre);
Y = Y(1,1:N_incre/2+1);
Y(1) = Y(1)*1/N_incre;
Y(end) = Y(end)*1/N_incre;
Y(2:end-1) = Y(2:end-1)*2/N_incre;
Y = abs(Y);
f = fnyq*(0:(N_incre/2))/N_incre;
plot(f/2/bw,func_dbv(Y),'r')
%calculate sqnr
NBW = (norm(w,2)/norm(w,1))^2;
psig = sum(Y(1,fbin-3:fbin+3).^2)/2/(NBW*N_incre);
ptotal = sum(Y(1,1:fbin_bw).^2)/2/(NBW*N_incre);
pnoise = ptotal - psig;
psig = func_dbp(psig);
pnoise = func_dbp(pnoise);
sqnr_l2min2 = psig - pnoise

ylabel('PSD (dB)','Color','k','FontName','Times New Roman','FontSize',12);
xlabel('f/f_{nyq}','Color','k','FontName','Times New Roman','FontSize',12);
set(gca,'GridLineStyle','--');
set(gca,'XTickLabelRotation', 0);
set(gca,'XTick',0:0.1:0.5,'FontName','Times New Roman','FontSize',12);
set(gca,'YTick',-240:40:0,'FontName','Times New Roman','FontSize',12);
xlim([0 0.5]);
ylim([-240 0]);
legend('CoI2','CoI3','sinc2','sinc3','L2min2s','L2min2')