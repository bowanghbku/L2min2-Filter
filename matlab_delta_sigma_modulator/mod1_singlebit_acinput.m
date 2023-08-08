% Title: simplified 1st-order single-bit incremental modulator
% Author Affiliation: Hamad Bin Khalifa University, 34110, Qatar
% Code function: basic matlab simulation of the output quantization error using
% different reconstruction filters for a single-bit 1st-order modulator, including CoI1,
% CoI2, sinc2, and the proposed L2min2, L2min2s filters
% Pure differential NTF is used for the modulator model
% Contact: Bo Wang, bwang@hbku.edu.qa

%% simulation parameter setup
vin             = 2/3;                           %input sin wave amplitude
inte_scaling    = 0.75;                             %integrator gain to limit its maximum output level
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

%% dither signal setup
dither_seed             = 2*rand(40000,1)-1;   %uniform dither seed

%% decimated output storage
deci_out_coi1           = zeros(1,N/OSR);      %decimated output for incremental operation
deci_out_coi2           = zeros(1,N/OSR);      %decimated output for incremental operation
deci_out_sinc2          = zeros(1,N/OSR);      %decimated output for incremental operation
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

coefficient         = ones(1,floor((OSR+1)/2));
coefficient_sinc2   = conv(coefficient,coefficient);
    
%% simulate the modulator with different input
for itotal  =  1:N/OSR               %total number of incremental conversion points
    q             = zeros(1,OSR);    %stores the quantizer output
    U1            = zeros(1,OSR);    %integrator output
    max_U1        = zeros(1,OSR);
    min_U1        = zeros(1,OSR);
    for i=1:OSR
        if (i==1)
            U1(1,i)  = inte_scaling*(vinsp(i+(itotal-1)*OSR));                     %no feedback in first cycle
        else
            U1(1,i)  = U1(1,i-1)+inte_scaling*(vinsp(i+(itotal-1)*OSR)-q(1,i-1));
        end
        %if U1(1,i) >= 0                                                            %uncomment if no dither required
        if (U1(1,i)+1/3*dither_seed(unidrnd(40000)))>=  0                    %uncomment by using +/-1/3Vref uniform dither
            q(1,i)   =   1;
        else
            q(1,i)   =  -1;
        end
    end
    deci_out_coi1(1,itotal)     = sum(q)/OSR;                                       %output reconstruction
    deci_out_coi2(1,itotal)     = sum(q.*[OSR:-1:1])/sum([OSR:-1:1]);                  
    deci_out_sinc2(1,itotal)    = sum(q(1:length(coefficient_sinc2)).*coefficient_sinc2)/sum(coefficient_sinc2);            
    deci_out_l2min2(1,itotal)   = sum(q(1:length(coefficient_l2min2)).*coefficient_l2min2)/sum(coefficient_l2min2);
    deci_out_l2min2s(1,itotal)  = sum(q(1:length(coefficient_l2min2s)).*coefficient_l2min2s)/sum(coefficient_l2min2s);
end

%% freq. domain analysis of the reconstructed output
N_incre = N/OSR;                            %number of output samples
fnyq = fs/OSR;                              %nyquist frequency

%% min4 term window function
w=0.35875-0.48829*cos(2*pi*(0:N_incre-1)/N_incre)+0.14128*cos(4*pi*(0:N_incre-1)/N_incre)-0.01168*cos(6*pi*(0:N_incre-1)/N_incre);

%% for coi1
Y = deci_out_coi1.*w;
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
sqnr_coi1 = psig - pnoise

%% for coi2
Y = deci_out_coi2.*w;
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
sqnr_coi2 = psig - pnoise

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
set(gca,'YTick',-160:40:0,'FontName','Times New Roman','FontSize',12);
xlim([0 0.5]);
ylim([-160 0]);
legend('CoI1','CoI2','sinc2','L2min2s','L2min2')