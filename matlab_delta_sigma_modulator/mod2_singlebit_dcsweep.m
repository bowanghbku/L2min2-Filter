% Title: simplified 2nd-order single-bit incremental modulator with CIFB topology
% Author Affiliation: Hamad Bin Khalifa University, 34110, Qatar
% Code function: basic matlab simulation of the output quantization error using
% different reconstruction filters for a single-bit 2nd-order modulator, including CoI2, CoI3, sinc2, sinc3, and the proposed L2min2, L2min2s filters
% Pure differential NTF is used for the modulator model
% Contact: Bo Wang, bwang@hbku.edu.qa

%% simulation parameter setup
clear;
format long;
OSR             = 1300;          %oversampling ratio of the converter
step_size       = 1.2715e-6;    %simulation steps, use small, irregular number
start_point     = -2/3;         %minimum input
end_point       = 2/3;          %maximum input
vinsp           = [start_point:step_size:end_point]; %input sweep vector
N               = length(vinsp);
full_scale      = vinsp(N)-vinsp(1);
q               = zeros(N,OSR); %stores the quantizer digital output bitstream
noise_flag      = 0;            %1-turn on thermal noise simulation; 0-turnoff noise
noise_seed      = normrnd(0,64e-6,[40000,1]);%example thermal noise seed

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

%% simulate the modulator with different input
U1_max = zeros(1,N);    %monitor the integrator output
U1_min = zeros(1,N);    %monitor the integrator output
U2_max = zeros(1,N);    %monitor the integrator output
U2_min = zeros(1,N);    %monitor the integrator output
Y_max  = zeros(1,N);    %monitor the quantizer input
Y_min  = zeros(1,N);    %monitor the quantizer input
for i = 1:N
    U1 = zeros(1,OSR);
    U2 = zeros(1,OSR);
    Y = zeros(1,OSR);
    for j=1:OSR
        noise_in = 0;%noise_flag*noise_seed(unidrnd(40000),1); %uncomment if thermal noise simulation is needed
        if (j==1)
            U1(1,j) = vinsp(i)*a1+ a1*noise_in;
            U2(1,j) = 0;
            Y(1,j) = signal_scaling*(vinsp(i)*a3+c1*U1(1,j)+U2(1,j));
            q(i,j) = 0;
        else
            U1(1,j) = U1(1,j-1)+a1*vinsp(i)-b1*q(i,j-1) + a1*noise_in;
            U2(1,j) = U2(1,j-1)+U1(1,j-1)+a2*vinsp(i)-b2*q(i,j-1)+a2*noise_in;
            Y(1,j) = a3*vinsp(i)+c1*U1(1,j)+U2(1,j)+a3*noise_in;
            if Y(1,j) >= 0
                q(i,j) = 1;
            else
                q(i,j) = -1;
            end
        end
    end
    U1_max(1,i) = max(U1); %maximum integrator output
    U1_min(1,i) = min(U1); %minimum integrator output
    U2_max(1,i) = max(U2);
    U2_min(1,i) = min(U2);
    Y_max(1,i)  = max(Y);
    Y_min(1,i)  = min(Y);
    if i == floor(N/2)
        fprintf('---50%%--- finish\n'); %to display the simulation progress
    end
end
fprintf('--100%%--- finish\n');

%% check the quantizer stability
if (max(Y_max)>1.5) %quantizer full scale is 1 + LSB/2 where LSB=1 for 1bit quantizer
    warning('quantizer might be out-of-range, please check Y_max for details. For implementable modulator design, please reduce the input range (decrease end_point)');
end
if (min(Y_min)<-1.5)
    warning('quantizer might be out-of-range, please check Y_min for details. For implementable modulator design, please reduce the input range (increase start_point)');
end

%% output decoding by different filters
vo_sinc2    = zeros(N,1);
vo_sinc3    = zeros(N,1);
vo_coi2     = zeros(N,1);
vo_coi3     = zeros(N,1);
vo_l2min2   = zeros(N,1);
vo_l2min2s  = zeros(N,1);

q_linear = q;
q_linear(:,1)=[];       %delete the first column with all zeros
OSR = OSR - 1;          %data length reduces by 1

%% proposed filter weighting function
coefficient_l2min2 = OSR*(OSR+1)/2 - ((1:1:OSR)-1).*(1:1:OSR)/2;
M = floor((OSR+1)/2);
wn = zeros(1,M);
for i=1:M
    wn(i) = M*(M+1)/2 - sum(0:1:i-1);
end
wn = [flip(wn),wn(2:end)];
coefficient_l2min2s = wn;

%% bitstream reconstruction
parfor i = 1:N
    temp = q_linear(i,:);
    %coi2 filter
    vo_coi2(i,1) = sum(temp.*(OSR:-1:1))/sum(1:1:OSR);
    %coi3 filter
    coefficient_tmp = conv(1:1:OSR,ones(1,OSR));
    coefficient_coi3 = flip(coefficient_tmp(1:OSR));
    vo_coi3(i,1) = sum(temp.*coefficient_coi3)/sum(coefficient_coi3(1:OSR));  
    %sinc2 filter
    coefficient = ones(1,floor((OSR+1)/2));
    coefficient_sinc2 = conv(coefficient,coefficient);
    vo_sinc2(i,1) =    sum(temp(1:length(coefficient_sinc2)).*coefficient_sinc2)/sum(coefficient_sinc2);
    %sinc3 filter
    coefficient = ones(1,floor((OSR+2)/3));
    coefficient_sinc3 = conv(conv(coefficient,coefficient),coefficient);
    vo_sinc3(i,1)   =    sum(temp(1:length(coefficient_sinc3)).*coefficient_sinc3)/sum(coefficient_sinc3);
    %l2min2 filter
    vo_l2min2(i,1) = sum(temp.*coefficient_l2min2)/sum(coefficient_l2min2(1:OSR));
    %l2min2s filter
    vo_l2min2s(i,1)= sum(temp(1:length(coefficient_l2min2s)).*coefficient_l2min2s)/sum(coefficient_l2min2s);
end

%% calculate the MSE of the conversion errors (change the input range to make the calculation below meaningful)
if (max(Y_max)>1.5 || min(Y_min)<-1.5)
    warning('the MSE calculation is not meaningful as the quantizer is out-of-range');
end
qpower_coi2     = mse(vinsp'-vo_coi2);
qpower_coi3     = mse(vinsp'-vo_coi3);
qpower_sinc2    = mse(vinsp'-vo_sinc2);
qpower_sinc3    = mse(vinsp'-vo_sinc3);
qpower_l2min2   = mse(vinsp'-vo_l2min2);
qpower_l2min2s  = mse(vinsp'-vo_l2min2s);

%% plot the output curve
nlsb = full_scale/2^16;   % error normalized to 16-bit lsb, can change based on requirement
hold on
plot(vinsp,(vinsp-vo_sinc2')/nlsb,'b')
plot(vinsp,(vinsp-vo_l2min2s')/nlsb,'c')
plot(vinsp,(vinsp-vo_l2min2')/nlsb,'r')
plot(vinsp,(vinsp-vo_coi2')/nlsb,'k')
plot(vinsp,(vinsp-vo_sinc3')/nlsb,'y')
plot(vinsp,(vinsp-vo_coi3')/nlsb,'g')
legend('sinc2','l2min2s','l2min2','CoI2','sinc3','coi3')