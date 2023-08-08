% Title: simplified 1st-order single-bit incremental modulator
% Author Affiliation: Hamad Bin Khalifa University, 34110, Qatar
% Code function: basic matlab simulation of the output quantization error using
% different reconstruction filters for a single-bit 1st-order modulator, including CoI1,
% CoI2, sinc2, and the proposed L2min2, L2min2s filters
% Pure differential NTF is used for the modulator model
% Contact: Bo Wang, bwang@hbku.edu.qa

%% simulation parameter setup
OSR             = 1024;         %oversampling ratio of the converter
step_size       = 19.031e-6;    %simulation steps, use small, irregular number
start_point     = -2/3;         %minimum input
end_point       = 2/3;          %maximum input
inte_scaling    = 0.75;         %integrator gain to limit its maximum output level
vinsp           = [start_point:step_size:end_point]; %input sweep vector
N               = length(vinsp);
full_scale      = vinsp(N)-vinsp(1);
q               = zeros(N,OSR); %stores the quantizer digital output bitstream
noise_flag      = 0;            %1-turn on thermal noise simulation; 0-turnoff noise
noise_seed      = normrnd(0,64e-6,[40000,1]);%example thermal noise seed

%% dither signal setup
dither_seed     = 2*rand(40000,1)-1;   %uniform dither seed

U1_max = zeros(1,N);
U1_min = zeros(1,N);
%% simulate the modulator for different input
for i = 1:N
    U1 = zeros(1,OSR);
    for j = 1:OSR
        noise_in = 0;%noise_flag*noise_seed(unidrnd(40000),1);      %uncomment if thermal noise simulation is needed
        if (j==1)
            U1(1,j) = inte_scaling*(vinsp(i)+noise_in);            %no feedback in the fist cycle
        else
            U1(1,j) = U1(1,j-1)+inte_scaling*(vinsp(i)+noise_in-q(i,j-1));
        end
        %if U1(1,j) >= 0                                            %uncomment if no dither required
        if (U1(1,j)+1/3*dither_seed(unidrnd(40000)))>=  0           %uncomment by using +/-1/3Vref uniform dither
            q(i,j) = 1;
        else
            q(i,j) = -1;
        end
    end
    U1_max(1,i) = max(U1);   %to store the maximum integrator output for stability check
    U1_min(1,i) = min(U1);   %to store the minimum integrator output for stability check
    if i == floor(N/2)       %to visualize the simulation progress
        fprintf('---50%%--- finish\n');
    end
end
fprintf('--100%%--- finish\n');

%% to store the output reconstructed by different filters
vo_coi1     = zeros(N,1);
vo_coi2     = zeros(N,1);
vo_sinc2    = zeros(N,1);
vo_l2min2   = zeros(N,1);
vo_l2min2s  = zeros(N,1);

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
    temp = q(i,:);
    %coi1 filter
    vo_coi1(i,1)    = sum(temp)/OSR;
    %coi2 filter
    vo_coi2(i,1)    = sum(temp.*(OSR:-1:1))/sum(1:1:OSR);
    %sinc2 filter
    coefficient     = ones(1,floor((OSR+1)/2));
    coefficient_sinc2 = conv(coefficient,coefficient);
    vo_sinc2(i,1)   = sum(temp(1:length(coefficient_sinc2)).*coefficient_sinc2)/sum(coefficient_sinc2);
    %l2min2 filter
    vo_l2min2(i,1)  = sum(temp(1:length(coefficient_l2min2)).*coefficient_l2min2)/sum(coefficient_l2min2);
    %l2min2s filter
    vo_l2min2s(i,1) = sum(temp(1:length(coefficient_l2min2s)).*coefficient_l2min2s)/sum(coefficient_l2min2s);   
end

%% reconstructed output MSE
qpower_coi1     = mse(vinsp' - vo_coi1);
qpower_coi2     = mse(vinsp' - vo_coi2);
qpower_sinc2    = mse(vinsp' - vo_sinc2);
qpower_l2min2   = mse(vinsp' - vo_l2min2);
qpower_l2min2s  = mse(vinsp' - vo_l2min2s);

%% plot the output curve
nlsb = full_scale/2^10;   % error normalized to 10-bit lsb here, can change based on requirement
plot(vinsp,(vinsp-vo_coi1')/nlsb,'k')
hold on
plot(vinsp,(vinsp-vo_sinc2')/nlsb,'b')
plot(vinsp,(vinsp-vo_l2min2s')/nlsb,'c')
plot(vinsp,(vinsp-vo_coi2')/nlsb,'g')
plot(vinsp,(vinsp-vo_l2min2')/nlsb,'r')
legend('CoI1','sinc2','L2min2s','CoI2','L2min2')