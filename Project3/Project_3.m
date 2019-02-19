clear all;
close all;
clc
%A.1
%create demanded 4n bit series for random N
N_bits = 300;
b = (sign(randn(4*N_bits,1))+1)/2;

%A.2
%create demanded function bits_to_4_PAM
A=1;%random value
X = bits_to_4_PAM(b, A);

%A.3
%design first 2N_bits at 4PAM symbols Xi and next 2N_bits at 4PAM symbols
Xi = X(1:N_bits);
Xq = X(N_bits+1:2*N_bits);

%A.4
%given values
T=1;
over=10;
Ts = T/over;
Fs = 1/Ts;

%some random values for a and A to create phi signal
A1 = 4;
a = 0.5;
N = 2048; %we use a large N

%Frequency vector
F= -Fs/2:Fs/N:Fs/2-Fs/N;

[phi,t] = srrc_pulse(T, Ts, A1 ,a);%Time will cover from 0 until time of Í symbol

%first create signals
Xi_delta = 1/Ts * upsample(Xi,over);
Xi_delta_conv = conv(Xi_delta, phi)*Ts;
Xq_delta = 1/Ts * upsample(Xq,over);
Xq_delta_conv = conv(Xq_delta, phi)*Ts;

%then define time vector
T_plot = 0:Ts:N_bits-Ts;
ti_conv = linspace(T_plot(1)+t(1), T_plot(end)+t(end),length(Xi_delta_conv));%generates length(X_delta_conv) points
tq_conv = linspace(T_plot(1)+t(1), T_plot(end)+t(end),length(Xq_delta_conv));%generates length(X_delta_conv) points

%show all those
figure;
subplot(2,1,1);
plot(ti_conv,Xi_delta_conv);
title('XI(t) waveform');
ylabel('XI');
xlabel('Time axis');
subplot(2,1,2);
plot(tq_conv,Xq_delta_conv);
title('XQ(t) waveform');
ylabel('XQ');
xlabel('Time axis');

%from previous project we know how to find the demanded PXF
Ti_total = length(ti_conv)*T;
Tq_total = length(tq_conv)*T;

PXFi = ((abs(fftshift(fft(Xi_delta_conv,N))).^2)*Ts)./Ti_total;
PXFq = ((abs(fftshift(fft(Xq_delta_conv,N))).^2)*Ts)./Tq_total;

%design all periodograms using plot
figure;
subplot(2,1,1);
plot(F, PXFi);
title('Periodogram of PFXI using plot');
xlabel('Frequency axis');
ylabel('PFXI');
subplot(2,1,2);
plot(F, PXFq);
title('Periodogram of PFXQ using plot');
xlabel('Frequency axis');
ylabel('PFXQ');

%design all periodograms using semilogy
figure;
subplot(2,1,1);
semilogy(F, PXFi);
title('Periodogram of PFXI using semilogy');
xlabel('Frequency axis');
ylabel('PFXI');
subplot(2,1,2);
semilogy(F, PXFq);
title('Periodogram of PFXQ using semilogy');
xlabel('Frequency axis');
ylabel('PFXQ');

%A.5
%data given from exercise
Fo = 2;
Ximod = 2*Xi_delta_conv.*cos(2*pi*Fo*ti_conv);
Xqmod = -2*Xq_delta_conv.*sin(2*pi*Fo*tq_conv);

%design Ximod and Xqmod
figure;
subplot(2,1,1);
plot(ti_conv,Ximod);
title('XImod waveform');
xlabel('XImod');
ylabel('Time axis');
subplot(2,1,2);
plot(tq_conv,Xqmod);
title('XQmod waveform');
xlabel('XQmod');
ylabel('Time axis');

%find periodogram of those signals
PXFimod = ((abs(fftshift(fft(Ximod,N))).^2)*Ts)./Ti_total;
PXFqmod = ((abs(fftshift(fft(Xqmod,N))).^2)*Ts)./Tq_total;

%design all periodograms using plot
figure;
subplot(2,1,1);
plot(F, PXFimod);
title('Periodogram of PFXImod using plot');
xlabel('Frequency axis');
ylabel('PFXImod');
subplot(2,1,2);
plot(F, PXFqmod);
title('Periodogram of PFXQmod using plot');
xlabel('Frequency axis');
ylabel('PFXQmod');

%design all periodograms using semilogy
figure;
subplot(2,1,1);
semilogy(F, PXFimod);
title('Periodogram of PFXImod using semilogy');
xlabel('Frequency axis');
ylabel('PFXImod');
subplot(2,1,2);
semilogy(F, PXFqmod);
title('Periodogram of PFXQmod using semilogy');
xlabel('Frequency axis');
ylabel('PFXQmod');

%A.6
XmodTotal = Ximod + Xqmod;

figure;
%ti_conv is the same as tq_conv so it does not matter which one we choose
plot(ti_conv ,XmodTotal);
title('Xmod = XImod + XQmod waveform');
ylabel('Xmod');
xlabel('Time axis');

PXFmodTotal = ((abs(fftshift(fft(XmodTotal,N))).^2)*Ts)./Ti_total;

%periodogram of XmodTotal
figure;
subplot(2,1,1);
plot(F, PXFmodTotal);
title('Periodogram of PFXmod using plot');
xlabel('Frequency axis');
ylabel('PFXmod');
subplot(2,1,2);
semilogy(F, PXFmodTotal);
title('Periodogram of PFXQmod using semilogy');
xlabel('Frequency axis');
ylabel('PFXmod');


%A.7 & A.8
%Add gaussian noise
SNR = 22;
var_w = (10*A^2)/(Ts*(10^(SNR/10)));
W_sig = sqrt(var_w)*randn(1, length(XmodTotal));
W = W_sig + XmodTotal;

%A.9
WI = W.*cos(2*pi*Fo*ti_conv);
WQ = W.*(-sin(2*pi*Fo*tq_conv));

%show both WI and WQ 
figure;
subplot(2,1,1);
plot(ti_conv ,WI);
title('WI after Gaussian noise');
ylabel('WI');
xlabel('Time axis');
subplot(2,1,2);
plot(ti_conv ,WQ);
title('WQ after Gaussian noise');
ylabel('WG');
xlabel('Time axis');

PXFWI = ((abs(fftshift(fft(WI,N))).^2)*Ts)./Ti_total;
PXFWQ = ((abs(fftshift(fft(WQ,N))).^2)*Ts)./Tq_total;

%design all periodograms using plot
figure;
subplot(2,1,1);
plot(F, PXFWI);
title('Periodogram of PFXWI using plot');
xlabel('Frequency axis');
ylabel('PFXImod');
subplot(2,1,2);
plot(F, PXFWQ);
title('Periodogram of PFXWQ using plot');
xlabel('Frequency axis');
ylabel('PFXQmod');

%design all periodograms using semilogy
figure;
subplot(2,1,1);
semilogy(F, PXFWI);
title('Periodogram of PFXWI using semilogy');
xlabel('Frequency axis');
ylabel('PFXImod');
subplot(2,1,2);
semilogy(F, PXFWQ);
title('Periodogram of PFXWQ using semilogy');
xlabel('Frequency axis');
ylabel('PFXQmod');

%A.10
sigWI = Ts*conv(WI,phi);
sigWQ = Ts*conv(WQ,phi);

ti_conv_sig = linspace(ti_conv(1)+t(1), ti_conv(end)+t(end),length(sigWI));
tq_conv_sig = linspace(tq_conv(1)+t(1), tq_conv(end)+t(end),length(sigWQ));

figure;
subplot(2,1,1);
plot(ti_conv_sig, sigWI);
title('YI waveform');
ylabel('YI(t)');
xlabel('Time axis');
subplot(2,1,2);
plot(tq_conv_sig, sigWQ);
title('YQ waveform');
ylabel('YQ(t)');
xlabel('Time axis');

%compute periodogramms
PXFYI = ((abs(fftshift(fft(sigWI,N))).^2)*Ts)./Ti_total;
PXFYQ = ((abs(fftshift(fft(sigWQ,N))).^2)*Ts)./Tq_total;

%design all periodograms using plot
figure;
subplot(2,1,1);
plot(F, PXFYI);
title('Periodogram of YI using plot');
xlabel('Frequency axis');
ylabel('YI');
subplot(2,1,2);
plot(F, PXFYQ);
title('Periodogram of YQ using plot');
xlabel('Frequency axis');
ylabel('YQ');

%design all periodograms using semilogy
figure;
subplot(2,1,1);
semilogy(F, PXFYI);
title('Periodogram of WI using semilogy');
xlabel('Frequency axis');
ylabel('WI');
subplot(2,1,2);
semilogy(F, PXFYQ);
title('Periodogram of WQ using semilogy');
xlabel('Frequency axis');
ylabel('WQ');

%A.11
%find all the negative values in time domain with the following loop
cutoffCounter=0;
for p=1:length(ti_conv_sig)
    if(ti_conv_sig(1,p)<0)
        cutoffCounter = cutoffCounter+1;
    else
        break;
    end
end

%take only the positive values
YIAfteCutoff = sigWI(cutoffCounter:(length(ti_conv_sig) - (cutoffCounter+1)));
YQAfterCutoff = sigWQ(cutoffCounter:(length(tq_conv_sig) - (cutoffCounter+1)));

%decreases sample rate by over=10
YIdecresed = downsample(YIAfteCutoff,over);
YQdecreased = downsample(YQAfterCutoff,over);

figure();
scatter(YIdecresed,YQdecreased);
title('Symbols of Y');
xlabel('t(s)');
ylabel('Y(t)')

%A.12

YIestimate = detect_4_PAM(YIdecresed,A);
YQestimate = detect_4_PAM(YQdecreased,A);

figure();
scatter(YIestimate,YQestimate);
title('Symbols of Y');
xlabel('t(s)');
ylabel('Y(t)')

%A.13
ErrorI = 0;
ErrorQ = 0;

%Xi and Xq have the same length so it does not really matter which one we
%choose
for i=1:length(Xi)
    %compute I and Q errors seperately
    %also use round function to scale the number into the nearest int and
    %the comparison can be done easier
    if(YIestimate(i) ~= round(Xi(i)))
        ErrorI = ErrorI + 1;
    end
    if(YQestimate(i) ~= round(Xq(i)))
        ErrorQ = ErrorQ + 1;
    end
end
%compute total Errors
totalError = ErrorI+ErrorQ;

disp(['Total Errors ',num2str(totalError),' with ',num2str(ErrorI),' I Errors and ',num2str(ErrorQ)  ,' Q Error']);

%A.14
est_bitI = PAM_4_to_bits(YIestimate, A);
est_bitQ = PAM_4_to_bits(YQestimate, A);
%est_bit = [est_bitI est_bitQ];

est_bit = zeros(1,4*N_bits);
est_bit(1:2*N_bits) = est_bitI;
est_bit(2*N_bits+1:4*N_bits) = est_bitQ;

%A.15
bit_error = 0;
for i=1:length(b)
    if(b(i) ~= est_bit(i))
        bit_error = bit_error + 1;
    end
end
disp(['Bit error: ',num2str(bit_error)]);

%Part B
%1
SNRdb = 0:2:16;
K = 200;
for n=1:length(SNRdb)
    totalError = 0;
    bit_error = 0;
    for j=1:K
        %the code below is the same as part A
        var_w = (10*A^2)/(Ts*(10^(SNRdb(n)/10)));
        W_sig = sqrt(var_w)*randn(1, length(XmodTotal));
        W = W_sig + XmodTotal;
    
        WI = W.*cos(2*pi*Fo*ti_conv);
        WQ = W.*(-sin(2*pi*Fo*tq_conv));
        sigWI = Ts*conv(WI,phi);
        sigWQ = Ts*conv(WQ,phi);
        
        cutoffCounter=0;
        for p=1:length(ti_conv_sig)
            if(ti_conv_sig(1,p)<0)
                cutoffCounter = cutoffCounter+1;
            else
                break;
            end
        end
        %take only the positive values
        YIAfteCutoff = sigWI(cutoffCounter:(length(ti_conv_sig) - (cutoffCounter+1)));
        YQAfterCutoff = sigWQ(cutoffCounter:(length(tq_conv_sig) - (cutoffCounter+1)));
        %decreases sample rate by over=10
        YIdecresed = downsample(YIAfteCutoff,over);
        YQdecreased = downsample(YQAfterCutoff,over);
        YIestimate = detect_4_PAM(YIdecresed,A);
        YQestimate = detect_4_PAM(YQdecreased,A);
        
        ErrorI = 0;
        ErrorQ = 0;
        %Xi and Xq have the same length so it does not really matter which one we
        %choose
        for i=1:length(Xi)
        %compute I and Q errors seperately
        %also use round function to scale the number into the nearest int and
        %the comparison can be done easier
            if(YIestimate(i) ~= round(Xi(i)))
                ErrorI = ErrorI + 1;
            end
            if(YQestimate(i) ~= round(Xq(i)))
                ErrorQ = ErrorQ + 1;
            end
        end
        %compute total Errors
        totalError = totalError+ ErrorI+ErrorQ;
          
        est_bitI = PAM_4_to_bits(YIestimate, A);
        est_bitQ = PAM_4_to_bits(YQestimate, A);
        
        est_bit = zeros(1,4*N_bits);
        est_bit(1:2*N_bits) = est_bitI;
        est_bit(2*N_bits+1:4*N_bits) = est_bitQ;
        for i=1:length(b)
            if(b(i) ~= est_bit(i))
                bit_error = bit_error + 1;
            end
        end
    end
      IQerr_exp(1,n) =totalError/(N_bits*K);
      ber_exp(1,n) = bit_error/(N_bits*K*4);
end

%sources for theoritical T_symbol and T_bit
%https://www.mathworks.com/matlabcentral/fileexchange/30074-ber-curve-for-qam-16-in-gaussian-environment?focused=5176925&tab=function
%http://www.academia.edu/5770384/Simulation_of_OFDM_and_BERvsSNR_plots_in_matlab
%Check if correct
T_symbol = 3/2.*erfc(sqrt(0.1*(10.^(SNRdb/10))))-(1/4)*3/2.*erfc(sqrt(0.1*(10.^(SNRdb/10))));
T_bit = (1/4)*3/2.*T_symbol;

%2
figure()
semilogy(SNRdb,IQerr_exp, 'DisplayName', 'Experimental SER')
title('SER Monte Carlo method');
xlabel('SNR in db')
ylabel('Symbol error rate for 16-QAM')
hold on
semilogy(SNRdb,T_symbol,'DisplayName', 'Theoritical SER')
legend('show');

%3
figure()
semilogy(SNRdb,ber_exp, 'DisplayName', 'Experimental BER');
title('BER Monte Carlo method');
xlabel('SNR in db')
ylabel('Bit error rate for 16-QAM')
hold on
semilogy(SNRdb,T_bit,'DisplayName', 'Theoritical BER');
legend('show');


