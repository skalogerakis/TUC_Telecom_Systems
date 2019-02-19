clear all;
close all;
clc;
%------------------------------------------------------
%Part A
%A.1
T = 10.^-2;
over = 10;
Ts = T/over;
A = 5;

%roll-off a for a = 0 0.5 1
[phiA,tA] = srrc_pulse(T, Ts, A, 0);
[phiB,tB] = srrc_pulse(T, Ts, A, 0.5);
[phiC,tC] = srrc_pulse(T, Ts, A, 1);

figure (1);
%property DisplayName can be used to show described labels
plot(tA, phiA,'DisplayName','a = 0');
hold on; %used for multiple  plots and not delete the old ones
plot(tB, phiB,'DisplayName','a = 0.5');
plot(tC, phiC,'DisplayName','a = 1');
legend('show');%visibility of legend

title('Time domain SRRC pulses of different roll-off');
xlabel('Time (sec)');
ylabel('SRRC pulses'); 
hold off; %used to clear old plots

%---------------------------------------------------------

%A.2
Fs = 1/Ts; %Ts is the sampling period so Fs is the sampling frequency
N = 4096; %length of signal(it is more effective when power of 2)

PHIA_f = fftshift(fft(phiA,N)*Ts); %Fourier for each function
PHIB_f = fftshift(fft(phiB,N)*Ts); %fftshift rearranges a Fourier transform X by shifting the zero-frequency component to the center of the array.
PHIC_f = fftshift(fft(phiC,N)*Ts);

%Frequency vector
F= -Fs/2:Fs/N:Fs/2-Fs/N;

spectrum_A_f = abs(PHIA_f).^2;
spectrum_B_f = abs(PHIB_f).^2;
spectrum_C_f = abs(PHIC_f).^2;

%design using plot
figure(2);
plot(F, spectrum_A_f,'DisplayName','a = 0');
hold on;
plot(F, spectrum_B_f,'DisplayName','a = 0.5');
plot(F, spectrum_C_f,'DisplayName','a = 1');
legend('show');

title('Frequency domain energy spectrum SRRC pulses for different roll-off');
xlabel('Frequency (Hz)');
ylabel('Energy spectrum of pulses'); 
hold off;

%desing using semilogy
%semilogy plots data with logarithmic scale
figure(3);
semilogy(F, spectrum_A_f,'DisplayName','a = 0');
hold on;
semilogy(F, spectrum_B_f,'DisplayName','a = 0.5');
semilogy(F, spectrum_C_f,'DisplayName','a = 1');
legend('show');

title('Semi logarithmic frequency domain energy spectrum SRRC pulses for different roll-off');
xlabel('Frequency (Hz)');
ylabel('Energy spectrum of pulses'); 
hold off;

%------------------------------------------------------
%A.3
C = T/(10.^3);
C1 = T/(10.^5);

figure(4);
semilogy(F, spectrum_A_f,'DisplayName','a = 0');
hold on;
semilogy(F, spectrum_B_f,'DisplayName','a = 0.5');
semilogy(F, spectrum_C_f,'DisplayName','a = 1');

semilogy(F, C*ones(length(F)),'HandleVisibility','off');
semilogy(F, C1*ones(length(F)),'HandleVisibility','off');
title('Semi logarithmic  frequency domain energy spectrum SRRC pulses for different roll-off');
xlabel('Frequency (Hz)');
ylabel('Energy spectrum of pulses'); 
hold off;
legend('show');

%-------------------------------------------------------
%Part B
%B.1
%k = 1;
for k=0:1
    % just design phi(t-kT) for k=0,1 and a=0,0.5,1
    figure
    plot(tA,phiA,'DisplayName','ö(t)');
    hold on;
    plot(tA+k*T,phiA, 'DisplayName','ö(t-kT)');% we compute tA+k*T because we want to shift the signal to the right
    title(['Roll-off a=0 and k =' num2str(k)]);
    hold off;
    legend('show');
    
    figure
    plot(tB,phiB,'DisplayName','ö(t)');
    hold on;
    plot(tB+k*T,phiB, 'DisplayName','ö(t-kT)');
    title(['Roll-off a=0.5 and k =' num2str(k)]);
    hold off;
    legend('show');
    
    figure
    plot(tC,phiC, 'DisplayName','ö(t)');
    hold on;
    plot(tC+k*T,phiC, 'DisplayName','ö(t-kT)');
    title(['Roll-off a=1 and k =' num2str(k)]);
    hold off;
    legend('show');
end

for k=0:1
    figure
    phiA_kT=[zeros(1,length(0:Ts:k*T)) phiA(1:end-length(0:Ts:k*T))];
    phiA_prod = phiA.*phiA_kT;
    plot(tA,phiA_prod, 'DisplayName','a=0');
    phiA_integ(k+1) = sum(phiA_prod)*Ts;
    
    %figure
    phiB_kT=[zeros(1,length(0:Ts:k*T)) phiB(1:end-length(0:Ts:k*T))];
    phiB_prod = phiB.*phiB_kT;
    hold on;
    plot(tB,phiB_prod, 'DisplayName','a=0.5');
    phiB_integ(k+1) = sum(phiB_prod)*Ts;
    
    %figure
    phiC_kT=[zeros(1,length(0:Ts:k*T)) phiC(1:end-length(0:Ts:k*T))];
    phiC_prod = phiC.*phiC_kT;
    plot(tC,phiC_prod, 'DisplayName','a=1');
    title(['Product of ö(t) and ö(t-kT) when k=' num2str(k)]);
    legend('show')
    hold off;
    phiC_integ(k+1) = sum(phiC_prod)*Ts;
end    

disp('Integral Á: '); disp(phiA_integ)
disp('Integral Â: '); disp(phiB_integ) 
disp('Integral C: '); disp(phiC_integ) 


%--------------------------------------------------------
%Part C
%from data given only the value of T is different from question A
TC = 1;
T_s_C = TC/over;%the T has changed so does the T_s

%C.1
%create demanded bits
N_bits = 100;
%Create n bits series
b = (sign(randn(N_bits,1))+1)/2;

%C.2.a
%call created function as requested
X = bits_to_2PAM(b);

%C.2.b
T_plot = 0:T_s_C:N_bits-T_s_C;%Time will cover from 0 until time of Í symbol
X_delta = 1/T_s_C * upsample(X,over);
figure(9);
plot(T_plot,X_delta);

%C.2.c
figure
[phiB_CC,tB_CC] = srrc_pulse(TC, T_s_C, A, 0.5);
X_delta_conv = conv(X_delta, phiB_CC)*T_s_C;
t_conv = linspace(T_plot(1)+tB_CC(1), T_plot(end)+tB_CC(end),length(X_delta_conv));%generates length(X_delta_conv) points
plot(t_conv,X_delta_conv);

%C.2.d
figure
phiB_C = phiB_CC(end:-1:1);%reverse signal
tB_CCC = -tB_CC(end):T_s_C:-tB_CC(1);
Z = conv(X_delta_conv, phiB_C)*T_s_C;
t_Z= linspace(t_conv(1)+tB_CCC(1), t_conv(end)+tB_CCC(end),length(Z));%generates length(X_delta_conv) points
plot(t_Z,Z);
hold on;
stem([0:N_bits-1]*TC,X);


