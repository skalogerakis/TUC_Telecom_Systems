clc
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Part A%%%%%%%%%%%%%%%%%%%%%%%%%%%

%A.1
T = 10^-3;
over = 10;
Ts = T/over;
A = 3;
a = 0.5;
N = 4096;%choose a large N as requested

[phi,t] = srrc_pulse(T, Ts, A, a);

PHI_f = fftshift(fft(phi,N)*Ts); %Fourier for phi function

Fs = 1/Ts; %Ts is the sampling period so Fs is the sampling frequency

%Frequency vector
F= -Fs/2:Fs/N:Fs/2-Fs/N;

spectrum_f = abs(PHI_f).^2;

figure;
semilogy(F,spectrum_f);


%A.2
%create demanded bits
N_bits = 100;
%Create n bits series
b = (sign(randn(N_bits,1))+1)/2;

X = bits_to_2PAM(b);% the demanded decoding

T_plot = 0:Ts:N_bits-Ts;%Time will cover from 0 until time of Í symbol
X_delta = 1/Ts * upsample(X,over);

X_delta_conv = conv(X_delta, phi)*Ts;
t_conv = linspace(T_plot(1)+t(1), T_plot(end)+t(end),length(X_delta_conv));%generates length(X_delta_conv) points
figure;
plot(t_conv,X_delta_conv);

Sx = (var(X_delta_conv)/Ts).*spectrum_f;%commm

%A.3.a
T_total = length(t_conv);
%T_total = (max(t_conv)-min(t_conv))*T;

FX = (abs(fftshift(fft(X_delta_conv,N))).^2)*Ts;
PXF = FX./T_total;
%F_Px = linspace(-Fs/2, Fs/2, length(PXF));
figure;
plot(F,PXF);

figure;
semilogy(F,PXF);

%A.3.b
reps = 1000;
for n = 1:reps
    b = (sign(randn(N_bits,1))+1)/2;    %generate different symbols

    X1 = bits_to_2PAM(b);% the demanded decoding

    %T_plot = 0:Ts:N_bits-Ts;%Time will cover from 0 until time of Í symbol
    X_delta1 = 1/Ts * upsample(X1,over);

    X_delta_conv1 = conv(X_delta1, phi)*Ts;
    %t_conv = linspace(T_plot(1)+t(1), T_plot(end)+t(end),length(X_delta_conv));%generates length(X_delta_conv) points
    PFX_keep(n,:) = ((abs(fftshift(fft(X_delta_conv1,N)))).^2)./T_total;      %save for different symbols
    
end

%Sx = (var(X_delta_conv)/T).*spectrum_f;

PSD_exp = sum(PFX_keep, 1)./reps;

figure;
semilogy(F, Sx);
hold on;
semilogy(F, PSD_exp);


%A.4
%N_bits = 100 and b = (sign(randn(N_bits,1))+1)/2 from A.2;

X2 = bits_to_4PAM(b);% the demanded decoding

TX2_plot = 0:Ts:N_bits/2-Ts;%Time will cover from 0 until time of Í symbol
X_delta = 1/Ts * upsample(X,over);

X_delta_conv = conv(X_delta, phi)*Ts;
t_conv = linspace(TX2_plot(1)+t(1), TX2_plot(end)+t(end),length(X_delta_conv));%generates length(X_delta_conv) points
figure;
plot(t_conv,X_delta_conv);





