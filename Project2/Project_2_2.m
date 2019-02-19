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
N = 2048;%choose a large N as requested

[phi,t] = srrc_pulse(T, Ts, A, a);
[phib,tb] = srrc_pulse(T, Ts, A, 0.8);
[phic,tc] = srrc_pulse(T, Ts, 8, a);
[phid,td] = srrc_pulse(T, Ts, 8, 0.8);

PHI_f = fftshift(fft(phi,N)*Ts); %Fourier for phi function
PHIB_f = fftshift(fft(phib,N)*Ts); %Fourier for phi function
PHIC_f = fftshift(fft(phic,N)*Ts); %Fourier for phi function
PHID_f = fftshift(fft(phid,N)*Ts); %Fourier for phi function

Fs = 1/Ts; %Ts is the sampling period so Fs is the sampling frequency

%Frequency vector
F= -Fs/2:Fs/N:Fs/2-Fs/N;

spectrum_f = abs(PHI_f).^2;
spectrumb_f = abs(PHIB_f).^2;
spectrumc_f = abs(PHIC_f).^2;
spectrumd_f = abs(PHID_f).^2;

figure;
semilogy(F,spectrum_f,'DisplayName','A=3 and a=0.5');
hold on;
semilogy(F,spectrumb_f,'DisplayName','A=3 and a=0.8');
hold off;
legend('show');
title('Energy Spectral Density for A=3');

figure;
semilogy(F,spectrumc_f,'DisplayName','A=8 and a=0.5');
hold on;
semilogy(F,spectrumd_f,'DisplayName','A=8 and a=0.8');
hold off;
legend('show');
title('Energy Spectral Density for A=8');


%A.2
%create demanded bits
N_bits = 100;
%Create n bits series
b = (sign(randn(N_bits,1))+1)/2;

X = bits_to_2PAM(b);% the demanded decoding

T_plot = 0:Ts:N_bits-Ts;%Time will cover from 0 until time of Í symbol
X_delta = 1/Ts * upsample(X,over);

X_delta_conv = conv(X_delta, phi)*Ts;
X_delta_convb = conv(X_delta, phib)*Ts;
X_delta_convc = conv(X_delta, phic)*Ts;
X_delta_convd = conv(X_delta, phid)*Ts;

t_conv = linspace(T_plot(1)+t(1), T_plot(end)+t(end),length(X_delta_conv));%generates length(X_delta_conv) points
t_convb = linspace(T_plot(1)+tb(1), T_plot(end)+tb(end),length(X_delta_convb));
t_convc = linspace(T_plot(1)+tc(1), T_plot(end)+tc(end),length(X_delta_convc));
t_convd = linspace(T_plot(1)+td(1), T_plot(end)+td(end),length(X_delta_convd));

Sx = (var(X_delta_conv)/T).*spectrum_f;
Sxb = (var(X_delta_convb)/T).*spectrumb_f;
Sxc = (var(X_delta_convc)/T).*spectrumc_f;
Sxd = (var(X_delta_convd)/T).*spectrumd_f;

figure;
subplot(2,1,1);
plot(t_conv,X_delta_conv,'DisplayName', 'a=0.5');
title('X(t) when A=3');
legend('show');
hold on;
subplot(2,1,2);
plot(t_convb,X_delta_convb,'DisplayName', 'a=0.8');
legend('show');
hold off;

figure;
subplot(2,1,1);
plot(t_convc,X_delta_convc,'DisplayName', 'a=0.5');
title('X(t) when A=8');
legend('show');
hold on;
subplot(2,1,2);
plot(t_convd,X_delta_convd,'DisplayName', 'a=0.8');
legend('show');
hold off;


%A.3.a
T_total = length(t_conv)*T;
T_totalb = length(t_convb)*T;
T_totalc = length(t_convc)*T;
T_totald = length(t_convd)*T;

PXF = ((abs(fftshift(fft(X_delta_conv,N))).^2)*Ts)./T_total;
PXFb = ((abs(fftshift(fft(X_delta_convb,N))).^2)*Ts)./T_totalb;
PXFc = ((abs(fftshift(fft(X_delta_convc,N))).^2)*Ts)./T_totalc;
PXFd = ((abs(fftshift(fft(X_delta_convd,N))).^2)*Ts)./T_totald;
%PXF = FX./T_total;

figure;
subplot(4,1,1)
plot(F,PXF);
hold on;
title('Px(F) A=3 and a=0.5');
subplot(4,1,2)
plot(F,PXFb);
hold on;
title('Px(F) A=3 and a=0.8');
subplot(4,1,3)
plot(F,PXFc);
hold on;
title('Px(F) A=8 and a=0.5');
subplot(4,1,4)
plot(F,PXFd);
hold on;
title('Px(F) A=8 and a=0.8');

figure;
subplot(2,2,1)
semilogy(F,PXF);
hold on;
title('Px(F) A=3 and a=0.5');
subplot(2,2,2)
semilogy(F,PXFb);
hold on;
title('Px(F) A=3 and a=0.8');
subplot(2,2,3)
semilogy(F,PXFc);
hold on;
title('Px(F) A=8 and a=0.5');
subplot(2,2,4)
semilogy(F,PXFd);
hold on;
title('Px(F) A=8 and a=0.8');

%i choose to use phi pulse that has the default values and not the ones i
%randomly used
%A.3.b
reps = 100;
for n = 1:reps
    b = (sign(randn(N_bits,1))+1)/2;    %generate different symbols
    X1 = bits_to_2PAM(b);% the demanded decoding
    X_delta1 = 1/Ts * upsample(X1,over);
    X_delta_conv1 = conv(X_delta1, phi)*Ts;
    PFX_keep(n,:) = ((abs(fftshift(fft(X_delta_conv1,N)))).^2)./T_total;      %save for different symbols
end

PSD_exp = sum(PFX_keep, 1)*Ts./reps;

figure;
semilogy(F, PSD_exp,'DisplayName', 'experimental values');
hold on;
semilogy(F, Sx,'DisplayName', 'theoritical values');
title('Power Spectral Density of X for 2PAM');
legend('show');

%A.4
%N_bits = 100 and b = (sign(randn(N_bits,1))+1)/2 from A.2;

X2 = bits_to_4PAM(b);% the demanded decoding
TX2_plot = 0:Ts:N_bits/2-Ts;%Time will cover from 0 until time of Í symbol
XA4_delta = 1/Ts * upsample(X2,over);

XA4_delta_conv = conv(XA4_delta, phi)*Ts;

tA4_conv = linspace(TX2_plot(1)+t(1), TX2_plot(end)+t(end),length(X_delta_conv));%generates length(X_delta_conv) points
%figure;
%plot(t_conv,X_delta_conv)
T_total_A4 = length(tA4_conv)*T;

FX4 = (abs(fftshift(fft(XA4_delta_conv,N))).^2)*Ts;
PXF4 = FX4./T_total_A4;

SxA4 = (var(XA4_delta_conv)/T).*spectrum_f;

for n = 1:reps
    b = (sign(randn(N_bits,1))+1)/2;    %generate different symbols

    XA4 = bits_to_4PAM(b);% the demanded decoding

    %T_plot = 0:Ts:N_bits-Ts;%Time will cover from 0 until time of Í symbol
    X_deltaA4 = 1/Ts * upsample(XA4,over);

    X_delta_convA4 = conv(X_deltaA4, phi)*Ts;
    %t_conv = linspace(T_plot(1)+t(1), T_plot(end)+t(end),length(X_delta_conv));%generates length(X_delta_conv) points
    PFX_keepA4(n,:) = ((abs(fftshift(fft(X_delta_convA4,N)))).^2)./T_total_A4;      %save for different symbols
    
end

PSD_expA4 = sum(PFX_keepA4, 1)*Ts./reps;

figure;
semilogy(F, PSD_expA4,'DisplayName', 'experimental values');
hold on;
semilogy(F, SxA4,'DisplayName', 'theoritical values');
title('Power Spectral Density of X for 4PAM');
legend('show');

%A.5
Tt = 2*T;   %as requested
%create new signals for new time
[phi,t] = srrc_pulse(Tt, Ts, A, a);
[phib,tb] = srrc_pulse(Tt, Ts, A, 0.8);
[phic,tc] = srrc_pulse(Tt, Ts, 8, a);
[phid,td] = srrc_pulse(Tt, Ts, 8, 0.8);

PHI_f = fftshift(fft(phi,N)*Ts); %Fourier for phi function
PHIB_f = fftshift(fft(phib,N)*Ts); %Fourier for phi function
PHIC_f = fftshift(fft(phic,N)*Ts); %Fourier for phi function
PHID_f = fftshift(fft(phid,N)*Ts); %Fourier for phi function

spectrum_f = abs(PHI_f).^2;
spectrumb_f = abs(PHIB_f).^2;
spectrumc_f = abs(PHIC_f).^2;
spectrumd_f = abs(PHID_f).^2;

N_bits = 100;
%Create n bits series
b = (sign(randn(N_bits,1))+1)/2;

X = bits_to_2PAM(b);% the demanded decoding

T_plot = 0:Ts:N_bits-Ts;%Time will cover from 0 until time of Í symbol
X_delta = 1/Ts * upsample(X,over);

X_delta_conv = conv(X_delta, phi)*Ts;
X_delta_convb = conv(X_delta, phib)*Ts;
X_delta_convc = conv(X_delta, phic)*Ts;
X_delta_convd = conv(X_delta, phid)*Ts;

Sx = (var(X_delta_conv)/T).*spectrum_f;
Sxb = (var(X_delta_convb)/T).*spectrumb_f;
Sxc = (var(X_delta_convc)/T).*spectrumc_f;
Sxd = (var(X_delta_convd)/T).*spectrumd_f;

t_conv = linspace(T_plot(1)+t(1), T_plot(end)+t(end),length(X_delta_conv));%generates length(X_delta_conv) points
t_convb = linspace(T_plot(1)+tb(1), T_plot(end)+tb(end),length(X_delta_convb));
t_convc = linspace(T_plot(1)+tc(1), T_plot(end)+tc(end),length(X_delta_convc));
t_convd = linspace(T_plot(1)+td(1), T_plot(end)+td(end),length(X_delta_convd));

T_total = length(t_conv)*Tt;
T_totalb = length(t_convb)*Tt;
T_totalc = length(t_convc)*Tt;
T_totald = length(t_convd)*Tt;

PXF = ((abs(fftshift(fft(X_delta_conv,N))).^2)*Ts)./T_total;
PXFb = ((abs(fftshift(fft(X_delta_convb,N))).^2)*Ts)./T_totalb;
PXFc = ((abs(fftshift(fft(X_delta_convc,N))).^2)*Ts)./T_totalc;
PXFd = ((abs(fftshift(fft(X_delta_convd,N))).^2)*Ts)./T_totald;

figure;
subplot(4,1,1)
plot(F,PXF);
hold on;
title('Px(F) A=3 and a=0.5');
subplot(4,1,2)
plot(F,PXFb);
hold on;
title('Px(F) A=3 and a=0.8');
subplot(4,1,3)
plot(F,PXFc);
hold on;
title('Px(F) A=8 and a=0.5');
subplot(4,1,4)
plot(F,PXFd);
hold on;
title('Px(F) A=8 and a=0.8');

figure;
subplot(2,2,1)
semilogy(F,PXF);
hold on;
title('Px(F) A=3 and a=0.5');
subplot(2,2,2)
semilogy(F,PXFb);
hold on;
title('Px(F) A=3 and a=0.8');
subplot(2,2,3)
semilogy(F,PXFc);
hold on;
title('Px(F) A=8 and a=0.5');
subplot(2,2,4)
semilogy(F,PXFd);
hold on;
title('Px(F) A=8 and a=0.8');

%i choose to use phi pulse that has the default values and not the ones i
%randomly used
%A.5.b
reps = 100;
for n = 1:reps
    b = (sign(randn(N_bits,1))+1)/2;    %generate different symbols
    X1 = bits_to_2PAM(b);% the demanded decoding
    X_delta1 = 1/Ts * upsample(X1,over);
    X_delta_conv1 = conv(X_delta1, phi)*Ts;
    PFX_keep(n,:) = ((abs(fftshift(fft(X_delta_conv1,N)))).^2)./T_total;      %save for different symbols
end

PSD_exp = sum(PFX_keep, 1)*Ts./reps;

figure;
semilogy(F, PSD_exp,'DisplayName', 'experimental values');
hold on;
semilogy(F, Sx/2,'DisplayName', 'theoritical values');
title('Power Spectral Density of X for 2PAM');
legend('show');


