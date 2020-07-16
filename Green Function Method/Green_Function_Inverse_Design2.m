clc;clear all;close all;
% HologramHelperClass.supperposition(conj(E_sample))

%% Input image (virtual plane)
global lambda k AmpImage N M dW z;
N = 3;

dW = 20e-6; %Spacing between elements
z = 100e-6; % Distance of Image plane from the resonator plane

dW = 100e-6; %Spacing between elements
z = 1000e-6; % Distance of Image plane from the resonator plane

dW = 200e-6; %Spacing between elements
z = 1000e-6; % Distance of Image plane from the resonator plane

dW = 50e-6; %Spacing between elements
z = 200e-6; % Distance of Image plane from the resonator plane

lambda = 1.550e-6;
k = 2 * pi / lambda;

A_virtual = zeros(N, N); % initial an array to store sample image for calculation
A_virtual((N+1)/2, (N+1)/2) = 1;
[M, N] = size(A_virtual);

figure;
imagesc(abs(A_virtual));
colormap gray
colorbar
title('Image Plane(desired pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

%% Sample image initial state (sample plane)
A_sample = ones(N, N); % initial an array to store sample image for calculation
A_sample(1, :) = 1;
A_sample(N, :) = 1;
A_sample(:, 1) = 1;
A_sample(:, N) = 1;

figure;
imagesc(abs(A_sample));
colormap gray
colorbar
title('Image Plane(desired pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size


    
    E_sample = HologramHelperClass.supperposition((A_virtual));
    
    E_holo = HologramHelperClass.supperposition(conj(E_sample));
    
    oriphase = angle(conj(E_sample));
    
%% virtual plane

figure;
imagesc(abs(E_holo));
colormap gray
colorbar
title('Image Plane(recovered pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

figure;
imagesc(angle(E_holo));
colormap gray
colorbar
title('Image Plane(recovered pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

%% sample profile

figure;
imagesc(abs(E_sample));
colormap jet
colorbar
title('Sample Plane Amplitude')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

figure;
imagesc(angle(E_sample),[-pi, pi]);
colormap jet
colorbar
title('Sample Plane Phase')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size



%% sample plane profile quantization
E_sample_nomalized = HologramHelperClass.normalize(E_sample);
E_quantized_sample = HologramHelperClass.quantizeAmp(E_sample_nomalized)...
    .*exp(1i*HologramHelperClass.quantizePhase3(angle(E_sample)));
E_quantized_holo = HologramHelperClass.supperposition(conj(E_quantized_sample));
am = abs(E_quantized_sample);
pm = angle(conj(E_quantized_sample));

%% sample profile

figure;
imagesc(abs(E_quantized_sample));
colormap gray
colorbar
title('Sample Plane Amplitude')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

figure;
imagesc(angle(E_quantized_sample), [-pi, pi]);
colormap jet
colorbar
title('Sample Plane Phase')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

%% virtual plane

figure;
imagesc(abs(E_quantized_holo));
colormap gray
colorbar
title('Image Plane(recovered pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

figure;
imagesc(angle(E_quantized_holo));
colormap gray
colorbar
title('Image Plane(recovered pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

%% AVG SNR
I_holo = abs(E_holo).*abs(E_holo);
SNRA = I_holo((N+1)/2, (N+1)/2)/ (sum(sum(I_holo)) - I_holo((N+1)/2, (N+1)/2)) * (N*N-1);

%% SNR
screen = ones(N,N);
screen((N+1)/2, (N+1)/2) = 0;
I_holo = abs(E_holo).*abs(E_holo);
SNRS = I_holo((N+1)/2, (N+1)/2)/ (max(max(I_holo.*screen)));

data = [3,20.9841962400000,10.5065537600000;5,23.2812373900000,8.10447850900000;7,52.9323449900000,18.7397361900000;9,88.3430482100000,16.6502890400000;11,154,25.0604257200000;13,190,28.4645057900000;15,299,29.1835535300000;17,409,44.2370793900000;19,461,53.9294961400000;33,1240,131;55,3190,176];
Nx = data(:, 1);
SNR1 = data(:,2);
SNR2 = data(:,3);
figure;
plot(Nx, SNR1, 'b-*');
hold on 
% plot(Nx, SNR2, 'b-*');
title('Signal to Noise Ratio')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size


