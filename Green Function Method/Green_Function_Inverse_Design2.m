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

lambda = 1.550e-6;
k = 2 * pi / lambda;

A_virtual = zeros(N, N); % initial an array to store sample image for calculation
A_virtual(2, 2) = 1;
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
imagesc(angle(E_sample));
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

%% sample profile

figure;
imagesc(abs(E_quantized_sample));
colormap gray
colorbar
title('Sample Plane Amplitude')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

figure;
imagesc(angle(E_quantized_sample));
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

