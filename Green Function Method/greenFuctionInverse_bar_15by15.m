clc;
clear all;
close all;

% Data of Amplitude on the Image Plane

global lambda k AmpImage N M dW z;
N = 15;
M = 15;

% Data of Amplitude on the Image Plane
%80 by 80 lasers

% PhaseImage = (rand(20, 20) - 0.5);
% AmpImage = AmpImage .* exp(1i * PhaseImage);
figure;
imagesc(abs(AmpImage));
colormap gray
colorbar
title('Image Plane(desired pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

dW = 20e-6; %Spacing between elements
z = 100e-6; % Distance of Image plane from the resonator plane

lambda = 1.550e-6;
k = 2 * pi / lambda;

E_Sample = zeros(N, M); % initial an array to store sample image for calculation
E_Sample(8, 6:10) = 1;
[M, N] = size(E_Sample);

figure;
imagesc(abs(E_Sample));
colormap gray
colorbar
title('Image Plane(desired pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

E_Holo = HologramHelperClass.supperposition(E_Sample);

E_Holo = HologramHelperClass.normalize(E_Holo);

figure;
imagesc(abs(E_Holo));
colormap jet
colorbar
title('Sample Plane Amplitude')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

figure;
imagesc(angle(E_Holo));
colormap jet
colorbar
title('Sample Plane Phase')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size
