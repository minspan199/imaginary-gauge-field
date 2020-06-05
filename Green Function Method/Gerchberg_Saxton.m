clc;clear all;close all;
% HologramHelperClass.supperposition(conj(E_sample))

%% Input image (virtual plane)
global lambda k AmpImage N M dW z;
N = 15;

dW = 20e-6; %Spacing between elements
z = 100e-6; % Distance of Image plane from the resonator plane

lambda = 1.550e-6;
k = 2 * pi / lambda;

A_virtual = zeros(N, N); % initial an array to store sample image for calculation
A_virtual(8, 6:10) = 1;
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

%% Gerchberg-Saxton iteration

E_sample = A_sample;
E_virtual = A_virtual.*exp(2*pi*1i*(rand(N,N)-0.5));

iter = 100;

while iter > 0
    
    E_sample2 = HologramHelperClass.supperposition((E_virtual));
    
    E_sample = A_sample*angle(E_sample2);
    
    E_holo = HologramHelperClass.supperposition(conj(E_sample));
    
    E_virtual = A_virtual*angle(E_holo);
    
    iter = iter - 1; 
    
    disp(iter)
   
end

%% After iteration
E_virtual = HologramHelperClass.supperposition(conj(E_sample));

figure;
imagesc(abs(E_virtual));
colormap gray
colorbar
title('Image Plane(recovered pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

figure;
imagesc(angle(E_virtual));
colormap gray
colorbar
title('Image Plane(recovered pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size


%Test
% 
% E1 = HologramHelperClass.supperposition((A_virtual));
%     
% E2 = HologramHelperClass.supperposition(conj(E1));
%     
% figure;
% subplot(1,2,1)
% imagesc(abs(E2));
% colormap gray
% colorbar
% title('Image Plane(recovered pattern)')
% set(gca,'FontSize', 12) % Font Size
% subplot(1,2,2)
% imagesc(abs(A_virtual));
% colormap gray
% colorbar
% title('Image Plane(desired pattern)')
% set(gca,'FontSize', 12) % Font Size
% 
% figure;
% subplot(1,2,1)
% imagesc(abs(E1));
% colormap gray
% colorbar
% title('Image Plane(recovered pattern)')
% set(gca,'FontSize', 12) % Font Size
% subplot(1,2,2)
% imagesc(angle(E1));
% colormap gray
% colorbar
% title('Image Plane(desired pattern)')
% set(gca,'FontSize', 12) % Font Size
% 
