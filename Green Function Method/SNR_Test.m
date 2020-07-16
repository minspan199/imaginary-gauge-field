clc;clear all;close all;
% HologramHelperClass.supperposition(conj(E_sample))

%% Input image (virtual plane)
global lambda k AmpImage N M dW z;
N = 7;

dW = 20e-6; %Spacing between elements
z = 100e-6; % Distance of Image plane from the resonator plane

dW = 100e-6; %Spacing between elements
z = 1000e-6; % Distance of Image plane from the resonator plane

dW = 200e-6; %Spacing between elements
z = 1000e-6; % Distance of Image plane from the resonator plane





lambda = 1.550e-6;
k = 2 * pi / lambda;

A_virtual = zeros(N, N); % initial an array to store sample image for calculation
A_virtual((N+1)/2, (N+1)/2) = 1;
[M, N] = size(A_virtual);


for ii = 1:100
    for jj = 1: 100
        
        SNR1(ii,jj)=calculate_SNR(20e-6 + ii*5e-6, 10e-6 + jj*10e-6);
        
    end
end

function [SNR1, SNR2] = calculate_SNR(dW, z)
% dW = 100e-6; %Spacing between elements
% z = 200e-6; % Distance of Image plane from the resonator plane
    E_sample = HologramHelperClass.supperposition((A_virtual));
    
    E_holo = HologramHelperClass.supperposition(conj(E_sample));

%% AVG SNR
I_holo = abs(E_holo).*abs(E_holo);
SNR1 = I_holo((N+1)/2, (N+1)/2)/ (sum(sum(I_holo)) - I_holo((N+1)/2, (N+1)/2)) * (N*N-1);

%% SNR
screen = ones(N,N);
screen((N+1)/2, (N+1)/2) = 0;
I_holo = abs(E_holo).*abs(E_holo);
SNR2 = I_holo((N+1)/2, (N+1)/2)/ (max(max(I_holo.*screen)));

end

