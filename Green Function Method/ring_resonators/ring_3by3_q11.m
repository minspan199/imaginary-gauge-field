clc;
clear all;
close all;

%% define a backbone
global lambda k AmpImage N M meshSize z;
N = 3; %3by3
space = 200e-6; %Spacing between elements
z = 1000e-6; % Distance of Image plane from the resonator plane

lambda = 1.550e-6;
k = 2 * pi / lambda;

Amp_sample = ones(N, N);
Phase_sample = [-2.09439510239320	0	0
            0	0	2.09439510239320
            0	2.09439510239320	2.09439510239320];

% sizeOfResonator = 10e-6; %10um.
spotSize = 20e-6;
meshSize = 10000E-9;

boundary = 00e-6;

nearFields{N, N} = 0;

totalMesh = (space + boundary + spotSize) * (N + 1) / meshSize;

[xsNear, ysNear] = meshgrid((0:totalMesh) .* meshSize); % defining near field grid + xN/2

nearField = 0;

%% Generate the sample plane profile
for col = 1:N

    for row = 1:N
        nearFields{row, col} = Amp_sample(row, col) * exp((-((xsNear - boundary * 2 - (col - 0.25) * space).^2 ...
            + (ysNear - boundary * 2 - (row - 0.25) * space).^2) / (spotSize^2))) * exp(1i * Phase_sample(row, col));
        nearField = nearField + nearFields{row, col};
    end

end

figure;
imagesc(abs(nearField));

figure;
imagesc(angle(nearField));
colorbar

%% calculate the hologram
E_holo = HologramHelperClass2.supperposition(conj(HologramHelperClass.normalize(nearField)));

figure;
imagesc(abs((E_holo .* E_holo)));

figure;
imagesc(angle(E_holo));
