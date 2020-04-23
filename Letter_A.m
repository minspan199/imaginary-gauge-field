clc;
clear all;
close all;

A = importdata("A.txt");

A = 1 ./ (1 + A);
A = floor(A);
figure;
imagesc(A)
x_pitch = 1;
y_pitch = 1;

[width, length] = size(A);
Phi = 2i*pi*(rand(width, length) - 0.5);

[xs, ys] = meshgrid(1:length, 1:width); % defining near field grid
U_ap = zeros(width, length);

i=0;
parfor ii = 1:width
    i=i+1
    for jj = 1:length

        U_ap = U_ap + A(ii, jj) * exp((-((xs - (jj - 1) * x_pitch).^2 ...
            + (ys - (ii - 1) * y_pitch).^2) / (5^2))) * exp(1i * Phi(ii, jj));
    end

end

figure;
imagesc(U_ap);

a = fftshift(ifft2(fftshift(Ap)));

figure;
imagesc(abs(a));