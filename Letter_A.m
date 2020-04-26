clc;
clear all;
close all;

spotSize = 5e-6;
dx = 100e-9; %width of each pixel ???

if exist('letterA.mat',  'file')
    load('letterA.mat')
else
    A = importdata("A.txt");
    A = 1 ./ (1 + A);
    A = floor(A);
end

[width, length] = size(A);

if width ~= length
    A = convert2Square(A);
end

[N, ~] = size(A);
screenSize = N * dx;
% parfor ii = 1:width
%
%     for jj = 1:length
%         if(A(ii, jj)>0.5)
%             Phi(ii, jj) = 0;
%         else
%             Phi(ii, jj) = 2*pi*(rand(1, 1) - 0.5);
%         end
%     end
% end

[xsNear, ysNear] = meshgrid((1:N) * dx); % defining near field grid
lambda = 1550e-9;
k = 2 * pi / lambda;

Phi = 2 * pi * cos(k * xsNear) .* cos(k * ysNear);

z = 10;

figure;
imagesc(abs(A))
figure;
imagesc(Phi)

%
if ~exist('U_ap',  'var')

    U_ap = zeros(width, length);

    % Phi = 2*pi*cos(0.05*ys);

%     parfor ii = 1:width
% 
%         for jj = 1:length
% 
%             U_ap = U_ap + A(ii, jj) * exp((-((xsNear - (jj - 1) * spotSize).^2 ...
%                 + (ysNear - (ii - 1) * spotSize).^2) / (20^2))) .* exp(1i * Phi);
%         end
% 
%     end

end

figure;
imagesc(abs(U_ap));

figure;
imagesc(angle(U_ap));

AmpIFFT = fftshift(ifft2(fftshift(U_ap))); %Near Field
nearFieldsIFFT = exp(1i * k * z) * exp(1i * k * (xsNear.^2 + ysNear.^2) / (2 * z)) .* AmpIFFT / (1i * lambda * z);

figure;
imagesc(abs(nearFieldsIFFT))

figure;
imagesc(angle(nearFieldsIFFT))

function SMatrix = convert2Square(A)

    [width, length] = size(A);
    N = width + length;
    left = floor(length / 2);
    right = width - left;
    B = [zeros(width, left) A zeros(width, right)];
    top = floor(width / 2);
    bottom = length - top;
    SMatrix = [zeros(top, N); B; zeros(bottom, N)];
end
