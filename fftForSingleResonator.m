clc;
clear all;
close all;

sizeOfResonator = 10e-6;
spotSize = 3e-6;
xNum = 6;
yNum = 6;
N = 2000; %number of pixels in each dimension (determines fidelity and processing time)

% Amp = zeros(xNum,yNum);
Amp = [0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 0 0 1 0 0 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    0 0 0 0 0 0 0
    ];

AmpPhase = [
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 pi / 2 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        ];

dx = 500e-9; %width of each pixel ???
xN = xNum * sizeOfResonator * 2;
yN = yNum * sizeOfResonator * 2;
screenSize = N * dx;

[xsNear, ysNear] = meshgrid((-N / 2:N / 2 - 1) .* dx); % defining near field grid + xN/2
lambda = 1550e-9;
k = 2 * pi / lambda;
z = 10e1;

% grid_size = N * dx; % number of pixels in one dimension times pixel width
% [xs ys] = meshgrid((-xN / 2:xN / 2 - 1)); % defining near field grid
dxFar = lambda * z / screenSize; % far field pixel width
[xdFar, ydFar] = meshgrid((-N / 2:N / 2 - 1) .* dxFar); % defining far field grid + yN/2

nearField = zeros(N, N);
nearFields{xNum, yNum} = nearField;
farField = zeros(N, N);
farFields{xNum, yNum} = farField;

N_LOOP = xNum * yNum;

parfor loopNum = 1:N_LOOP

    [ii, jj] = getInd(loopNum - 1, xNum);
    nearFields{loopNum} = Amp(ii, jj) * exp((-((xsNear - (jj - 1) * sizeOfResonator).^2 ...
        + (ysNear - (ii - 1) * sizeOfResonator).^2) / (spotSize^2))) * exp(1i * AmpPhase(ii, jj));
    nearField = nearField + nearFields{loopNum};

end

figure;
imagesc(abs(nearField));
figure;
imagesc(angle(nearField));
Angular_Spectra = zeros(N, N);

parfor loopNum = 1:N_LOOP

    [ii, jj] = getInd(loopNum - 1, xNum);

    Angular_Spectrum{loopNum} = fftshift(fft2(fftshift(nearFields{loopNum}))); %Near Field
    farFields{loopNum} = exp(1i * k * z) * exp(1i * k * (xdFar.^2 + ydFar.^2) / (2 * z)) .* Angular_Spectrum{loopNum} / (1i * lambda * z);
    farField = farField + farFields{loopNum};
    Angular_Spectra = Angular_Spectra + Angular_Spectrum{loopNum};
    
end

theta_sim = 10;

% PLOTTING FAR FIELDS
angle_x = atan(xdFar / z) * 360 / (2 * pi); %translating to angular values
angle_y = atan(ydFar / z) * 360 / (2 * pi);

angleX_Min = min(min(angle_x)); %translating to angular values
angleX_Max = max(max(angle_x));
angleY_Min = min(min(angle_y)); %translating to angular values
angleY_Max = max(max(angle_y));

farFieldIntensity = farField.^2;
% normalizing intensity
farFieldNormalizedIntensity = farFieldIntensity / max(max(farFieldIntensity));

figure;
surf(angle_x, -angle_y, abs(farFieldNormalizedIntensity), ...
    'LineStyle',  'none',  'FaceColor',  'interp',  'FaceLighting',  'phong', ...
    'AmbientStrength', 0.3), shading flat;
axis([-theta_sim theta_sim -theta_sim theta_sim 0 1]); % setting axis limits
set(gca,  'Visible',  'off',  'plotboxaspectratio', [1, 1, 3]);
camva(3);
grid off;
view([0 90]);
axis on

figure;
imagesc(angle(farField));
% axis([-theta_sim theta_sim -theta_sim theta_sim]);

AmpIFFT = fftshift(ifft2(fftshift(farField))); %Near Field
nearFieldsIFFT = exp(1i * k * z) * exp(1i * k * (xsNear.^2 + ysNear.^2) / (2 * z)) .* AmpIFFT / (1i * lambda * z);

figure;
imagesc(abs(nearFieldsIFFT))

figure;
imagesc(angle(nearFieldsIFFT))

function [ind1, ind2] = getInd(loopNum, size)
    ind1 = floor(loopNum / size) + 1;
    ind2 = rem(loopNum, size) + 1;
end
