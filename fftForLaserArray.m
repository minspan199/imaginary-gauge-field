
clc;
clear all;
close all;

radiumOfPixel = 10e-6;
spotSize = 5e-6;
xNum = 6;
yNum = 6;
N = 2500; %number of pixels in each dimension (determines fidelity and processing time)


% Amp = zeros(xNum,yNum);
Amp = [0.9911505	0.982379315	0.973685749	0.965069118	0.956528739	0.948063938
1	9.07E-01	0.914947229	0.923116346	0.931358402	0.939674047
0.931358402	0.973685749	0.965069118	0.956528739	0.948063938	0.939674047
0.939674047	0.898825231	0.906850404	9.15E-01	0.923116346	0.931358402
0.875173319	0.965069118	0.956528739	0.948063938	0.939674047	0.931358402
0.882987315	0.890871078	0.898825231	0.906850404	0.914947229	0.923116346
];
AmpPhase = [1.570796327	3.141592654	-1.570796327	-4.21E-15	1.570796327	3.141592654
0	-1.570796327	3.14E+00	1.570796327	-4.71E-15	-1.570796327
-1.570796327	-3.37E-16	1.570796327	3.14E+00	-1.570796327	-3.69E-15
3.141592654	1.570796327	-4.09E-15	-1.570796327	-3.141592654	1.570796327
1.570796327	3.141592654	-1.570796327	-1.86E-16	1.570796327	3.141592654
2.15E-16	-1.570796327	3.141592654	1.570796327	-2.92E-15	-1.570796327
];

dx = 300e-9; %width of each pixel ???
xN = xNum*radiumOfPixel*2;
yN = yNum*radiumOfPixel*2;
screenSize = N*dx;


[xsNear, ysNear] = meshgrid((-N / 2:N / 2 - 1) .* dx); % defining near field grid + xN/2
lambda = 1550e-9;
k = 2*pi/lambda;
z = 10;

% grid_size = N * dx; % number of pixels in one dimension times pixel width
% [xs ys] = meshgrid((-xN / 2:xN / 2 - 1)); % defining near field grid
dxFar = lambda * z / screenSize; % far field pixel width
[xdFar, ydFar] = meshgrid((-N / 2:N / 2 - 1) .* dxFar); % defining far field grid + yN/2

nearField = zeros(N,N);
nearFields{xNum, yNum} = nearField;
farField = zeros(N,N);
farFields{xNum, yNum} = farField;

N_LOOP = xNum*yNum;

parfor loopNum = 1:N_LOOP

    [ii, jj] = getInd(loopNum-1, xNum);
    nearFields{loopNum} =  Amp(ii, jj) * exp((-((xsNear - (jj - 1) * radiumOfPixel).^2 ...
            + (ysNear - (ii - 1) * radiumOfPixel).^2) / (spotSize^2))) * exp(1i * AmpPhase(ii, jj));
    nearField = nearField + nearFields{loopNum};
end

figure;
imagesc(abs(nearField));
figure;
imagesc(imag(nearField));

parfor loopNum = 1:N_LOOP
    
   [ii, jj] = getInd(loopNum-1, xNum);

       Angular_Spectrum = fftshift(fft2(fftshift(nearFields{loopNum}))); %Near Field
       farFields{loopNum} = exp(1i*k*z)*exp(1i*k*(xdFar.^2+ydFar.^2)/(2*z)).*Angular_Spectrum/(1i*lambda*z);
       farField = farField + farFields{loopNum};
   
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
        'LineStyle', 'none', 'FaceColor', 'interp', 'FaceLighting', 'phong', ...
        'AmbientStrength', 0.3), shading flat;
axis([-theta_sim theta_sim -theta_sim theta_sim 0 1]); % setting axis limits
set(gca, 'Visible', 'off', 'plotboxaspectratio', [1, 1, 3]);
camva(3);
grid off;
view([0 90]);
axis on

figure;
imagesc(angleX_Min, -angleY_Min, angle(farFieldNormalizedIntensity));
axis([-theta_sim theta_sim -theta_sim theta_sim]);

AmpIFFT = fftshift(ifft2(fftshift(farField))); %Near Field
bearFieldsIFFT = exp(1i*k*z)*exp(1i*k*(xsNear.^2+ysNear.^2)/(2*z)).*AmpIFFT/(1i*lambda*z);

figure;
imagesc(abs(AmpIFFT))

figure;
imagesc(angle(AmpIFFT))


function [ind1, ind2] = getInd(loopNum, size)
    ind1 = floor(loopNum/size) + 1;
    ind2 = rem(loopNum, size) + 1;
end


