clc;
clear all;
close all;

global lambda k AmpImage N M dW z;

% Data of Amplitude on the Image Plane
%20 by 20 lasers
AmpImage = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0
        0 0 0 0 0 0 1 1 0 0 0 0 1 1 0 0 0 0 0 0
        0 0 0 0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
        0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
        0 0 0 0 1 1 0 0 0 0 0 0 0 0 1 1 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

        ];

PhaseImage = 2 * pi * (rand(20, 20) - 0.5);
AmpImage = AmpImage .* exp(1i * PhaseImage);
figure;
imagesc(AmpImage);
colormap gray
colorbar
title('Image Plane(desired pattern)')

dW = 10e-6; %Spacing between elements
z = 50e-6; % Distance of Image plane from the resonator plane

lambda = 1.550e-6;
k = 2 * pi / lambda;
[N, M] = size(AmpImage); % Image plane discretization size
E_Sample = zeros(N, M); % initial an array to store sample image for calculation

E_Sample = supperposition(E_Sample);

E_Sample = normalize(E_Sample);

figure;
imagesc(abs(E_Sample), [0 1]);
colormap jet
colorbar
title('Sample Plane Amplitude')

figure;
imagesc(angle(E_Sample));
colormap jet
colorbar
title('Sample Plane Phase')

E_Sample_Intensity = (abs(E_Sample)).^2;
E_Sample_Intensity = Binarize(E_Sample_Intensity);
figure;
imagesc(E_Sample_Intensity);
colormap jet
colorbar
title('Sample Plane Intensity')

E_Sample_Phase = quantizePhase(angle(E_Sample));
% E_Sample_Phase = (angle(E_Sample));
figure;
imagesc(E_Sample_Phase);
colormap jet;
colorbar;
title('Sample Plane Quantized Phase');

E_Sample_recovered = E_Sample_Intensity .* exp(-1i * E_Sample_Phase);
E_Image_recovered = supperposition(E_Sample_recovered);
E_Image_recovered = normalize(E_Image_recovered);
figure;
imagesc(abs(E_Image_recovered), [0.1 1]);
colormap jet;
colorbar;
title('Recovered Image Plane Intensity');

function E_Sample = supperposition(E_Sample)

    [M, N] = size(E_Sample);

    for ii = 1:M

        for jj = 1:N

            E_Sample = E_Sample + Green(ii, jj);
        end

    end

end

function quantizedPhase = quantizePhase(Array)

    [M, N] = size(Array);
    quantized = [-pi, -pi / 2, 0, pi / 2, pi]; % four level quantization

    for ii = 1:M

        for jj = 1:N
            [~, idx] = min(abs(quantized - Array(ii, jj)));
            Array(ii, jj) = quantized(idx);

        end

    end

    quantizedPhase = Array;
end

function BinarizedArray = Binarize(Array)

    [M, N] = size(Array);

    for ii = 1:M

        for jj = 1:N

            if (Array(ii, jj) >= 0.25)% half amplitude
                Array(ii, jj) = 1;
            else
                Array(ii, jj) = 0;
            end

        end

    end

    BinarizedArray = Array;

end

function E_Sample = Green(coordImageX, coordImageY)

    global lambda k AmpImage N M dW z;
    [N, M] = size(AmpImage);
    amplitude = AmpImage(coordImageX, coordImageY); % inversely solved electric field at the sample plane

    SamplePlane = zeros(N, M); % electric field at sample plane.

    [coordSampleX, coordSampleY] = meshgrid((1:N) * dW, (1:M) * dW);

    for ii = 1:N

        for jj = 1:M

            distance = (((coordImageX - ii) * dW)^2 + ((coordImageY - jj) * dW)^2 + z^2)^(-0.5);
            SamplePlane(ii, jj) = amplitude * distance * exp(1i * k * distance);

        end

    end

    E_Sample = SamplePlane;
end

function normalizedIntensity = normalize(Original)
    [M, N] = size(Original);
    Emax = max(max(abs(Original)));

    for ii = 1:M

        for jj = 1:N

            Original(ii, jj) = Original(ii, jj) / Emax;

        end

    end

    normalizedIntensity = Original;
end
