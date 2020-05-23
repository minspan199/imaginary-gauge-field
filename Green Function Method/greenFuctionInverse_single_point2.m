clc;
clear all;
close all;

% Data of Amplitude on the Image Plane

global lambda k AmpImage N M dW z;

% Data of Amplitude on the Image Plane
%80 by 80 lasers
AmpImage = zeros(100, 100);
AmpImage(50, 50) = 1;

% PhaseImage = (rand(20, 20) - 0.5);
% AmpImage = AmpImage .* exp(1i * PhaseImage);
figure;
imagesc(abs(AmpImage));
colormap gray
colorbar
title('Image Plane(desired pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

dW = 0.01e-6; %Spacing between elements
z = 80e-6; % Distance of Image plane from the resonator plane

lambda = 1.550e-6;
k = 2 * pi / lambda;
[N, M] = size(AmpImage); % Image plane discretization size
E_Sample = zeros(N, M); % initial an array to store sample image for calculation

E_Sample = supperposition(E_Sample);

E_Sample = normalize(E_Sample);

figure;
imagesc(abs(E_Sample));
colormap jet
colorbar
title('Sample Plane Amplitude')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

figure;
imagesc(angle(E_Sample));
colormap jet
colorbar
title('Sample Plane Phase')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

function E_Sample = supperposition(E_Sample)

    [M, N] = size(E_Sample);

    for ii = 1:M

        for jj = 1:N

            E_Sample = E_Sample + Green(ii, jj);
            
        end

    end

end

function E_Sample = Green(coordImageX, coordImageY)

    global lambda k AmpImage N M dW z;
    [N, M] = size(AmpImage);
    amplitude = AmpImage(coordImageX, coordImageY); % inversely solved electric field at the sample plane

    SamplePlane = zeros(N, M); % electric field at sample plane.

%     [coordSampleX, coordSampleY] = meshgrid((1:N) * dW, (1:M) * dW);

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
