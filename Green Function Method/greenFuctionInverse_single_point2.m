clc;
clear all;
close all;

% Data of Amplitude on the Image Plane

global lambda k AmpImage N M dW z;
N = 100;
M = 100;

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

dW = 0.01e-6; %Spacing between elements
z = 80e-6; % Distance of Image plane from the resonator plane

lambda = 1.550e-6;
k = 2 * pi / lambda;

E_Sample = zeros(N, M); % initial an array to store sample image for calculation
E_Sample(50, 50) = 1;
[M, N] = size(E_Sample);


E_Holo = supperposition(E_Sample);

E_Holo = normalize(E_Holo);

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

function E_Sample2 = supperposition(E_Sample)

    [M, N] = size(E_Sample);
    E_Sample2 = zeros(M, N);
    
    for ii = 1:M

        for jj = 1:N

            E_Sample2 = E_Sample2 + Green(ii, jj, E_Sample);
            
        end

    end

end

function E_Sample = Green(coordImageX, coordImageY, E_Sample)

    global lambda k AmpImage N M dW z;
    [N, M] = size(E_Sample);
    amplitude = E_Sample(coordImageX, coordImageY); % inversely solved electric field at the sample plane

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
