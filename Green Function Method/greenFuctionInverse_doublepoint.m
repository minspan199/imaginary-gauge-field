clc;
clear all;
close all;

global lambda k AmpImage N M dW z;

% Data of Amplitude on the Image Plane
AmpImage = [0 0 0 0 0 0 0 0 0 0
        0 0 0 0 1 0 0 0 0 0
        0 0 0 1 1 1 0 0 0 0
        0 0 0 0 1 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0
        0 0 0 0 1 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0
        ];
AmpImage = [0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0
        0 0 0 0 1 1 0 0 0 0
        0 0 0 0 1 1 0 0 0 0
        0 0 1 1 1 1 1 1 0 0
        0 0 1 1 1 1 1 1 0 0
        0 0 0 0 1 1 0 0 0 0
        0 0 0 0 1 1 0 0 0 0
        0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0
        ];
AmpImage = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
        ];

figure;
imagesc(AmpImage);
colormap gray
colorbar
title('Image Plane(desired pattern)')

dW = 10e-6; %Spacing between elements
z = 100e-6; % Distance of Image plane from the resonator plane

lambda = 1.550e-6;
k = 2 * pi / lambda;
[N, M] = size(AmpImage); % Image plane discretization size
E_Sample = zeros(N, M); % initial an array to store sample image for calculation

for ii = 1:N

    for jj = 1:M

        E_Sample = E_Sample + Green(ii, jj);
    end

end

figure;
imagesc(abs(E_Sample));
colormap jet
colorbar
title('Sample Plane Intensity')

figure;
imagesc(angle(E_Sample));
colormap jet
colorbar
title('Sample Plane Phase')

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
