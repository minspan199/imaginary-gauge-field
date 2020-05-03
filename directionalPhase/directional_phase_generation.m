clc;

global dx dy N
N = 20; %20 by 20 lasers
Phase = zeros(N, N);
dx = [1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0
0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0
0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0
0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0
0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0
0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0
0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0
0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0
0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0
0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0
0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
    ];
dy = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
1	-1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	-1
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
    ];

[coordX, coordY] = meshgrid(1:N, 1:N);

%% Generation of one Euler Travel
CurrentPhase = 0;
indX = 1;
indY = 1;
CoordCollec{N * N} = 0;

for ii = 1:N * N
    
    CoordCollec{ii} = [indX, indY];
    temptY = indY + dx(indX, indY);
    temptX = indX - dy(indX, indY);
    indY = temptY;
    indX = temptX;
    
        
end


for ii = 1:N * N
    
    indX = CoordCollec{ii}(1);
    indY = CoordCollec{ii}(2);
    Phase(indX, indY) = CurrentPhase; 
    CurrentPhase = CurrentPhase + pi / 2;
    if (CurrentPhase > pi)
        CurrentPhase = CurrentPhase - 2*pi;
    end
    
end


figure
quiver(coordX, flipud(coordY), dx, dy);

figure
imagesc(Phase, [-pi, pi])
colorbar

save('Phase.mat','Phase')

nearField = abs(ones(N, N)) .* (exp(1i * Phase));
Angular_Spectrum = (fft2(fftshift(nearField))); %Near Field
figure
imagesc(abs(Angular_Spectrum).^2);
set(gcf,  'Position', [00, 00, 350, 300])
set(gca,  'FontSize', 10)% Font Size
colormap jet
colorbar
title('Recovered Far Field Image')

function [indX, indY] = getInd(ii)
    global dx dy N

    if (ii <= N)
        indX = mod(ii, N);

    end

end
