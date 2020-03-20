function plotring(E)
    %E: target intensity array
    Ny = length(E(1, :)); % ring number in y
    Nx = length(E(:, 1)); % ring number in x

    tt = 50; % 50*50 pixels per ring, must be changed together with file edgeplot1.m
    f = zeros(tt * Nx, tt * Ny);

    for i = 1:Nx

        for j = 1:Ny
            f = edgeplot1(i, j, E(i, j) + 0.005, f); % 0.005 is an offset to plot the ring shape even at 0 intensity site
        end

    end

    %start to plot
    figure;
    ha = bar3(f); 
    zlim([0 1.2]);
    box on;
    set(gca,'AmbientLightColor',...
    [0.941176470588235 0.941176470588235 0.941176470588235],'Color',...
    [0.941176470588235 0.941176470588235 0.941176470588235],...
    'PlotBoxAspectRatio',[1.60837408433945 1.61803398874989 1])
    colormap('hot');

    for n = 1:numel(ha)
        cdata = get(ha(n), 'zdata'); set(ha(n), 'cdata', cdata, 'facecolor', 'interp');
    end

    set(ha, 'edgecolor', 'none');
