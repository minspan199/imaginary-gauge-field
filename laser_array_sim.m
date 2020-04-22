%Program:   VCSEL array model
%Author:    Matthew Johnson
%Date:      June 2013
%Modeling far field from array of Gaussian sources defined by parameters in excel file

% DEFINING SOME INITIAL VALUES
clear all;
close all;
ext = '';
column = 2;
filename_data = 'variables.xlsx'; % name of file containing simulation parameters
num_var = 26; % number of variable parameters in excel document, have to update this number if you add another variable
sheet = 3; % sheet the excel data is on
plot_style_vec = {'b'; ':r'; '--c'; '-.m'}; % linestyles to plot
legend_vec = {'4'; '3'; '2'; '1'}; % strings to plot in legend

% READING IN PARAMETERS FROM EXCEL FILE
[num, txt, raw] = xlsread(filename_data, sheet, ['A2:B', num2str(num_var + 1)]); % pulling parameters from excel
vars = genvarname(raw(1:num_var, 1)); % assigning values to variable names defined in excel file

for ii = 1:length(vars)
    str = char(vars(ii));
    eval([str '=num(ii,1)'])
end

% OBTAINING INTENSITY AND PHASE MAP
[num1, txt1, raw1] = xlsread(filename_data, map_sheet, 'A1:H25'); %intensity/phase info in different sheet

% DEFINING/INITIALIZING OTHER VARIABLES
N_theta = ceil(theta_sim * N / 10);
FWHMs = zeros(1, map_N);
Intensity_max = zeros(1, map_N);

% INITIATING LOOP TO GO THROUGH MULTIPLE INTENSITY/PHASE MAPS
for NN = 1:map_N
    Amp = num1((NN - 1) * 10 + 1:(NN - 1) * 10 + map_dim, 1:map_dim); % intensity and phase maps offset by 10 rows
    Phi = num1((NN - 1) * 10 + map_dim + 2:(NN - 1) * 10 + 2 * map_dim + 1, 1:map_dim);
    plot_style = char(plot_style_vec(NN));
    fignum = (NN - 1) * 10; % different figure numbers for different maps

    if NN > 1
        cleared = 0;
    end

    % GAUSSIAN VALUES %
    k = 2 * pi / lambda; % wave vector

    % SAMPLE AND DETECTOR GRIDS %
    grid_size = N * dx; % number of pixels in one dimension times pixel width
    [xs ys] = meshgrid((-N / 2:N / 2 - 1) .* dx); % defining near field grid
    dx_det = lambda * z / grid_size; % far field pixel width
    [x_det y_det] = meshgrid((-N / 2:N / 2 - 1) .* dx_det); % defining far field grid

    % INITIALIZING APERTURES %
    U_ap_total = zeros(N);
    U_total_Phi = zeros(N);
    Intensity_total_incoh = zeros(N);

    for ii = 1:map_dim

        for jj = 1:map_dim
            U_ap(:, :, ii, jj) = zeros(N); %setting near field matrices to zero
            U(:, :, ii, jj) = zeros(N); %setting far field matrices to zero
        end

    end

    % Assigning Gaussian values to apertures with corresponding amplitude/phase
    %%% if hexagonal near-field pattern
    if Hex == 1
        x_pitch = y_pitch * sin(pi / 3);

        for ii = 1:map_dim

            for jj = 1:map_dim

                if mod(jj, 2) == 0    % if ii is even / xs - near field grids
                    U_ap(:, :, ii, jj) = Amp(ii, jj) * exp(-(((xs - (jj - 1) * x_pitch) * 0).^2 ...
                        + ((ys - (ii * 2 - 1) * y_pitch / 2) * 0).^2) / (w^2)) * exp(1i * Phi(ii, jj));
                    else                 % if ii is odd
                    U_ap(:, :, ii, jj) = Amp(ii, jj) * exp(-(((xs - (jj - 1) * x_pitch) * 0).^2 ...
                        + ((ys - (ii - 1) * y_pitch) * 0).^2) / (w^2)) * exp(1i * Phi(ii, jj));
                end

            end

        end

        %%% if rectangular near-field pattern
    elseif Hex == 0
        x_pitch = y_pitch;

        for ii = 1:map_dim

            for jj = 1:map_dim

                U_ap(:, :, ii, jj) = Amp(ii, jj) * exp((-((xs - (jj - 1) * x_pitch).^2 ...
                    + (ys - (ii - 1) * y_pitch).^2) / (w^2))) * exp(1i * Phi(ii, jj));
            end

        end

    end

    for ii = 1:map_dim

        for jj = 1:map_dim
            U(:, :, ii, jj) = Fraunhoffer(U_ap(:, :, ii, jj), k, z, x_det, y_det, lambda);
            U_ap_total = U_ap_total + U_ap(:, :, ii, jj);
            U_total_Phi = U_total_Phi + U(:, :, ii, jj);
            Intensity_total_incoh = Intensity_total_incoh + (abs(U(:, :, ii, jj)).^2);
        end

    end

    %ADDING COMPLEX FIELDS AND SOLVING FOR FAR FIELD INTENSITY
    % from Equations 2.2-2.5
    Intensity_total_coh = U_total_Phi.^2;
    Intensity_total_Phi = (1 - coh) * Intensity_total_incoh + coh * ...
        Intensity_total_coh;
    % normalizing intensity
    Intensity_normalized_Phi = Intensity_total_Phi / max(max(Intensity_total_Phi));

    % PLOTTING NEAR AND FAR FIELDS
    angle_x = atan(x_det / z) * 360 / (2 * pi); %translating to angular values
    angle_y = atan(y_det / z) * 360 / (2 * pi);

    % PLOT NEAR FIELD
    Fig2 = figure(2 + fignum);

    if cleared == 1
        clf;
    end

    set(Fig2, 'Position', [0 0 350 300]);
    colormap(jet);
    surf(xs * 10^6, -ys * 10^6, abs(U_ap_total), ...
        'LineStyle', 'none', 'FaceColor', 'interp', ...
        'FaceLighting', 'phong', 'AmbientStrength', 0.3);
    shading flat;
    axis([x_ap_min x_ap_max y_ap_min y_ap_max]); % setting axis limits

    if sum(sum(abs(Phi))) > 1
        caxis([-1 1]);
    else
        caxis([0 1]);
    end

    set(gca, 'dataaspectratio', [1, 1, 1]);
    grid off; view([0 90]);
    xlabel('\mum', 'FontSize', fs, 'FontName', 'Calibri');
    ylabel('\mum', 'FontSize', fs, 'FontName', 'Calibri');
    set(gcf, 'PaperPositionMode', 'auto');

    % PLOT FAR FIELD PHASE
    Fig3 = figure(8 + column + fignum); clf;
    set(Fig3, 'Position', [300 + 50 * (column - 1) 400 350 300]);
    surf(angle_x, -angle_y, angle(U_total_Phi), ...
        'LineStyle', 'none', 'FaceColor', 'interp', 'FaceLighting', 'phong', ...
        'AmbientStrength', 0.3), shading flat;
    axis([-theta_sim theta_sim -theta_sim theta_sim 0 1]); % setting axis limits
    set(gca, 'Visible', 'off', 'plotboxaspectratio', [1, 1, 3]);
    %     camva(3);
    grid off;
    view([0 90]);
    figure; imagesc(angle(U_total_Phi))

    % PLOT FAR FIELD
    Fig3 = figure(8 + column + fignum + 10); clf;
    set(Fig3, 'Position', [300 + 50 * (column - 1) 400 350 300]);
    surf(angle_x, -angle_y, abs(Intensity_normalized_Phi), ...
        'LineStyle', 'none', 'FaceColor', 'interp', 'FaceLighting', 'phong', ...
        'AmbientStrength', 0.3), shading flat;
    axis([-theta_sim theta_sim -theta_sim theta_sim 0 1]); % setting axis limits
    set(gca, 'Visible', 'off', 'plotboxaspectratio', [1, 1, 3]);
    camva(3);
    grid off;
    view([0 90]);

    % adding 10 degree circle and cut line
    hold on;
    x1 = 0:pi / 500:10;
    t = ones(1, length(x1));
    %     plot3(theta_circle * cos(x1), theta_circle * sin(x1), t, 'r', 'LineWidth', 2);
    %     text(0.8 * theta_circle, 0.8 * theta_circle, 1, [num2str(theta_circle), 'ÃƒÆ’Ã‚Â¯Ãƒâ€šÃ‚Â¿Ãƒâ€šÃ‚Â½'], 'Color', 'r', 'FontSize', fs);

    if plot_cutline == 1   % plots line along slice angle
        plot3(0, 0, 1, '.k');
        plot3(x1, x1 * tan(angcut * 2 * pi / 360), t, 'k', 'Linewidth', 2);
        plot3(-x1, -x1 * tan(angcut * 2 * pi / 360), t, 'k', 'Linewidth', 2);
    end

    hold off;

    % Calculating power within 'theta_center_lobe' degrees
    power_total = sum(sum(abs(Intensity_normalized_Phi)));
    power_center_lobe = 0;

    for ii = 1:length(xs)

        for jj = 1:length(ys)

            if sqrt(angle_x(ii, jj)^2 + angle_y(ii, jj)^2) <= theta_center_lobe
                power_center_lobe = power_center_lobe + abs(Intensity_normalized_Phi(ii, jj));
            end

        end

    end

    percent_center_lobe = power_center_lobe / power_total;

    % PLOT FAR FIELD INTENSITY SLICE
    % Plots intensity along slice angle
    Fig8 = figure(8);
    set(Fig8, 'Position', [10 10 500 300]);

    if cleared == 1
        clf;
    end

    hold on;
    % interpolating to find values along cut
    point1 = [0; 0];
    point2 = [theta_circle; theta_circle * tan(-angcut * 2 * pi / 360)];
    t = 0:1 / (N_theta):1;
    xi = point2(1) * t + point1(1) * (1 - t);
    yi = point2(2) * t + point1(2) * (1 - t);
    Z_interp1 = interp2(angle_x, angle_y, abs(Intensity_normalized_Phi), xi, yi);
    Z_interp2 = interp2(angle_x, angle_y, abs(Intensity_normalized_Phi), -xi, -yi);
    xyi_both(N_theta + 1:-1:1) = -sqrt(xi.^2 + yi.^2);
    xyi_both(N_theta + 2:2 * (N_theta + 1) - 1) = sqrt(xi(2:N_theta + 1).^2 + yi(2:N_theta + 1).^2);
    Z_both(N_theta + 1:-1:1) = Z_interp2;
    Z_both(N_theta + 2:2 * (N_theta + 1) - 1) = Z_interp1(2:N_theta + 1);
    plot(xyi_both, Z_both, plot_style, 'linewidth', 2);
    axis([-theta_sim * sim_trim theta_sim * sim_trim 0 1]);
    figure(8);
    xlabel('angle(\circ)', 'FontSize', fs, 'FontName', 'Calibri');
    ylabel('Intensity (a.u.)', 'FontSize', fs, 'FontName', 'Calibri');

    if plot_legend == 1
        legend(legend_vec);
    end

end
