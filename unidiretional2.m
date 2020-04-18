clc
clear all
close all
Nd = 15;
t = 1;
h = -0.01;
k1 = exp(-h);
k2 = exp(h);
% k_link_1 = 5;
ind = 1;
g = 0.01;
% k_link_1 = 5;
k_link_2 = 0.00;

for k_link_1 = 0:0.01:5

    %     k_link_2 = 1 / k_link_1;
    % positive imaginary part is gain

    Nr = (Nd + 1) / 4;
    F0 = sparse(zeros(Nr, Nr));
    E1 = sparse(fliplr(diag([k2 0 0 k1])));
    E2 = sparse(fliplr(diag([0 k1 k2 0])));
    E2_ = sparse(fliplr(diag([0 k1 k2 0])));
    E3 = sparse(fliplr(diag([k2 0 0 k1])));

    D1 = diag(-g * 1i * ones(1, Nr)) + diag(k1 * ones(1, Nr - 1), 1) + diag(k2 * ones(1, Nr - 1), -1);
    D2 = diag(-g * 1i * ones(1, Nr)) + diag([k1 0 k2], 1) + diag([k2 0 k1], -1);
    D3 = diag(-g * 1i * ones(1, Nr)) + diag([k1 0 k2], 1) + diag([k2 0 k1], -1);
    D4 = diag(-g * 1i * ones(1, Nr)) + diag(k1 * ones(1, Nr - 1), 1) + diag(k2 * ones(1, Nr - 1), -1);

    H = sparse([D1 E1 F0 F0;
            E1 D2 E2 F0;
            F0 E2_ D3 E3;
            F0 F0 E3 D4]);
    H(Nd + 2, Nd + 2) = 0;
    H(Nd + 1, Nd + 2) = k_link_1;
    H(Nd + 2, Nd + 1) = k_link_2;
    [V, A] = eig(full(H));
    lam = diag(A);
    Lasing = find(abs(real(lam)) < 1e-6);
    data(ind, :) = [lam(Lasing)];
    ind = ind + 1
end

g = 0.08;

Lasing = 9

bn = 9
figure
bar(angle(V(:, bn)), 'b')
hold on
plot(abs(V(:, bn)) / max(abs(V(:, bn))), 'r*')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca, 'FontSize', 14) % Font Size
axis([0 Nd + 1 -pi pi])

Intensity = abs(V(:, Lasing)) .* abs(V(:, Lasing));
Intensity = Intensity ./ max(Intensity);
Phase = angle(V(:, Lasing));
Nall = sqrt(Nd + 1);

for k = 1:1:Nall

    if mod(k, 2)
        RT(k, :) = Intensity((k - 1) * Nall + (1:Nall));
    else
        RT(k, :) = fliplr(Intensity((k - 1) * Nall + (1:Nall))');
    end

end

for k = 1:1:Nall

    if mod(k, 2)
        RA(k, :) = Phase((k - 1) * Nall + (1:Nall));
    else
        RA(k, :) = fliplr(Phase((k - 1) * Nall + (1:Nall))');
    end

end

plotring(RT, [0 1.2]);
set(gcf, 'Position', [00, 00, 400, 300]);
set(gca, 'FontSize', 14);

figure;
imagesc(RA);
colormap([0.5 0.5 0.5; 0 1 1; 0 1 0; 1 0.5 0; 0.5 0.5 0.5]);
colorbar;
set(gcf, 'Position', [00, 00, 400, 300]);
set(gca, 'FontSize', 14);
