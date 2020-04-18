clc
clear all
close all
Nd = 15; 
t = 1; 
h = - 0.02;
k1 = exp(4*h); 
k2 = exp(-0*h); 
% positive imaginary part is gain
g = k1 - k2;
g = 0.08;
Nr = (Nd + 1)/4;
F0 = sparse(zeros(Nr, Nr));
E1 = sparse(fliplr(diag([k2 0 0 k1])));
E2 = sparse(fliplr(diag([0 k1 k2 0])));
E2_ = sparse(fliplr(diag([0 k1 k2 0])));
E3 = sparse(fliplr(diag([k2 0 0 k1])));

D1 = diag(-g*1i*ones(1,Nr)) + diag(k1*ones(1,Nr-1), 1) + diag(k2*ones(1,Nr-1), -1);
D2 = diag(-g*1i*ones(1,Nr)) + diag([k1 0 k2], 1) + diag([k2 0 k1], -1);
D3 = diag(-g*1i*ones(1,Nr)) + diag([k1 0 k2], 1) + diag([k2 0 k1], -1);
D4 = diag(-g*1i*ones(1,Nr)) + diag(k1*ones(1,Nr-1), 1) + diag(k2*ones(1,Nr-1), -1);

H = sparse([D1 E1 F0 F0;
            E1 D2 E2 F0;
            F0 E2_ D3 E3;
            F0 F0 E3 D4]);

% H = sparse([D1 E1;
%             E1 D1;
%             ]);

[V, A] = eig(full(H));
lam = diag(A);
% [~, idx] = sort(real(lam));
% lam1 = lam(idx);
figure
plot(real(lam),'b*')
hold on
plot(imag(lam)*50,'r*')
set(gcf, 'Position', [00, 00, 350, 300])
% axis([0 Nd + 1 -2.5 2.5])
set(gca,'FontSize', 14) % Font Size
Lasing = find(imag(diag(V)) > 0);
Lasing=9

bn = 9
figure
bar(angle(V(:,bn)),'b')
hold on
plot(abs(V(:,bn))/max(abs(V(:,bn))),'r*')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 14) % Font Size
axis([0 Nd + 1 -pi pi])


Intensity = abs(V(:, Lasing)) .* abs(V(:, Lasing));
Intensity = Intensity ./ max(Intensity);
Phase = angle(V(:, Lasing));
Nall = sqrt(Nd+1);
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
