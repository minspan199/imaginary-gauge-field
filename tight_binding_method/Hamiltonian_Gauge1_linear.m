clc
clear all
close all

Nd = 8; t1 = 1; t2 = 1; g = 0.2; h = 0.01;loss = -0.;

% h = -0.3466; h = 0.2466;

H = diag([repmat([1i*loss -1i*g],[1 Nd]) 1i*g]) +...
    diag(repmat([t1*exp(h) t2*exp(h)],[1 Nd]),1) +...
    diag(repmat([t1*exp(-h) t2*exp(-h)],[1 Nd]),-1);  % Assembling of Hamiltonian matrix
% H(1,1) = 0.62i; H(2*Nd + 1,2*Nd + 1) = H(2*Nd + 1,2*Nd + 1) - 0.1i;
% H(2*Nd + 1,2*Nd + 1) = 0.81i; H(1,1) = H(1,1) - 0.5i;
[V, A] = eig(H);
lam = diag(A);
[~, idx] = sort(real(lam));
lam1 = lam(idx);

figure
% plot(real(lam1),'b*')
hold on
plot(imag(lam1),'r*')
set(gcf, 'Position', [00, 00, 350, 300])
% axis([0 2*Nd + 2 -4 4])
set(gca,'FontSize', 14) % Font Size
box on

bn = 9;
figure
bar((V(:,bn))/max(abs(V(:,bn))))
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 14) % Font Size
axis([0 2*Nd + 2 -1 1])



