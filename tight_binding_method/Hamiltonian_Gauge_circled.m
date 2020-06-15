clc
clear all
close all
Nd = 7; 
t = 1; 
h = - 1;
k1 = exp(h*0.5); 
k2 = exp(-h*0.5); 
g = k1 - k2;

% h = -0.3466; h = 0.2466;

H = diag(g*1i*ones(1,Nd + 1)) + ...
    diag(t*k1*ones(1,Nd),1) + ...
    diag(t*k2*ones(1,Nd),-1);  % Assembling of Hamiltonian matrix
H(1,Nd + 1) = t*k2; H(Nd + 1,1) = t*k1;
% H((Nd + 1)/2 + 1,(Nd + 1)/2 + 1) = H((Nd + 1)/2 + 1,(Nd + 1)/2 + 1) - 0.1i;
% H(Nd + 1,Nd + 1) = 0.25i;
H(1,1) = g*1i;
[V, A] = eig(H);
lam = diag(A);
[~, idx] = sort(real(lam));
lam1 = lam(idx);

figure
% plot(real(lam1),'b*')
hold on
plot(imag(lam1),'r*')
set(gcf, 'Position', [00, 00, 350, 300])
axis([0 Nd + 1 -2.5 2.5])
set(gca,'FontSize', 14) % Font Size

bn = 5;
figure
bar(angle(V(:,bn)),'b')
hold on
plot(abs(V(:,bn)),'r*')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 14) % Font Size
axis([0 Nd + 1 -pi pi])

bn = 5;
figure
hold on
bar(abs(V(:,bn)),'b')
set(gcf, 'Position', [00, 00, 250, 200])
set(gca,'FontSize', 14) % Font Size
box on

phi = linspace(2*pi/(Nd + 1),2*pi,Nd + 1);
phi2 = circshift(phi,1);
xRing = [cos(phi);cos(phi)*1.5;cos(phi2)*1.5;cos(phi2)];
yRing = [sin(phi);sin(phi)*1.5;sin(phi2)*1.5;sin(phi2)];
ringData = angle(V(:,bn));

figure
patch(xRing,yRing,ringData, 'Edgecolor','none');
set(gca,'cLim',[-pi pi]);
axis square
axis off
set(gcf,'color','w');
colormap('hsv')
colorbar