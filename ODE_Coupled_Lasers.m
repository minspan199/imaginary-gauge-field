clc 
clear all 
close all 
global Nd t0 t1 t2 g h N0 tau_p pA pB dB tau_s alpha H
Nd = 15; t1 = 1; t2 = 3; t0 = 1; g = 0.5; h = -0.5493; N0 = 1e24;
tau_p = 0.02;  pB = g*tau_p; dB = 0; tau_s = 2000*tau_p; alpha = 3; 

k1 = exp(h); k2 = exp(-h); g0 = k2 - k1;pA = g*tau_p;
H = - diag(g0*ones(1,Nd + 1)) ...
    - 1i*diag(t0*k1*ones(1,Nd),1) ...
    - 1i*diag(t0*k2*ones(1,Nd),-1);  % Assembling of Hamiltonian matrix
H(1,Nd + 1) = -1i*t0*k2; H(Nd + 1,1) = -1i*t0*k1;

T_span = [0 30000]; 
y0 = [rand(Nd + 1,1);
    rand(Nd,1); 
    zeros(Nd + 1,1)]; 
[T,Y] = ode45(@rate_equation_SSH,T_span,y0); 

% h = -0.3466; h = 0.2466;

figure
N_span = length(T);
for x0 = 1:1:Nd + 1
    n0 = 1:10:N_span;
    Yi = angle(Y(n0,x0));
    p = plot3(T(n0),x0*ones(length(n0),1),Yi);
    hold on
end
set(gcf, 'Position', [00, 00, 350, 300])
zlim([-3.5,3.5])
xlabel('Time')
ylabel('Site number')
zlabel('Phase')
% subplot(2,1,2)

figure
% subplot(2,1,1) 
plot(T,abs(Y(:,1:Nd + 1))) 
ylim([-0 0.6])
title('Density carriers') % carriers density in high laser level 
xlabel('Time')
ylabel('Density')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 14) % Font Size

figure
Hf = surf(abs(Y(:,1:Nd + 1)));
set(gcf, 'Position', [00, 00, 350, 300])
set(Hf,'LineStyle','none');
zlim([0,0.5])
ylabel('Time')
xlabel('Site number')
zlabel('Carrier densities')

colormap(jet)

figure
N_span = length(T);
for x0 = 1:1:Nd + 1
    n0 = 1:10:N_span;
    Yi = abs(Y(n0,x0));
    plot3(T(n0),x0*ones(length(n0),1),Yi,'LineWidth',1);
    hold on
end
set(gcf, 'Position', [00, 00, 350, 300])
zlim([0,0.5])
xlabel('Time')
ylabel('Site number')
zlabel('Carrier densities')
% subplot(2,1,2)

figure
plot(T,angle(Y(:,1:Nd + 1))) 
title('Density of the photons in the active medium') % photons density in activer region And this is the function. 
set(gca,'FontSize', 12) % Font Size
set(gcf, 'Position', [00, 00, 350, 300])
xlabel('Time')
ylabel('Phase')

figure
subplot(2,1,1) 
plot(T,abs(Y(:,Nd + 2:2*Nd + 1))) 
ylim([-1 2])
title('Densità dei portatori nel livello laser superiore') % carriers density in high laser level 
subplot(2,1,2) 
plot(T,angle(Y(:,Nd + 2:2*Nd + 1))) 
title('Densità dei fotoni nel mezzo attivo') % photons density in activer region And this is the function. 

clearvars -except Y T

function dy = rate_equation_SSH(t,y) 

global Nd t0 t1 t2 g h N0 tau_p pA pB dB tau_s alpha H
y1 = y(1:Nd + 1);
y2 = y(Nd + 2:2*Nd + 1);
y3 = y(2*Nd + 2:3*Nd + 2);
dy1 = (1 - 1i*alpha)*y3.*y1 - 1i*t1*tau_p*exp(h)*[y2;0] - 1i*t2*tau_p*exp(-h)*[0;y2];
dy2 = -(pB - 1i*dB*tau_p)*y2 - 1i*t1*tau_p*exp(-h)*y1(1:Nd) - 1i*t2*tau_p*exp(h)*y1(2:Nd + 1);
dy3 = pA*tau_p/tau_s - tau_p/tau_s*y3 - tau_p/tau_s*(1 + 2*y3).*abs(y1).*abs(y1);
dy = [dy1;dy2;dy3];
end

% this is a very simple (dimensionless) and efficient version of the rate equations (note you % can add another equation, dy(3)(not coupled), that gives you the phase(t) of your laser!
function dy = rate_equation_circular(t,y) 

global Nd t0 t1 t2 g h N0 tau_p pA pB dB tau_s alpha H
y1 = y(1:Nd + 1);
y2 = y(Nd + 2:2*Nd + 2);

dy1 = (1 - 1i*alpha)*y2.*y1 + H*y1*tau_p;
dy2 = pA*tau_p/tau_s - tau_p/tau_s*y2 - tau_p/tau_s*(1 + 2*y2).*abs(y1).*abs(y1);
dy = [dy1;dy2];
t
end
