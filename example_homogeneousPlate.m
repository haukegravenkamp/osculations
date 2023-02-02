%% example: homogeneous isotropic layer 
% This is the example in section 3.5.1 in
% Gravenkamp, H., Plestenjak, B., & Kiefer, D. A. (2023). 

%% input

load('matrices_homogeneousPlate.mat');

ka = 1;                                                                     % first wavenumber for testing decomposability
kb = 2;                                                                     % second wavenumber for testing decomposability
kC = linspace(0,40,200);                                                    % wavenumbers for computing dispersion curves
thB = 1e-6;                                                                 % threshold for determining block structure

omB=eigencurves_withRepeatedEV(E0,E1,E2,M,ka,kb,kC,thB);                    % call routines for blockdiagonalization and computing eigencurves

%% plot dispersion curves
% See Figure 5 in the paper
figure
set(gcf,'defaulttextinterpreter','latex')
h1=plot(kC,real(omB{1}),'Linewidth',2,'Color',[0.004 0.23 0.4],'DisplayName','block 2');
hold on
h2=plot(kC,real(omB{2}),'Linewidth',2,'Color',[0.38 0.65 0.75],'DisplayName','block 1');
xlabel('$k$','FontSize',14)
ylabel('$\omega$','FontSize',14)
legend([h1(1),h2(1)],'Location','southeast','FontSize',12,'Interpreter','latex')
ylim([0 35])
