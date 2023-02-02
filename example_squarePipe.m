%% example: square pipe
% This is the example in section 3.5.2 in
% Gravenkamp, H., Plestenjak, B., & Kiefer, D. A. (2023).

%% input

load('matrices_squarePipe.mat');

ka = 1;                                                                     % first wavenumber for testing decomposability
kb = 2;                                                                     % second wavenumber for testing decomposability
kC = linspace(0,10,200);                                                    % wavenumbers for computing dispersion curves
th = 4;                                                                     % number of significant digits for testing block-structure

omB=eigencurves_withRepeatedEV(E0,E1,E2,M,ka,kb,kC,th);                     % call routines for blockdiagonalization and computing eigencurves

%% plot dispersion curves
% See Figure 10 in the paper
figure
set(gcf,'defaulttextinterpreter','latex')
for i=1:numel(omB)
    subplot(3,2,i)
    h=plot(kC,real(omB{i}),'Linewidth',1,'Color',[0.004 0.23 0.4],'DisplayName',['block',num2str(i)]);
    legend(h(1),'Location','southeast','FontSize',12,'Interpreter','latex')
end

subplot(3,2,6)
hold all
for i=1:numel(omB)
    h=plot(kC,real(omB{i}),'Linewidth',1,'Color',[0.004 0.23 0.4],'DisplayName','all blocks');
    legend(h(1),'Location','southeast','FontSize',12,'Interpreter','latex')
end

for i=1:6
    subplot(3,2,i)
    xlabel('$k$','FontSize',14)
    ylabel('$\omega$','FontSize',14)
    ylim([0 10])
end
box on

