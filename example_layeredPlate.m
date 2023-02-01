%% example: layered plate
% This is the example in section 3.5.2 in
% Gravenkamp, H., Plestenjak, B., & Kiefer, D. A. (2023).

%% input

load('matrices_layeredPlate.mat');

ka = 1;                                                                     % first wavenumber for testing decomposability
kb = 2;                                                                     % second wavenumber for testing decomposability
kC = linspace(0,4,200);                                                     % wavenumbers for computing dispersion curves
th = 6;                                                                     % number of significant digits for testing block-structure

omB=eigencurves_withRepeatedEV(E0,E1,E2,M,ka,kb,kC,th);                     % call routines for blockdiagonalization and computing eigencurves

%% plot dispersion curves
% See Figure 7 in the paper
unitBlocks = (cellfun(@(x)size(x,2),omB)==1);
SHblocks   = find(unitBlocks);
LambBlocks = find(~unitBlocks);
figure
set(gcf,'defaulttextinterpreter','latex')
subplot(2,2,1)
hold all
for i=1:numel(SHblocks)
    hS=plot(kC,real(omB{SHblocks(i)}),'Linewidth',1,'Color',[0.004 0.23 0.4],'DisplayName','SH modes');
end
legend(hS(1),'Location','southeast','FontSize',12)

for i=1:numel(LambBlocks)
    subplot(2,2,i+1)
    hL=plot(kC,real(omB{LambBlocks(i)}),'Linewidth',1,'Color',[0.004 0.23 0.4],'DisplayName',['block',num2str(i)]);
    legend(hL(1),'Location','southeast','FontSize',12)
end
for i=1:4
    subplot(2,2,i)
    xlabel('$k$','FontSize',14)
    ylabel('$\omega$','FontSize',14)
    ylim([0 2*pi])
end

figure
set(gcf,'defaulttextinterpreter','latex')
hold all
for i=1:numel(omB)
    hS=plot(kC,real(omB{i}),'Linewidth',1,'Color',[0.004 0.23 0.4],'DisplayName','all blocks');
end
xlabel('$k$','FontSize',14)
ylabel('$\omega$','FontSize',14)
ylim([0 2*pi])

