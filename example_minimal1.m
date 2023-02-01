%% Minimal example 1
% This is the first example in section 3.2 in
% Gravenkamp, H., Plestenjak, B., & Kiefer, D. A. (2023). 
%
% This example uses the minimalistic code 'eigencurves' which is the same
% as printed in the paper. 

%% input
M  = [2 1; 1 2];                                                            % coefficient matrices of matrix flow
E0 = M/3;
E1 = zeros(2);
E2 = 3/2*[1 -1; -1 1];

ka = 1;                                                                     % first wavenumber for testing decomposability
kb = 2;                                                                     % second wavenumber for testing decomposability
kC = linspace(0,5,100);                                                     % wavenumbers for computing dispersion curves
th = 6;                                                                     % number of significant digits for testing block-structure

omB=eigencurves(E0,E1,E2,M,ka,kb,kC,th);                                    % call routines for blockdiagonalization and computing eigencurves

%% plot dispersion curves
% See Figure 2(a) in the paper
figure
set(gcf,'defaulttextinterpreter','latex')
plot(kC,real(omB{1}),'Linewidth',2,'Color',[0.004 0.23 0.4],'DisplayName','$\omega_1$')
hold on
plot(kC,real(omB{2}),'Linewidth',2,'Color',[0.38 0.65 0.75],'DisplayName','$\omega_2$')
xlabel('$k$','FontSize',14)
ylabel('$\omega$','FontSize',14)
legend('Location','southeast','FontSize',12)
ylim([0 4])
