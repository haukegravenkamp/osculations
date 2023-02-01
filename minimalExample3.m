%% Minimal example 3
% This is the example in section 3.3 in
% Gravenkamp, H., Plestenjak, B., & Kiefer, D. A. (2023). 
%
% In this example, we choose ka at the crossing point where a
% block-diagonalization is not possible. In those cases, the function 
% eigencurves_withRepeatedEV keeps testing for decomposability at every new
% k-value until all repeated values are removed or all values in kC have
% been tested.

%% input
M  = eye(2);                                                                % coefficient matrices of matrix flow
E0 = [3 2; 2 3]/4;
E1 = zeros(2);
E2 = [1 -1; -1 1]/2;

ka = 1;                                                                     % first wavenumber for testing decomposability
kb = 2;                                                                     % second wavenumber for testing decomposability
kC = linspace(0,4,100);                                                     % wavenumbers for computing dispersion curves
th = 6;                                                                     % number of significant digits for testing block-structure

omB=eigencurves_withRepeatedEV(E0,E1,E2,M,ka,kb,kC,th);                     % call routines for blockdiagonalization and computing eigencurves

%% plot dispersion curves
% See Figure 3 in the paper
figure
set(gcf,'defaulttextinterpreter','latex')
plot(kC,omB{1},'Linewidth',2,'Color',[0.004 0.23 0.4],'DisplayName','$\omega_1$')
hold on
plot(kC,omB{2},'Linewidth',2,'Color',[0.38 0.65 0.75],'DisplayName','$\omega_2$')
xlabel('$k$','FontSize',14)
ylabel('$\omega$','FontSize',14)
legend('Location','southeast','FontSize',12)
ylim([0 4.5])