%% Minimal example 2
% This is the second example in section 3.2 in
% Gravenkamp, H., Plestenjak, B., & Kiefer, D. A. (2023). 
%
% This example uses the minimalistic code 'eigencurves' which is the same
% as printed in the paper. 

%% input
M = [...                                                                    % coefficient matrices of matrix flow
    2, 0, 1, 0; ...
    0, 2, 0, 1; ...
    1, 0, 2, 0; ...
    0, 1, 0, 2];
E0 = [...
    2,   0, 1,   0; ...
    0, 2/3, 0, 1/3; ...
    1,   0, 2,   0; ...
    0, 1/3, 0, 2/3];
E1 = diag([-1,-1,1,1]); E1=1i*E1(end:-1:1,:);

E2 = 1/2*[...
     1,  0, -1,  0; ...
     0,  3,  0, -3; ...
    -1,  0,  1,  0; ...
     0, -3,  0,  3];

ka = 1;                                                                     % first wavenumber for testing decomposability
kb = 2;                                                                     % second wavenumber for testing decomposability
kC = linspace(0,5,100);                                                     % wavenumbers for computing dispersion curves
thB = 1e-6;                                                                  % threshold for determining block structure

omB=eigencurves(E0,E1,E2,M,ka,kb,kC,thB);                                    % call routines for blockdiagonalization and computing eigencurves

%% plot dispersion curves
% See Figure 2(b) in the paper
figure
set(gcf,'defaulttextinterpreter','latex')
h1=plot(kC,real(omB{1}),'Linewidth',2,'Color',[0.004 0.23 0.4],'DisplayName','block 1');
hold on
h2=plot(kC,real(omB{2}),'Linewidth',2,'Color',[0.38 0.65 0.75],'DisplayName','block 2');
xlabel('$k$','FontSize',14)
ylabel('$\omega$','FontSize',14)
legend([h1(1),h2(1)],'Location','southeast','FontSize',12,'Interpreter','latex')
ylim([0 3.5])
