%% computation of eigencurves utilizing uniform block-diagonalization
% input:
% E0,  E1, E2, M: coefficient matrices of matrix flow
% ka:  k-value for first decomposition
% kb:  k-value for second decomposition
% kC:  k-values for computing eigencurves
% thB: threshold for determining block structure

%% main function
function omB=eigencurves(E0,E1,E2,M,ka,kb,kC,thB)
Ea=matrixFlow(E0,E1,E2,ka);                     % evaluate matrix flow at ka
[Phi,~]=eig(Ea,M);                              % eigenvalue decomposition
Eb=matrixFlow(E0,E1,E2,kb);                     % evaluate matrix flow at kb
[ind,nBl]=decompose(Eb,Phi,thB);                % block decomposition of Eb
omB=solveBlocks(Phi,E0,E1,E2,M,ind,nBl,kC);     % compute eigencurves of blocks
end

%% definition of matrix flow
function E=matrixFlow(E0,E1,E2,k)
E = k^2*E0 - k*E1 + E2;                         % evaluate matrix flow
end

%% decomposition of matrix flow E using eigenvectors Phi with accuracy thB
function [ind,nBl]=decompose(E,Phi,thB)
B = Phi'*E*Phi;                                 % apply transformation
B(abs(B)/norm(B)<thB)=0;                        % neglect small values
[p,~,r,~,~,~] = dmperm(B);                      % permutation
nBl = numel(r)-1;                               % number of blocks
ind=cellfun(@(i)p(r(i):r(i+1)-1),...
    num2cell(1:nBl),'UniformOutput',false);     % store block indices
end

%% compute eigencurves for each block
function omB=solveBlocks(Phi,E0,E1,E2,M,ind,nBl,kC)
omB{nBl} = [];                                  % allocate frequencies
for i = 1:nBl                                   % loop blocks
    Phic=Phi(:,ind{i});                         % eigenvectors of current block
    Mb = Phic'*M*Phic;                          % decomposed mass matrix
    omB{i} = zeros(numel(kC),numel(ind{i}));    % allocate frequencies
    for j=1:numel(kC)                           % loop k-values
        E=matrixFlow(E0,E1,E2,kC(j));           % evaluate matrix flow at k
        Eb = Phic'*E*Phic;                      % current block
        omB{i}(j,:)=sort(sqrt(eig(Eb,Mb)));     % compute eigenvalues
    end
end
end