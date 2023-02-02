%% omB=eigencurves_withRepeatedEV(E0,E1,E2,M,ka,kb,kC,th,nAttempt)
% input:
% E0, E1, E2, M:            coefficient matrices of matrix flow
% ka (default 1):           k-value for first  decomposition
% kb (default 2):           k-value for second decomposition
% kC (default 0:0.1:10):    k-values for computing eigencurves
% thB (default 1e-6):       threshold for determining block structure
% thR (default 1e-6):       threshold for determining whether two eigenvalues are the same
% nAttempt (default 5):     number of attempts to decompose in case of repeated EVs

%% main function
function omB=eigencurves_withRepeatedEV(E0,E1,E2,M,varargin)

[ka,kb,kC,thB,thR,flagInv,nAttempt] = checkInput(varargin{:});              % check input for inconsistencies
[E0,E1,E2] = generalizedToStandardEVP(E0,E1,E2,M);                          % transform to standard eigenvalue problem

Ea = matrixFlow(E0,E1,E2,ka);                                               % evaluate matrix flow at ka
Eb = matrixFlow(E0,E1,E2,kb);                                               % evaluate matrix flow at kb
[Phi,Wa] = eig(Ea,'vector');                                                % eigenvalue decomposition
[~,rFlag] = findRepeatedEV(Wa,thR);                                         % check for repeated eigenvalues

[ind,nBl] = decompose(Eb,Phi,thB);                                          % block decomposition of Eb
omB = solveBlocks(Phi,E0,E1,E2,ind,nBl,kC,rFlag,thB,thR,flagInv,nAttempt);  % compute eigencurves of blocks

end


%% any conditions we want to require for the input
function [ka,kb,kC,thB,thR,flagInverse,nAttempt]=checkInput(varargin)

nVar = numel(varargin);                                                     % number of optional inputs

% set default values
if (nVar>0) && ~isempty(varargin{1})
    ka = varargin{1};
else
    ka = 1;
end
if (nVar>1) && ~isempty(varargin{2})
    kb = varargin{2};
else
    kb = 2;
end
if (nVar>2) && ~isempty(varargin{3})
    kC = varargin{3};
else
    kC = 0:0.1:10;
end
if (nVar>3) && ~isempty(varargin{4})
    thB = varargin{4};
else
    thB = 6;
end
if (nVar>4) && ~isempty(varargin{5})
    thR = varargin{5};
else
    thR = 1e-6;
end
if (nVar>5) && ~isempty(varargin{6})
    nAttempt = varargin{6};
else
    nAttempt = 5;
end

% check if order of k-values should be inverted
if abs(ka-kC(1)) < thR                                                      % if ka is equal to the first requested k-value
    kC = kC(end:-1:1);                                                      % invert order
    flagInverse = true;                                                     % set a flag so we can re-order later
    % otherwise we would apply the decomposition at ka at the first k in
    % kC, finding that the matrix flow can be diagonalized
else
    flagInverse = false;
end

end

%% transformation to standard EVP
function [E0,E1,E2]=generalizedToStandardEVP(E0,E1,E2,M)

iM = M^(-0.5);                                                              % M^(-1/2)

E0 = iM*E0*iM';                                                             % transform E0, E1, E2 accordingly
E1 = iM*E1*iM';
E2 = iM*E2*iM';

E0 = (E0'+E0)/2;                                                            % enforce Hermitian structure
E1 = (E1'+E1)/2;                                                            
E2 = (E2'+E2)/2;                                                            

end

%% definition of matrix flow
function E=matrixFlow(E0,E1,E2,k)
E = k^2*E0 - k*E1 + E2;                                                     % evaluate matrix flow
end

%% decomposition of matrix flow E using eigenvectors Phi with accuracy th
function [ind,nBl]=decompose(E,Phi,thB)
B = Phi'*E*Phi;                                                             % apply transformation
B(abs(B)/norm(B)<thB)=0;                                                    % neglect values below threshold

[p,~,r,~,~,~] = dmperm(B);                                                  % permutation
nBl = numel(r)-1;                                                           % number of blocks
ind = cellfun(@(i)p(r(i):r(i+1)-1),num2cell(1:nBl),'UniformOutput',false);  % store block indices
indEmpty = cellfun(@isempty,ind);                                           % in rare cases, this approach returns empty blocks
ind = ind(~indEmpty);                                                       % remove indices of empty blocks
nBl = sum(~indEmpty);                                                       % correct number of blocks if necessary

end

%% compute eigencurves for each block
% this version considers repeated eigenvalues
function omB=solveBlocks(Phi,E0,E1,E2,ind,nBl,kC,rFlag,thB,thR,flagInv,nAttempt)

omB{nBl} = [];                                                              % allocate frequencies of all blocks
iBl = 1;                                                                    % index of current block

while iBl <= nBl                                                            % loop blocks
    if rFlag                                                                % if there were repeated eigenvalues in the global matrix
        rFlagBlock = true;                                                  % assume there are repeated EVs in current block
    else                                                                    % otherwise no repeated EVs had been found from the beginning
        rFlagBlock = false;                                                 % then don't worry about repeated EVs
    end
    Phic = Phi(:,ind{iBl});                                                 % eigenvectors of current block
    omB{iBl} = nan(numel(kC),numel(ind{iBl}));                              % allocate frequencies
    E0b =  Phic'*E0*Phic;                                                   % current block
    E1b =  Phic'*E1*Phic;                                                   % current block
    E2b =  Phic'*E2*Phic;                                                   % current block
    for j=1:numel(kC)                                                       % loop k-values    
        Eb = matrixFlow(E0b,E1b,E2b,kC(j));                                 % evaluate matrix flow at k
        if rFlagBlock ...                                                   % if the block previously had repeated EVs
                && ~(j==1) ...                                              % and it's not the first wavenumber
                && j<=nAttempt                                              % and we haven't exceeded the allowed number of attempts
            [indSub,nBlSub] = decompose(Eb,PhiB,thB);                       % try to decompose block with previous eigenvectors
            if nBlSub>1                                                     % if block can be further decomposed
                Phi(:,ind{iBl}) = Phic*PhiB;                                % update eigenvectors in global matrix
                Phic = Phic*PhiB(:,indSub{1});                              % update eigenvectors of current block

                for i=1:nBlSub                                              % loop subblocks
                    indSub{i} = ind{iBl}(indSub{i});                        % translate to numbering of global matrix
                end

                ind = [ind,indSub{2:end}];                                  % append indices of subblocks
                ind{iBl} = indSub{1};                                       % update indices of current block
                nBl = nBl+nBlSub-1;                                         % add nBlSub-1 new blocks
                E0b = Phic'*E0*Phic;                                        % current block
                E1b = Phic'*E1*Phic;                                        % current block
                E2b = Phic'*E2*Phic;                                        % current block
                Eb = matrixFlow(E0b,E1b,E2b,kC(j));                         % evaluate matrix flow at k
                omB{iBl} = omB{iBl}(:,1:numel(ind{iBl}));                   % remove obsolete columns after block size is reduced
            end


        end
        [PhiB,Om] = eig(Eb,'vector');                                       % compute eigenvalues and eigenvectors
        if rFlagBlock&&~(j==1) && (nBlSub>1) 
            [~,rFlagBlock] = findRepeatedEV(Om,thR);                        % check if there are still repeated eigenvalues
        end
        omB{iBl}(j,1:numel(Om)) = sort(sqrt(Om));                           % sort and store square root of eigenvalues
    end                                                                     % end k-loop
    if flagInv                                                              % if order of k was reversed
        omB{iBl} = omB{iBl}(end:-1:1,:);                                    % correct order of omega
    end
    iBl = iBl +1;                                                           % go to next block
end                                                                         % end block-loop

end

%% check for repeated eigenvalues
function [rep,rFlag]=findRepeatedEV(w,thR)

[~,IW,IU] = uniquetol(w,thR);                                               % find unique EVs
h = histcounts(IU,numel(IW));                                               % multiplicity of EVs
iR = find(h>1);                                                             % indices of repeated EVs
mulMax = max(h);                                                            % maximum multiplicity
if mulMax>1                                                                 % if there are repeated EVs
    rFlag = true;                                                           % flag for repeated EVs
    rep = nan(numel(iR),mulMax);                                            % allocate array of repeated indices
    for i=1:numel(iR)                                                       % loop repeated EVs
        dupI = find(IU==iR(i));                                             % find indices of repeated EVs
        rep(i,1:numel(dupI)) = dupI;                                        % store in array
    end
else
    rFlag = false;                                                          % flag for no repeated EVs
    rep=[];
end

end
