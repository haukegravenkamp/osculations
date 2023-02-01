%% input:
% E0, E1, E2, M: coefficient matrices of matrix flow
% ka: k-value for first  decomposition
% kb: k-value for second decomposition
% kC: k-values for computing eigencurves
% th: threshold for determining block structure (significant digits)

%% main function
function omB=eigencurves_withRepeatedEV(E0,E1,E2,M,ka,kb,kC,th)

[E0,E1,E2,M,ka,kb,kC,th,flagInverse] = checkInput(E0,E1,E2,M,ka,kb,kC,th);  % check input for inconsistencies
[E0,E1,E2] = generalizedToStandardEVP(E0,E1,E2,M);                          % transform to standard eigenvalue problem

Ea = matrixFlow(E0,E1,E2,ka);                                               % evaluate matrix flow at ka
Eb = matrixFlow(E0,E1,E2,kb);                                               % evaluate matrix flow at kb
[Phi,Wa] = eig(Ea,'vector');                                                % eigenvalue decomposition
[~,rFlag] = findRepeatedEV(Wa,th);                                          % check for repeated eigenvalues

[ind,nBl] = decompose(Eb,Phi,th);                                           % block decomposition of Eb
omB = solveBlocks(Phi,E0,E1,E2,ind,nBl,kC,rFlag,th,flagInverse);            % compute eigencurves of blocks

end


%% any conditions we want to require for the input
function [E0,E1,E2,M,ka,kb,kC,th,flagInverse]=checkInput(E0,E1,E2,M,ka,kb,kC,th)

if isempty(ka)                                                              % k-value for first decomposition  not given
    ka = kC(end);                                                           % use last k-value
end
if isempty(kb)                                                              % k-value for second decomposition  not given
    ka = mean(kC);                                                          % use average k-value
end
if (nargin<8)||isempty(th)                                                  % threshold not given
    th = 6;                                                                 % default six significant digits
end
if abs(ka-kC(1)) < 10^(-th)                                                 % if ka is equal to the first requested k-value
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
E0=1/2*(E0'+E0);        % make hermitian
E1=1/2*(E1'+E1);        % make hermitian
E2=1/2*(E2'+E2);        % make hermitian



end

%% definition of matrix flow
function E=matrixFlow(E0,E1,E2,k)
E = k^2*E0 - k*E1 + E2;                                                     % evaluate matrix flow
end

%% decomposition of matrix flow E using eigenvectors Phi with accuracy th
function [ind,nBl]=decompose(E,Phi,th)

B = round(Phi'*E*Phi,th);                                                   % decompose E
[p,~,r,~,~,~] = dmperm(B);                                                  % permutation
nBl = numel(r)-1;                                                           % number of blocks
ind = cellfun(@(i)p(r(i):r(i+1)-1),num2cell(1:nBl),'UniformOutput',false);  % store block indices
indEmpty = cellfun(@isempty,ind);                                           % in rare cases, this approach returns empty blocks
ind = ind(~indEmpty);                                                       % remove indices of empty blocks
nBl = sum(~indEmpty);                                                       % correct number of blocks if necessary

end

%% compute eigencurves for each block
% this version considers repeated eigenvalues
function omB=solveBlocks(Phi,E0,E1,E2,ind,nBl,kC,rFlag,th,flagInverse)

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
    for j=1:numel(kC)                                                       % loop k-values
        E = matrixFlow(E0,E1,E2,kC(j));                                     % evaluate matrix flow at k
        Eb =  Phic'*E*Phic;                                                 % current block
        if rFlagBlock&&~(j==1)                                                       % if the block previously had repeated EVs
            [indSub,nBlSub] = decompose(Eb,PhiB,th);                        % try to decompose block with previous eigenvectors
            if nBlSub>1                                                     % if block can be further decomposed
                Phi(:,ind{iBl}) = Phic*PhiB;                                % update eigenvectors in global matrix
                Phic = Phic*PhiB(:,indSub{1});                              % update eigenvectors of current block

                for i=1:nBlSub                                              % loop subblocks
                    indSub{i} = ind{iBl}(indSub{i});                        % translate to numbering of global matrix
                end

                ind = [ind,indSub{2:end}];                                  % append indices of subblocks
                ind{iBl} = indSub{1};                                       % update indices of current block
                nBl = nBl+nBlSub-1;                                         % add nBlSub-1 new blocks

                Eb =  Phic'*E*Phic;                                         % update matrix flow of current block
                %                 Om = eig(Eb,'vector');                                      % compute eigenvalues of first subblock
                omB{iBl} = omB{iBl}(:,1:numel(ind{iBl}));                   % remove obsolete columns after block size is reduced
            end


        end
        [PhiB,Om] = eig(Eb,'vector');                                       % compute eigenvalues and eigenvectors
        if rFlagBlock&&~(j==1) && (nBlSub>1) 
            [~,rFlagBlock] = findRepeatedEV(Om,th);                     % check if there are still repeated eigenvalues
        end
        omB{iBl}(j,1:numel(Om)) = sort(sqrt(Om));                           % sort and store square root of eigenvalues
    end                                                                     % end k-loop
    if flagInverse                                                          % if order of k was reversed
        omB{iBl} = omB{iBl}(end:-1:1,:);                                    % correct order of omega
    end
    iBl = iBl +1;                                                           % go to next block
end                                                                         % end block-loop

end

%% check for repeated eigenvalues
function [rep,rFlag]=findRepeatedEV(w,th)

[~,IW,IU] = (unique(round(w,th)));                                          % find unique EVs
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
