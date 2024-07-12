% Sandpile model
number=1;
nrsteps = 25000;
% preallocate output  
nrgrainstest100 = zeros(nrsteps,number);   % nr grains in the sandpile
avsizetest100 = zeros(nrsteps,number);   % avalanche size
nrouttest100 = zeros(nrsteps,number);   % grains falling from the grid
for b=1:number;
rand('seed',sum(100*clock())) ;
siz = [64 64];
 
if numel(siz)~=2;
    error('siz must be a two-element vector')
end
 
 
% if siz(1)*siz(2) > 900;
%     warning('Caution! I suggest to take a smaller grid (e.g. siz = [20 20]).')
% end
 
nrsteps = ceil(nrsteps);
if ~isscalar(nrsteps) || any(nrsteps<1);
    error('nrsteps must be a positive integer')
end
 
% Ceil siz in case non-integers are provided
siz = ceil(siz);
 
% Nr of cells
nrc = siz(1)*siz(2);
 
% Grid containing sand
Z   = zeros(siz);
 
% critical height
Zcr = 4;                     
 
 
 
% define adjacency
ixc             = reshape((1:nrc)',siz);
ixnu            = nan(siz);
ixnu(2:end,:)   = ixc(1:end-1,:); %upper neighbor;
ixnb            = nan(siz);
ixnb(1:end-1,:) = ixc(2:end,:); %lower neighbor;
ixnl            = nan(siz);
ixnl(:,2:end,:) = ixc(:,1:end-1); %right neighbor;
ixnr            = nan(siz);
ixnr(:,1:end-1) = ixc(:,2:end); %left neighbor;
 
% add one more cell that captures grains falling from the grid
ixg             = repmat(ixc(:),4,1);
ixr             = [ixnu(:); ixnb(:); ixnl(:); ixnr(:)];
ixr(isnan(ixr)) = nrc+1;
 
% create sparse distribution matrix (M)
M = sparse(ixg,ixr,ones(length(ixg),1)*0.25,nrc+1,nrc+1);
 
% add one more cell also to Z
Z = [Z(:); 0]; 
 
% waitbar
h = waitbar(0,'Grains fall...');
 
for r1 = 1:nrsteps;
    waitbar(r1/nrsteps)
    % Put a grain at a random location of the grid
    ixdrop = ceil(nrc*rand(1));
    Z(ixdrop) = Z(ixdrop)+1;
    
    % Check if any location exceeds the critical height Zcr
    I = Z>=Zcr;
    
    % If any, distribute
    while any(I);
        Z = (Z-4*I) + M'*(4*I);
        nrouttest100(r1,b) = nrouttest100(r1,b)+Z(end);
        Z(end) = 0;
        avsizetest100(r1,b) = avsizetest100(r1,b)+sum(I);
        I = Z>=Zcr;
 
    end    
   
    nrgrainstest100(r1,b) = sum(Z);
%     G(:,r1)=Z;
    Ratio3grainstest100(r1,b)=(sum(Z(1:end)==3))/nrc;   %ratio of 3-grains%
%     AvsizeNrout25A10=avsize64+nrout64;
 
end
end    
close(h)

