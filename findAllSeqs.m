function [ordseq,ordreg,wT] = findAllSeqs(mz,ppm,aminoLabels,minlength,maxlength,seqsize,nseeds)
    %oi
[~, wT, ~] = calculateDistances(mz,ppm);

[seedi,seedj] = find(wT~=0);

%seedi = flipud(seedi); % flip seed points
%seedj = flipud(seedj);

% skip very long outputs
% if length(seedi) > 200
%     ordseq{1} = '';
%     ordreg{1} = [];
%     return
% end

S = [];
R = [];
if nargin < 7
    nseeds = length(seedi);
end

for m=1:min(nseeds,length(seedi))  
    idx = 1;
    i = seedi(m);
    j = seedj(m);
    seq = repmat(aminoLabels(wT(i,j)),seqsize,1);
    reg = repmat({[i j]},seqsize,1);
    [idx,seq,reg] = findSeqs(i,j,wT,aminoLabels,idx,seq,reg);
    %[nxt,seq,idx] = findNext(i,j,wT,aminoLabels,idx,seq);

    reg(cellfun(@length,seq)<minlength) = [];
    %reg(cellfun(@length,seq)>maxlength) = [];
    seq(cellfun(@length,seq)<minlength) = [];    
    %seq(cellfun(@length,seq)>maxlength) = [];
    
    S = [S; seq];
    R = [R; reg];
    disp([m length(seedi)])
end

%S = unique(S);
if ~isempty(S)
    [~,ord] = sort(cellfun(@length,S));
    ordseq = flipud(S(ord));
    ordreg = flipud(R(ord));    
else
    ordseq = {''};
    ordreg = {[0 0]};
end