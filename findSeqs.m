function [idx,seq,reg] = findSeqs(i,j,wT,aminoLabels,idx0,seq,reg)

idx = idx0;
i = j;
J = find(wT(i,:) ~= 0); % get candidates

if ~isempty(J)
    for n=1:length(J)
        %disp([i j idx])   
        j = J(n);
        seq{idx} = [seq{idx} aminoLabels{wT(i,j)}];
        reg{idx} = [reg{idx} ; i j];
        [idx,seq,reg] = findSeqs(i,j,wT,aminoLabels,idx,seq,reg);
        %if length(seq{idx}) > 5
            idx = idx + 1;  
        %end
    end
end
    


end