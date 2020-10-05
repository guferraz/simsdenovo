function [fD, T, D] = calculateDistances(mz,ppm)
% calculate pair-wise distances and assign possible losses (ppm based)

% load residues
R = importdata('aminoResidues1.txt');
aminoLabels(:) = R.textdata;
aminoLabels{20} = 'x';
aminoResidues = R.data;

D = zeros(length(mz));
T = D;
%clear fD T 
fD{length(mz),length(mz)} = '';

k = 1;
for i=1:length(mz)
    for j=i:length(mz)        
        ppmD = 10^6*(abs(abs(mz(j) - mz(i))-aminoResidues)./aminoResidues);
        D(i,j) = min(ppmD); % to be used as weights
        w = find(ppmD == min(ppmD),1);
        %f = find(ppmD < ppm*mz(i)/100,1); %ppm WEIGHED BY MASS%_____
        f = find(ppmD < ppm); %ppm
        if ~isempty(f)
            %D(i,j) = ppmD(f); % to be used as weights
            fD{i,j} = aminoLabels{f};
            T(i,j) = f;            
            k = k + 1;
        end
    end
    disp(i)
end