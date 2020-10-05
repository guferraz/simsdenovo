function [mz,massIntensity] = loadPeakListAuto(file,crop)

% load peak list
%datalist = dir('protein_1.TXT');
M  = importdata(file);
%M  = importdata('protein_1.TXT');
if isstruct(M)
    M = M.data;
    mz  = M(:,1);
    massIntensity = M(:,3);
else
    mz  = M(:,1);
    massIntensity = M(:,3);
end
%massIntensity = massIntensity./max(massIntensity);

if crop
    figure, 
    subplot(211), plot(mz,massIntensity)
    text(mz(round(end/2)),max(massIntensity),['max: ' num2str(max(massIntensity)) ' min: ' num2str(min(massIntensity))])
    text(mz(round(end/2)),max(massIntensity),['max: ' num2str(max(massIntensity)) ' mean: ' num2str(mean(massIntensity))])
    text(mz(round(end/2)),max(massIntensity),['max: ' num2str(max(massIntensity)) ' mean: ' num2str(mean(massIntensity))])
    text(mz(massIntensity==max(massIntensity)),mean(massIntensity),['mean: ' num2str(mean(massIntensity))])

    threshold = inputdlg({['Set MAX threshold. (max: ' num2str(ceil(max(massIntensity))) ')'],...
                          ['Set MIN threshold. (min: ' num2str(floor(min(massIntensity))) ')'],...
                          ['Set MIN mass. (min: ' num2str(min(mz)) ')'],...
                          ['Set MAX mass. (max: ' num2str(max(mz)) ')']},...
                          'filter',[1 40], {num2str(round(max(massIntensity))),num2str(0),num2str(min(mz)),num2str(max(mz))});
    [tmax,tmin,mmin,mmax] = threshold{:};
    tmax = str2num(tmax);
    tmin = str2double(tmin);
    mmin = str2double(mmin);
    mmax = str2double(mmax);

    mz(massIntensity>tmax | massIntensity<tmin) = [];
    massIntensity(massIntensity>tmax | massIntensity<tmin) = [];

    massIntensity(mz<mmin | mz>mmax) = [];
    mz(mz<mmin | mz>mmax) = [];
    subplot(212), plot(mz,massIntensity), drawnow
end