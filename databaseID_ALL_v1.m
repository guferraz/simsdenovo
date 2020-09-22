function [hits,results,seqhits,uniquehits] = databaseID_ALL_v1(database,ordseq,taxid,minlength,parent)

% load database
taxonomy = database.taxonomy;
refString = database.string;
protein = database.protein;

% **TEST** load trio frequency (ONLY 9913 for now) **TEST**
% load('singles_9913.mat','tabletters')
% patt1  = tabletters{:,1};
% freq1  = tabletters{:,2};
% freq1  = freq1/max(freq1);
% 
% load('pairs_9913.mat','tabpairs')
% patt2  = tabpairs{:,1};
% freq2  = tabpairs{:,2};
% freq2  = freq2/max(freq2);
% 
% load('trios_9606.mat','tabtrios')
% patt3  = tabtrios{:,1};
% freq3  = tabtrios{:,2};
% freq3  = 100*freq3/sum(freq3);
% 
% 
% freq = [freq1 ; freq2 ; freq3];
% patt = [patt1 ; patt2 ; patt3];
% 
% freq = freq3;
% patt = patt3;

Lweight = [0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

type = find(taxonomy == taxid);
refString = refString(type);
protein = protein(type);

% load residues
R = importdata('aminoResidues1.txt');
aminoLabels(:) = R.textdata;
aminoLabels{20} = 'x';
aminoResidues = R.data;

shortrefString = refString;

for i=1:length(refString)
    shortrefString{i} = shortrefString{i}([1:min(25,end) max(1,end-25):end]);
    %shortrefString{i} = shortrefString{i}([1:min(15,end)]);
    %hits(i) = hits(i) + 1*contains(refString{i}(max(1,end-15):end),ordseq{s});
end

%f = figure('color','w','Position',[150 400 300 200]);

% id from database
hits = zeros(1,length(refString));
seqhits = zeros(1,length(ordseq));
uniquehits = seqhits;
for s = 1:length(ordseq)
    
    
    disp([s length(ordseq)])
    
    if length(ordseq{s}) >= minlength && ~isempty(shortrefString) %&& testCandidate(ordseq{s},refString,aminoLabels,aminoResidues)
    %fac = 1/sum(contains(refString,ordseq{s}));
%     clear pattfac
%     for i=1:length(patt)
%         %pattfac(i) = contains(ordseq{s},patt{i})/(1+freq(i));
%         pattfac(i) = contains(ordseq{s},patt{i})*freq(i);
%     end
% 
%     fac = (sum(pattfac));
%     fac = 1./fac;
    fac = 1;
    refString = database.string;
    refString = refString(type);
    uniquehits(s) = sum(contains(refString,ordseq{s}));
    
    for i=1:length(shortrefString)
               
        hits(i) = hits(i) + Lweight(length(ordseq{s})).*fac*contains(shortrefString{i},ordseq{s});
        %hits(i) = hits(i) + Lweight(length(ordseq{s}))*contains(refString{i},ordseq{s}) + fac;
        seqhits(s) = seqhits(s) + contains(shortrefString{i},ordseq{s});
    end    
  
     if sum(hits==0)>0
         refString(hits==0) = {''};
         %hits = hits(hits>0);    
     end
    
     for k=1:length(ordseq)
         
%          clear pattfac
%          for i=1:length(patt)
%             %pattfac(i) = contains(ordseq{s},patt{i})/(1+freq(i));
%             pattfac(i) = contains(ordseq{s},patt{i})*freq(i);
%          end
% 
%          fac = (sum(pattfac));
%          fac = 1./fac; 
         fac = 1;
         
         if length(ordseq{k}) >= minlength && ~isempty(refString) %&& testCandidate(ordseq{k},refString,aminoLabels,aminoResidues)
             sj = find(hits>0);
             for j=sj             
                %hits(i) = hits(j) + length(ordseq{k})^1*contains(refString{j},ordseq{k});
                %disp([s k j]) 
                 if contains(refString{j},ordseq{k})
                    hits(j) = hits(j) + Lweight(length(ordseq{k}))*fac;% + 0*contains(refString{j},ordseq{k});
                    seqhits(k) = seqhits(k) + 1;
                 end
             end
         end
     end
    
       
    end 
    %figure(parent)
    doplot = 0;
    if doplot
        cla(parent), plot(parent,hits), hold(parent,'on')
        leads = find(hits >= 1*max(hits)/4);
        if length(leads) < 50
            text(parent,leads,hits(leads),protein(leads),'Rotation',45,'Interpreter','none')        
        end
        %ylim(parent,[0 1.5*max(hits)])
        top = find(seqhits==max(seqhits),1);
        if isempty(top)
            topseq = 'none';
        else
            topseq = ordseq{top};
        end
        title(parent,topseq)
        drawnow
    end
end
%
%hits = hits/sum(hits);
% top = find(seqhits==max(seqhits),1);
% if isempty(top)
%     topseq = 'none';
% else
%     topseq = ordseq{top};
% end
topseq = 'none';
results = protein(hits>=max(hits)/2);
[hits, ordhits]  = sort(hits,'descend'); 
results = protein(ordhits);
%results = protein(hits>1);
