% 1. Define a high threshold above which sites are within a compartment
% and a low threshold below which sites are in different compartments
% 2. For each site, prepare a list of high threshold similar sites
% 3. Intersect the low threshold similarities of all sites in the list
% 4. Compare between the intersection and the list.
% If not all the sites in the list appear in the intersection, this means
% that the site being tested is similar to at least two sites that are not
% supposed to be in the same compartment and thus is marked as a
% non-compartment site.


close all;
clear all;

% threshold (max 100%)
th = 45;    % high threshold (within compartment). ideal: 45
tl = 25;    % low threshold (between compartments). ideal: 25
c = 90;     % criterion for built in hierarchical clustering. ideal: 75

% simulation parameters
RIPNUM = 2;
RDPNUM = 0;
intrvl = 300;
cutoff = 300;

% load correlation data
crrdir = sprintf('../main/extras/cross_correlations/data/RIP%d-RDP%d-intrvl%d/matrix', RIPNUM, RDPNUM, intrvl);
crrfilename = sprintf('%s/interval%d-cutoff%d',crrdir, intrvl, cutoff);
R = dlmread(crrfilename,'\t');
NLOGSYNS = length(R);

% cut threshold
Mh = R >= th/100;
Ml = R >= tl/100;

T = ones(NLOGSYNS,1);
for i = 1:NLOGSYNS
    Il = ones(NLOGSYNS,1);   % low threshold intersection
    Lh = find(Mh(:,i));         % high threshold list
    Ls = length(Lh);
    for j = 1:Ls
        Il = Il.*Ml(:,Lh(j));
    end
    Chl = Il.*Mh(:,i);       % compare high low
    if (sum(Chl) < Ls)
        T(i) = 0;
    end
end
sprintf('number of remaining compartment sites %d', sum(T))

% reduced similarity matrix for remaining sites
SCi = find(T);   % sites within compartments
Mr = R(SCi,SCi);
M = 1 - Mr;   % for clustering distance is used: high correlation => short distance
z = find(M < 1e-3);
M(z) = 0;
Y = squareform(M);
Z = linkage(Y,'single'); % nearest neighbour is the most relevant it captures the problem 
TC = cluster(Z,'cutoff',c/100); % compartment clusters
nclust = max(TC);
sprintf('%d clusters formed', nclust)

T(SCi) = TC;

S = ones(NLOGSYNS,1)*-1;
[SC,h] = silhouette([], TC, Y);
S(SCi) = SC;

% save
clstfilename = sprintf('%s/cluster-interval%d-cutoff%d-th%d-tl%d-c%d',crrdir, intrvl, cutoff, th, tl, c);%, i-1);
dlmwrite(clstfilename,[T S],'\t');

% =========================================================================
% demonstration
    
    
figure(1)            
imagesc(R);
colormap(jet);
axis square
colorbar

figure(2)            
imagesc(Ml);
colormap(gray);
axis square
colorbar

figure(3)            
imagesc(Mh);
colormap(gray);
axis square
colorbar

figure(4)            
imagesc(Mr);
colormap(jet);
axis square
colorbar

if 0
i = 179; % a trunk site
Il = ones(NLOGSYNS,1);   % low threshold intersection
Lh = find(Mh(:,i));         % high threshold list
Ls = length(Lh);
for j = 1:Ls
    Il = Il.*Ml(:,Lh(j));
end
Chl = Il.*Mh(:,i);       % compare high low
% if (sum(Chl) < Ls)
%     T(i) = 0;
% end
%j1 = 
exmpfilename = sprintf('%s/example-Lh-interval%d-cutoff%d-th%d-tl%d-c%d-logsyn%d',crrdir, intrvl, cutoff, th, tl, c);
dlmwrite(exmpfilename,Lh,'\t');

exmpfilename = sprintf('%s/example-Ll%d-interval%d-cutoff%d-th%d-tl%d-c%d-logsyn%d',crrdir, j1, intrvl, cutoff, th, tl, c);
dlmwrite(exmpfilename,Ml(:,Lh(j1)),'\t');

exmpfilename = sprintf('%s/example-Ll%d-interval%d-cutoff%d-th%d-tl%d-c%d-logsyn%d',crrdir, j2, intrvl, cutoff, th, tl, c);
dlmwrite(exmpfilename,Ml(:,Lh(j1)),'\t');

end