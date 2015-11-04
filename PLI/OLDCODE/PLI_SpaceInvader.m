% PLI Space Invader Plot


g = defgrids();

g(2).layout


colMap = jet(128);

%% Delta Layout

layout = g(1).layout;

% realPLI = PLI(:,:,1);

chans = 1:16;

% Define lowerleft corner for each channel

[~,chanLLIdx] = intersect(layout,chans);

[ii,jj] = ind2sub([4,4], chanLLIdx);



