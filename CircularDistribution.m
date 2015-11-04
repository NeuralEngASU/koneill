figure (2)
clf
subplot(2,2,1)
theta = linspace(-pi, pi, 1000);

rho = ones(1,length(theta));

t = 0 : .01 : 2 * pi;
P = polar(t, 2 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
p = polar(theta,rho);
title(sprintf('No Coupling, Flat Phase Distribution of Phase\n PLI = 0. R = 0'))
%%
T = findall(gca,'type','text');
%%
delete(T([1:2, 4:end]))
%%
h = findall(gca,'type','line');

% remove the handle for the polar plot line from the array
h(h == p) = [];
% delete all other lines
%%
delete(h([3:4, 6:end]));
%%
subplot(2,2,2)
theta = linspace(-pi, pi, 1000);

rho = (cos(theta)+1)/2 + 0.5;

t = 0 : .01 : 2 * pi;
P = polar(t, 2 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
p = polar(theta,rho);
title(sprintf('Coupling around 0 mod pi\n PLI = 0. 0 < R < 1'))
%%
T = findall(gca,'type','text');
%%
delete(T([1:2, 4:end]))
%%
h = findall(gca,'type','line');

% remove the handle for the polar plot line from the array
h(h == p) = [];
% delete all other lines
%%
delete(h([3:4, 6:end]));

%%
subplot(2,2,3)
theta = linspace(-pi, pi, 1000);

rho = (cos(theta+(pi/4))+1)/2 + 0.5;

t = 0 : .01 : 2 * pi;
P = polar(t, 2 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
p = polar(theta,rho);
title(sprintf('Coupling around pi/4 away from 0 mod pi\n 0 < PLI < 1. 0 < R < 1'))
%%
T = findall(gca,'type','text');
%%
delete(T([1:2, 4:end]))
%%
h = findall(gca,'type','line');

% remove the handle for the polar plot line from the array
h(h == p) = [];
% delete all other lines
%%
delete(h([3:4, 6:end]));


%%
subplot(2,2,4)
theta = linspace(-pi, pi, 1000);

rho = 2*(cos(theta+(pi/2)));

t = 0 : .01 : 2 * pi;
P = polar(t, 2 * ones(size(t)));
set(P, 'Visible', 'off')
hold on
p = polar(theta,rho);
title(sprintf('Coupling around pi/2 away from 0 mod pi\n PLI ~ 1. 0 << R < 1'))
%%
T = findall(gca,'type','text');
%%
delete(T([1:2, 4:end]))
%%
h = findall(gca,'type','line');

% remove the handle for the polar plot line from the array
h(h == p) = [];
% delete all other lines
%%
delete(h([3:4, 6:end]));