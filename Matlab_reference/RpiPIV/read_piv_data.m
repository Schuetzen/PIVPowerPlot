clear all
close all
clc

load velocity_field_all_after_clean.mat

[~,~,N] = size(u);

% reverse direction of v, positive to be up direction
v = -v;

% obtain time-averaged velocity field
um = nanmean(u,3);
vm = nanmean(v,3);

figure
subplot(121)
imagesc(um);
colorbar
subplot(122);
imagesc(vm);
colorbar;

Ur = mean(um,1); % radial velocity
Uz = mean(vm,1); % axial velocity

% find the location of the maximal velocity profile
I = find(Uz == max(Uz));

% convert yphy to r
r = yphy - yphy(I);

figure
h1 = plot(r,mean(vm,1),'o-');
xx=ylabel('$U_z$ (m/s)'); set(xx,'interpreter','latex');
yyaxis right
h2 = plot(r,mean(um,1),'s-');
xx=ylabel('$U_r$ (m/s)'); set(xx,'interpreter','latex');
ylim([-0.01 0.01])
h = legend('$U_z$','$U_r$');
set(h,'interpreter','latex');
xx = xlabel('$r$ (m)'); set(xx,'interpreter','latex');
xlim([-0.04 0.04])
set(gca,'fontsize',15)


%% calculate turbulence parameters



% save RpiPIV data
rpiv.r = r;
rpiv.Ur = Ur;
rpiv.Uz = Uz;

