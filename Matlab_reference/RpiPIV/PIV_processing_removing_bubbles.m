% Clearing the environment, closing all figures, and clearing the command window
clear all;
close all;
clc;

% Loading resolution data
load resolution;

% Setting up resolution and other parameters
res = RpiPIV_res; % unit: meter/pixel
threshold = 40000; % particle threshold for binarization
minBubbleSize = 100; % minimum size of bubbles to be filtered
mm = 800;
nn = 1600;

dt = 3e-3; % Time step: 3 ms 
grid_pixel_size = 20;
[xi, yi] = meshgrid(350:grid_pixel_size:nn-400,200:grid_pixel_size:500);
[m,n] = size(xi);
ue = zeros(m,n); % Initializing u-component of velocity
ve = zeros(m,n); % Initializing v-component of velocity

xphy = [1:n]*res*grid_pixel_size;
yphy = [1:m]*res*grid_pixel_size;
um = zeros(m,n);
vm = zeros(m,n);

% Setting PIV parameters for first pass
pp.method = 0;  
pp.nx = 64;
pp.ny = 64;
pp.Vel_range = [-20 20 -30 30];
pp.sub_pixel_method = 0;

% Setting PIV parameters for second pass
p.method = 0;
p.nx = 24;
p.ny = 24;
p.Vel_range = [-20 20 -30 30];
p.sub_pixel_method = 1;

% Preparing the figure for displaying results
figure(1)
set(gcf,'position',[50 50 400 500]);

% Range of image indices for PIV analysis
iRange = 0:1:1198;

k = 0;
for i = iRange
    k = k+1;
    % Reading consecutive images
    AID = sprintf('G:/PIV_compare/Bubble_in_chain/12012023_Pylon_RaspPi_Cam/2023-12-01_11-45/image%07d.tif',i);
    BID = sprintf('G:/PIV_compare/Bubble_in_chain/12012023_Pylon_RaspPi_Cam/2023-12-01_11-45/image%07d.tif',i+1);

    A = imread(AID);
    B = imread(BID);
    
    % Converting images to double for processing
    A = double(A);
    B = double(B);

    % Binarizing and filtering bubbles from images
    BWa = imbinarize(A,threshold);
    filterA = bwareaopen(BWa, minBubbleSize);
    A(filterA) = NaN;

    BWb = imbinarize(B,threshold);
    filterB = bwareaopen(BWb, minBubbleSize);
    B(filterB) = NaN;

    % PIV analysis: First pass
    imagesc(B-A);
    colormap gray;
    hold on;
    [u,v,c] = mat_piv(A,B,xi,yi,ue,ve,pp);

    % Median filtering and further processing
    umed = medfilt2(u,[3,3],'symmetric');
    vmed = medfilt2(v,[3,3],'symmetric');
    I = isnan(umed); II = I==1; umed(II) = 0;
    I = isnan(vmed); II = I==1; vmed(II) = 0;
    [u,v,c] = mat_piv(A,B,xi,yi,umed,vmed,p);
    c(c>1) = NaN;

    % Post-processing velocity fields
    umed = medfilt2(u,[3,3],'symmetric');
    vmed = medfilt2(v,[3,3],'symmetric');
    flag1 = abs(u-umed)<3;
    flag2 = abs(v-vmed)<3;
    flag = flag1.*flag2.*(c>0.7);
    u1 = u.*flag.*res./dt + umed.*(1-flag).*res./dt;
    v1 = v.*flag.*res./dt + vmed.*(1-flag).*res./dt;
    u = u1;
    v = -v1;
    
    % Displaying and saving velocity vectors
    if nanmean(c(:))>0.6
        vel(k).u = u;
        vel(k).v = v;
        vel(k).c = c;
        quiver(xi, yi, u, v, 6, 'y');
        drawnow;
    end
    hold off;
end

% Preparing for saving results
[m,n] = size(um);
xphy = (1:m)*res*grid_pixel_size;
yphy = (1:n)*res*grid_pixel_size; 

% Saving the velocity field and related parameters
save velocity_field
