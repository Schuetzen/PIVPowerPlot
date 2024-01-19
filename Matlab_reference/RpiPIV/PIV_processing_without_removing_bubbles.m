clear all
close all
clc

load resolution

res = RpiPIV_res; % unit: meter/pixel
threshold = 20000; % particle 
minBubbleSize = 15;
mm = 800;
nn = 1600;  

dt = 3e-3; % 3 ms 
grid_pixel_size = 20;
[xi, yi] = meshgrid(350:grid_pixel_size:nn-400,200:grid_pixel_size:500);
[m,n] = size(xi);
ue = zeros(m,n);  
ve = zeros(m,n);

xphy = [1:n]*res*grid_pixel_size;
yphy = [1:m]*res*grid_pixel_size;
um = zeros(m,n);
vm = zeros(m,n);


pp.method = 0;  
pp.nx = 64;
pp.ny = 64;
pp.Vel_range = [-20 20 -30 30];
pp.sub_pixel_method = 0;

p.method = 0;
p.nx = 24;
p.ny = 24;
p.Vel_range = [-20 20 -30 30];
p.sub_pixel_method = 1;

figure(1)
set(gcf,'position',[50 50 400 500]);
% figure(2)
% set(gcf,'position',[450 50 400 400]);

iRange = 0:1:1198;

k = 0;
for i =iRange
    k = k+1;
    AID = sprintf('../RpiPIV/2023-12-01_11-32/image%07d.tif',i);
    BID = sprintf('../RpiPIV/2023-12-01_11-32/image%07d.tif',i+1);

    A = imread(AID);
    B = imread(BID);
    
    A = double(A);
    B = double(B);

%     BWa = imbinarize(A,threshold);
%     filterA = bwareaopen(BWa, minBubbleSize);
%     A(filterA) = NaN;
% 
%     BWb = imbinarize(B,threshold);
%     filterB = bwareaopen(BWb, minBubbleSize);
%     B(filterB) = NaN;



%     figure
%     subplot(121)
%     imagesc(A);
%     colormap gray
%     axis equal
%     axis([1200 1500 800 1000])
%     
%     subplot(122)
%     imagesc(BWa);
%     colormap gray
%     axis equal
%     axis([1200 1500 800 1000])

   
    imagesc(B-A);
    colormap gray
    hold on

    [u,v,c] = mat_piv(A,B,xi,yi,ue,ve,pp);

%     quiver(xi,yi,u,-v,6);
    

    umed = medfilt2(u,[3,3],'symmetric');
    vmed = medfilt2(v,[3,3],'symmetric');
    I = isnan(umed); II = I==1; umed(II) = 0;
 
    I = isnan(vmed); II = I==1; vmed(II) = 0;

    [u,v,c] = mat_piv(A,B,xi,yi,umed,vmed,p);
    I = (c>0.7); I = sum(I(:)); validity = I/m/n;
    title([i mean(c(:))]);
    umed = medfilt2(u,[3,3],'symmetric');
    vmed = medfilt2(v,[3,3],'symmetric');
    flag1 = abs(u-umed)<3;
    flag2 = abs(v-vmed)<3;
    flag = flag1.*flag2.*(c>0.7);
    u1 = u.*flag.*res./dt + umed.*(1-flag).*res./dt;
    v1 = v.*flag.*res./dt + vmed.*(1-flag).*res./dt;
%     u = medfilt2(u,[2 2],'symmetric');
%     v = medfilt2(v,[2 2],'symmetric');

    u = u1;
    v = -v1;
    
    
    if mean(c(:))>0.6
        vel(k).u = u;
        vel(k).v = v;
        vel(k).c = c;
        quiver (xi,yi,u,v,6,'y');
        
        drawnow;
    end

    hold off

%     max(v(:).^2)
%     pause
%     figure(2)
%     subplot(221)
%     imagesc(u)
%     title(i)
%     colorbar
%     subplot(222)
%     imagesc(v)
%     title(max(u(:).^2+v(:).^2))
%     colorbar
%     subplot(223)
%     hist(um(:),20);
%     subplot(224)
%     hist(vm(:),20);
end

[m,n] = size(um);
xphy = (1:m)*res*grid_pixel_size;
yphy = (1:n)*res*grid_pixel_size; 

save velocity_field_wo_removing_bubble_2023_12_01_11_32 vel res xphy yphy dt 
