% Clear workspace, close all figures, and clear command window
clear all;
close all;
clc;

% Load resolution parameters

% Set parameters
%res = ; % Resolution: meters per pixel
%threshold = 40000; % Particle detection threshold
%minBubbleSize = 3000; % Minimum size of bubbles
%mm = 1000; % Image dimension
%nn = 1600; % Image dimension
dt = 5e-3; % Time step: 3 ms
grid_pixel_size = 20; % Size of grid in pixels

% Define the directory where the folders are located
rootDir = 'G:/PIV_compare/Bubble_in_chain/12012023_Phantom_High_Speed_Cam/'; % Replace with your directory path
folders = dir(fullfile(rootDir, 'Test*')); % Adjust the pattern as needed

resolutionPath = fullfile(rootDir, 'resolution.mat');
load(resolutionPath); % Ensure 'resolution' file exists with 'RpiPIV_res'

% Check folders
if isempty(folders)
    disp('No folders found.');
    return;
end

% Specify the folder name
folderName = 'data';

% Check if the folder exists
if ~exist(folderName, 'dir')
    % Folder does not exist, so create it
    mkdir(folderName);
    disp(['Folder "', folderName, '" was created.']);
else
    % Folder already exists
    disp(['Folder "', folderName, '" already exists.']);
end

% Iterate over each folder
for folderIdx = 1:length(folders)
    folderName = folders(folderIdx).name;
    folderPath = fullfile(rootDir, folderName);

     % Find the first image file in the folder
     imageFiles = dir(fullfile(folderPath, '*.tif')); % Adjust file extension if necessary
     if isempty(imageFiles)
         disp(['No image files found in ', folderName]);
         continue;
     end
     firstImageFile = fullfile(folderPath, imageFiles(1).name);
 
     % Read the dimensions of the first image
     firstImage = imread(firstImageFile);
     [mm, nn] = size(firstImage);
    
    % Count the number of pictures
    files = dir(folderPath);
    files = files(~[files.isdir]);
    numFiles = numel(files);
    numTiff = numFiles - 10;

    % Initialize variables for each folder
    %[xi, yi] = meshgrid(340:grid_pixel_size:1600, 200:grid_pixel_size:100);
    [xi, yi] = meshgrid(600:grid_pixel_size:nn-400,200:grid_pixel_size:mm-700);

    [m, n] = size(xi);
    ue = zeros(m, n);  
    ve = zeros(m, n);

    xphy = [1:n]*res*grid_pixel_size;
    yphy = [1:m]*res*grid_pixel_size;
    um = zeros(m,n);
    vm = zeros(m,n);

    % vel = struct; % Structure to store velocity data

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

    % Image processing range
    iRange = 0:1:numTiff;
    k = 0;

    for i = iRange
        k = k + 1;
        AID = fullfile(folderPath, sprintf('Img%06d.tif', i));
        BID = fullfile(folderPath, sprintf('Img%06d.tif', i + 1));

        % Load and process images
        A = imread(AID);
        B = imread(BID);

        A = double(A);
        B = double(B);

        % Binarizing and filtering bubbles from images
        %BWa = imbinarize(A,threshold);
        %filterA = bwareaopen(BWa, minBubbleSize);
        %A(filterA) = NaN;

        %BWb = imbinarize(B,threshold);
        %filterB = bwareaopen(BWb, minBubbleSize);
        %B(filterB) = NaN;

        % PIV analysis: First pass
        imagesc(B-A);
        colormap gray
        hold on

        [u,v,c] = mat_piv(A, B, xi, yi, ue, ve, pp);

        
        % Median filtering and further processing
        umed = medfilt2(u, [3, 3], 'symmetric');
        vmed = medfilt2(v, [3, 3], 'symmetric');
        I = isnan(umed); II = I==1; umed(II) = 0;
        I = isnan(vmed); II = I==1; vmed(II) = 0;
        [u,v,c] = mat_piv(A,B,xi,yi,umed,vmed,p);
        % I = (c > 0.7); I = sum(I(:)); validity = I / (m * n);
        c(c>1) = NaN;


        % Post-processing velocity fields
        umed = medfilt2(u, [3, 3], 'symmetric');
        vmed = medfilt2(v, [3, 3], 'symmetric');
        flag1 = abs(u - umed) < 3;
        flag2 = abs(v - vmed) < 3;
        flag = flag1 .* flag2 .* (c > 0.7);
        u1 = u .* flag .* res ./ dt + umed .* (1 - flag) .* res ./ dt;
        v1 = v .* flag .* res ./ dt + vmed .* (1 - flag) .* res ./ dt;

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

    % Saving the results in a file named after the folder
    saveFileName = sprintf('data/velocity_field_removing_bubble_%s.mat', strrep(folderName, '_', ' '));
    
    % Preparing for saving results
    [m,n] = size(um);
    xphy = (1:m)*res*grid_pixel_size;
    yphy = (1:n)*res*grid_pixel_size; 

    % Saving the velocity field and related parameters
    save(fullfile(rootDir, saveFileName), 'vel', 'res', 'xphy', 'yphy', 'dt');
end
