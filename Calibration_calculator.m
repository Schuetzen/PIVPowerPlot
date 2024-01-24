rootDir = 'G:/PIV_compare/Bubble_in_chain/12012023_Pylon_RaspPi_Cam/';
im = imread('G:/PIV_compare/Bubble_in_chain/12012023_Pylon_RaspPi_Cam/Calibration.bmp');
imshow(im)
[x, y] = ginput(2); % Graphical input from the user
pixelDistance = sqrt((x(2) - x(1))^2 + (y(2) - y(1))^2);

% Step 5: Calculate Pixel Length
realWorldDistance = 0.05; % Replace with the actual distance in your chosen units 5 cm
res = realWorldDistance / pixelDistance;

fprintf('Each pixel represents %f real-world units\n', res);

% Saving the pixel length to a file named 'resolution.mat'
save(fullfile(rootDir,'resolution.mat'), 'res');