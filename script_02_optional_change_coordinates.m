clearvars
clc
% Prompt user to select TIFF image file
[fileName, filePath] = uigetfile('*.tif', 'Select TIFF Image File');

% Check if user canceled file selection
if isequal(fileName, 0)
    disp('User canceled file selection.');
    return;
end

% Build full file path
fullFilePath = fullfile(filePath, fileName);

% Read TIFF image
tiffImage = imread(fullFilePath);

% Resize the image for better quality
scaleFactor = 0.5; % You can adjust this scale factor as needed
resizedImage = imresize(tiffImage, scaleFactor);

% Display the image with improved quality
figure('Name', 'TIFF Image', 'Position', [100, 100, size(resizedImage, 2), size(resizedImage, 1)]);
imshow(resizedImage);
% title('TIFF Image');

clearvars
clc

data_table = readtable("CO_to_PHB.csv");

reactions = data_table.ReactionName;
message = reactions;
message(end+1)={'DOUBLE-CLICK!! in the last selection'};
h = msgbox(message);
set(h,'Position',[50 50 300 400]);%[x_begin y_begin length height]

[axis_x,axis_y] = getpts;
    if axis_x & axis_y
       delete(h);
    end
clear data_table h message
save reactions_and_positions
