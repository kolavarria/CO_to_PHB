clearvars
clc

load CO_to_PHB_model

% Read TIFF image
tiffImage = imread('metabolic_network_PK.tif');
% Resize the image for better quality
scaleFactor = 0.5; % You can adjust this scale factor as needed
resizedImage = imresize(tiffImage, scaleFactor);
% Display the image with improved quality
figure('Name', 'TIFF Image', 'Position', [100, 100, size(resizedImage, 2), size(resizedImage, 1)]);
imshow(resizedImage);
%loading the .mat file containing the list of reactions and the coordinates
load reactions_and_positions
rxnID = findRxnIDs(model,reactions);

    for i=1:length(rxnID)    
        fluxes(i)=FBAsolution.x(rxnID(i));
        labels{i} = num2str(sprintf('%g',round(fluxes(i)*100)/100));
        text(axis_x(i),axis_y(i),labels(i),'FontSize',8,'Color','red','FontWeight','bold');
    end
    
% Save the resulting figure as a TIFF file with a resolution of 300 dpi
print('flux_distribution', '-dtiff', '-r300');