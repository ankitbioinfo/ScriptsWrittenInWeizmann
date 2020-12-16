

clear all 

filename='cell_pos_14_MorphologicalWatershedFilter_Out1_1500_0.4.tif';

I.info = imfinfo(filename);
I.size = [I.info(1).Height, I.info(1).Width, length(I.info)];

 if I.info(1).MaxSampleValue == intmax('uint16') && I.info(1).MinSampleValue == intmin('uint16')
        I.data_type = 'uint16';
        I.bits = 16;
    elseif I.info(1).MaxSampleValue == intmax('uint8') && I.info(1).MinSampleValue == intmin('uint8')
        I.data_type = 'uint8';
        I.bits = 8;
    elseif I.info(1).MaxSampleValue == intmax('int16') && I.info(1).MinSampleValue == intmin('int16')
        I.data_type = 'int16';
        I.bits = 16;
    elseif I.info(1).MaxSampleValue == intmax('int8') && I.info(1).MinSampleValue == intmin('int8')
        I.data_type = 'int8';
        I.bits = 8;
    else
        disp('Data type could not be identified - contact Tomer.');
        return;
    end

% I.img = zeros(I.size, I.data_type);
% fprintf('Loading image... ');
% for i = 1 : I.size(3)
%      str = [num2str(i), '/', num2str(I.size(3))];
%       fprintf(str);
%     I.img(:,:,i) = imread(filename, 'Index', i);
%     fprintf(repmat('\b', 1, length(str)));
% end






predicted_cent=load('Dcut10/column_degree30_1.dat');
characteristics = zeros(I.size, I.data_type);


load('c_n_pos14 (Characteristics).mat')
index=1;
for i=1: size(C.coords,1)
    pos=C.coords(i,:);
    centroids = calc_centroids(C.masks(i), C.origins(i,:));
%     for j=1:size(predicted_cent,1)
%         dist=pdist([predicted_cent(j,1:3);centroids]);
%         if dist<0.1
%             matching(index)=i;
%             index=index+1;
            characteristics(pos(1):pos(2),pos(3):pos(4),pos(5):pos(6))= characteristics(pos(1):pos(2),pos(3):pos(4),pos(5):pos(6)) | C.masks{i};
%         end
%     end
    
end

 characteristics=255*characteristics;
    
%matching 

rt=Tiff('cell_pos_14_MorphologicalWatershedFilter_Out1_1500_0.4.tif','r+');
t = Tiff('myColumnPrediction.tif','w');
tagstruct.ImageLength = I.size(1);
tagstruct.ImageWidth =  I.size(2);
tagstruct.SubFileType= Tiff.SubFileType.Page;
tagstruct.SampleFormat = Tiff.SampleFormat.UInt; % uint
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.PackBits;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
for ii=1:I.size(3)
   setTag(t,tagstruct);
   write(t,characteristics(:,:,ii));
   writeDirectory(t);
end
close(t)


% t.setTag('ImageLength',I.size(1));
% t.setTag('ImageWidth', I.size(2));
% t.setTag('Photometric', getTag(rt,'Photometric'));
% t.setTag('BitsPerSample',  getTag(rt,'BitsPerSample'));
% t.setTag('SamplesPerPixel', getTag(rt,'SamplesPerPixel'));
% t.setTag('TileWidth', 128);
% t.setTag('TileLength',128);
% t.setTag('Compression', getTag(rt,'Compression'));
% t.setTag('PlanarConfiguration', getTag(rt,'PlanarConfiguration'));
% t.setTag('Software', 'MATLAB');
% t.write(characteristics);
% t.close();


while 0
    for i = 1:I.size(3)
%         subplot(1,2,1)
%         imagesc(I.img(:,:,i));
%         title(strcat('my id', num2str(i)))
        
        subplot(1,2,2)
        imagesc(characteristics(:,:,i));
        title(strcat('my id', num2str(i)))
        
        pause(0.15)
    end
end


    

function centroids = calc_centroids(masks, origins)
centroids = nan(length(masks),3);
for i = 1 : length(masks)
    [x,y,z] = ind2sub(size(masks{i}), find(masks{i}));
    centroids(i,:) = (mean([x,y,z],1) + double(origins(i,:) - 1));
end
end


