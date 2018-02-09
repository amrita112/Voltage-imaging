clear all;
%%

    dr = uigetdir;
    D = dir([dr filesep '*.tif'])';
    [fnSave, drSave] = uiputfile([dr filesep '..' filesep 'dataset.mat']);
    
    fsize = [D.bytes];
    select = abs(fsize-median(fsize))<0.1*median(fsize);
    D = D(select);
    fNames = sort_nat({D.name});
    clear D
    
    %make a preview file
    fnums = round(linspace(1, length(fNames), 100)); 
    for f_ix = length(fnums):-1:1 
        A = tiffread2([dr filesep fNames{fnums(f_ix)}]);
        preview(:,:,f_ix) =  A.data;
    end
    figure('name', 'doubleclick when done selecting subregion to load');imagesc(mean(preview,3));
    %     ROI = imrect;
    %     pos = wait(ROI);
    %     boundX = max(1,round(pos(1))):min(size(preview,2),round(pos(1)+pos(3)));
    %     boundY = max(1,round(pos(2))):min(size(preview,1),round(pos(2)+pos(4)));
    
    x = inputdlg('Enter the video frame rate in Hz (e.g. 500)',...
        'Sample Freq', 1, {'500'});
    sampleRate = str2double(x{:});
    
    boundX = 1:size(preview,2);
    boundY = 1:size(preview,1);

    blocklength = min(length(fNames), ceil(length(fNames)/ ceil(length(fNames)/40000)));
    blockstarts = 1:blocklength:length(fNames);    
    for blockN = 1:length(blockstarts)
        data = zeros([length(boundY) length(boundX) blocklength], 'single');
        for fnum = blockstarts(blockN):min(length(fNames), blockstarts(blockN)+blocklength-1)
            if ~mod(fnum,100)
                disp(['reading file:' int2str(fnum)])
            end
            A = tiffread2([dr filesep fNames{fnum}]);
            data(:,:,fnum-blockstarts(blockN)+1) = single(A.data(boundY,boundX));
        end
        
        %save data
        savefast([drSave fnSave(1:end-4) 'block' int2str(blockN) '.mat'], 'data', 'sampleRate')
    end

