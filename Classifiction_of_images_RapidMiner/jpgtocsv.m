val = 'E:\Lectures\2ndSem\IntroToDS\proj4\test' ;
dirPath = [val,'\Class'] ;
outDirPath = [val,'\csvfiles\' ] ;
nClasses = 6 ;
for i = 1 : nClasses
    dirPath1 = [dirPath,' ', num2str(i), '\'] ;
    outDirPath1 = [outDirPath, num2str(i)] ;
    mkdir(outDirPath) ;
    F = dir(dirPath1) ;
    imgReshaped = [];
    for j =  1 : length(F)-2
    	imgReshapedInt = [];
        imgName = F(j+2).name ;
        fullImgName = strcat(dirPath1, imgName) ;
        img = imread(fullImgName) ;
        [r, c, d] = size(img) ;
        for k = 1 : r
        	for l = 1 : c
        		rvalue = floor(double(img(l,k,1))/16);
        		gvalue = floor(double(img(l,k,2))/16);      		
        		bvalue = floor(double(img(l,k,3))/16);     		
        		%classNum = ['Class 1'];
        		imgReshapedInt = [imgReshapedInt rvalue gvalue bvalue];
        	end
        end
        imgReshaped = [imgReshaped; imgReshapedInt];
        %imgReshaped = [img(1,:,1) ; img(1,:,2) ; img(1,:,3)] ;
        %splitImgName = regexp(imgName, '\.', 'split') ;
        %imgNameWithoutExtension = splitImgName{1} ;
    end
    fileName = ['.csv'] ;
    fullFileName = [outDirPath1, fileName] ;
    %fid = fopen(fullFileName,'w');
    csvwrite(fullFileName, imgReshaped);    
end

