%%%%%%%%%%%%%%%%%%%%%%%%%%%Bayesian Matting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The following program is an implementation of the paper: A Bayesian
%approach to digital matting by Chuang et al.
%N.B. The total computation may take several minutes according to the size
%of the unknown area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%%%%%%%%%%%%%%%%Constants declaration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIGMA      -> sigma value for the spatial Gaussian falloff. It is set = 8  
%              as suggested in the paper
%WINDOWSIZE -> size of the sliding window The wider is the unknown area the 
%              bigger should be the size of the window. In case the window
%              is not big enough the algorithm COULD NOT END because the 
%              window, at some point, will only include unknown values. If
%              it is too big instead the quality of the result may be
%              reduced.
%SIGMAC     -> tunable parameter
%G          -> spatial Gaussian falloff
%CLUSTERMVAR-> limit used for computing the clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIGMA = 8;
WINDOWSIZE = 45;
SIGMAC = 1;
CLUSTERMVAR = 0.05;
G = fspecial('gaussian', WINDOWSIZE, SIGMA);


%load the sample image
img = im2double(imread('SampleImages2\GT03.png'));
%load its trimap
trimap = im2double(imread('SampleImages2\Trimap1\GT03.png'));


%Create B and F masks according to the trimap
B_mask = zeros(size(trimap));
B_mask(trimap == 0) = 1;
B = bsxfun(@times, img, B_mask);
F_mask = ones(size(trimap));
F_mask(trimap < 1) = 0;
F = bsxfun(@times, img, F_mask);

%Alpha is initialised with the values of the trimap even though the unknown
%values in the trimap are represented with NaN
alpha = trimap;
alpha(alpha > 0 & alpha < 1) = NaN;
%Store the position of the unknown values (NaN). This list is used in the 
%following for loops to understand on which pixels the window should slide
[row, column] = find(isnan(alpha));
unknownCoords = [row, column];

%The source image and the provided trimap are displayed
figure('Name','Bayesian Matting','NumberTitle','off')
subplot(2,2,1)
imshow(img);
title('Source image');
subplot(2,2,2)
imshow(trimap);
title('Trimap');

%temp variables used to store the pixels that are still unknown
newUnknownCoords = unknownCoords;
todelete = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (size(unknownCoords,1)>=1)
    %temp alpha value is displayed after every loop so that it is possible 
    %to veryfy the improvements
    subplot(2,2,3)
    imshow(alpha);
    title('Alpha channel');
    refreshdata
    drawnow
    
    for k = 1 : size(unknownCoords,1)
        i = unknownCoords(k,1);
        j = unknownCoords(k,2);
        %C = current pixel RGB colour
        C = reshape(img(i,j,:),3,1);    
        
        %compute the top, bottom, left and right coordinates boundaries for
        %the sliding window. This is done especially to manage special
        %cases in which the distance of the current pixel from the image 
        %borders is smaller than the preset window size (i.e. I(1,1). In 
        %that case we  apply an offset to slightly move the window towards 
        %the center of  the image. This way the window size is respected 
        %and the target pixel will still part of the window.
        [top, bottom, left, right] = getWindowCoordinates(WINDOWSIZE,i,j,img);

        %By exploiting the previously computed values, the window is
        %created separately for the foreground(Wf), Background (Wb) and
        %alpha (Walpha) then it is filled with the value of the neighbour 
        %for the current pixel. 
        Wb = B(top:bottom,left:right,:);
        Wf = F(top:bottom,left:right,:);
        Walpha = alpha(top:bottom,left:right,:);

        %compute the weight for each pixel in the windows. We only store the
        %values > 0 since these will be used for computing mean and covariance
        %The function used for computing the clustes requires the windows
        %to be shpaed as a nx3 matrix where the 3 columns represent the R,G
        %and B values for that pixel.
        %Formulas (8) and (9) - report
        
        %Foreground. 
        weightF = (Walpha.^2) .* G;
        Wf = reshape(Wf, WINDOWSIZE * WINDOWSIZE, 3);
        Wf = Wf(weightF > 0,:);
        weightF = weightF(weightF>0);

        %Background
        weightB = ((1 - Walpha).^2).*G;
        Wb = reshape(Wb, WINDOWSIZE * WINDOWSIZE, 3);
        Wb = Wb(weightB > 0, :);
        weightB = weightB(weightB>0);

        if(~isempty(Wf) && ~isempty(Wb))
            
            %cluster creation performed with Orchard and Bouman method
            [muF,SigmaF]= getClusters(Wf,weightF,CLUSTERMVAR);
            [muB,SigmaB]= getClusters(Wb,weightB,CLUSTERMVAR);  

            %temp value of alpha set = to the window mean (as described in
            %the original paper)
            tempAlpha = nanmean(Walpha(:));                    

            %solve the log likelihood for the current pixel and update the
            %values
            [f,b,a] = computeAlpha(muF,SigmaF,muB,SigmaB,SIGMAC,C,tempAlpha);
            F(i,j,:)=f;
            B(i,j,:)=b;
            alpha(i,j) = a;
            
            %store the coordinates of the new known values so that we can
            %remove them from the list of the unknown coordinates
            todelete=[todelete;k];
        end
    end
    %remove the new known pixels coordinates from the list
    unknownCoords(todelete,:)=[];
    todelete = [];
    disp(['Unknown pixels:', num2str(size(unknownCoords,1))]);
end

U = bsxfun(@times, img, alpha);
subplot(2,2,4)
imshow(U);
title('Foreground');