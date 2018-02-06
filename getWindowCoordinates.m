function [top, bottom, left, right] = getWindowCoordinates(WINDOWSIZE,i,j,img)
    %size of the original image
    [r, c] = size(img(:,:,1));

    %check if the size of the window can fit into the boundary of the image
    %if not, the position of the windows is adapted to respect the border
    %of the image but keeping the correct size
    if(i-floor(WINDOWSIZE/2)<=1)
        top = 1;
        bottom = WINDOWSIZE;
    else
        if(i+floor(WINDOWSIZE/2)>r)
            bottom = r;
            top = r-WINDOWSIZE+1;
        else
            top = i-floor(WINDOWSIZE/2);
            bottom = i+floor(WINDOWSIZE/2);
        end
    end

    if(j-floor(WINDOWSIZE/2)<=1)
        left = 1;
        right = WINDOWSIZE;
    else
        if(j+floor(WINDOWSIZE/2)>c)
            right = c;
            left = c-WINDOWSIZE+1;
        else
            left = j-floor(WINDOWSIZE/2);
            right = j+floor(WINDOWSIZE/2);
        end
    end
end

