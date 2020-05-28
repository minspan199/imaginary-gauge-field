
img=imread('Capture.PNG'); 

img=rgb2gray(img);

Nt = 100;

[M, N] = size(img);
m = (M/Nt);
n = (N/Nt);

img_binarized = zeros(Nt, Nt);

for ii = 1:1:Nt
    
    for jj = 1:1:Nt
        
        img_binarized(ii, jj) = img(round(ii*m), round(jj*n));
        
    end
    
end

for ii = 1:1:Nt
    
    for jj = 1:1:Nt
        
        if img_binarized(ii, jj)>200
            img_binarized(ii, jj) = 0;
        else
            img_binarized(ii, jj) = 1;
        end
        
    end
    
end

figure;
imagesc(abs(img));
colormap gray
colorbar
title('Original Image')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

figure;
imagesc(abs(img_binarized));
colormap gray
colorbar
title('Binarized Image(desired pattern)')
set(gcf, 'Position', [00, 00, 350, 300])
set(gca,'FontSize', 12) % Font Size

