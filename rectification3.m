imdir1="/Users/tianqitang/CV_Research/images/imagePair1/Resizerect_001_0_r5000.png";
load('/Users/tianqitang/CV_Research/images/imagePair1/Data12.mat')
im1=imread(imdir1);
ydiv=40;
imsize1=size(im1);
% imsize2=size(im2);
gridwidth1=imsize1(1)/ydiv;
% gridwidth2=imsize2(1)/ydiv;
for m=1:gridwidth1:imsize1(1)
    im1(m,:,:)=0;
end
for m=1:gridwidth1:imsize1(2)
    im1(:,m,:)=0;
end
% for m=1:gridwidth2:imsize2(1)
%     im2(m,:,:)=0;
% end
% for m=1:gridwidth1:imsize2(2)
%     im2(:,m,:)=0;
% end
[recim,epi,rmin,rmax,thmin,thmax,dth,imsize]=rectify2(im1,F,1);
imshow(recim)
hold on
plot(round(size(recim,2)/2)*ones(1,size(recim,1)+1),linspace(1,size(recim,1),size(recim,1)+1),'r')
hold off