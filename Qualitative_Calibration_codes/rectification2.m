imdir1="/Users/tianqitang/CV_Research/images/imagePair1/Resizerect_001_0_r5000.png";
imdir2="/Users/tianqitang/CV_Research/images/imagePair1/Resizerect_002_0_r5000.png";
load('/Users/tianqitang/CV_Research/images/imagePair1/Data12.mat')
im1=imread(imdir1);
im2=imread(imdir2);
ydiv=40;
imsize1=size(im1);
imsize2=size(im2);
gridwidth1=imsize1(1)/ydiv;
gridwidth2=imsize2(1)/ydiv;
for m=1:gridwidth1:imsize1(1)
    im1(m,:,:)=0;
end
for m=1:gridwidth1:imsize1(2)
    im1(:,m,:)=0;
end
for m=1:gridwidth2:imsize2(1)
    im2(m,:,:)=0;
end
for m=1:gridwidth1:imsize2(2)
    im2(:,m,:)=0;
end
[recim1,epi1,rmin1,rmax1,thmin1,thmax1,dth1,imsize1,invert1]=rectify(im1,F,1);
[recim2,epi2,rmin2,rmax2,thmin2,thmax2,dth2,imsize2,invert2]=rectify(im2,F,2);
recim1_deci = recim1(290:300,290:300,:);
recimsize1=size(recim1);
recimsize2=size(recim2);
thetadiv = 9;
% figure
% imshow(recim1)
% figure
% imshow(recim2)
fig1 = figure;
imshow(recim1)
% hold on
fig2 = figure;
imshow(recim2)
% hold on
saveas(fig1,'recim1.jpg');
saveas(fig2,'recim2.jpg');
% colorset = ['y','m','c','r','g','b','w','k'];
% for n = 1:thetadiv-1
%     theta1=thmin1+n*(thmax1-thmin1)/thetadiv;
%     gamma1=[imsize1(2);tan(theta1)*(imsize1(2)-epi1(2))+epi1(1);1];
%     perp2=gamma1'*F';
%     theta2=atan(-perp2(1)/perp2(2));
%     if theta2<thmin2
%         theta2=theta2+pi;
%     end
%     if theta2>thmax2
%         theta2=theta2-pi;
%     end
%     if invert2 ~= invert1
%         theta2 = (thmax2+thmin2)-theta2;
%     end
%     slope=-perp2(1)/perp2(2);
%     intersect=-perp2(3)/perp2(2);
%     intersect2=slope*(0-epi2(1))+epi2(2);
%     figure(fig1)
%     plot(linspace(1,recimsize1(2),recimsize1(2)),ones(1,recimsize1(2))*round((theta1-thmin1)/dth1),colorset(n))
%     if theta2>=thmin2 && theta2 <= thmax2
%     figure(fig2)
%     hold on
%     plot(linspace(1,recimsize2(2),recimsize2(2)),ones(1,recimsize2(2))*round((theta2-thmin2)/dth2),colorset(n))
%     end
% end