%% Determine epipoles
imdir="/Users/tianqitang/CV_Research/images/imagePair1/Resizerect_001_0_r5000.png";
load('/Users/tianqitang/CV_Research/images/imagePair1/Data12.mat')
im=imread(imdir);
im=rgb2gray(im);
imsizes=size(im);
if length(imsizes) == 3
    imsize=imsizes(1:2);
    colors = imsizes(3);
else
    imsize=imsizes;
    colors = 1;
end
epi=null(F,'r');
epi=epi/epi(3);
temp=epi(2);
epi(2)=epi(1);
epi(1)=temp;
%% Calculate rmax,rmin,thetamax,thetamin
cornerdiff=[epi(1),epi(1)-imsize(1),epi(2),epi(2)-imsize(2)];
yup = epi(1)>imsize(1)/2;
xup = epi(2)>imsize(2)/2;
determine=[cornerdiff(1)>0,cornerdiff(2)>0,cornerdiff(3)>0,cornerdiff(4)>0];
invert=false;
if cornerdiff(1)<0 && cornerdiff(2)<0 && cornerdiff(3)<0 && cornerdiff(4)<0
    rmin=abs(norm(cornerdiff(1),cornerdiff(3)));
    rmax=abs(norm(cornerdiff(2),cornerdiff(4)));
    thmin=atan(cornerdiff(1)/cornerdiff(4));
    thmax=atan(cornerdiff(2)/cornerdiff(3));
end
if cornerdiff(1)>=0 && cornerdiff(2)<0 && cornerdiff(3)<0 && cornerdiff(4)<0
    rmin=abs(cornerdiff(3));
    if yup
        rmax=abs(norm([cornerdiff(1),cornerdiff(4)]));
    else
        rmax=abs(norm([cornerdiff(2),cornerdiff(4)]));
    end
    thmin=atan(cornerdiff(1)/cornerdiff(3));
    thmax=atan(cornerdiff(2)/cornerdiff(3));
end
if cornerdiff(1)>=0 && cornerdiff(2)>=0 && cornerdiff(3)<0 && cornerdiff(4)<0
    rmin=abs(norm([cornerdiff(2),cornerdiff(3)]));
    rmax=abs(norm([cornerdiff(1),cornerdiff(4)]));
    thmin=atan(cornerdiff(1)/cornerdiff(3));
    thmax=atan(cornerdiff(2)/cornerdiff(4));
end
if cornerdiff(1)>=0 && cornerdiff(2)>=0 && cornerdiff(3)>=0 && cornerdiff(4)<0
    rmin=abs(cornerdiff(2));
    if xup
        rmax=abs(norm([cornerdiff(1),cornerdiff(3)]));
    else
        rmax=abs(norm([cornerdiff(1),cornerdiff(4)]));
    end
    thmin=atan(cornerdiff(2)/cornerdiff(3))+pi;
    thmax=atan(cornerdiff(2)/cornerdif(4))+2*pi;
end
if cornerdiff(1)>=0 && cornerdiff(2)>=0 && cornerdiff(3)>=0 && cornerdiff(4)>=0
    rmin=abs(norm([cornerdiff(2),cornerdiff(4)]));
    rmax=abs(norm([cornerdiff(1),cornerdiff(3)]));
    thmin=atan(cornerdiff(2)/cornerdiff(3))+pi;
    thmax=atan(cornerdiff(1)/cornerdiff(4))+pi;
    thinvert=true;
    rinvert=true;
end
if cornerdiff(1)>=0 && cornerdiff(2)<0 && cornerdiff(3)>=0 && cornerdiff(4)>=0
    rmin=abs(cornerdiff(4));
    if yup
        rmax=abs(norm([cornerdiff(1),cornerdiff(3)]));
    else
        rmax=abs(norm([cornerdiff(2),cornerdiff(3)]));
    end
    thmin=atan(cornerdiff(2)/cornerdiff(4));
    thmax=atan(cornerdiff(1)/cornerdiff(4));
    thinvert=true;
    rinvert=true;
end
if cornerdiff(1)<0 && cornerdiff(2)<0 && cornerdiff(3)>=0 && cornerdiff(4)>=0
    rmin=abs(norm([cornerdiff(1),cornerdiff(4)]));
    rmax=abs(norm([cornerdiff(2),cornerdiff(3)]));
    thmin=atan(cornerdiff(2)/cornerdiff(4))+pi;
    thmax=atan(cornerdiff(1)/cornerdiff(3))+pi;
    thinvert=true;
    rinvert=true;
end
if cornerdiff(1)<0 && cornerdiff(2)<0 && cornerdiff(3)>=0 && cornerdiff(4)<0
    rmin=abs(cornerdiff(1));
    if xup
        rmax=abs(norm([cornerdiff(2),cornerdiff(3)]));
    else
        rmax=abs(norm([cornerdiff(2),cornerdiff(4)]));
    end
    thmin=atan(cornerdiff(1)/cornerdiff(4));
    thmax=atan(cornerdiff(1)/cornerdiff(3))+pi;        
end
dth=abs(atan(1/rmax));
thsize=round(abs((thmax-thmin)/dth))+1;
rsize=round(rmax-rmin+1)+1;
recim=zeros(thsize,rsize,colors);
for recth = 1:thsize
    th=recth*dth+thmin;
    for recr = 1:rsize
        r=recr+rmin;
        y=round(r*sin(th)+epi(1));
        x=round(r*cos(th)+epi(2));
        if y>1 && y<imsize(1) && x>1 && x<imsize(2)
            if thinvert == false && rinvert == false
                for color = 1:colors
                    val=im(y,x);
                    recim(recth,recr,color)=val;
                end
            else if thinvert == true && rinvert == true
                    for color = 1:colors
                    val=im(y,x);
                    recim(thsize-recth+1,rsize-recr+1,color)=val;
                    end
                end
            end
        end
    end
end
recim=uint8(recim);
figure
imshow(recim)
% xlabel('\nu')
% ylabel('\theta')
% xlim([rmin,rmax])
% ylim([thmin,thmax])
%% Matching theta
