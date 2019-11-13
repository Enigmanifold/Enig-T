function [recim,epi,rmin,rmax,thmin,thmax,dth,imsize] = rectify2(im,F,imorder)
imsizes=size(im);
if length(imsizes) == 3
        imsize=imsizes(1:2);
        colors = imsizes(3);
    else
        imsize=imsizes;
        colors = 1;
end
if imorder == 1
        epi=null(F,'r');
    else
        epi=null(F','r');
end
rmin=0;
rmax=norm(imsize(1),imsize(2));
thmin=0;
dth=abs(atan(1/rmax));
thmax=pi-dth;
thsize=round(pi/dth);
rsize=round(4*rmax+1);
recim=zeros(thsize,rsize,colors);
recimsize=size(recim);
recim_epi=[round(thsize/2);0];
for theta=1:thsize
    th=theta*dth;
    for r=-2*rmax:1:2*rmax
        x=round(r*cos(th)+epi(2));
        y=round(r*sin(th)+epi(1));
        if x>=1 && x<=imsize(2) && y>=1 && y<=imsize(1)
            for color=1:colors
                val=im(y,x,color);
                recim(theta,r+2*rmax+1,color)=val;
            end
        else
            continue
        end
    end
end
recim=uint8(recim);