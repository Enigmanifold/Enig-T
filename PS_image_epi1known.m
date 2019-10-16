function [PS_x,sorted_xstd]=PS_image_epi1known(epi1_sph,p1,p2,ptsnum,psruns)
    diff_amplifier=1000000;
    nvars=2;
    lb=[0,-pi];
    ub=[pi/2-0.01,pi];
    fun=@(x)calc_phi3_epi1known(x,epi1_sph,p1,p2,ptsnum,diff_amplifier);
    options = optimoptions('particleswarm','SwarmSize',10000,'FunctionTolerance',1e-11);
    flag=1;
    xstd=[];
    counter=1;
    while flag
        [x,minc0std]=particleswarm(fun,nvars,lb,ub,options);
        xstd(counter,:)=[x,minc0std];
        if minc0std < -10
            flag=0;
        end
        counter=counter+1;
        if counter > psruns
            flag=0;
        end
    end
    sorted_xstd=sortrows(xstd,3);
    PS_x=sort(sorted_xstd(1:5,1:2));
    PS_x=PS_x(2,:)+PS_x(3,:)+PS_x(4,:);
    PS_x=PS_x./3;
end