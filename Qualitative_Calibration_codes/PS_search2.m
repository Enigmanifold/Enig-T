run('generate_RT.m');
e1Theta=acos(1/norm(e1));
e2Theta=acos(1/norm(e2));
e1cosTheta=cos(e1Theta);
e1sinTheta=sin(e1Theta);
e2cosTheta=cos(e2Theta);
e2sinTheta=sin(e2Theta);
e1Phi=atan(e1(2)/e1(1));
if e1(1)<0
    e1Phi=e1Phi+pi;
    if e1Phi>pi
        e1Phi=e1Phi-2*pi;
    end
end
e2Phi=atan(e2(2)/e2(1));
if e2(1)<0
    e2Phi=e2Phi+pi;
    if e2Phi>pi
        e2Phi=e2Phi-2*pi;
    end
end
e12=[e1Theta,e1Phi,e2Theta,e2Phi];
epi1_sph=e12(1:2)';
ptsnum=30;
diff_amplifier=1000000;
nvars=2;
lb=[0,-pi];
ub=[pi/2-0.01,pi];
fun_PSstdphi=@(x)calc_phi3_epi1known(x,epi1_sph,p1,p2,ptsnum,diff_amplifier);
% options = optimoptions('particleswarm','InitialSwarmMatrix',e12s_sph,'SwarmSize',size(e12s_sph,1));
flag=1;
xstd=[];
counter=1;
while flag
    [x,minc0std]=particleswarm(fun_PSstdphi,nvars,lb,ub);
    xstd(counter,:)=[x.*180./pi,minc0std];
    if minc0std < -10
        flag=0;
    end
    counter=counter+1;
end
sorted_xstd=sortrows(xstd,3);
sx_array=sorted_xstd';
sx_array=sx_array(1:2,1:12);
sx_array=[epi2_sph.*180./pi,sx_array];
% total_array=reshape(cell2mat(total_cell(1:12,2)'),[4,12]);