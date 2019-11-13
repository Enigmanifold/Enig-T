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
e12_degree=e12.*180./pi;
ptsnum=30;
diff_amplifier=1000000;
nvars=4;
lb_PSRE=[0,-pi,0,-pi];
ub_PSRE=[pi/2-0.01,pi,pi/2-0.01,pi];
fun_PSRE=@(x)epipole_corrs_to_RE(x,p1,p2,ptsnum,'s');
flag=1;
ps_RE_results=[];
counter_PSRE=1;
e12s_sph_counter=1;
for m=1:size(e1s_sph,2)
    e1_sph=e1s_sph(:,m);
    for n=1:size(e2s_sph,2)
        e2_sph=e2s_sph(:,n);
        x=[e1_sph',e2_sph'];
        min_val_initial=norm(fun_PSRE(x)).*K(1);
        e12s_sph_ordered(e12s_sph_counter,1:5)=[x,min_val_initial];
        e12s_sph_counter=e12s_sph_counter+1;
    end
end
e12s_sph_ordered=sortrows(e12s_sph_ordered,5);
e12s_sph_ordered=e12s_sph_ordered(1:round(size(e12s_sph_ordered,1)*top_frac),:);
options = optimoptions('particleswarm','InitialSwarmMatrix',e12s_sph_ordered(:,1:4),'SwarmSize',size(e12s_sph_ordered,1));
while flag
    [x,min_avg_reproj_err]=particleswarm(fun_PSRE,nvars,lb_PSRE,ub_PSRE,options);
    ps_RE_results(counter_PSRE,:)=[x.*180./pi,calculate_ang_dist(e12(1:2),x(1:2))*180/pi,calculate_ang_dist(e12(3:4),x(3:4))*180/pi,min_avg_reproj_err*K(1,1)];
    counter_PSRE=counter_PSRE+1;
    if counter_PSRE > 10
        flag=0;
    end
end
ps_RE_results=sortrows(ps_RE_results,5);
epi_gt_row=[e12_degree,0];
ps_RE_results=[e12_degree;ps_RE_results(:,1:4)]';
