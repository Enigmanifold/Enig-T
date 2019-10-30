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
e12=[e1Theta,e1Phi,e2Theta,e2Phi].*180./pi;
epi1_sph=[e1Theta;e1Phi];
ptsnum=30;
diff_amplifier=1000000;
nvars=2;
lb_PSRE_epi1known=[0,-pi];
ub_PSRE_epi1known=[pi/2-0.01,pi];
fun_PSRE_epi1known=@(x)epipole_corrs_to_RE_epi1known(x,epi1_sph,p1,p2,ptsnum,'s');
flag=1;
ps_RE=[];
counter_PSRE_epi1known=1;
while flag
    [x,min_avg_reproj_err]=particleswarm(fun_PSRE_epi1known,nvars,lb_PSRE_epi1known,ub_PSRE_epi1known);
    ps_RE(counter_PSRE_epi1known,:)=[x.*180./pi,min_avg_reproj_err*K(1,1)];
    counter_PSRE_epi1known=counter_PSRE_epi1known+1;
    if counter_PSRE_epi1known > 12
        flag=0;
    end
end
ps_RE=sortrows(ps_RE,3)';
ps_RE=ps_RE(1:2,1:12);
ps_RE=[epi2_sph.*180./pi,ps_RE];