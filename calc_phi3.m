function [minc0std,minc0stdind,chosen_phis]=calc_phi3(epi12_sph,p1,p2,ptsnum,diff_amplifier)
    epi1_plane=[tan(epi12_sph(1))*cos(epi12_sph(2));tan(epi12_sph(1))*sin(epi12_sph(2));1];
    epi2_plane=[tan(epi12_sph(3))*cos(epi12_sph(4));tan(epi12_sph(3))*sin(epi12_sph(4));1];
    [phi1,phi2,phi1_alt,phi2_alt]=calc_phi(p1(:,1:ptsnum),p2(:,1:ptsnum),epi1_plane,epi2_plane,ptsnum);
    phis=[phi1;phi2;phi1_alt;phi2_alt];
    phis_c0=-ones(4,ptsnum);
    for c=1:4               
        for d=1:ptsnum
            [~,phi2_c0_ind]=min([abs(phis(c,d)),abs(phis(c,d)-pi),abs(phis(c,d)-2*pi)]);
            switch phi2_c0_ind
                case 1
                    phis_c0(c,d)=phis(c,d);
                case 2
                    phis_c0(c,d)=phis(c,d)-pi;
                case 3
                    phis_c0(c,d)=phis(c,d)-2*pi;
            end
        end
    end
    mean1=mean(phis_c0(1,:));
    mean2=mean(phis_c0(2,:));
    mean3=mean(phis_c0(3,:));
    mean4=mean(phis_c0(4,:));
    phi1_c0std=std((phis_c0(1,:)-ones(1,ptsnum).*mean1).*diff_amplifier+ones(1,ptsnum).*mean1.*180./pi);
    phi2_c0std=std((phis_c0(2,:)-ones(1,ptsnum).*mean2).*diff_amplifier+ones(1,ptsnum).*mean2.*180./pi);
    phi3_c0std=std((phis_c0(3,:)-ones(1,ptsnum).*mean3).*diff_amplifier+ones(1,ptsnum).*mean3.*180./pi);
    phi4_c0std=std((phis_c0(4,:)-ones(1,ptsnum).*mean4).*diff_amplifier+ones(1,ptsnum).*mean4.*180./pi);
%     phi1_c0std=log2(std((phis_c0(1,:)-ones(1,ptsnum).*mean1).*diff_amplifier+ones(1,ptsnum).*mean1.*180./pi));
%     phi2_c0std=log2(std((phis_c0(2,:)-ones(1,ptsnum).*mean2).*diff_amplifier+ones(1,ptsnum).*mean2.*180./pi));
%     phi3_c0std=log2(std((phis_c0(3,:)-ones(1,ptsnum).*mean3).*diff_amplifier+ones(1,ptsnum).*mean3.*180./pi));
%     phi4_c0std=log2(std((phis_c0(4,:)-ones(1,ptsnum).*mean4).*diff_amplifier+ones(1,ptsnum).*mean4.*180./pi));
    [minc0std,minc0stdind]=min([phi1_c0std,phi2_c0std,phi3_c0std,phi4_c0std]);
    chosen_phis=phis(minc0stdind,:);
end