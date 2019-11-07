function epiplane_angle=calc_epiplane_angle(epi_sph,pt1,pt2)
    epi=imsph2implane(epi_sph);
%     v=epi./norm(epi);
%     epiline_seg1=pt1-epi;
%     epiline_seg2=pt2-epi;
% %     u1=epiline_seg1./norm(epiline_seg1);
% %     u2=epiline_seg2./norm(epiline_seg2);
%     
%     epiplane_vec1=cross(pt1,v);
%     epiplane_vec2=cross(pt2,v);
%     epiplane_normvec1=epiplane_vec1./norm(epiplane_vec1);
%     epiplane_normvec2=epiplane_vec2./norm(epiplane_vec2);
%     
%     intersection1=cross(epiplane_normvec1,v);
%     intersection2=cross(epiplane_normvec2,v);
%     intersection1=epi'*pt1*epi-epi'*epi*pt1;
%     intersection2=epi'*pt2*epi-epi'*epi*pt2;
%     intersection1=intersection1./norm(intersection1);
%     intersection2=intersection2./norm(intersection2);
%     epiplane_angle=compute_T_diff(intersection1,intersection2);
    vec1=cross(epi,pt1);
    vec2=cross(epi,pt2);
    vec1=vec1./norm(vec1);
    vec2=vec2./norm(vec2);
    epiplane_angle=compute_T_diff(vec1,vec2);
end