run_num=1;
e1th_max=pi/2;
e2th_max=pi/2;
eta2_results=[];
for m=1:run_num
    [R,T,e1,e2,phi,K,p1,p2,p1ImageOnImage,p2ImageOnImage]=gen_RT(e1th_max,e2th_max);
    xi2=e2(1);
    eta2=e2(2);
    pt1=p1(:,1);
    pt2=p1(:,2);
    pt3=p2(:,1);
    pt4=p2(:,2);
%     lhsrhs=calc_eta2s_lhsrhs(e1,e2,xi2,pt1,pt2,pt3,pt4);
    [et2s,C4,C3,C2,C1,C0]=calc_eta2s(e1,xi2,pt1,pt2,pt3,pt4);
%     [~,C4,C3,C2,C1,C0]=calc_eta2s_coeffs(e1,xi2,pt1,pt2,pt3,pt4);
    quartic_result=C4*eta2^4+C3*eta2^3+C2*eta2^2+C1*eta2+C0;
end