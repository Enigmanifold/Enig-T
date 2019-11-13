function [C0,C1,C2,C3,C4]=calc_eta2s_coeffs(e1,xi2,p1,p2,p3,p4)
% Knowing the first epipole and one component of the second epipole, using
% 2 pairs of image points, calculate the coefficients of the quartic
% equation of the second component of the second epipole.
%     [xi1,eta1]=deal(e1(1:2));
%     [x1,y1]=deal(p1(1:2));
%     [x2,y2]=deal(p2(1:2));
    [x3,y3]=deal(p3(1),p3(2));
    [x4,y4]=deal(p4(1),p4(2));
    e1g1=cross(e1,p1);
    e1g2=cross(e1,p2);
    e1g1De1g2=e1g1'*e1g2;
    e1g1De1g2squared=e1g1De1g2^2;
    e1g1squared=e1g1'*e1g1;
    e1g2squared=e1g2'*e1g2;
    lhsC0=((x3-xi2)^2*(x4-xi2)^2+(x4-xi2)^2*y3^2+(x4-xi2)^2*xi2^2*y3^2+(x3-xi2)^2*y4^2+(x3-xi2)^2*xi2^2*y4^2+y3^2*y4^2+2*xi2^2*y3^2*y4^2+xi2^4*y3^2*y4^2)*e1g1De1g2squared;
    lhsC1=(-2*(x4-xi2)^2*y3-2*x3*(x4-xi2)^2*xi2*y3-2*(x3-xi2)^2*y4-2*x4*(x3-xi2)^2*xi2*y4-2*y3^2*y4-2*x4*xi2*y3^2*y4-2*xi2^2*y3^2*y4-2*x4*xi2^3*y3^2*y4-2*y3*y4^2-2*x3*xi2*y3*y4^2-2*xi2^2*y3*y4^2-2*x3*xi2^3*y3*y4^2)*e1g1De1g2squared;
    lhsC2=((x3-xi2)^2+x4^2*(x3-xi2)^2+(x4-xi2)^2+x3^2*(x4-xi2)^2+y3^2+x4^2*y3^2+xi2^2*y3^2+x4^2*xi2^2*y3^2+4*y3*y4+4*x3*xi2*y3*y4+4*x4*xi2*y3*y4+4*x3*x4*xi2^2*y3*y4+y4^2+x3^2*y4^2+xi2^2*y4^2+x3^2*xi2^2*y4^2)*e1g1De1g2squared;
    lhsC3=(-2*y3-2*x4^2*y3-2*x3*xi2*y3-2*x3*x4^2*xi2*y3-2*y4-2*x3^2*y4-2*x4*xi2*y4-2*x3^2*x4*xi2*y4)*e1g1De1g2squared;
    lhsC4=(1+x3^2+x4^2+x3^2*x4^2)*e1g1De1g2squared;
    rhsC0=((x3-xi2)^2*(x4-xi2)^2+2*(x3-xi2)*(x4-xi2)*y3*y4+2*(x3-xi2)*(x4-xi2)*xi2^2*y3*y4+y3^2*y4^2+2*xi2^2*y3^2*y4^2+xi2^4*y3^2*y4^2)*e1g1squared*e1g2squared;
    rhsC1=(-2*(x3-xi2)*(x4-xi2)*y3-2*x4*(x3-xi2)*(x4-xi2)*xi2*y3-2*(x3-xi2)*(x4-xi2)*y4-2*x3*(x3-xi2)*(x4-xi2)*xi2*y4-2*y3^2*y4-2*x4*xi2*y3^2*y4-2*xi2^2*y3^2*y4-2*x4*xi2^3*y3^2*y4-2*y3*y4^2-2*x3*xi2*y3*y4^2-2*xi2^2*y3*y4^2-2*x3*xi2^3*y3*y4^2)*e1g1squared*e1g2squared;
    rhsC2=(2*(x3-xi2)*(x4-xi2)+2*x3*x4*(x3-xi2)*(x4-xi2)+y3^2+2*x4*xi2*y3^2+x4^2*xi2^2*y3^2+4*y3*y4+2*x3*x4*y3*y4+2*x3*xi2*y3*y4+2*x4*xi2*y3*y4+2*xi2^2*y3*y4+4*x3*x4*xi2^2*y3*y4+y4^2+2*x3*xi2*y4^2+x3^2*xi2^2*y4^2)*e1g1squared*e1g2squared;
    rhsC3=(-2*y3-2*x3*x4*y3-2*x4*xi2*y3-2*x3*x4^2*xi2*y3-2*y4-2*x3*x4*y4-2*x3*xi2*y4-2*x3^2*x4*xi2*y4)*e1g1squared*e1g2squared;
    rhsC4=(1+2*x3*x4+x3^2*x4^2)*e1g1squared*e1g2squared;
    C0=lhsC0-rhsC0;
    C1=lhsC1-rhsC1;
    C2=lhsC2-rhsC2;
    C3=lhsC3-rhsC3;
    C4=lhsC4-rhsC4;
end