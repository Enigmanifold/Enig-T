function [theta,phi]=implane2imsphere(vec)
    theta=acos(1/norm(vec));
    phi=atan(vec(2)/vec(1));
    if vec(1)<0
        phi=phi+pi;
    end
    if phi>pi
        phi=phi-2*pi;
    end
end