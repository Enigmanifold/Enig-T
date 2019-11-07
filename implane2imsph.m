function vec_sph=implane2imsph(vec)
% Same as implane2imsph, but easier to use.
    vec_sph=-ones(2,1);
    theta=acos(1/norm(vec));
    phi=atan(vec(2)/vec(1));
    if vec(1)<0
        phi=phi+pi;
    end
    if phi>pi
        phi=phi-2*pi;
    end
    vec_sph(1)=theta;
    vec_sph(2)=phi;
end