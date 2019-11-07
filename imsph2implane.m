function plane_coords=imsph2implane(sph_coords)
    plane_coords=ones(3,1);
    plane_coords(1)=tan(sph_coords(1))*cos(sph_coords(2));
    plane_coords(2)=tan(sph_coords(1))*sin(sph_coords(2));
end