function rot_matrix=rot2sph_pt(sph_pt)
% Input: Point in spherical coordinates. Output: Rotation marix that
% rotates theta=0 point to sph_pt, assuming sph_pt has radial coordinate 1.
    rot_axis=cross([0,0,1],[sin(sph_pt(1))*cos(sph_pt(2)),sin(sph_pt(1))*sin(sph_pt(2)),cos(sph_pt(1))]);
    rot_axis=rot_axis./norm(rot_axis);
    %rot_angle=phi;
    rot_axis_cross=cross_matrix(rot_axis);
    rot_axis_cross2=rot_axis_cross*rot_axis_cross;
    rot_matrix=eye(3)+sin(sph_pt(1)).*rot_axis_cross+(1-cos(sph_pt(1))).*rot_axis_cross2;
end