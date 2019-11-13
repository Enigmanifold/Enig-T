function angle_diff=compute_R_diff(R1,R2)
    angle_diff=acos((trace(R1'*R2)-1)/2);
end