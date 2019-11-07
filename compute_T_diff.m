function angle_diff=compute_T_diff(T1,T2)
    T1=T1./norm(T1);
    T2=T2./norm(T2);
    angle_diff=acos(T1'*T2);
end