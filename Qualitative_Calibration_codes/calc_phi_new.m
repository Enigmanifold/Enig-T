%% Inputs: Image points, epipoles, number of points. Outputs: Four sets of phi's corresponding to four ambiguities.
function [phi1_best,phi2_best,phi1_alt_best,phi2_alt_best]=calc_phi_new(p1,p2,e1,e2,ptsnum)
    if ptsnum<1
        error('Need more points.')
    end
    % Calculate phi's.
    gamma_cross_e=zeros(3,ptsnum);
    e_norm=norm(e1);
    ebar_norm=norm(e2);
    vbar=e2./ebar_norm;
    ebar_dot_gamma=zeros(1,ptsnum);
    e_dot_ebar=e1'*e2;
    eebar_norm=e_norm*ebar_norm;
    eebar_ebare=e_norm.*e2+ebar_norm.*e1;
    ebar_cross_e=cross(e2,e1);
    eeee=eebar_norm+e_dot_ebar;
    eeeee=ebar_norm*eeee;
    R2_additional_axis_y=-(e2(1)+1)/e2(2);
    R2_additional_cross=[0,-1,R2_additional_axis_y;1,0,-1;-R2_additional_axis_y,1,0];
    R2_additional_cross2=R2_additional_cross*R2_additional_cross;
    R2_additional=eye(3)+2*e2(2)^2/(2*e2(2)^2+e2(1)^2+2*e2(1)+1)*R2_additional_cross2;
    phi1=zeros(1,ptsnum);
    phi1_alt=zeros(1,ptsnum);
    for m=1:ptsnum
        gamma_cross_e(:,m)=cross(p1(:,m),e1);
        ebar_dot_gamma(m)=e2'*p1(:,m);
    end
    for m=1:ptsnum
        N2=eeeee.*gamma_cross_e(:,m)+p1(:,m)'*ebar_cross_e.*eebar_ebare;
        N2_alt=R2_additional*N2;
        tanphi1=-(p2(:,m)'*N2)/(p2(:,m)'*cross(vbar,N2));
        tanphi1_alt=-(p2(:,m)'*N2_alt)/(p2(:,m)'*cross(vbar,N2_alt));
        phi1_new=atan(tanphi1);
        phi1_alt_new=atan(tanphi1_alt);
        phi1(m)=phi1_new;
        phi1_alt(m)=phi1_alt_new;
    end
    % Progressively add pi to phi's to ensure that std(phi) is minimized.
    phi1_sorted=sort(phi1);
    phi1_alt_sorted=sort(phi1_alt);
    phi1_sorted_array=repmat(phi1_sorted,ptsnum,1);
    phi1_alt_sorted_array=repmat(phi1_alt_sorted,ptsnum,1);
    for m=1:ptsnum-1
        phi1_sorted_array(m,1:m)=phi1_sorted_array(m,1:m)+ones(1,m).*pi;
        phi1_alt_sorted_array(m,1:m)=phi1_alt_sorted_array(m,1:m)+ones(1,m).*pi;
    end
    % Select the phi's that yields the minimal std.
    [~,min_std_ind]=min(std(phi1_sorted_array,0,2));
    [~,min_std_ind_alt]=min(std(phi1_alt_sorted_array,0,2));
    phi1_best=phi1_sorted_array(min_std_ind,:);
    phi1_alt_best=phi1_alt_sorted_array(min_std_ind_alt,:);
    phi2_best=phi1_best-ones(1,ptsnum).*pi;
    phi2_best=phi2_best+(phi2_best<-pi).*2.*pi;
    phi2_alt_best=phi1_alt_best-ones(1,ptsnum).*pi;
    phi2_alt_best=phi2_alt_best+(phi2_alt_best<-pi).*2.*pi;
end