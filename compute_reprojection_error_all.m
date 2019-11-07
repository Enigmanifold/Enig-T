function errs=compute_reprojection_error_all(R,T,p1,p2,ptsnum,varargin)
% Compute reprojection errors using R,T and correspondences, either in
% focal length or pixel regime. If in pixel regime, should specify K in
% varargin.
    narginchk(5,6);
    p1=p1(:,1:ptsnum);
    p2=p2(:,1:ptsnum);
    T_cross=cross_matrix(T);
    E=T_cross*R;
    errs=zeros(1,ptsnum);
    p2_trans=p2';
    if nargin==5 % Focal length regime.
        abc=E*p1;
        ab=abc(1:2,:);
        norm_const = vecnorm(ab,2,1);
        for n = 1:ptsnum
            err = abs(p2_trans(n,:)*abc(:,n));
            errs(1,n) = err;
        end
        errs = errs./norm_const;
    end
    if nargin==6 % Pixel regime.
        K=varargin{1};
        invK=inv(K);
        F=invK'*E*invK;
        abc=F*p1;
        ab=abc(1:2,:);
        norm_const=vecnorm(ab,2,1);
        for n=1:ptsnum
            err=abs(p2_trans(n,:)*abc(:,n));
            errs(1,n)=err;
        end
        errs = errs./norm_const;
    end
end