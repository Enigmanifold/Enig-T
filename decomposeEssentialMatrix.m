function [R1,T1,R2,T2,R3,T3,R4,T4] = decomposeEssentialMatrix(E)
    TT = 0.5 * trace(E * E') * eye(3) - E*E';
    
    [~,pos] = max(sum(TT,2));
    T = (TT(pos,:) ./ norm(TT(pos,:)));
    T = TT(pos,:) ./ sqrt(TT(pos,pos));
    Tx = [0 -T(3) T(2);
      T(3) 0  -T(1);
      -T(2) T(1) 0];
    cofactorE = [cross(E(:,2),E(:,3)), cross(E(:,3),E(:,1)),cross(E(:,1),E(:,2))]';
    Ra = (transpose(cofactorE) - Tx * E) ./ (T*T');
    Rb = (transpose(cofactorE) + Tx * E) ./ (T*T');
    
    R1 = Ra; T1 = T;
    R2 = Ra; T2 = -T;
    R3 = Rb; T3 = T;
    R4 = Rb; T4 = -T;
    
    T1 = T1 ./ norm(T1);
    T2 = T2 ./ norm(T2);
    T3 = T3 ./ norm(T3);
    T4 = T4 ./ norm(T4);
end