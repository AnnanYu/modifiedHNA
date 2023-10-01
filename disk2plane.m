function [Ad,Bd,Cd,Dd] = disk2plane(A,B,C,D)
    Ad = (eye(length(A)) + A) \ (A - eye(length(A)));
    Bd = sqrt(2) * ((eye(length(A)) + A) \ B);
    Cd = sqrt(2) * ((eye(length(A)) + A') \ C')';
    Dd = D - C * ((eye(length(A)) + A) \ B);
end