function [Ap,Bp,Cp,Dp] = plane2disk(A,B,C,D)
    Ap = ((eye(length(A))-A') \ (A' + eye(length(A))))';
    Bp = ((eye(length(A))+Ap) * B) / sqrt(2);
    Cp = (C * (eye(length(A))+Ap)) / sqrt(2);
    Dp = D + Cp * ((eye(length(A))+Ap) \ Bp);
end