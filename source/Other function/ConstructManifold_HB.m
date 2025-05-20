function [yn, dyn] = ConstructManifold_HB(a, phi, ndof, X_HB)
Qh = a*IVec(X_HB(1:end-2),ndof);
w = X_HB(end-1);
H = (size(Qh,2) - 1)/2;

Hamonic = zeros(2*H+1,size(phi,2));
Hamonic(1,:) = ones(size(phi));
for i = 1:H
    Hamonic(2*i,:) = cos(i*phi);
    Hamonic(2*i+1,:) = sin(i*phi);
end
yn = (Qh*Hamonic)';
dyn = (w*Qh*Delta(H)*Hamonic)';
end