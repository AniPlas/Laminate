function [SE, P] = eshelby_tensor_aniso (C_tensor, a, b, c, tol)
% compute Eshelby tensor in anisotropic elasticity (using numerical integration)
% C_tensor: elastic stiffness tensor
% a, b, c: semi-axes of ellipsoid
%

% numerical integration - anisotropic elasticity
D_tensor_aniso = int_Gijkl_on_ellipse_aniso([a b c],C_tensor,tol);
SE_tensor_aniso = zeros(size(C_tensor)); P_tensor_aniso = zeros(size(C_tensor));

for i=1:3, for j=1:3, for k=1:3, for l=1:3
    P_tensor_aniso(i,j,k,l) = -0.25*(D_tensor_aniso(i,k,l,j)+D_tensor_aniso(j,l,k,i)+D_tensor_aniso(i,l,k,j)+D_tensor_aniso(j,k,l,i));
end; end; end; end

for i=1:3, for j=1:3, for m=1:3, for n=1:3, for k=1:3, for l=1:3
    %SE_tensor_aniso(i,j,m,n) = SE_tensor_aniso(i,j,m,n)-0.5*C_tensor(l,k,m,n)*(D_tensor_aniso(i,k,l,j)+D_tensor_aniso(j,k,l,i));
    SE_tensor_aniso(i,j,m,n) = SE_tensor_aniso(i,j,m,n)+P_tensor_aniso(i,j,k,l)*C_tensor(l,k,m,n);
end; end; end; end; end; end

voigt_ind = [ 1 6 5
              6 2 4
              5 4 3 ];

SE = zeros(6,6); D = zeros(6,6); C = zeros(6,6);
for i=1:3, for j=1:3, for k=1:3, for l=1:3,
       I = voigt_ind(i,j);
       J = voigt_ind(k,l);
       SE(I,J) = SE_tensor_aniso(i,j,k,l)*(1+(I>=4));
       P (I,J) = P_tensor_aniso (i,j,k,l)*(1+(I>=4))*(1+(J>=4));
       C (I,J) = C_tensor(i,j,k,l);
end; end; end; end
    
end
