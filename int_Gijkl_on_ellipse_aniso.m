function [int_G] = int_Gijkl_on_ellipse_aniso(ellipse,C_tensor,tol)
% compute the auxiliary tensor D_ijkl needed to form Eshelby tensor S_ijkl
% integral is performed on the surface of ellipse in real-space
% applicable only to points inside ellipsoidal inclusion
% but the elastic medium can be anisotropic

if (~exist('ellipse'))
    ellipse = [5 4 7];
end

if (~exist('tol'))
    tol = 1e-6;
end

a=ellipse(1);
b=ellipse(2);
c=ellipse(3);

int_theta_phi = quadv(@integrand_phi,0,pi,tol);

int_G = zeros(3,3,3,3);
n=1;
for i=1:3, for j=1:3, for k=1:3, for l=1:3
    int_G(i,j,k,l) = int_theta_phi(n);
    n = n+1;
end; end; end; end

int_G = - (a*b*c)/(4*pi) * int_G;

    function int_theta = integrand_phi(phi)        
    int_theta = quadv(@integrand_theta,0,2*pi,tol);

            function data = integrand_theta(theta)                
            % integrand is (inv(zz))_ij z_k z_l sin(phi) beta^(-3)

            % area of surface patch dA_k = n_k * dS
            a = ellipse(1); b = ellipse(2); c = ellipse(3);
            beta = sqrt( (a^2*cos(theta)^2 + b^2*sin(theta)^2)*sin(phi)^2 + c^2*cos(phi)^2 );

            z = [cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi)];
            zz = zeros(3,3);
            for i=1:3, for j=1:3, for k=1:3, for l=1:3,
                            zz(j,k) = zz(j,k) + z(i)*C_tensor(i,j,k,l)*z(l);
            end; end; end; end
            invzz = inv(zz);

            data=zeros(3^4,1); n=1;
            for i=1:3, for j=1:3, for k=1:3, for l=1:3,
                data(n) = invzz(i,j)*z(k)*z(l)*sin(phi)./beta.^3;
                n = n+1;
            end; end; end; end

            end 

    end

end