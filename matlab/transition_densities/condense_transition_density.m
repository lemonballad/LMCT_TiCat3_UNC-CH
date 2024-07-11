function [cond_rho,err_term]=condense_transition_density(center,radius,X,Y,Z,rho_)
% rho_=trans_rho(:,:,:,8);
dV=abs(X(1)-X(2))*abs(Y(1)-Y(2))*abs(Z(1)-Z(2));
NX=length(X);
NY=length(Y);
NZ=length(Z);
NCenters=length(radius);
dist=zeros(NZ,NY,NX,'double');
I_mat=ones(NZ,NY,NX,'double');
cond_rho(1:4)=0;
remains=rho_;
for ic=1:NCenters
    [Y_mat,Z_mat,X_mat]=...
        meshgrid((Y-center(2,ic)).^2,(Z-center(3,ic)).^2,(X-center(1,ic)).^2);
    if ic==1 % Account for infinite cylindrilical regions
            dist(:,:,:)=sqrt(X_mat+Y_mat)-radius(ic)*I_mat;
    else % Account for spherical regions
            dist(:,:,:)=sqrt(X_mat+Y_mat+Z_mat)-radius(ic)*I_mat;
    end
    cond_rho(ic)=trapz(rho_(dist<0))*dV;
    remains((dist<0))=0;
end
err_term=trapz(trapz(trapz(remains)))*dV;
end