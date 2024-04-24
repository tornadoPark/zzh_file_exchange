function uLid = getuLid(i,ncellr,dr,Q_evap,l_cond,rho_l,hfg_vapor)

global rMin

% x = (i-0.5)*dr;
% xMax = (ncellr)*dr;
% uLidy = -x*(x-xMax) ;
% uLidy = -2e-3;
uLidy = Q_evap/(2*rMin*pi*l_cond)/rho_l/hfg_vapor;

uLid = [0,uLidy];
end

