function [vE,vN] = FNC_EV2siteSV(EV,SITES)


% [vE,vN] = FNC_EV2siteSV(EV,SITES)
%
% This function uses an Euler vector in input to calculate a set of surface 
% velocities, and outputs arrays usable as txt-saved files in plots. 
% 
% EV is a 1x11 array containing the Euler vector and its covariance. 
% Format is the usual:
% DRorEV(1) = initial time, or NaN
% DRorEV(2) = final time, or NaN
% DRorEV(3) = Euler pole longitude (deg-E)
% DRorEV(4) = Euler pole latitude (deg-N)
% DRorEV(5) = angular velocity (deg/Myr)
% DRorEV(6) = covariance XX (rad/Myr)^2
% DRorEV(7) = covariance XY (rad/Myr)^2
% DRorEV(8) = covariance XZ (rad/Myr)^2
% DRorEV(9) = covariance YY (rad/Myr)^2
% DRorEV(10) = covariance YZ (rad/Myr)^2
% DRorEV(11) = covariance ZZ (rad/Myr)^2
%
% SITES [n x 2] contains longitude and latitude (deg) of sites where to
% calculate velocities.
%
% vE/vN [n x 1 each] are arrays (mm/yr) containing the East/North components of the
% velocity predicted by EV at SITES.

%--------------- PARAMETERS -----------------------------------------------
Re_m = 6371e3;
Re_mm = 1e3*Re_m;
Nsamples = 2;

SITES_LL = SITES;
S_LLr = SITES.*(pi/180);
NST = length(SITES_LL(:,1));

Lo_MTX = repmat(SITES_LL(:,1),[1 Nsamples]);
La_MTX = repmat(SITES_LL(:,2),[1 Nsamples]);
%--------------------------------------------------------------------------




%--------------- GENERATES EV ENSEMBLE ------------------------------------
[wx,wy,wz] = sph2cart(EV(3).*(pi/180),EV(4).*(pi/180),EV(5).*(pi/180));%deg/Myr > rad/Myr;
EVc = [wx,wy,wz];

EVcvr = EV(1,6:11);

dx = zeros(1,Nsamples);
dy = dx;
dz = dx;

EV_ENSEMBLE(:,1) = [EVc(1)+real(dx)'].*1e-6;%rad/yr
EV_ENSEMBLE(:,2) = [EVc(2)+real(dy)'].*1e-6;
EV_ENSEMBLE(:,3) = [EVc(3)+real(dz)'].*1e-6;
%--------------------------------------------------------------------------


    
%--------------- VELOCITIES FROM EV ENSEMBLE ------------------------------
[x,y,z] = FNC_pLoLaHt2pXYZ_ELLIPSOID(SITES_LL(:,1),SITES_LL(:,2),zeros(size(SITES_LL(:,1))));%m
x=1e3.*x;%mm
y=1e3.*y;%mm
z=1e3.*z;%mm

%mm/yr
vx_MTX = repmat(EV_ENSEMBLE(:,2),[1 length(NST)]).*repmat(z',[Nsamples 1]) - repmat(EV_ENSEMBLE(:,3),[1 length(NST)]).*repmat(y',[Nsamples 1]);
vy_MTX = repmat(EV_ENSEMBLE(:,3),[1 length(NST)]).*repmat(x',[Nsamples 1]) - repmat(EV_ENSEMBLE(:,1),[1 length(NST)]).*repmat(z',[Nsamples 1]);
vz_MTX = repmat(EV_ENSEMBLE(:,1),[1 length(NST)]).*repmat(y',[Nsamples 1]) - repmat(EV_ENSEMBLE(:,2),[1 length(NST)]).*repmat(x',[Nsamples 1]);

[vE_MTX,vN_MTX,vU_MTX] = FNC_vXYZ2vEN_CONVERSION(Lo_MTX',La_MTX',vx_MTX,vy_MTX,vz_MTX);
%--------------------------------------------------------------------------


vE = vE_MTX(1,:)';
vN = vN_MTX(1,:)';


