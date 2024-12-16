function OVL_2EVs = FNC_OVL_2EVs(EV_A,EV_B,N_SAMPLES,LoLa_RES,N_NORM_BINS)

% OVL_2EVs = FNC_OVL_2EVs(EV_A,EV_B,N_SAMPLES,LoLa_RES,N_NORM_BINS)
%
% This function calculates the overlap between the probability density
% functions (PDFs) of two Euler vectors in input.
%
% The two Euler vectors in input (EV_A and EV_B) are formatted as follows:
% Column 1: initial age of the stage [Ma], or NaN.
% Column 2: final age of the stage [Ma], or NaN.
% Column 3: longitude (deg East) of the stage Euler pole.
% Column 4: latitude (deg North) of the stage Euler pole.
% Column 5: stage angular velocity (deg/Myr).
% Columns 6-11: Covariance entries (xx xy xz yy yz zz) in (rad/Myr)^2.
% Note that stages do not have to necessarly represent a time
% series, hence the option of setting not-a-numbers (NaN) for the initial
% and final ages of the stage. This is useful, for instance, if one wants
% to simply compare two or more Euler vectors that do not necessarly 
% represent a time series.
%
% N_SAMPLES is the number of samples to draw in order to reconstruct the
% PDF of each Euler vector. Set this to no less than 1e5 for sufficiently
% accurate sampling.
%
% The calculation is performed in spherical coordinates. Therefore,
% LoLa_RES (1 x 1, deg) is the step to use when binning the 2D Lon/Lat 
% space of Euler poles. Instead, N_NORM_BINS is the number of bins to use
% when binning angular velocities, between minimum and maximum values 
% observed across the two ensembles of EV_ and EV_B.
%
% OVL_2EVs is the overlap value obtained from the two ensembles.


%- Sampling 
EVs = [EV_A;EV_B];
Nsampled = 0;
    
[x,y,z] = sph2cart(EVs(:,3).*(pi/180),EVs(:,4).*(pi/180),EVs(:,5));
EVc = [x y z];
        
EVcvr = EVs(:,6:11).*(180/pi).^2;

for i2 = 1:length(EVc(:,1))
        
    cov_mtx_3D = [EVcvr(i2,[1 2 3]);EVcvr(i2,[2 4 5]);EVcvr(i2,[3 5 6])];
    ndata = CORRELATED_ENSEMBLE_nD(cov_mtx_3D,N_SAMPLES);
    dx = ndata(1,:);
    dy = ndata(2,:);
    dz = ndata(3,:);
    

    EV2add(1,:) = EVc(i2,1)+real(dx)';
    EV2add(2,:) = EVc(i2,2)+real(dy)';
    EV2add(3,:) = EVc(i2,3)+real(dz)';
       
    eval(['EVc_STG' num2str(i2) ' = EV2add'';';])
        
    clear EV2add* cov_mtx_3D dx dy dz
        
end
%-
  

ENSEMBLE_A = EVc_STG1;
ENSEMBLE_B = EVc_STG2;


%- XYZ > SPH
[th,ph,N_A] = cart2sph(ENSEMBLE_A(:,1),ENSEMBLE_A(:,2),ENSEMBLE_A(:,3));
Lo_A = th.*(180/pi);
La_A = ph.*(180/pi);

[th,ph,N_B] = cart2sph(ENSEMBLE_B(:,1),ENSEMBLE_B(:,2),ENSEMBLE_B(:,3));
Lo_B = th.*(180/pi);
La_B = ph.*(180/pi);
%- 


%- Pole
N1 = 180/LoLa_RES;
N_LoLa = 180/LoLa_RES*360/LoLa_RES;

IDX_LoLa_A = N1*floor((Lo_A+180)./LoLa_RES)+ceil((La_A+90)./LoLa_RES);
IDX_LoLa_B = N1*floor((Lo_B+180)./LoLa_RES)+ceil((La_B+90)./LoLa_RES);
%- 

%- Norm
Mn = 0.99.*min([N_A;N_B]);
Mx = 1.01.*max([N_A;N_B]);
dN = (Mx-Mn)/(N_NORM_BINS);

IDX_NORM_A = ceil((N_A-Mn)./dN);
IDX_NORM_B = ceil((N_B-Mn)./dN);
%-

%- ENTIRE ENSEMBLE
IDX_A = IDX_LoLa_A + N_LoLa.*(IDX_NORM_A-1);
IDX_B = IDX_LoLa_B + N_LoLa.*(IDX_NORM_B-1);

N_LLN = N_LoLa*(N_NORM_BINS);

EDGES = [1:N_LLN+1]-0.5;
[HC_A,EDGES_A] = histcounts(IDX_A,EDGES);
[HC_B,EDGES_B] = histcounts(IDX_B,EDGES);

PDF_A = HC_A./length(ENSEMBLE_A(:,1));
PDF_B = HC_B./length(ENSEMBLE_B(:,1));

OVL_2EVs = sum(min([PDF_A;PDF_B]));
%-






