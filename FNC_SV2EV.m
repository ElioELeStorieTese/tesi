function FNC_SV2EV(DR,INPUT_FILENAME,Nsamples,FILES_LABEL,dV_V_mmyr2deg,PLATE_CONTOUR_FILENAME,varargin)

% -
% SYNOPSIS
% FNC_SV2EV(DR,INPUT_FILENAME,Nsamples,FILES_LABEL,dV_V_mmyr2deg,PLATE_CONTOUR_FILENAME,varargin)
%
%
% DESCRIPTION 
% This function uses a set of input surface velocities to calculate an 
% Euler vector via linear regression (least square method). 
% 
% The problem is solved through a linear system w = A**-1 * b.
%
% DR is a string with the path to the directory where input data are 
% stored, and where results are also to be stored. DR is expected to 
% contain a file named as the string INPUT_FILENAME.

% INPUT_FILENAME is a string with the name of the txt file containing input
% velocities, located within DR. INPUT_FILENAME shall be formatted as 
% follows:
% Col. 1: longitude (degE) of the site.
% Col. 2: latitude (degN) of the site.
% Col. 3: height above sea level (m) of the site.
% Col. 4: East-velocity (mm/yr) at the site.
% Col. 5: North-velocity (mm/yr) at the site.
% Col. 6: Standard deviation of East-velocity (mm/yr) at the site.
% Col. 7: Standard deviation of North-velocity (mm/yr) at the site.
% Col. 8: Covariance between East and North velocities (mm/yr)^2.
%
% Nsamples is the number of samples to draw for each surface velocities.
% This means that Nsamples realisations of the Euler vectors are
% calculated, and then used to calculate a nominal (average) Euler vector
% and associated covariances.
%
% FILES_LABEL is a string to append to saved results/figure files, within DR.
%
% dV_V_mmyr2deg (1 x 2) contains two scaling factors for maps generated
% by the function. They are the number of longitudinal degrees 
% corresponding to 1 mm/yr in the maps of residuals and velocities,
% respectively.
%
% PLATE_CONTOUR_FILENAME is a string containing the name of the .txt file
% that contains the coordinates (lon-deg lat-deg) of the plate contour (or 
% any contour) to plot on the maps.
%
% varargin{1} is an optional input (1 x n) that contains the progressive 
% numbers of the sites (rows of SV_LLHENssc.txt) to be excluded from the 
% Euler vector calculation. 
%
% The resulting Euler vector is stored in DR in a file named 
% EV_'FILES_LABEL'.txt. Its format is:
% NaN NaN lon-degE lat-degN AV(deg/Myr) covariances (rad/Myr)^2 Nsamples
%
% A file named REPORT_'FILES_LABEL'.txt is also stored in DR. 
% Each row in it contains:
% 01. site-# 
% 02. NaN 
% 03. site-lon(degE)
% 04. site-lat(degN)
% 05. Observed East component of velocity vE (m/yr)
% 06. Observed North component of velocity vN (m/yr)
% 07. standard-deviation of vE (m/yr)
% 08. standard-deviation of vN (m/yr)
% 09. Residual (observation minus prediction) of vE, dvE(m/yr)
% 10. Residual of vN, dvN(m/yr)
% 11. standard-deviation of dvE (m/yr)
% 12. standard-deviation of dvN (m/yr)
% 13. covariance vET-vNT(m/yr)^2
% 14. covariance dvE-dvN(m/yr)^2
% If some sites have been excluded, there are also listed, but separated 
% from the used ones by row of NaNs.
% -

%--------------- PARAMETERS -----------------------------------------------

FILES_LBL = FILES_LABEL;

SWTC_MAP = 1;
deg_dv = dV_V_mmyr2deg(1);
deg_v  = dV_V_mmyr2deg(2);

eval(['cd ' DR])
eval(['data = load(''./' INPUT_FILENAME ''');'])
data = [[1:length(data(:,1))]' data];

eval(['BDR = load(''./' PLATE_CONTOUR_FILENAME ''');'])

RMVD = [];
USD = data;
if nargin == 7
    
    S2RMV = varargin{1};
    
    Col4Sel = 1;
    for i2 = 1:length(S2RMV)
    
        [r,c] = find(USD(:,Col4Sel)==S2RMV(i2));
        
        if ~isempty(r)
        
            RMVD = [RMVD;USD(r(1),:)];
            USD(r(1),:) = [];
            
        end
        
    end
        
    
end

if ~isempty(RMVD)
    SITES_REPORT = [ USD(:,1) NaN.*ones(size(USD(:,1))) USD(:,2:3) zeros(length(USD(:,1)),10) ;NaN.*ones(1,14);RMVD(:,1) NaN.*ones(size(RMVD(:,1))) RMVD(:,2:3)  zeros(length(RMVD(:,1)),10)];
    REPORT_LoLaHt = [USD(:,2:4) ; NaN.*ones(1,3) ; RMVD(:,2:4) ];
else
    SITES_REPORT = [ USD(:,1) NaN.*ones(size(USD(:,1))) USD(:,2:3) zeros(length(USD(:,1)),10)];
    REPORT_LoLaHt = data(:,2:4);
end
clear REPORT_LoLaHt

SITES_LLH = data(:,2:4);

SITES_SV = [data(:,5:8).*1e-3 data(:,9).*1e-6];


Re = 6371e3;
if length(data(1,:))==10
    NST = data(:,[1 10]);
else
    NST = data(:,1);
end
[x,y,z] = FNC_pLoLaHt2pXYZ_ELLIPSOID(SITES_LLH(:,1),SITES_LLH(:,2),SITES_LLH(:,3));%m

Lo_MTX = repmat(SITES_LLH(:,1),[1 Nsamples]);
La_MTX = repmat(SITES_LLH(:,2),[1 Nsamples]);

dvEN_res_vct = [round(log10(Nsamples))-1].*[10 10];

cmap = [[1 0 0];%  dSV/SV and 68% confidence
        [1 0 0];%  95% confidence
        0 0 0;%         Site number
        [0 0 0].*0;%  plate contour
        [0 0.7 0];%  removed site number
        [0 0.7 0];%  removed site dSV/SV
        [0 0 1];%       predicted SV
        [1 1 1].*0.8;%       coastline
        ];
%cmap = jet(8);
    
lwdt = [1;    %main 
        1;      %68% confidence
        0.5;    %95% confidence
        0.5;      %plate contour
        0.5;      %coastline
        0.5;      %Euler pole contours
        ];
    
mksz = [4];%main     
WHTL = 0.2;  
switch_Xnorm_XY_PDD = 0;
%--------------------------------------------------------------------------


%--------------- GENERATE ENSEMBLES OF VELOCITIES -------------------------
vEN = SITES_SV(:,1);
vNN = SITES_SV(:,2);
    
for i1 = 1:length(SITES_LLH(:,1))
    
    
    c11 = SITES_SV(i1,3).^2;
    c12 = SITES_SV(i1,5).*[SITES_SV(i1,3)*SITES_SV(i1,4)];
    c13 = 0;
    c22 = SITES_SV(i1,4).^2;
    c23 = 0;
    c33 = 1;
    
    CMTX = [ c11 c12;c12 c22];
    dv_ensemble = CORRELATED_ENSEMBLE_nD(CMTX,Nsamples);
    dvx = dv_ensemble(1,:);
    dvy = dv_ensemble(2,:);
    
    vE(i1,:) = vEN(i1,1) + dvx;
    vN(i1,:) = vNN(i1,1) + dvy;
        
end

[vx,vy,vz] = FNC_vEN2vXYZ_CONVERSION(Lo_MTX,La_MTX,vE,vN,zeros(size(vN)));

vEN_MTX = vE';
vNN_MTX = vN';    
%--------------------------------------------------------------------------




%--------------- INVERSION FOR EV ENSEMBLE FROM VELOCITIES ----------------
A11 = 1 * sum( [y(USD(:,1),:).^2 + z(USD(:,1),:).^2] );
A22 = 1 * sum( [x(USD(:,1),:).^2 + z(USD(:,1),:).^2] );
A33 = 1 * sum( [x(USD(:,1),:).^2 + y(USD(:,1),:).^2] );

A12 = -1 * sum( [x(USD(:,1),:).*y(USD(:,1),:)] );
A21 = A12;

A13 = -1 * sum( [x(USD(:,1),:).*z(USD(:,1),:)] );
A31 = A13;

A23 = -1 * sum( [y(USD(:,1),:).*z(USD(:,1),:)] );
A32 = A23;

A = [A11 A12 A13;A21 A22 A23;A31 A32 A33];
iA = inv(A);


bx = 1 .* sum( repmat(y(USD(:,1),:),[1 Nsamples]).*vz(USD(:,1),:) - repmat(z(USD(:,1),:),[1 Nsamples]).*vy(USD(:,1),:) );
by = 1 .* sum( repmat(z(USD(:,1),:),[1 Nsamples]).*vx(USD(:,1),:) - repmat(x(USD(:,1),:),[1 Nsamples]).*vz(USD(:,1),:) );
bz = 1 .* sum( repmat(x(USD(:,1),:),[1 Nsamples]).*vy(USD(:,1),:) - repmat(y(USD(:,1),:),[1 Nsamples]).*vx(USD(:,1),:) );

EV_ENSEMBLE(:,1) = [ iA(1,1).*bx + iA(1,2).*by + iA(1,3).*bz ]';%rad/yr
EV_ENSEMBLE(:,2) = [ iA(2,1).*bx + iA(2,2).*by + iA(2,3).*bz ]';
EV_ENSEMBLE(:,3) = [ iA(3,1).*bx + iA(3,2).*by + iA(3,3).*bz ]';
%--------------------------------------------------------------------------


%--------------- COVARIANCES AND SAVING TO TXT ----------------------------
xEV = EV_ENSEMBLE(:,1).*1e6;%rad/Myr
yEV = EV_ENSEMBLE(:,2).*1e6;
zEV = EV_ENSEMBLE(:,3).*1e6;
N = length(EV_ENSEMBLE(:,3));

[th,ph,av] = cart2sph(mean(xEV),mean(yEV),mean(zEV));
EVtxt(1:5) = [NaN NaN [th ph av].*(180/pi)];%deg/Myr

cMTX(1,1) = sum(xEV.*xEV)/N - sum(xEV)/N*sum(xEV)/N; 
cMTX(2,2) = sum(yEV.*yEV)/N - sum(yEV)/N*sum(yEV)/N; 
cMTX(3,3) = sum(zEV.*zEV)/N - sum(zEV)/N*sum(zEV)/N; 
%-XY    
cMTX(1,2) = sum(xEV.*yEV)/N - sum(xEV)/N*sum(yEV)/N;   
cMTX(2,1) = cMTX(1,2);
%-XZ    
cMTX(1,3) = sum(xEV.*zEV)/N - sum(xEV)/N*sum(zEV)/N;   
cMTX(3,1) = cMTX(1,3);
%-YZ
cMTX(2,3) = sum(yEV.*zEV)/N - sum(yEV)/N*sum(zEV)/N; 
cMTX(3,2) = cMTX(2,3);
    
EVtxt(6:12) = [cMTX(1,1) cMTX(1,2) cMTX(1,3) cMTX(2,2) cMTX(2,3) cMTX(3,3) Nsamples];
    
eval(['save -ascii EV_' FILES_LBL '.txt EVtxt'])    
%--------------------------------------------------------------------------



    
%--------------- RESIDUALS FROM EV ENSEMBLE -------------------------------
vxT_MTX = repmat(EV_ENSEMBLE(:,2),[1 length(NST(:,1))]).*repmat(z',[Nsamples 1]) - repmat(EV_ENSEMBLE(:,3),[1 length(NST(:,1))]).*repmat(y',[Nsamples 1]);
vyT_MTX = repmat(EV_ENSEMBLE(:,3),[1 length(NST(:,1))]).*repmat(x',[Nsamples 1]) - repmat(EV_ENSEMBLE(:,1),[1 length(NST(:,1))]).*repmat(z',[Nsamples 1]);
vzT_MTX = repmat(EV_ENSEMBLE(:,1),[1 length(NST(:,1))]).*repmat(y',[Nsamples 1]) - repmat(EV_ENSEMBLE(:,2),[1 length(NST(:,1))]).*repmat(x',[Nsamples 1]);

[vET_MTX,vNT_MTX,vUT_MTX] = FNC_vXYZ2vEN_CONVERSION(Lo_MTX',La_MTX',vxT_MTX,vyT_MTX,vzT_MTX);

% Here instead I calculate residuals as the difference between the ensemble 
% of theoretical velocities and the ensemble of observed velocities (i.e., 
% the one actually used in each single inversion for EV).
dvN_MTX = vNN_MTX - vNT_MTX;
dvE_MTX = vEN_MTX - vET_MTX;

%- Averages for report
for i1 = 1:length(USD(:,1))
    
    [r,c] = find(SITES_REPORT(:,1) == USD(i1,1));
    SITES_REPORT(r(1),5) = mean(vET_MTX(:,USD(i1,1)));
    SITES_REPORT(r(1),6) = mean(vNT_MTX(:,USD(i1,1)));
    SITES_REPORT(r(1),7) = std(vET_MTX(:,USD(i1,1)));
    SITES_REPORT(r(1),8) = std(vNT_MTX(:,USD(i1,1)));
    
    SITES_REPORT(r(1),9) = mean(dvE_MTX(:,USD(i1,1)));
    SITES_REPORT(r(1),10) = mean(dvN_MTX(:,USD(i1,1)));
    SITES_REPORT(r(1),11) = std(dvE_MTX(:,USD(i1,1)));
    SITES_REPORT(r(1),12) = std(dvN_MTX(:,USD(i1,1)));
    
    %Covariance vET - vNT
    SITES_REPORT(r(1),13) = sum(vET_MTX(:,USD(i1,1)).*vNT_MTX(:,USD(i1,1)))/N - sum(vET_MTX(:,USD(i1,1)))/N*sum(vNT_MTX(:,USD(i1,1)))/N;   
    
    %Covariance dvE - dvN
    SITES_REPORT(r(1),14) = sum(dvE_MTX(:,USD(i1,1)).*dvN_MTX(:,USD(i1,1)))/N - sum(dvE_MTX(:,USD(i1,1)))/N*sum(dvN_MTX(:,USD(i1,1)))/N;   
    
end

if ~isempty(RMVD)
    for i1 = 1:length(RMVD(:,1))
    
        [r,c] = find(SITES_REPORT(:,1) == RMVD(i1,1));
        SITES_REPORT(r(1),5) = mean(vET_MTX(:,RMVD(i1,1)));
        SITES_REPORT(r(1),6) = mean(vNT_MTX(:,RMVD(i1,1)));
        SITES_REPORT(r(1),7) = std(vET_MTX(:,RMVD(i1,1)));
        SITES_REPORT(r(1),8) = std(vNT_MTX(:,RMVD(i1,1)));
        
        SITES_REPORT(r(1),9) = mean(dvE_MTX(:,RMVD(i1,1)));
        SITES_REPORT(r(1),10) = mean(dvN_MTX(:,RMVD(i1,1)));
        SITES_REPORT(r(1),11) = std(dvE_MTX(:,RMVD(i1,1)));
        SITES_REPORT(r(1),12) = std(dvN_MTX(:,RMVD(i1,1)));
    
        %Covariance vET - vNT
        SITES_REPORT(r(1),13) = sum(vET_MTX(:,RMVD(i1,1)).*vNT_MTX(:,RMVD(i1,1)))/N - sum(vET_MTX(:,RMVD(i1,1)))/N*sum(vNT_MTX(:,RMVD(i1,1)))/N;   
    
        %Covariance dvE - dvN
        SITES_REPORT(r(1),14) = sum(dvE_MTX(:,RMVD(i1,1)).*dvN_MTX(:,RMVD(i1,1)))/N - sum(dvE_MTX(:,RMVD(i1,1)))/N*sum(dvN_MTX(:,RMVD(i1,1)))/N;   
    
    end
end
%-
eval(['save -ascii REPORT_' FILES_LBL '.txt SITES_REPORT']) 
%--------------------------------------------------------------------------




%------------ PLOTS -------------------------------------------------------
if SWTC_MAP==0
    
    %------------ PLOTTING RESIDUALS --------------------------------------
    figure
    hold on

    if ~isempty(BDR)
        plot(BDR(:,1),BDR(:,2),'color',cmap(4,:),'linewidth',lwdt(4))
    end
    
    ax=axis;
    dax = 0.05.*min(abs(diff(ax(1:2))),diff(ax(3:4)));
    axis([ax+dax.*[-1 1 -1 1]])
    ax=axis;
    pbaspect([abs(diff(ax(1:2)))/abs(diff(ax(3:4))) 1 1])

    map_f = deg_dv/1e-3;

    for i1 = 1:length(USD(:,1))
        
        xp = SITES_LLH(USD(i1,1),1)+map_f.*[0 mean(dvE_MTX(:,USD(i1,1)))];
        yp = SITES_LLH(USD(i1,1),2)+map_f.*[0 mean(dvN_MTX(:,USD(i1,1)))];
        plot(xp,yp,'color',cmap(1,:),'linewidth',lwdt(1))
        plot(xp(1),yp(1),'.','color',cmap(1,:))
    
        [mtx_x,mtx_y,mtx_r] = XY_ARRAY_PDD([dvE_MTX(:,USD(i1,1)) dvN_MTX(:,USD(i1,1))],dvEN_res_vct,switch_Xnorm_XY_PDD);
        ctr68 = XY_ARRAY_PDD_contour(mtx_x,mtx_y,mtx_r,68);
        ctr95 = XY_ARRAY_PDD_contour(mtx_x,mtx_y,mtx_r,95);
        ctr99 = XY_ARRAY_PDD_contour(mtx_x,mtx_y,mtx_r,99);
    
        %xp = SITES_LLH(USD(i1,1),1)+map_f.*ctr68(:,1);
        %yp = SITES_LLH(USD(i1,1),2)+map_f.*ctr68(:,2);
        %plot(xp,yp,'color',cmap(1,:),'linewidth',lwdt(2))
    
        xp = SITES_LLH(USD(i1,1),1)+map_f.*ctr95(:,1);
        yp = SITES_LLH(USD(i1,1),2)+map_f.*ctr95(:,2);
        plot(xp,yp,'color',cmap(2,:),'linewidth',lwdt(2))

        xp = SITES_LLH(USD(i1,1),1)+map_f.*ctr99(:,1);
        yp = SITES_LLH(USD(i1,1),2)+map_f.*ctr99(:,2);
        plot(xp,yp,'-','color',cmap(2,:),'linewidth',lwdt(3))
    
        if length(NST(1,:)) == 1
            LBL = [num2str(NST(USD(i1,1),1))];
        elseif length(NST(1,:)) == 2
            LBL = [num2str(NST(USD(i1,1),1)) ',t' num2str(NST(USD(i1,1),2))];
        end
        xp = SITES_LLH(USD(i1,1),1)-map_f.*0.15*mean(dvE_MTX(:,USD(i1,1)));
        yp = SITES_LLH(USD(i1,1),2)-map_f.*0.15*mean(dvN_MTX(:,USD(i1,1)));
        text(xp,yp,LBL,'color',cmap(3,:),'horizontalalignment','right')
        
    end

    box on
    grid off
    
    if ~isempty(RMVD)
        for i2=1:length(RMVD(:,1))
            plot3(RMVD(i2,2),RMVD(i2,3),-1,'.','color',cmap(5,:))
            
            if length(NST(1,:)) == 1
                LBL = [ num2str(NST(RMVD(i2,1),1)) ];
            else
                LBL = [ num2str(NST(RMVD(i2,1),1)) ',t' num2str(NST(RMVD(i2,1),2)) ];
            end
            text(RMVD(i2,2),RMVD(i2,3),-1,LBL,'color',cmap(5,:),'horizontalalignment','right')
        end
    end
        
        
    xSCL=ax(1)+2*dax;%0.1*diff(ax(1:2));
    ySCL=ax(3)+1*dax;%+0.1*diff(ax(3:4));
    ySCL2=ySCL+1.3*dax;%0.03*diff(ax(3:4));
    plot(xSCL+[0 deg_dv],ySCL.*[1 1],'color','k','linewidth',lwdt(1))
    LBL = ['1 mm/yr'];
    text(mean(xSCL+[0 deg_dv]),ySCL2,LBL,'color','k','horizontalalignment','center')

    ax=axis;
    axis(ax)
    
    load coastline
    plot(coastline(:,1),coastline(:,2),'-','color',cmap(8,:),'linewidth',lwdt(5))

    LBL = ['Site residuals (obs. minus pred.) (r: used - 95/99% confidence, g: unused)'];
    title(LBL)

    
    %- Add sites not used in the inversion
    if ~isempty(RMVD)
        for i1 = 1:length(RMVD(:,1))
    
            xp = SITES_LLH(RMVD(i1,1),1)+map_f.*[0 mean(dvE_MTX(:,RMVD(i1,1)))];
            yp = SITES_LLH(RMVD(i1,1),2)+map_f.*[0 mean(dvN_MTX(:,RMVD(i1,1)))];
            plot3(xp,yp,-1.*ones(size(xp)),'color',cmap(6,:),'linewidth',lwdt(3))
        
            [mtx_x,mtx_y,mtx_r] = XY_ARRAY_PDD([dvE_MTX(:,RMVD(i1,1)) dvN_MTX(:,RMVD(i1,1))],dvEN_res_vct,switch_Xnorm_XY_PDD);
            ctr68 = XY_ARRAY_PDD_contour(mtx_x,mtx_y,mtx_r,68);
            ctr95 = XY_ARRAY_PDD_contour(mtx_x,mtx_y,mtx_r,95);
    
            xp = SITES_LLH(RMVD(i1,1),1)+map_f.*ctr68(:,1);
            yp = SITES_LLH(RMVD(i1,1),2)+map_f.*ctr68(:,2);
            plot3(xp,yp,-1.*ones(size(xp)),'color',cmap(6,:),'linewidth',lwdt(2))
    
            xp = SITES_LLH(RMVD(i1,1),1)+map_f.*ctr95(:,1);
            yp = SITES_LLH(RMVD(i1,1),2)+map_f.*ctr95(:,2);
            plot3(xp,yp,-1.*ones(size(xp)),'color',cmap(6,:),'linewidth',lwdt(3))
    
        end
    end
    eval(['print -r200 -dpng SITE_RESIDUALS_' FILES_LBL ])
    %-
    
    %----------------------------------------------------------------------
    
    
    %------------ PLOTTING VELOCITIES --------------------------------------
    figure
    hold on

    if ~isempty(BDR)
        plot(BDR(:,1),BDR(:,2),'color',cmap(4,:),'linewidth',lwdt(4))
    end
    
    ax=axis;
    dax = 0.05.*min(abs(diff(ax(1:2))),diff(ax(3:4)));
    axis([ax+dax.*[-1 1 -1 1]])
    ax=axis;
    pbaspect([abs(diff(ax(1:2)))/abs(diff(ax(3:4))) 1 1])

    map_f = deg_v/1e-3;

    for i1 = 1:length(USD(:,1))
    
        xp = SITES_LLH(USD(i1,1),1)+map_f.*[0 mean(vET_MTX(:,USD(i1,1)))];
        yp = SITES_LLH(USD(i1,1),2)+map_f.*[0 mean(vNT_MTX(:,USD(i1,1)))];
        plot3(xp,yp,ones(size(yp)),'color',cmap(7,:),'linewidth',lwdt(1))
        
        xp = SITES_LLH(USD(i1,1),1)+map_f.*[0 mean(vEN_MTX(:,USD(i1,1)))];
        yp = SITES_LLH(USD(i1,1),2)+map_f.*[0 mean(vNN_MTX(:,USD(i1,1)))];
        plot(xp,yp,'color',cmap(1,:),'linewidth',lwdt(1))
        plot(xp(1),yp(1),'.','color',cmap(1,:))
    
        [mtx_x,mtx_y,mtx_r] = XY_ARRAY_PDD([vEN_MTX(:,USD(i1,1)) vNN_MTX(:,USD(i1,1))],dvEN_res_vct,switch_Xnorm_XY_PDD);
        ctr68 = XY_ARRAY_PDD_contour(mtx_x,mtx_y,mtx_r,68);
        ctr95 = XY_ARRAY_PDD_contour(mtx_x,mtx_y,mtx_r,95);
    
        xp = SITES_LLH(USD(i1,1),1)+map_f.*ctr68(:,1);
        yp = SITES_LLH(USD(i1,1),2)+map_f.*ctr68(:,2);
        plot(xp,yp,'color',cmap(1,:),'linewidth',lwdt(2))
    
        xp = SITES_LLH(USD(i1,1),1)+map_f.*ctr95(:,1);
        yp = SITES_LLH(USD(i1,1),2)+map_f.*ctr95(:,2);
        plot(xp,yp,'color',cmap(2,:),'linewidth',lwdt(3))
    
    
        if length(NST(1,:)) == 1
            LBL = [num2str(NST(USD(i1,1),1))];
        elseif length(NST(1,:)) == 2
            LBL = [num2str(NST(USD(i1,1),1)) ',t' num2str(NST(USD(i1,1),2))];
        end
        xp = SITES_LLH(USD(i1,1),1)-map_f.*0.1*mean(vEN_MTX(:,USD(i1,1)));
        yp = SITES_LLH(USD(i1,1),2)-map_f.*0.1*mean(vNN_MTX(:,USD(i1,1)));
        text(xp,yp,LBL,'color',cmap(3,:),'horizontalalignment','right')
    
    end

    box on
    grid off
    
    if ~isempty(RMVD)
        for i2=1:length(RMVD(:,1))
            plot3(RMVD(i2,2),RMVD(i2,3),-1,'.','color',cmap(5,:))
            
            if length(NST(1,:)) == 1
                LBL = [ num2str(NST(RMVD(i2,1),1)) ];
            else
                LBL = [ num2str(NST(RMVD(i2,1),1)) ',t' num2str(NST(RMVD(i2,1),2)) ];
            end
            text(RMVD(i2,2),RMVD(i2,3),-1,LBL,'color',cmap(5,:),'horizontalalignment','right')
        end
    end
        
        
    
    xSCL=ax(1)+2*dax;%0.1*diff(ax(1:2));
    ySCL=ax(3)+1*dax;%+0.1*diff(ax(3:4));
    ySCL2=ySCL+1.3*dax;%0.03*diff(ax(3:4));
    plot(xSCL+[0 5*deg_v],ySCL.*[1 1],'color','k','linewidth',lwdt(1))
    LBL = ['5 mm/yr'];
    text(mean(xSCL+[0 deg_v]),ySCL2,LBL,'color','k','horizontalalignment','center')

    ax=axis;
    axis(ax)
    
    load coastline
    plot(coastline(:,1),coastline(:,2),'-','color',cmap(8,:),'linewidth',lwdt(5))

    LBL = ['Site velocities (b: predicted, r: observed with 68/95% confidence, g: unused)'];
    title(LBL)

    
    %- Add sites not used in the inversion
    if ~isempty(RMVD)
        for i1 = 1:length(RMVD(:,1))%length(dvE_MTX(1,:))
    
            xp = SITES_LLH(RMVD(i1,1),1)+map_f.*[0 mean(vET_MTX(:,RMVD(i1,1)))];
            yp = SITES_LLH(RMVD(i1,1),2)+map_f.*[0 mean(vNT_MTX(:,RMVD(i1,1)))];
            plot3(xp,yp,-1.*ones(size(xp)),'--','color',cmap(6,:),'linewidth',lwdt(3))
        
            
            xp = SITES_LLH(RMVD(i1,1),1)+map_f.*[0 mean(vEN_MTX(:,RMVD(i1,1)))];
            yp = SITES_LLH(RMVD(i1,1),2)+map_f.*[0 mean(vNN_MTX(:,RMVD(i1,1)))];
            plot3(xp,yp,-1.*ones(size(xp)),'color',cmap(6,:),'linewidth',lwdt(3))
        
            [mtx_x,mtx_y,mtx_r] = XY_ARRAY_PDD([vEN_MTX(:,RMVD(i1,1)) vNN_MTX(:,RMVD(i1,1))],dvEN_res_vct,switch_Xnorm_XY_PDD);
            ctr68 = XY_ARRAY_PDD_contour(mtx_x,mtx_y,mtx_r,68);
            ctr95 = XY_ARRAY_PDD_contour(mtx_x,mtx_y,mtx_r,95);
    
            xp = SITES_LLH(RMVD(i1,1),1)+map_f.*ctr68(:,1);
            yp = SITES_LLH(RMVD(i1,1),2)+map_f.*ctr68(:,2);
            plot3(xp,yp,-1.*ones(size(xp)),'color',cmap(6,:),'linewidth',lwdt(2))
    
            xp = SITES_LLH(RMVD(i1,1),1)+map_f.*ctr95(:,1);
            yp = SITES_LLH(RMVD(i1,1),2)+map_f.*ctr95(:,2);
            plot3(xp,yp,-1.*ones(size(xp)),'color',cmap(6,:),'linewidth',lwdt(3))
    
        end
    end
    eval(['print -r200 -dpng SITE_VELOCITIES_' FILES_LBL ])
    %-
    
    close all
    %----------------------------------------------------------------------
    





    %------------ PLOTTING EULER VECTOR -----------------------------------
    [th,ph,av] = cart2sph(mean(EV_ENSEMBLE(:,1)),mean(EV_ENSEMBLE(:,2)),mean(EV_ENSEMBLE(:,3)));
    xpM = th.*(180/pi);
    ypM = ph.*(180/pi);
    avM = av.*(180/pi*1e6);%rad/yr > deg/Myr

    [thE,phE,av] = cart2sph(EV_ENSEMBLE(:,1),EV_ENSEMBLE(:,2),EV_ENSEMBLE(:,3));
    xp = thE.*(180/pi);
    yp = phE.*(180/pi);
    avE = av.*(180/pi*1e6);%rad/yr > deg/Myr

    cov_MTX(1,1) = std(xp)^2;
    cov_MTX(2,2) = std(yp)^2;
    
    cov_EN = sum(xp.*yp)/Nsamples - sum(xp)/Nsamples*sum(yp)/Nsamples;   
    cov_MTX(1,2) = cov_EN;
    cov_MTX(2,1) = cov_EN;
    
    [E_1S,E_2S,E_3S] = COV2ELLIPSE_2D(cov_MTX);
    
    xp = xpM + E_1S(:,1);
    yp = ypM + E_1S(:,2);
    CTR_68p = [xp yp];
    
    xp = xpM + E_2S(:,1);
    yp = ypM + E_2S(:,2);
    CTR_95p = [xp yp];
    
    xp = xpM + E_3S(:,1);
    yp = ypM + E_3S(:,2);
    CTR_99p = [xp yp];
    
    
    figure
    
    subplot(2,1,1)
    LBL = ['Distributions of Euler pole (top, 68/95/99% confidence are shown) and angular velocity (bottom)'];
    title(LBL)
    hold on
    box on
    plot(CTR_68p(:,1),CTR_68p(:,2),'-','color',cmap(7,:),'linewidth',lwdt(6))
    plot(CTR_95p(:,1),CTR_95p(:,2),'--','color',cmap(7,:),'linewidth',lwdt(6))
    plot(CTR_99p(:,1),CTR_99p(:,2),':','color',cmap(7,:),'linewidth',lwdt(6))
    %plot(xpM,ypM,'o','color',cmap(7,:),'markerface',cmap(7,:))
    
    if ~isempty(BDR)
        plot(BDR(:,1),BDR(:,2),'color',cmap(4,:),'linewidth',lwdt(4))
    end
    load coastline
    plot(coastline(:,1),coastline(:,2),'-','color',cmap(8,:),'linewidth',lwdt(5))
    
    ax(1) = min(min(BDR(:,1)),min(CTR_99p(:,1)))-5;
    ax(2) = max(max(BDR(:,1)),max(CTR_99p(:,1)))+5;
    ax(3) = min(min(BDR(:,2)),min(CTR_99p(:,2)))-5;
    ax(4) = max(max(BDR(:,2)),max(CTR_99p(:,2)))+5;
    axis(ax);
    pbaspect([abs(diff(ax(1:2)))/abs(diff(ax(3:4))) 1 1])

    subplot(2,1,2)
    hold on
    box on
    histogram(avE)
    xlabel('[deg/Myr]')
    ylabel('Recurrence')
    pbaspect([abs(diff(ax(1:2)))/abs(diff(ax(3:4))) 1 1])
    eval(['print -r200 -dpng EULER_VECTOR_' FILES_LBL ])
    close all
    %----------------------------------------------------------------------
    

end

