function FNC_SV_PLOT(DR,INPUT_FILENAME,V_mmyr2deg,PLATE_CONTOUR_FILENAME)

% -
% SYNOPSIS
% FNC_SV_PLOT(DR,INPUT_FILENAME,V_mmyr2deg,PLATE_CONTOUR_FILENAME)
%
%
% DESCRIPTION 
% This function generates a plot of a set of input surface velocities. 
% 
% DR is a string with the path to the directory where input data are 
% stored. 
%
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
% V_mmyr2deg (scalar) is a scaling factor for the map generated
% by the function. It is the number of longitudinal degrees 
% corresponding to 1 mm/yr in the map of site velocities.
%
% PLATE_CONTOUR_FILENAME is a string containing the name of the .txt file
% that contains the coordinates (lon-deg lat-deg) of the plate contour (or 
% any contour) to plot on the maps.
%
% The generated plot is also saved within DR as SITE_VELOCITIES.png
% -

%--------------- PARAMETERS -----------------------------------------------


Nsamples = 1e5;

SWTC_MAP = 1;
deg_v  = V_mmyr2deg(1);


eval(['cd ' DR])
eval(['data = load(''./' INPUT_FILENAME ''');'])
data = [[1:length(data(:,1))]' data];

eval(['BDR = load(''./' PLATE_CONTOUR_FILENAME ''');'])

RMVD = [];
USD = data;

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

cmap = [[0 0 1];%  SV and 68% confidence
        [0 0 1];%  95% confidence
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




%------------ PLOTS -------------------------------------------------------
if SWTC_MAP==1
    
    %------------ PLOTTING VELOCITIES --------------------------------------
    figure
    hold on

    box on
    grid off
    
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

    
    xSCL=ax(1)+2*dax;%0.1*diff(ax(1:2));
    ySCL=ax(3)+1*dax;%+0.1*diff(ax(3:4));
    ySCL2=ySCL+1.3*dax;%0.03*diff(ax(3:4));
    plot(xSCL+[0 5*deg_v],ySCL.*[1 1],'color','k','linewidth',lwdt(1))
    LBL = ['5 mm/yr'];
    text(mean(xSCL+[0 deg_v]),ySCL2,LBL,'color','k','horizontalalignment','center')

    load coastline
    plot(coastline(:,1),coastline(:,2),'-','color',cmap(8,:),'linewidth',lwdt(5))

    LBL = ['Site velocities (68/95% confidence are also shown)'];
    title(LBL)

    
    eval(['print -r200 -dpng SITE_VELOCITIES'])
    %-
    
    %----------------------------------------------------------------------
    

end

