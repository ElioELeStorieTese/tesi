clear, clc, close all;

if 1<10
% % POOL PER CALCOLO PARALLELO
% if isempty(gcp('nocreate'))
%         parpool; % Lavoratori paralleli
%     end
% \\\\\\\\\\ CARICAMENTO DATI //////////
earthRadius = 6371;
% PARAMETRI D'INDAGINE: stazioni, rumore (sensibilità dei dati), deviazione standard (varianza dei dati)
% Stazioni
RN_STAZIONI = [ 5 20;
                20 100;
                100 500;
                500 1000;
                ];
% Deviazione standard
stdDevCoeff = 0.6:0.1:1.5;
% Rumore
noiseRange = 0.2:0.2:2; %da 0.2 a 2, passi di 0.2mm/yr
% Vettori di Eulero
N_samples_vettore = 500;
N_samples_overlap = 1e5;
LoLa_RES = 0.5;         % Risoluzione coordinate polari
N_NORM_BINS = 20;       % Numero di bin per la normalizzazione
DR='C:/Users/utente/Documents/MATLAB_tesi';
% DR='~/Documents/WORK/TEACHING/TEACHING_AT_UNIPR/00_TESI_LAUREA/2024_ELIA_ROSSETTI/MATLAB/';
% Matrici vettori stazioni fittizie
outputDir = DR;
results_AD = cell(length(noiseRange), length(stdDevCoeff), length(RN_STAZIONI(:,1)));
results_AP = cell(length(noiseRange), length(stdDevCoeff), length(RN_STAZIONI(:,2)));
% Matrici overlap
overlapResults_AD = cell(size(results_AD));
overlapResults_AP = cell(size(results_AP));

% Caricamento vettori Eulero reali
EV_AD = load('EV_AD_0601_0812.txt');
EV_AP = load('EV_AP_1701_1909.txt');
% Caricamento dati placche
AD = load('PLATE_CONTOUR_ADRIA.txt');
AP = load('PLATE_CONTOUR_APULIA.txt');
lonAD = AD(:, 1);
latAD = AD(:, 2);
lonAP = AP(:, 1);
latAP = AP(:, 2);

% Caricamento dati città IN ADRIA
cityfileID_AD = fopen('Città Italiane Adria.txt', 'r');
data = textscan(cityfileID_AD, '%s %f %f %d %f %f', 'Delimiter', {' ', '\t'}, 'HeaderLines', 0);
fclose(cityfileID_AD);
cityID_AD = data{1};
cityLon_AD = data{2};
cityLat_AD = data{3};
cityPopul_AD = data{4};
cityArea_AD = data{5};
cityAltitude_AD = data{6};
% Caricamento dati città IN APULIA
cityfileID_AP = fopen('Città Italiane Apulia.txt', 'r');
data = textscan(cityfileID_AP, '%s %f %f %d %f %f', 'Delimiter', {' ', '\t'}, 'HeaderLines', 0);
fclose(cityfileID_AP);
cityID_AP = data{1};
cityLon_AP = data{2};
cityLat_AP = data{3};
cityPopul_AP = data{4};
cityArea_AP = data{5};
cityAltitude_AP = data{6};
% Caricamento dati città FUORI DALLE PLACCHE
cityfileID_fuori = fopen('Città Italiane esterne.txt', 'r');
data = textscan(cityfileID_fuori, '%s %f %f %d %f %f', 'Delimiter', {' ', '\t'}, 'HeaderLines', 0);
fclose(cityfileID_fuori);
cityID_fuori = data{1};
cityLon_fuori = data{2};
cityLat_fuori = data{3};
cityPopul_fuori = data{4};
cityArea_fuori = data{5};
cityAltitude_fuori = data{6};
% Caricamento dati coste
load coastlines;


%---------------------------
if 1<0
% PREGENERAZIONE DI STAZIONI CASUALI
numCities_AD = length(cityID_AD);
numCities_AP = length(cityID_AP);
numCities_fuori = length(cityID_fuori);
% Numero di iterazioni per stationIdx
StationIdx = length(numStazRange);
% Matrice per i numeri di stazioni (city x stationIdx)
numStazioniMatrix_AD = zeros(numCities_AD, StationIdx);
numStazioniMatrix_AP = zeros(numCities_AP, StationIdx);
numStazioniMatrix_fuori = zeros(numCities_fuori, StationIdx);
% Generazione numeri di stazioni unici per ogni città e stationIdx
for i = 1:numCities_AD
    numStazioniMatrix_AD(i, :) = randperm(maxStazioni - minStazioni + 1, StationIdx) + minStazioni - 1;  % Numeri casuali ma univoci unici tra min e max
end
for i = 1:numCities_AP
    numStazioniMatrix_AP(i, :) = randperm(maxStazioni - minStazioni + 1, StationIdx) + minStazioni - 1; 
end
for i = 1:numCities_fuori
    numStazioniMatrix_fuori(i, :) = randperm(maxStazioni - minStazioni + 1, StationIdx) + minStazioni - 1; 
end
end
%---------------------------

% checkpoint ciclo
totalIterations = length(noiseRange) * length(stdDevCoeff) * length(RN_STAZIONI); % Numero totale di iterazioni
checkpoint = 0.1; % ogni % di progresso sul totale
nextCheckpoint = checkpoint * totalIterations;

% \\\\\\\\\\ SEZIONE CALCOLI //////////
for noiseIdx = 1:length(noiseRange)
    noiseI = noiseRange(noiseIdx);
    for stdDevIdx = 1:length(stdDevCoeff)
        stdDevJ = stdDevCoeff(stdDevIdx);
        for stationIdx = 1:length(RN_STAZIONI)
            fprintf('Elaborazione sequenziale: Noise = %.1f, StdDev = %.1f, StationIdx = %d\n', noiseI, stdDevJ, stationIdx);
            stazLat_AD = [];
            stazLon_AD = [];
            cityIndex_AD = [];
            stazLat_AP = [];
            stazLon_AP = [];
            cityIndex_AP = [];
            stazLat_fuori = [];
            stazLon_fuori = [];
            cityIndex_fuori = [];
            %----- GENERAZIONE STAZIONI GNSS -----
            % GENERAZIONE STAZIONI GNSS ADRIA
            numStazioni_AD = randi(RN_STAZIONI(stationIdx,:),[length(cityID_AD) 1]);
            for i = 1:length(cityID_AD)
                %numStazioni_AD = numStazioniMatrix_AD(i, stationIdx);
                radius_AD = sqrt(cityArea_AD(i) / pi); % Raggio
                angles_AD = 2 * pi * rand(numStazioni_AD(i), 1); % Angoli casuali
                puntoStaz_AD = radius_AD * sqrt(rand(numStazioni_AD(i), 1)); % Distanza casuale entro il raggio
                deltaLat_AD = (puntoStaz_AD / earthRadius) * (180 / pi); % Variazione latitudine dalle coordinate d'origine
                deltaLon_AD = (puntoStaz_AD / earthRadius) .* (180 / pi) ./ cosd(cityLat_AD(i)); % Variazione longitudine dalle coordinate d'origine
                puntoStazLat_AD = cityLat_AD(i) + deltaLat_AD .* cos(angles_AD); % Latitudine stazioni
                puntoStazLon_AD = cityLon_AD(i) + deltaLon_AD .* sin(angles_AD); % Longitudine stazioni
                % Accumulo risultati
                stazLat_AD = [stazLat_AD; puntoStazLat_AD];
                stazLon_AD = [stazLon_AD; puntoStazLon_AD];
                cityIndex_AD = [cityIndex_AD; repmat(i, numStazioni_AD(i), 1)]; % ID città associato
            end
            % GENERAZIONE STAZIONI GNSS APULIA
            numStazioni_AP = randi(RN_STAZIONI(stationIdx,:),[length(cityID_AP) 1]);
            for i = 1:length(cityID_AP)
                %numStazioni_AP = numStazioniMatrix_AP(i, stationIdx);
                radius_AP = sqrt(cityArea_AP(i) / pi); % Raggio
                angles_AP = 2 * pi * rand(numStazioni_AP(i), 1); % Angoli casuali
                puntoStaz_AP = radius_AP * sqrt(rand(numStazioni_AP(i), 1)); % Distanza casuale entro il raggio
                deltaLat_AP = (puntoStaz_AP / earthRadius) * (180 / pi); % Variazione latitudine dalle coordinate d'origine
                deltaLon_AP = (puntoStaz_AP / earthRadius) .* (180 / pi) ./ cosd(cityLat_AP(i)); % Variazione longitudine dalle coordinate d'origine
                puntoStazLat_AP = cityLat_AP(i) + deltaLat_AP .* cos(angles_AP); % Latitudine stazioni
                puntoStazLon_AP = cityLon_AP(i) + deltaLon_AP .* sin(angles_AP); % Longitudine stazioni
                % Accumulo risultati
                stazLat_AP = [stazLat_AP; puntoStazLat_AP];
                stazLon_AP = [stazLon_AP; puntoStazLon_AP];
                cityIndex_AP = [cityIndex_AP; repmat(i, numStazioni_AP(i), 1)]; % ID città associato
            end
            % GENERAZIONE STAZIONI GNSS ESTERNE
            numStazioni_fuori = randi(RN_STAZIONI(stationIdx,:),[length(cityID_fuori) 1]);
            for i = 1:length(cityID_fuori)
                %numStazioni_fuori = numStazioniMatrix_fuori(i, stationIdx);
                radius_fuori = sqrt(cityArea_fuori(i) / pi); % Raggio
                angles_fuori = 2 * pi * rand(numStazioni_fuori(i), 1); % Angoli casuali
                puntoStaz_fuori = radius_fuori * sqrt(rand(numStazioni_fuori(i), 1)); % Distanza casuale entro il raggio
                deltaLat_fuori = (puntoStaz_fuori / earthRadius) * (180 / pi); % Variazione latitudine dalle coordinate d'origine
                deltaLon_fuori = (puntoStaz_fuori / earthRadius) .* (180 / pi) ./ cosd(cityLat_fuori(i)); % Variazione longitudine dalle coordinate d'origine
                puntoStazLat_fuori = cityLat_fuori(i) + deltaLat_fuori .* cos(angles_fuori); % Latitudine stazioni
                puntoStazLon_fuori = cityLon_fuori(i) + deltaLon_fuori .* sin(angles_fuori); % Longitudine stazioni
                % Accumulo risultati
                stazLat_fuori = [stazLat_fuori; puntoStazLat_fuori];
                stazLon_fuori = [stazLon_fuori; puntoStazLon_fuori];
                cityIndex_fuori = [cityIndex_fuori; repmat(i, numStazioni_fuori(i), 1)]; % ID città associato
            end

            fprintf('Stazioni totali: %d\n', length(stazLat_fuori)+length(stazLat_AP)+length(stazLat_AD));
            fprintf('Stazioni in Adria: %d\n', length(stazLat_AD));
            fprintf('Stazioni in Apulia: %d\n', length(stazLat_AP));
            
            % ----- CALCOLO DEVIAZIONE STANDARD, VELOCITà, RUMORE -----
            [vE_AD, vN_AD] = FNC_EV2siteSV(EV_AD, [stazLon_AD, stazLat_AD]);
            [vE_AP, vN_AP] = FNC_EV2siteSV(EV_AP, [stazLon_AP, stazLat_AP]);
            stdDev_AD = (0.5 + (1.5 - 0.5) * rand(size(vE_AD))) * stdDevJ;
            stdDev_AP = (0.5 + (1.5 - 0.5) * rand(size(vE_AP))) * stdDevJ;
            stdDevangle_AD = 2 * pi * rand(size(vE_AD));
            stdDevangle_AP = 2 * pi * rand(size(vE_AP));
            stdDevE_AD = stdDev_AD .* cos(stdDevangle_AD);
            stdDevN_AD = stdDev_AD .* sin(stdDevangle_AD);
            stdDevE_AP = stdDev_AP .* cos(stdDevangle_AP);
            stdDevN_AP = stdDev_AP .* sin(stdDevangle_AP);

            vE_AD_noise = vE_AD + (noiseI * rand(size(vE_AD)));
            vN_AD_noise = vN_AD + (noiseI * rand(size(vE_AD)));
            vE_AP_noise = vE_AP + (noiseI * rand(size(vE_AP)));
            vN_AP_noise = vN_AP + (noiseI * rand(size(vE_AP)));
            file_path_vel_AD = sprintf('SV_AD_Noise%.1f_StDev%.1f_Staz%d.txt', noiseI, stdDevJ, length(stazLat_AD));
            file_path_vel_AP = sprintf('SV_AP_Noise%.1f_StDev%.1f_Staz%d.txt', noiseI, stdDevJ, length(stazLat_AP));

            fileID_AD = fopen(file_path_vel_AD, 'w');
            for i = 1:size(stazLat_AD, 1)
                cityID_AD_VEL = cityIndex_AD(i);
                currentcityAltitude_AD = cityAltitude_AD(cityID_AD_VEL);
                fprintf(fileID_AD, '%.4f %.4f %.1f %1.7e %1.7e %1.7e %1.7e %1.7e\n', stazLon_AD(i), stazLat_AD(i), currentcityAltitude_AD, vE_AD_noise(i), vN_AD_noise(i), stdDevE_AD(i), stdDevN_AD(i), 0.0);
            end
            fclose(fileID_AD);

            fileID_AP = fopen(file_path_vel_AP, 'w');
            for i = 1:size(stazLat_AP, 1)
                cityID_AP_VEL = cityIndex_AP(i);
                currentcityAltitude_AP = cityAltitude_AP(cityID_AP_VEL);
                fprintf(fileID_AP, '%.4f %.4f %.1f %1.7e %1.7e %1.7e %1.7e %1.7e\n', stazLon_AP(i), stazLat_AP(i), currentcityAltitude_AP, vE_AP_noise(i), vN_AP_noise(i), stdDevE_AP(i), stdDevN_AP(i), 0.0);
            end
            fclose(fileID_AP);

            % ----- GENERAZIONE VETTORE DI EULERO -----
            EV_AD_result = FNC_SV2EV_C(DR,file_path_vel_AD,N_samples_vettore);
            EV_AP_result = FNC_SV2EV_C(DR,file_path_vel_AP,N_samples_vettore);
            
            % Salvataggio matrici risultati intermedi
            results_AD{noiseIdx, stdDevIdx, stationIdx} = struct('EV_AD', EV_AD_result, 'Noise', noiseI, 'StdDev', stdDevJ, 'Stations', length(stazLat_AD));
            results_AP{noiseIdx, stdDevIdx, stationIdx} = struct('EV_AP', EV_AP_result, 'Noise', noiseI, 'StdDev', stdDevJ, 'Stations', length(stazLat_AP));

            % Checkpoint
            currentIteration = (noiseIdx - 1) * length(stdDevCoeff) * length(RN_STAZIONI) + (stdDevIdx - 1) * length(RN_STAZIONI) + stationIdx;
            if currentIteration >= nextCheckpoint
                fprintf('Progress: %.1f%% completato\n', 100 * currentIteration / totalIterations);
                nextCheckpoint = nextCheckpoint + checkpoint * totalIterations;
                % Salvataggio dati intermedi
                fprintf('Salvataggio dati intermedi...\n');
                save(fullfile(outputDir, 'results_AD_intermediate.mat'), 'results_AD');
                save(fullfile(outputDir, 'results_AP_intermediate.mat'), 'results_AP');
            end
        end
    end
end

% Salvataggio matrici risultati finale
fprintf('Salvataggio finale dei dati\n');
save(fullfile(outputDir, 'results_AD_finale.mat'), 'results_AD', '-v7.3'); %
save(fullfile(outputDir, 'results_AP_finale.mat'), 'results_AP', '-v7.3');
fprintf('Salvataggio matrici risultati completato.\n');

% \\\\\\\\\\ CALCOLO OVERLAP //////////
fprintf('Calcolo degli overlap per tutti i vettori di Eulero...\n');
% Matrici per grafico
overlap_AD_4graph = zeros(length(noiseRange), length(stdDevCoeff), length(RN_STAZIONI(:,1)));
overlap_AP_4graph = zeros(length(noiseRange), length(stdDevCoeff), length(RN_STAZIONI(:,1)));

for noiseIdx = 1:length(noiseRange)
    for stdDevIdx = 1:length(stdDevCoeff)
        for stationIdx = 1:length(RN_STAZIONI(:,1))
            overlapResults_AD{noiseIdx, stdDevIdx, stationIdx} = FNC_OVL_2EVs(EV_AD, results_AD{noiseIdx, stdDevIdx, stationIdx}.EV_AD,N_samples_overlap,LoLa_RES,N_NORM_BINS);
            overlapResults_AP{noiseIdx, stdDevIdx, stationIdx} = FNC_OVL_2EVs(EV_AP, results_AP{noiseIdx, stdDevIdx, stationIdx}.EV_AP,N_samples_overlap,LoLa_RES,N_NORM_BINS);
            overlap_AD_4graph(noiseIdx, stdDevIdx, stationIdx) = overlapResults_AD{noiseIdx, stdDevIdx, stationIdx};
            overlap_AP_4graph(noiseIdx, stdDevIdx, stationIdx) = overlapResults_AP{noiseIdx, stdDevIdx, stationIdx};
        end
    end
end
% Salvataggio overlap finale
save(fullfile(outputDir, 'overlapResults_AD.mat'), 'overlapResults_AD');
save(fullfile(outputDir, 'overlapResults_AP.mat'), 'overlapResults_AP');
save(fullfile(outputDir, 'WORKSPACE.mat'));
else
    load('WORKSPACE.mat');
end

% SEZIONE RAPPRESENTAZIONE GRAFICA
% Creazione del grafico 3D per la sovrapposizione AD
[X, Y, Z] = ndgrid(noiseRange, stdDevCoeff, numStazRange);
overlap_AD_vector = overlap_AD_4graph(:);
overlap_AP_vector = overlap_AP_4graph(:);

% Grafico di overlap per Adria
figure;
scatter3(X(:), Y(:), Z(:), 50, overlap_AD_vector, 'filled');
xlabel('Noise');
ylabel('StdDev');
zlabel('Number of Stations');
colorbar;
title('Overlap tra i vettori di Eulero (Adria)');

% Grafico di overlap per Apulia
figure;
scatter3(X(:), Y(:), Z(:), 50, overlap_AP_vector, 'filled');
xlabel('Noise');
ylabel('StdDev');
zlabel('Number of Stations');
colorbar;
title('Overlap tra i vettori di Eulero (Apulia)');