clear, clc, close all;

% \\\\\\\\\\ CARICAMENTO DATI //////////
earthRadius = 6371;
% PARAMETRI D'INDAGINE: stazioni, rumore (sensibilità dei dati), deviazione standard (varianza dei dati)
% Stazioni
minStazioni = 5;
maxStazioni = 20;
n=10;
numStazRange = randi([minStazioni, maxStazioni], [1, n]);
% Deviazione standard
stdDevCoeff = 0.6:0.1:1.5;
% Rumore
noiseRange = 0.2:0.2:2; %da 0.2 a 2, passi di 0.2mm/yr
% Vettori di Eulero
DR='C:/Users/utente/Documents/MATLAB_tesi';
N_samples_vettore = 1e4;
N_samples_overlap = 1e5;

outputDir = 'C:/Users/utente/Documents/MATLAB_tesi';
% 1) Matrice risultati singola
results = cell(length(noiseRange), length(stdDevCoeff), length(numStazRange));
% 2) Matrice risultati doppia: una per placca
% results_AD = cell(length(noiseRange), length(stdDevCoeff), length(numStazRange));
% results_AP = cell(length(noiseRange), length(stdDevCoeff), length(numStazRange));

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

% Caricamento dati città
cityfileID = fopen('Città Italiane.txt', 'r');
data = textscan(cityfileID, '%s %f %f %d %f %f', 'Delimiter', {' ', '\t'}, 'HeaderLines', 0);
fclose(cityfileID);
cityID = data{1};
cityLon = data{2};
cityLat = data{3};
cityPopul = data{4};
cityArea = data{5};
cityAltitude = data{6};
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

% PREGENERAZIONE DI STAZIONI CASUALI
numCities_AD = length(cityID_AD);
numCities_AP = length(cityID_AP);
numCities_fuori = length(cityID_fuori);
% Numero di iterazioni per stationIdx
numStationIdx = length(numStazRange);
% Matrice per i numeri di stazioni (city x stationIdx)
numStazioniMatrix_AD = zeros(numCities_AD, numStationIdx);
numStazioniMatrix_AP = zeros(numCities_AP, numStationIdx);
numStazioniMatrix_fuori = zeros(numCities_fuori, numStationIdx);
% Generazione numeri di stazioni unici per ogni città e stationIdx
for i = 1:numCities_AD
    numStazioniMatrix_AD(i, :) = randperm(maxStazioni - minStazioni + 1, numStationIdx) + 4; % Numeri unici tra 5 e 20
end
for i = 1:numCities_AP
    numStazioniMatrix_AP(i, :) = randperm(maxStazioni - minStazioni + 1, numStationIdx) + 4;
end
for i = 1:numCities_fuori
    numStazioniMatrix_fuori(i, :) = randperm(maxStazioni - minStazioni + 1, numStationIdx) + 4;
end

% \\\\\\\\\\ SEZIONE CALCOLI //////////
for noiseIdx = 1:length(noiseRange)
    noiseI = noiseRange(noiseIdx);
    for stdDevIdx = 1:length(stdDevCoeff)
        stdDevJ = stdDevCoeff(stdDevIdx);
        for stationIdx = 1:length(numStazRange)
            % stationK = numStazRange(stationIdx); vecchio metodo
            stationK = numStazRange(stationIdx);

            fprintf('Elaborazione: Noise = %.1f, StdDev = %.1f, Stations = %d\n', noiseI, stdDevJ, stationK);

            %----- GENERAZIONE STAZIONI GNSS -----
            % GENERAZIONE STAZIONI GNSS ADRIA
            for i = 1:length(cityID_AD)
                numStazioni_AD = numStazioniMatrix_AD(i, stationIdx);
                radius_AD = sqrt(cityArea_AD(i) / pi); % Raggio
                angles_AD = 2 * pi * rand(numStazioni_AD, 1); % Angoli casuali
                puntoStaz_AD = radius_AD * sqrt(rand(numStazioni_AD, 1)); % Distanza casuale entro il raggio
                deltaLat_AD = (puntoStaz_AD / earthRadius) * (180 / pi); % Variazione latitudine dalle coordinate d'origine
                deltaLon_AD = (puntoStaz_AD / earthRadius) .* (180 / pi) ./ cosd(cityLat_AD(i)); % Variazione longitudine dalle coordinate d'origine
                puntoStazLat_AD = cityLat_AD(i) + deltaLat_AD .* cos(angles_AD); % Latitudine stazioni
                puntoStazLon_AD = cityLon_AD(i) + deltaLon_AD .* sin(angles_AD); % Longitudine stazioni
                % Accumulo risultati
                stazLat_AD = [stazLat_AD; puntoStazLat_AD];
                stazLon_AD = [stazLon_AD; puntoStazLon_AD];
                cityIndex_AD = [cityIndex_AD; repmat(i, numStazioni_AD, 1)]; % ID città associato
            end

            % stazLat_AD = [];
            % stazLon_AD = [];
            % cityIndex_AD = [];
            % for i = 1:length(cityID_AD)
            %     % Parametri città in Adria
            %     centroLat_AD = cityLat_AD(i);
            %     centroLon_AD = cityLon_AD(i);
            %     stationArea_AD = cityArea_AD(i);
            %     numStazioni_AD = randi([minStazioni, maxStazioni]);
            %     % Generazione stazioni GNSS
            %     radius_AD = sqrt(stationArea_AD / pi); % Raggio
            %     angles_AD = 2 * pi * rand(numStazioni_AD, 1); % Angoli casuali
            %     puntoStaz_AD = radius_AD * sqrt(rand(numStazioni_AD, 1)); % Distanza casuale entro il raggio
            %     deltaLat_AD = (puntoStaz_AD / earthRadius) * (180 / pi); % Variazione latitudine dalle coordinate d'origine
            %     deltaLon_AD = (puntoStaz_AD / earthRadius) .* (180 / pi) ./ cosd(centroLat_AD); % Variazione longitudine dalle coordinate d'origine
            %     puntoStazLat_AD = centroLat_AD + deltaLat_AD .* cos(angles_AD); % Latitudine stazioni
            %     puntoStazLon_AD = centroLon_AD + deltaLon_AD .* sin(angles_AD); % Longitudine stazioni
            %     % Accumulo risultati
            %     stazLat_AD = [stazLat_AD; puntoStazLat_AD];
            %     stazLon_AD = [stazLon_AD; puntoStazLon_AD];
            %     cityIndex_AD = [cityIndex_AD; repmat(i, numStazioni_AD, 1)]; % ID città associato
            % end

            % GENERAZIONE STAZIONI GNSS APULIA
            for i = 1:length(cityID_AP)
                numStazioni_AP = numStazioniMatrix_AP(i, stationIdx);
                radius_AP = sqrt(cityArea_AP(i) / pi); % Raggio
                angles_AP = 2 * pi * rand(numStazioni_AP, 1); % Angoli casuali
                puntoStaz_AP = radius_AP * sqrt(rand(numStazioni_AP, 1)); % Distanza casuale entro il raggio
                deltaLat_AP = (puntoStaz_AP / earthRadius) * (180 / pi); % Variazione latitudine dalle coordinate d'origine
                deltaLon_AP = (puntoStaz_AP / earthRadius) .* (180 / pi) ./ cosd(cityLat_AP(i)); % Variazione longitudine dalle coordinate d'origine
                puntoStazLat_AP = cityLat_AP(i) + deltaLat_AP .* cos(angles_AP); % Latitudine stazioni
                puntoStazLon_AP = cityLon_AP(i) + deltaLon_AP .* sin(angles_AP); % Longitudine stazioni
                % Accumulo risultati
                stazLat_AP = [stazLat_AP; puntoStazLat_AP];
                stazLon_AP = [stazLon_AP; puntoStazLon_AP];
                cityIndex_AP = [cityIndex_AP; repmat(i, numStazioni_AP, 1)]; % ID città associato
            end

            % stazLat_AP = [];
            % stazLon_AP = [];
            % cityIndex_AP = [];
            % for i = 1:length(cityID_AP)
            %     % Parametri città in Apulia
            %     centroLat_AP = cityLat_AP(i);
            %     centroLon_AP = cityLon_AP(i);
            %     stationArea_AP = cityArea_AP(i);
            %     numStazioni_AP = randi([minStazioni, maxStazioni]);
            %     % Generazione stazioni GNSS
            %     radius_AP = sqrt(stationArea_AP / pi); % Raggio
            %     angles_AP = 2 * pi * rand(numStazioni_AP, 1); % Angoli casuali
            %     puntoStaz_AP = radius_AP * sqrt(rand(numStazioni_AP, 1)); % Distanza casuale entro il raggio
            %     deltaLat_AP = (puntoStaz_AP / earthRadius) * (180 / pi); % Variazione latitudine
            %     deltaLon_AP = (puntoStaz_AP / earthRadius) .* (180 / pi) ./ cosd(centroLat_AP); % Variazione longitudine
            %     puntoStazLat_AP = centroLat_AP + deltaLat_AP .* cos(angles_AP); % Latitudine stazioni
            %     puntoStazLon_AP = centroLon_AP + deltaLon_AP .* sin(angles_AP); % Longitudine stazioni
            %     % Accumulo risultati
            %     stazLat_AP = [stazLat_AP; puntoStazLat_AP];
            %     stazLon_AP = [stazLon_AP; puntoStazLon_AP];
            %     cityIndex_AP = [cityIndex_AP; repmat(i, numStazioni_AP, 1)]; % ID città associato
            % end

            % GENERAZIONE STAZIONI GNSS ESTERNE
            for i = 1:length(cityID_fuori)
                numStazioni_fuori = numStazioniMatrix_fuori(i, stationIdx);
                radius_fuori = sqrt(cityArea_fuori(i) / pi); % Raggio
                angles_fuori = 2 * pi * rand(numStazioni_fuori, 1); % Angoli casuali
                puntoStaz_fuori = radius_fuori * sqrt(rand(numStazioni_fuori, 1)); % Distanza casuale entro il raggio
                deltaLat_fuori = (puntoStaz_fuori / earthRadius) * (180 / pi); % Variazione latitudine dalle coordinate d'origine
                deltaLon_fuori = (puntoStaz_fuori / earthRadius) .* (180 / pi) ./ cosd(cityLat_fuori(i)); % Variazione longitudine dalle coordinate d'origine
                puntoStazLat_fuori = cityLat_fuori(i) + deltaLat_fuori .* cos(angles_fuori); % Latitudine stazioni
                puntoStazLon_fuori = cityLon_fuori(i) + deltaLon_fuori .* sin(angles_fuori); % Longitudine stazioni
                % Accumulo risultati
                stazLat_fuori = [stazLat_fuori; puntoStazLat_fuori];
                stazLon_fuori = [stazLon_fuori; puntoStazLon_fuori];
                cityIndex_fuori = [cityIndex_fuori; repmat(i, numStazioni_fuori, 1)]; % ID città associato
            end

            % stazLat_fuori = [];
            % stazLon_fuori = [];
            % for i = 1:length(cityID_fuori)
            %     % Parametri della città Esterna
            %     centroLat_fuori = cityLat_fuori(i);
            %     centroLon_fuori = cityLon_fuori(i);
            %     stationArea_fuori = cityArea_fuori(i);
            %     numStazioni_fuori = randi([minStazioni, maxStazioni]);
            %     % Generazione stazioni GNSS
            %     radius = sqrt(stationArea_fuori / pi); % Raggio
            %     angles = 2 * pi * rand(numStazioni_fuori, 1); % Angoli casuali
            %     puntoStaz = radius * sqrt(rand(numStazioni_fuori, 1)); % Distanza casuale entro il raggio
            %     deltaLat = (puntoStaz / earthRadius) * (180 / pi); % Variazione latitudine
            %     deltaLon = (puntoStaz / earthRadius) .* (180 / pi) ./ cosd(centroLat_fuori); % Variazione longitudine
            %     puntoStazLat = centroLat_fuori + deltaLat .* cos(angles); % Latitudine stazioni
            %     puntoStazLon = centroLon_fuori + deltaLon .* sin(angles); % Longitudine stazioni
            %     % Accumulo risultati
            %     stazLat_fuori = [stazLat_fuori; puntoStazLat];
            %     stazLon_fuori = [stazLon_fuori; puntoStazLon];
            % end

            fprintf('Stazioni totali: %d\n', length(stazLat_fuori)+length(stazLat_AP)+length(stazLat_AP));
            fprintf('Stazioni in Adria: %d\n', length(stazLat_AP));
            fprintf('Stazioni in Apulia: %d\n', length(stazLat_AP));

            % ----- CALCOLO DEVIAZIONE STANDARD, VELOCITà, RUMORE -----
            [vE_AD, vN_AD] = FNC_EV2siteSV(EV_AD, [stazLon_AP, stazLat_AP]);
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
            file_path_vel_AD = sprintf('SV_AD_Noise%.1f_StDev%.1f_#Staz%d.txt', noiseI, stdDevJ, stationK);
            file_path_vel_AP = sprintf('SV_AP_Noise%.1f_StDev%.1f_#Staz%d.txt', noiseI, stdDevJ, stationK);

            fileID_AD = fopen(file_path_vel_AD, 'w');
            for i = 1:size(stazLat_AP, 1)
                cityID_AD = cityIndex_AP(i);
                cityAltitude_AD = cityAltitude(cityID_AD);
                fprintf(fileID_AD, '%.4f %.4f %.1f %1.7e %1.7e %1.7e %1.7e %1.7e\n', stazLon_AP(i), stazLat_AP(i), cityAltitude_AD, vE_AD_noise(i), vN_AD_noise(i), stdDevE_AD(i), stdDevN_AD(i), 0.0);
            end
            fclose(fileID_AD);

            fileID_AP = fopen(file_path_vel_AP, 'w');
            for i = 1:size(stazLat_AP, 1)
                cityID_AP = cityIndex_AP(i);
                cityAltitude_AP = cityAltitude(cityID_AP);
                fprintf(fileID_AP, '%.4f %.4f %.1f %1.7e %1.7e %1.7e %1.7e %1.7e\n', stazLon_AP(i), stazLat_AP(i), cityAltitude_AP, vE_AP_noise(i), vN_AP_noise(i), stdDevE_AP(i), stdDevN_AP(i), 0.0);
            end
            fclose(fileID_AP);

            % figure;
            % geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'k', 'HandleVisibility', 'off'); % Mappa costiera in nero
            % xlim([5 20]); % Intervallo longitudine
            % ylim([37 47]); % Intervallo latitudine
            % hold on;
            % plot(lonAD, latAD, '-b', 'LineWidth', 1.5, 'DisplayName', 'Microplacca Adria');
            % plot(lonAP, latAP, '-r', 'LineWidth', 1.5, 'DisplayName', 'Microplacca Apulia');
            % plot(allStazLon, allStazLat, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'm', 'DisplayName', 'Stazioni GNSS (altre)');
            % plot(stazLonAD, stazLatAD, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'DisplayName', 'Stazioni GNSS (Adria)');
            % plot(stazLonAP, stazLatAP, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', 'Stazioni GNSS (Apulia)');
            % xlabel('Longitudine (Deg E)');
            % ylabel('Latitudine (Deg N)');
            % title('PLACEHOLDER');
            % legend;
            % grid on;
            % hold off;

            % ----- GENERAZIONE VETTORE DI EULERO -----
            EV_AD_result = FNC_SV2EV_C(DR,file_path_vel_AD,N_samples_vettore);
            EV_AP_result = FNC_SV2EV_C(DR,file_path_vel_AP,N_samples_vettore);
                % 1a) Mat. singola: struct contenente vettore AD, vettore AP, e tre valori d'indice noiseI, stdDevJ, stationK
            results{noiseIdx, stdDevIdx, stationIdx} = struct('EV_AD',EV_AD_result, 'EV_AP',EV_AP_result, 'Noise',noiseI, 'StdDev',stdDevJ, 'Stations',stationK);
                % 1b) Mat. singola: struct contenente solo vettore AD e vettore AP, i valori d'indice sono inclusi nel label del vettore
            % label_AD = sprintf('EV_AD_Noise%.1f_StDev%.1f_Staz%d', noiseI, stdDevJ, stationK);
            % label_AP = sprintf('EV_AP_Noise%.1f_StDev%.1f_Staz%d', noiseI, stdDevJ, stationK);
            % results{noiseIdx, stdDevIdx, stationIdx} = struct(label_AD, EV_AD_result, label_AP, EV_AP_result);
                % 2a) Mat. doppia: per ciascuna, struct contenente il rispettivo vettore e i tre valori d'indice noiseI, stdDevJ, stationK
            % results_AD{noiseIdx, stdDevIdx, stationIdx} = struct();
            % results_AP{noiseIdx, stdDevIdx, stationIdx} = struct();
                % 2b) Mat. doppia: per ciascuna, struct contenente il
            % rispettivo vettore, i tre valori d'indice sono inclusi nel label del vettore
            % label_AD = sprintf('EV_AD_Noise%.1f_StDev%.1f_Staz%d', noiseI, stdDevJ, stationK);
            % label_AP = sprintf('EV_AP_Noise%.1f_StDev%.1f_Staz%d', noiseI, stdDevJ, stationK);
            % results_AD{noiseIdx, stdDevIdx, stationIdx} = struct();
            % results_AP{noiseIdx, stdDevIdx, stationIdx} = struct();

        end
    end
end


% % \\\\\\\\\\ RAPPRESENTAZIONE //////////
% % Elementi principali
% close all;
% figure;
% geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'k', 'HandleVisibility', 'off'); % Mappa costiera in nero
% xlim([5 20]); % Intervallo longitudine
% ylim([37 47]); % Intervallo latitudine
% hold on;
% plot(lonAD, latAD, '-b', 'LineWidth', 1.5, 'DisplayName', 'Microplacca Adria');
% plot(lonAP, latAP, '-r', 'LineWidth', 1.5, 'DisplayName', 'Microplacca Apulia');
% plot(allStazLon, allStazLat, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'm', 'DisplayName', 'Stazioni GNSS (altre)');
% plot(stazLonAdria, stazLatAdria, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'DisplayName', 'Stazioni GNSS (Adria)');
% plot(stazLonApulia, stazLatApulia, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'r', 'DisplayName', 'Stazioni GNSS (Apulia)');
% xlabel('Longitudine (Deg E)');
% ylabel('Latitudine (Deg N)');
% title('PLACEHOLDER');
% legend;
% grid on;
% hold off;
% 
% % Velocità fittizie
% figure;
% hold on;
% histogram(vTot_AD, 'FaceColor', 'b', 'EdgeColor', 'k', 'DisplayName', 'Dati velocità stazioni placca Adria');
% histogram(vTot_AP, 'FaceColor', 'r', 'EdgeColor', 'k', 'DisplayName', 'Dati velocità stazioni placca Apulia');
% xlabel('Velocità (mm/anno)');
% ylabel('Frequenza');
% title('Distribuzione delle velocità fittizie');
% legend;
% grid on;
% hold off;
% 
% fprintf('Numero totale di stazioni fittizie: %d\n', length(allStazLat));
% fprintf('Numero di stazioni in Adria: %d\n', length(stazLatAdria));
% fprintf('Numero di stazioni in Apulia: %d\n', length(stazLatApulia));
% fprintf('Numero di stazioni esterne alle microplacche: %d\n', length(allStazLat) - (length(stazLatApulia) + length(stazLatAdria)));
