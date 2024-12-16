% Carica il workspace
load('WORKSPACE.mat');
whos

% Prepara le variabili
noiseRange = 0.2:0.2:2; % range di rumore
stdDevCoeff = 0.6:0.1:1.5; % range di deviazione standard
RN_STAZIONI = [5 20; 20 100; 100 500; 500 1000]; % intervallo delle stazioni

% Preallocazione matrici per grafici
overlap_AD_2d = zeros(length(RN_STAZIONI(:,1)), length(noiseRange), length(stdDevCoeff));
overlap_AP_2d = zeros(length(RN_STAZIONI(:,1)), length(noiseRange), length(stdDevCoeff));
% Intervalli di stazioni (etichette)
stationLabels = ["5-20", "20-100", "100-500", "500-1000"];
stations_intervals = [5 20; 20 100; 100 500; 500 1000];
stationValues = mean(RN_STAZIONI, 2); % Media intervalli per rappresentazione
% Estrazione dati di overlap
overlap_AD_4graph = squeeze(mean(overlap_AD_4graph, 3)); % Aggregazione lungo l'asse delle stazioni
overlap_AP_4graph = squeeze(mean(overlap_AP_4graph, 3)); 
% Calcolo overlap medio rispetto al numero di stazioni
mean_overlap_AD = zeros(1, 4);  % Media overlap per Adria (5-20, 20-100, 100-500, 500-1000)
mean_overlap_AP = zeros(1, 4);  % Media overlap per Apulia (5-20, 20-100, 100-500, 500-1000)
% medie overlap per Adria
if iscell(overlapResults_AD) && size(overlapResults_AD, 3) == 4
    % Calcolo delle medie degli overlap per Adria
    mean_overlap_AD = zeros(1, 4);
    for i = 1:4
        % Estrae i dati dalla terza dimensione e calcola la media
        mean_overlap_AD(i) = mean(mean(cell2mat(overlapResults_AD(:,:,i)), 1), 2); % Media lungo la 1° e 2° dimensione
    end
else
    error('overlapResults_AD non è una cella con 4 intervalli di stazioni');
end
if iscell(overlapResults_AP) && size(overlapResults_AP, 3) == 4
    % Calcolo delle medie degli overlap per Apulia
    mean_overlap_AP = zeros(1, 4);
    for i = 1:4
        % Estrae i dati dalla terza dimensione e calcola la media
        mean_overlap_AP(i) = mean(mean(cell2mat(overlapResults_AP(:,:,i)), 1), 2); % Media lungo la 1° e 2° dimensione
    end
else
    error('overlapResults_AP non è una cella con 4 intervalli di stazioni');
end


% \\\\\ MICROPLACCA ADRIA ///// 
% Grafico 1: doppio grafico, Intervallo overlap 0-1
figure;
subplot(2, 1, 1);
plot(noiseRange, mean(overlap_AD_4graph, 2), 'o-', 'LineWidth', 1.5);
xlabel('Rumore (mm/yr)');
ylabel('Overlap');
ylim([0 1]);
title('Adria - Overlap vs Rumore (mm/yr)');
grid on;
subplot(2, 1, 2);
plot(stdDevCoeff, mean(overlap_AD_4graph, 1), 's-', 'LineWidth', 1.5);
xlabel('Deviazione Standard (mm/yr)');
ylabel('Overlap');
ylim([0 1]);
title('Adria - Overlap vs Deviazione Standard (mm/yr)');
grid on;
% Grafico 2: doppio grafico, intervallo overlap ristretto
figure;
subplot(2, 1, 1);
plot(noiseRange, mean(overlap_AD_4graph, 2), 'o-', 'LineWidth', 1.5);
xlabel('Rumore (mm/yr)');
ylabel('Overlap');
title('Adria - Overlap vs Rumore (mm/yr)');
grid on;
subplot(2, 1, 2);
plot(stdDevCoeff, mean(overlap_AD_4graph, 1), 's-', 'LineWidth', 1.5);
xlabel('Deviazione Standard (mm/yr)');
ylabel('Overlap');
title('Adria - Overlap vs Deviazione Standard (mm/yr)');
grid on;
% Grafico 3: grafico singolo, Overlap vs Rumore, intervallo 0-1
figure;
plot(noiseRange, mean(overlap_AD_4graph, 2), 'o-', 'LineWidth', 1.5);
xlabel('Rumore (mm/yr)');
ylabel('Overlap');
ylim([0 1]);
title('Adria - Overlap vs Rumore (mm/yr)');
grid on;
% Grafico 4: grafico singolo, Overlap vs Deviazione Standard, intervallo 0-1
figure;
plot(stdDevCoeff, mean(overlap_AD_4graph, 1), 's-', 'LineWidth', 1.5);
xlabel('Deviazione Standard (mm/yr)');
ylabel('Overlap');
ylim([0 1]);
title('Adria - Overlap vs Deviazione Standard (mm/yr)');
grid on;

% Grafico 5: grafico singolo, Overlap vs Rumore, intervallo ristretto
figure;
plot(noiseRange, mean(overlap_AD_4graph, 2), 'o-', 'LineWidth', 1.5);
xlabel('Rumore (mm/yr)');
ylabel('Overlap');
title('Adria - Overlap vs Rumore (mm/yr)');
grid on;
% Grafico 6: grafico singolo, Overlap vs Deviazione Standard, intervallo ristretto
figure;
plot(stdDevCoeff, mean(overlap_AD_4graph, 1), 's-', 'LineWidth', 1.5);
xlabel('Deviazione Standard (mm/yr)');
ylabel('Overlap');
title('Adria - Overlap vs Deviazione Standard (mm/yr)');
grid on;
% Grafico 7: grafico singolo, Overlap vs numero di stazioni
figure;
hold on;
plot(stations_intervals(:,1), mean_overlap_AD, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Adria');
xlabel('Intervallo di Stazioni');
ylabel('Overlap Medio');
xticks([[5 20 100 500]]);
xticklabels({'5-20', '20-100', '100-500', '500-1000'});
title('Adria - Overlap vs Intervallo di Stazioni');
grid on;

% figure;
% plot(stationValues, meanoverlap_AD_stations, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% xticks(stationValues);
% xticklabels(stationLabels);
% xlabel('Numero di Stazioni');
% ylabel('Overlap');
% ylim([0 1]);
% title('Adria: Overlap vs Numero di Stazioni');
% grid on;
% % Aggiunta etichette ai punti
% for i = 1:min(length(stationValues), length(meanoverlap_AD_stations))
%     text(stationValues(i), meanoverlap_AD_stations(i), sprintf('%.2f', meanoverlap_AD_stations(i)), ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10);
% end

% \\\\\ MICROPLACCA APULIA ///// 
% Grafico 1: doppio grafico, Intervallo overlap 0-1
figure;
subplot(2, 1, 1);
plot(noiseRange, mean(overlap_AP_4graph, 2), 'o-', 'LineWidth', 1.5);
xlabel('Rumore (mm/yr)');
ylabel('Overlap');
ylim([0 1]);
title('Apulia - Overlap vs Rumore (mm/yr)');
grid on;
subplot(2, 1, 2);
plot(stdDevCoeff, mean(overlap_AP_4graph, 1), 's-', 'LineWidth', 1.5);
xlabel('Deviazione Standard (mm/yr)');
ylabel('Overlap');
ylim([0 1]);
title('Apulia - Overlap vs Deviazione Standard (mm/yr)');
grid on;
% Grafico 2: doppio grafico, intervallo overlap ristretto
figure;
subplot(2, 1, 1);
plot(noiseRange, mean(overlap_AP_4graph, 2), 'o-', 'LineWidth', 1.5);
xlabel('Rumore (mm/yr)');
ylabel('Overlap');
title('Apulia - Overlap vs Rumore (mm/yr)');
grid on;
subplot(2, 1, 2);
plot(stdDevCoeff, mean(overlap_AP_4graph, 1), 's-', 'LineWidth', 1.5);
xlabel('Deviazione Standard (mm/yr)');
ylabel('Overlap');
title('Apulia - Overlap vs Deviazione Standard (mm/yr)');
grid on;
% Grafico 3: grafico singolo, Overlap vs Rumore, intervallo 0-1
figure;
plot(noiseRange, mean(overlap_AP_4graph, 2), 'o-', 'LineWidth', 1.5);
xlabel('Rumore (mm/yr)');
ylabel('Overlap');
ylim([0 1]);
title('Apulia - Overlap vs Rumore (mm/yr)');
grid on;
% Grafico 4: grafico singolo, Overlap vs Deviazione Standard, intervallo 0-1
figure;
plot(stdDevCoeff, mean(overlap_AP_4graph, 1), 's-', 'LineWidth', 1.5);
xlabel('Deviazione Standard (mm/yr)');
ylabel('Overlap');
ylim([0 1]);
title('Apulia - Overlap vs Deviazione Standard (mm/yr)');
grid on;
% Grafico 5: grafico singolo, Overlap vs Rumore, intervallo ristretto
figure;
plot(noiseRange, mean(overlap_AP_4graph, 2), 'o-', 'LineWidth', 1.5);
xlabel('Rumore (mm/yr)');
ylabel('Overlap');
title('Apulia - Overlap vs Rumore (mm/yr)');
grid on;
% Grafico 6: grafico singolo, Overlap vs Deviazione Standard, intervallo ristretto
figure;
plot(stdDevCoeff, mean(overlap_AP_4graph, 1), 's-', 'LineWidth', 1.5);
xlabel('Deviazione Standard (mm/yr)');
ylabel('Overlap');
title('Apulia - Overlap vs Deviazione Standard (mm/yr)');
grid on;
% Grafico 7: grafico singolo, Overlap vs numero di stazioni
figure;
hold on;
plot(stations_intervals(:,1), mean_overlap_AP, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Apulia');
xlabel('Intervallo di Stazioni');
ylabel('Overlap Medio');
xticks([5 20 100 500]);
xticklabels({'5-20', '20-100', '100-500', '500-1000'});
title('Apulia - Overlap vs Intervallo di Stazioni');
grid on;

% figure;
% plot(stationValues, meanoverlap_AP_stations, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% xticks(stationValues);
% xticklabels(stationLabels);
% xlabel('Numero di Stazioni');
% ylabel('Overlap');
% ylim([0 1]);
% title('Apulia: Overlap vs Numero di Stazioni');
% grid on;
% 
% % Aggiunta etichette ai punti
% for i = 1:min(length(stationValues), length(meanoverlap_AP_stations))
%     text(stationValues(i), meanoverlap_AP_stations(i), sprintf('%.2f', meanoverlap_AP_stations(i)), ...
%         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 10);
% end