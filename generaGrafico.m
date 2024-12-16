function generaGrafico(citta_AD, citta_AP, citta_esterne)
figure;

% Imposta la figura a schermo intero per una miglior visualizzazione
set(gcf, 'Position', get(0, 'Screensize'));

% Caricamento dati coste
load coastlines coastlat coastlon;

% Caricamento dati placche
AD = load('PLATE_CONTOUR_ADRIA.txt');
AP = load('PLATE_CONTOUR_APULIA.txt');
lonAD = AD(:, 1);
latAD = AD(:, 2);
lonAP = AP(:, 1);
latAP = AP(:, 2);

% Mostra la mappa dell'Italia
geoshow(coastlat, coastlon, 'DisplayType', 'line', 'Color', 'k', 'HandleVisibility', 'off'); % Mappa costiera in nero
xlim([5 20]);
ylim([37 47]);
hold on;

% Aggiunge al grafico le placche
plot(lonAD, latAD, '-b', 'LineWidth', 1.5, 'DisplayName', 'Microplacca Adria');
plot(lonAP, latAP, '-r', 'LineWidth', 1.5, 'DisplayName', 'Microplacca Apulia');

% Aggiunge i capoluoghi dell'Adria
firstCap = true; % Flag per aggiungere la legenda solo una volta
for i = 1:height(citta_AD)
    latCitta = citta_AD.Latitudine(i);
    lonCitta = citta_AD.Longitudine(i);
    if firstCap
        plot(lonCitta, latCitta, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'DisplayName', 'Capoluoghi');
        firstCap = false;
    else
        plot(lonCitta, latCitta, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
    end
end

% Aggiunge i capoluoghi dell'Apulia
for i = 1:height(citta_AP)
    latCitta = citta_AP.Latitudine(i);
    lonCitta = citta_AP.Longitudine(i);
    plot(lonCitta, latCitta, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
end

% Aggiunge i capoluoghi esterni alle placche
for i = 1:height(citta_esterne)
    latCitta = citta_esterne.Latitudine(i);
    lonCitta = citta_esterne.Longitudine(i);
    plot(lonCitta, latCitta, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
end

% Aggiunge le stazioni GNSS dell'Adria
firstStazioni = true; % Flag per aggiungere la legenda solo una volta
for i = 1:height(citta_AD)
    stazioni = citta_AD.Stazioni{i};
    latStazioni = stazioni(:, 1);
    lonStazioni = stazioni(:, 2);
    if firstStazioni
        plot(lonStazioni, latStazioni, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3, 'DisplayName', 'Stazioni GNSS');
        firstStazioni = false;
    else
        plot(lonStazioni, latStazioni, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3, 'HandleVisibility', 'off');
    end
end

% Aggiunge le stazioni GNSS dell'Apulia
for i = 1:height(citta_AP)
    stazioni = citta_AP.Stazioni{i};
    latStazioni = stazioni(:, 1);
    lonStazioni = stazioni(:, 2);
    plot(lonStazioni, latStazioni, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3, 'HandleVisibility', 'off');

end

% Aggiunge dettagli al grafico
legend;
xlabel('Longitudine (Deg E)');
ylabel('Latitudine (Deg N)');
title("Mappa dell'Italia con città e stazioni GNSS");
grid on;
hold off;
end
