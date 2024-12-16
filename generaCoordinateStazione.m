function [stazioneLat, stazioneLon] = generaCoordinateStazione(centroLat, centroLon, superficie)
    % Genera le coordinate delle stazioni GNSS fittizie 
    % distribuite casualmente in un'area circolare.
    %
    % Input:
    %   centroLat - Latitudine del centro città
    %   centroLon - Longitudine del centro città
    %   superficie - Area della città
    %
    % Output:
    %   stazioniLat - Latitudine della stazione generata
    %   stazioniLon - Longitudine della staziona generata
    
    RAGGIO_TERRA = 6371;

    % Calcola il raggio dalla superficie = pi * raggio^2
    raggioCitta = sqrt(superficie / pi);
    
    % % Genera un angolo casuale e una distanza casuale in coordinate polari
    % angolo_casuale = rand() * 2 * pi; % Angolo casuale (in radianti)
    % raggio_casuale = rand() * raggioCitta; % Distanza scalata al raggio massimo
    % 
    % % Converte le coordinate polari in variazioni relative
    % deltaX = raggio_casuale * cos(angolo_casuale); % Spostamento lungo l'asse X
    % deltaY = raggio_casuale * sin(angolo_casuale); % Spostamento lungo l'asse Y
    % 
    % % Converte le variazioni relative in latitudine e longitudine
    % deltaLat = (deltaY / RAGGIO_TERRA) * (180 / pi);
    % deltaLon = (deltaX / RAGGIO_TERRA) * (180 / pi) / cosd(centroLat);
    % 
    % % Calcola le coordinate finali
    % stazioneLat = centroLat + deltaLat;
    % stazioneLon = centroLon + deltaLon;

    % Alternativa con distribuzione migliore
    % Genera un angolo casuale per la posizione
    angolo = 2 * pi * rand();

    % Genera una distanza casuale all'interno del raggio
    distanza = raggioCitta * sqrt(rand());

    % Calcola le variazioni di latitudine e longitudine
    deltaLat = (distanza / RAGGIO_TERRA) * (180 / pi);
    deltaLon = (distanza / RAGGIO_TERRA) * (180 / pi) / cosd(centroLat);

    % Calcola le coordinate della stazione
    stazioneLat = centroLat + deltaLat * cos(angolo); % Latitudine
    stazioneLon = centroLon + deltaLon * sin(angolo); % Longitudine
end
