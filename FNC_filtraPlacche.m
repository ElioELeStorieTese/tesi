function [listaAD, listaAP, fuori] = FNC_filtraPlacche(data, lonAD, latAD, lonAP, latAP)
%% [listaAD, listaAP, fuori] = FNC_filtraPlacche(data, lonAD, latAD, lonAP, latAP);
% Separazione delle città nei gruppi
cityNames = data.Var1;       % Nomi delle città
cityLon = data.Var2;         % Longitudine
cityLat = data.Var3;         % Latitudine
cityPop = data.Var4;         % Popolazione
cityArea = data.Var5;        % Superficie
cityAlt = data.Var6;         % Altitudine

% Separazione delle città nei gruppi
listaAD = "";
listaAP = "";
fuori = "";

% Numero di città
numeroCitta = height(data);

for citta = 1:numeroCitta
    % Creazione della riga come stringa
    riga = sprintf('%s %.4f %.4f %.0f %.1f %.1f', cityNames{citta}, cityLon(citta), cityLat(citta), ...
        cityPop(citta), cityArea(citta), cityAlt(citta));

    if inpolygon(cityLon(citta), cityLat(citta), lonAD, latAD)
        listaAD = listaAD + riga + newline;
    elseif inpolygon(cityLon(citta), cityLat(citta), lonAP, latAP)
        listaAP = listaAP + riga + newline;
    else
        fuori = fuori + riga + newline;
    end
end
end
