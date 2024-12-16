function [stazioni] = creaStazioni(citta)

stazioni = cell(height(citta), 1);

for i = 1:height(citta)
    latCitta = citta.Latitudine(i);
    lonCitta = citta.Longitudine(i);
    supCitta = citta.Superficie(i);
    numStazCitta = citta.NumStazioni(i);
    % Prealloca una matrice per [Latitudine, Longitudine]
    stazioniCitta = zeros(numStazCitta, 2);
    
    for numStazione = 1:numStazCitta
        [stazioniCitta(numStazione, 1), stazioniCitta(numStazione, 2)] = generaCoordinateStazione(latCitta, lonCitta, supCitta);
        % disp(numStazione + ") Questa stazione a " + dataset.Nome{i} + " ha latitudine: " + stazioniCitta(numStazione, 1) + " e longitudine: " + stazioniCitta(numStazione, 2));
    end
    stazioni{i} = stazioniCitta;
    % disp(newline);
end
end
