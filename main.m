clear, clc, close all;

% TODO: 
% 1) Cicli for annidati per variare i deviazioneStandard/velocità/Eulero ecc..

% Usa la funzione custom per leggere i file e generare delle tabelle
citta_AD = leggiDatiCitta('citta/citta_dentro_AD.txt');
citta_AP = leggiDatiCitta('citta/citta_dentro_AP.txt');
citta_esterne = leggiDatiCitta('citta/citta_fuori.txt');

minStazioni = 5;
maxStazioni = 20;

% Genera un numero casuale di stazioni per ogni città
numStazioni_AD = randi([minStazioni, maxStazioni], height(citta_AD), 1);
numStazioni_AP = randi([minStazioni, maxStazioni], height(citta_AP), 1);

% Aggiunge la colonna stazioni alle rispettive tabelle
citta_AD.NumStazioni = numStazioni_AD;
citta_AP.NumStazioni = numStazioni_AP;

noiseRange = 0.2:0.2:2; % rumori da 0.2 a 2, passi di 0.2mm/yr
stdDevCoeff = 0.6:0.1:1.5; % deviazioni standard

% Caricamento vettori Eulero reali
EV_AD = load('EV_AD_0601_0812.txt');
EV_AP = load('EV_AP_1701_1909.txt');

% Creazione delle stazioni
stazioni_AD = creaStazioni(citta_AD);
stazioni_AP = creaStazioni(citta_AP);

% Aggiunge le stazioni alla tabella
citta_AD.Stazioni = stazioni_AD;
citta_AP.Stazioni = stazioni_AP;

% generaGrafico(citta_AD, citta_AP, citta_esterne);
