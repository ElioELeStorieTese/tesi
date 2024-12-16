% Definizione dei dati dei capoluoghi di provincia nel formato
% NOME/ID LONG(degE) LAT(degN) POPOLAZIONE SUPERFICIE UTILE(km2) ALTITUDINE(m)
dati_centri_abitati = {
% dati_centri_abitati_AD = {
% dati_centri_abitati_AP = {
% dati_centri_abitati_fuori = {
    'Alessandria', 8.6110, 44.9074, 406831, 17, 95;
    'Ancona', 13.5057, 43.6087, 100924, 24, 16;
    'Aosta', 7.3139, 45.7350, 34082, 8, 583;
    'Arezzo', 11.8824, 43.4633, 99425, 34, 296;
    'Ascoli_Piceno', 13.5719, 42.8540, 48195, 11, 154;
    'Asti', 8.2064, 44.9003, 207785, 21, 123;
    'Avellino', 14.7896, 40.9174, 52593, 11, 348;
    'Bari', 16.8724, 41.1136, 317241, 50, 5;
    'Barletta', 16.2776, 41.3117, 94616, 10.4, 15;
    'Belluno', 12.2176, 46.1392, 35372, 9, 389;
    'Benevento', 14.7815, 41.1306, 57444, 10.5, 130;
    'Bergamo', 9.6773, 45.6983, 120728, 33.5, 249;
    'Biella', 8.0538, 45.5613, 168707, 18.5, 420;
    'Bologna', 11.3426, 44.4949, 394463, 86.4, 54;
    'Brescia', 10.2118, 45.5416, 199724, 56, 149;
    'Brindisi', 17.9371, 40.6327, 86927, 13.5, 15;
    'Campobasso', 14.6664, 41.5618, 48462, 9.4, 701;
    'Caserta', 14.3332, 41.0734, 74954, 62, 68;
    'Chieti', 14.1620, 42.3512, 51316, 14.3, 330;
    'Como', 9.0852, 45.8081, 83651, 30, 201;
    'Cosenza', 16.2520, 39.3003, 67317, 11.2, 238;
    'Cremona', 10.0227, 45.1333, 71414, 12.2, 45;
    'Cuneo', 7.5508, 44.3842, 582194, 25, 534;
    'Fermo', 13.7211, 43.1583, 37400, 7, 319;
    'Ferrara', 11.6196, 44.8354, 132009, 33, 9;
    'Firenze', 11.2558, 43.7696, 372038, 63.6, 50;
    'Foggia', 15.5574, 41.4622, 147036, 15, 76;
    'Forlì', 12.0404, 44.2227, 118358, 40, 34;
    'Frosinone', 13.3426, 41.6410, 44847, 21.3, 291;
    'Genova', 8.9463, 44.4056, 562836, 60, 19;
    'Gorizia', 13.6205, 45.9411, 34168, 10.6, 84;
    'Grosseto', 11.1107, 42.7603, 81270, 12, 10;
    'Imperia', 8.0269, 43.8885, 42446, 15.8, 10;
    'Isernia', 14.2331, 41.5953, 21189, 3, 423;
    'L_Aquila', 13.3984, 42.3498, 69716, 15.6, 714;
    'La_Spezia', 9.8376, 44.1025, 92838, 17.5, 10;
    'Latina', 12.9050, 41.4676, 126470, 21, 21;
    'Lecce', 18.1709, 40.3515, 94896, 19.5, 49;
    'Lecco', 9.3904, 45.8566, 47198, 9.6, 214;
    'Livorno', 10.3116, 43.5485, 156998, 33.4, 3;
    'Lodi', 9.5037, 45.3145, 45270, 7, 87;
    'Lucca', 10.5034, 43.8442, 89746, 46, 19;
    'Macerata', 13.4537, 43.2997, 41492, 7, 315;
    'Mantova', 10.7914, 45.1564, 49547, 10.2, 19;
    'Massa', 10.1208, 44.0256, 69096, 34.4, 65;
    'Matera', 16.6069, 40.6672, 60316, 7.3, 401;
    'Milano', 9.1900, 45.4642, 1370996, 148, 122;
    'Modena', 10.9252, 44.6471, 186307, 40.7, 34;
    'Monza', 9.2769, 45.5845, 123240, 20.8, 162;
    'Napoli', 14.2681, 40.8518, 962589, 106, 17;
    'Novara', 8.6234, 45.4469, 364046, 22.2, 162;
    'Padova', 11.8768, 45.4064, 209829, 72, 12;
    'Parma', 10.3279, 44.8015, 198292, 38.2, 57;
    'Pavia', 9.1525, 45.1860, 72867, 15, 77;
    'Perugia', 12.3888, 43.1122, 166676, 60, 493;
    'Pesaro', 12.9119, 43.9125, 94345, 22.7, 11;
    'Pescara', 14.2161, 42.4618, 119217, 26.4, 4;
    'Piacenza', 9.6929, 45.0522, 103936, 27.8, 61;
    'Pisa', 10.4017, 43.7160, 89608, 23.7, 4;
    'Pistoia', 10.9252, 43.9335, 90034, 22, 67;
    'Pordenone', 12.6552, 45.9569, 51191, 39, 24;
    'Potenza', 15.8086, 40.6401, 67351, 19.6, 819;
    'Prato', 11.1022, 43.8804, 194590, 48, 65;
    'Ravenna', 12.2035, 44.4184, 159057, 18.5, 6;
    'Reggio_Emilia', 10.6315, 44.6983, 171944, 47, 58;
    'Rieti', 12.8580, 42.4041, 47287, 12.4, 405;
    'Rimini', 12.5686, 44.0678, 150576, 30.6, 5;
    'Roma', 12.4964, 41.9028, 2872800, 320, 21;
    'Rovigo', 11.7895, 45.0692, 51232, 12.6, 5;
    'Salerno', 14.7872, 40.6824, 133970, 11, 4;
    'Savona', 8.4791, 44.3084, 58712, 8.7, 4;
    'Siena', 11.3308, 43.3188, 53626, 16, 322;
    'Sondrio', 9.8713, 46.1698, 21784, 4, 307;
    'Taranto', 17.2665, 40.4596, 196702, 7.4, 15;
    'Teramo', 13.6997, 42.6589, 53934, 5.7, 265;
    'Terni', 12.6465, 42.5636, 111501, 24, 130;
    'Torino', 7.6869, 45.0703, 2203353, 128, 239;
    'Trento', 11.1217, 46.0669, 118902, 6.6, 194;
    'Treviso', 12.2459, 45.6665, 84198, 22.6, 15;
    'Trieste', 13.7768, 45.6495, 204345, 33, 2;
    'Udine', 13.2354, 46.0637, 100514, 44.2, 113;
    'Varese', 8.8236, 45.8206, 80994, 27, 382;
    'Venezia', 12.3155, 45.4408, 258685, 170, 1;
    'Verbania', 8.5567, 45.9299, 153844, 4.7, 197;
    'Vercelli', 8.4187, 45.3229, 165821, 11.6, 130;
    'Verona', 10.9916, 45.4384, 258031, 45, 59;
    'Vicenza', 11.5407, 45.5455, 111260, 30, 39;
    'Viterbo', 12.1077, 42.4207, 67432, 8.6, 326;
};

% Nome file CODICE ORIGINALE
file_path_centri_abitati = 'Città Italiane.txt';
% Abilita scrittura
fileID = fopen(file_path_centri_abitati, 'w');
% Scrittura & salvataggio file
for i = 1:size(dati_centri_abitati, 1)
    fprintf(fileID, '%s %.4f %.4f %d %.1f %.1f\n', dati_centri_abitati{i,1}, dati_centri_abitati{i,2}, dati_centri_abitati{i,3}, dati_centri_abitati{i,4}, dati_centri_abitati{i,5}, dati_centri_abitati{i,6});
end
% Chiudi il file
fclose(fileID);

% % salvataggio dataset città Adria
% file_path_centri_abitati_AD = 'Città Italiane Adria.txt';
% fileID_AD = fopen(file_path_centri_abitati_AD, 'w');
% for i = 1:size(dati_centri_abitati_AD, 1)
%     fprintf(fileID_AD, '%s %.4f %.4f %d %.1f %.1f\n', dati_centri_abitati_AD{i,1}, dati_centri_abitati_AD{i,2}, dati_centri_abitati_AD{i,3}, dati_centri_abitati_AD{i,4}, dati_centri_abitati_AD{i,5}, dati_centri_abitati_AD{i,6});
% end
% fclose(fileID_AD);

% % salvataggio dataset città Apulia
% file_path_centri_abitati_AP = 'Città Italiane Apulia.txt';
% fileID_AP = fopen(file_path_centri_abitati_AP, 'w');
% for i = 1:size(dati_centri_abitati_AP, 1)
%     fprintf(fileID_AP, '%s %.4f %.4f %d %.1f %.1f\n', dati_centri_abitati_AP{i,1}, dati_centri_abitati_AP{i,2}, dati_centri_abitati_AP{i,3}, dati_centri_abitati_AP{i,4}, dati_centri_abitati_AP{i,5}, dati_centri_abitati_AP{i,6});
% end
% fclose(fileID_AP);

% % salvataggio dataset città esterne alle placche
% file_path_centri_abitati_fuori = 'Città Italiane esterne.txt';
% fileID_fuori = fopen(file_path_centri_abitati_fuori, 'w');
% for i = 1:size(dati_centri_abitati_fuori, 1)
%     fprintf(fileID_fuori, '%s %.4f %.4f %d %.1f %.1f\n', dati_centri_abitati_fuori{i,1}, dati_centri_abitati_fuori{i,2}, dati_centri_abitati_fuori{i,3}, dati_centri_abitati_fuori{i,4}, dati_centri_abitati_fuori{i,5}, dati_centri_abitati_fuori{i,6});
% end
% fclose(fileID_fuori);
