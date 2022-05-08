% ОПИСАНИЕ:
%   isSatVisible = isVisible(constellation, epochIdx)
%   Функция isVisible возвращает столбец признаков нахождения КА в зоне
%   видимости для 1 выбранной эпохи.
% ВХОДНЫЕ ПАРАМЕТРЫ:
%   constellation - объект типа Constellation, содержащий данные об
%   орбитальной конфигурации и координтах каждого КА.
%   epochIdx - индекс выбранной эпохи в массиве точек расчётного времени
% ВЫХОДНЫЕ ЗНАЧЕНИЯ:
%   isSatVisible - вызвращаемый столбец признаков нахождения КА в зоне
%   видимости. Признаки имеют логический тип данных.

function isSatVisible = isVisible(constellation, epochIdx)
    % Чтение данных из файла gatewaysTest.json
    fileName = 'gatewaysTest.json';
    str = fileread(fileName);   % Получение данных в виде набора символов
    beginData = regexp(str, '{'); % Вычисление индексов начала полезных данных
    endData = regexp(str, '}'); % Вычисление индексов окончания полезных данных
    
    % Декодирование полезных данных в структуру
    data = cell2struct(cell(length(beginData),3),["lat","lon","altitude"],2);
    for dataIdx = 1:length(beginData)
        strData = str(beginData(dataIdx) : endData(dataIdx));
        strData = regexprep(strData, ',}', '}'); % Удаление лишних запятых перед '}'
        data(dataIdx) = jsondecode(strData);
    end
    clear str beginData endData strData dataIdx
    
    stationsCount = length(data);            % Количество шлюзов 
    coordStations = zeros(stationsCount, 3); % Инициализация матрицы координат шлюзов
        
    % Вычисление координат шлюзов в прямогульной системе координат с
    % началом отсчёта в центре Земли
    for stationsIdx = 1:stationsCount
        coordStations(stationsIdx, :) = [(constellation.earthRadius + data(stationsIdx).altitude) * ...
                                        cos(data(stationsIdx).lat) * cos(data(stationsIdx).lon), ...
                                        (constellation.earthRadius + data(stationsIdx).altitude) * ...
                                        cos(data(stationsIdx).lat) * sin(data(stationsIdx).lon), ...
                                        (constellation.earthRadius + data(stationsIdx).altitude) * ...
                                        sin(data(stationsIdx).lat)];
    end
  
    % Инициализация матрицы углов места
    elevation = zeros(constellation.totalSatCount, stationsCount);

    % Вычисление матрицы углов места для каждого спутника
    % Угол вычисляется как угол между векторами OG и SG (см. справочный
    % рисунок) за вычитом 90 градусов
    for satIdx = 1:constellation.totalSatCount
        for stationsIdx = 1:stationsCount
            vectorSat = coordStations(stationsIdx) - constellation.state.eci(satIdx, : , epochIdx);
            elevation(satIdx, stationsIdx) = acos(coordStations(stationsIdx, :) * vectorSat'  / ...
                                            (norm(vectorSat) * norm(coordStations(stationsIdx, :)))) - pi / 2;
        end % Окончание цикла по всем спутникам
    end % Окончание цикла по всем шлюзам
    
    % Условие видимости КА: угол места >= 25 градусов
    isSatVisible = any(elevation >= deg2rad(25), 2);
end