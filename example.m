clc
clear

tic
%создание объекта типа Constellation, инициализация параметрами группировки Stalink из конфига
constellation = Constellation('Starlink');

%вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState();

% определение точек на оси времени, в которые будут проихзводиться расчёты
epochs = (0: 1000: 6000);

% расчёт положений всех КА в заданные моменты времени
constellation.propagateJ2(epochs);
toc