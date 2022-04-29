#!/usr/bin/env python
# coding: utf-8

# In[7]:


from constellation import *

from random import random
from math import ceil 

#создание объекта типа Constellation, инициализация параметрами группировки Stalink из конфига
constellation = Constellation('Starlink')

#вычисление элементов орбиты для всех КА в начальный момент
constellation.getInitialState()

# определение точек на оси времени, в которые будут проихзводиться расчёты
epochs = (np.linspace(0, 1000, 1001)).tolist()

# расчёт положений всех КА в заданные моменты времени
constellation.propagateJ2(epochs)


# Координаты случайного КА (в инерциальных осях) после этого можно прочитать из constellation.stateEci
satIdx = ceil(constellation.totalSatCount * random())
epochIdx = ceil(len(epochs) * random())
print("Положение КА-" + str(satIdx) + " на эпоху " + str(epochs[epochIdx]) + ":")
print(constellation.stateEci[satIdx, :, epochIdx])


# In[ ]:




