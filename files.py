# -*- coding: utf-8 -*-
"""
Created on Sun Jan  6 15:50:16 2019

@author: Quentin
"""

files = [
"ECMA_Instances_2018-2019/20_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/40_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/60_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/80_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/100_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/120_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/140_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/160_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/180_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/200_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/250_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/300_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/350_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/400_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/450_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/500_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/550_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/600_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/650_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/700_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/750_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/800_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/850_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/900_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/950_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/1000_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/1100_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/1200_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/1300_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/1400_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/1500_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/1600_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/1700_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/1800_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/1900_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/2000_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/2100_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/2200_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/2300_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/2400_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/2500_USA-road-d.BAY.gr",
"ECMA_Instances_2018-2019/20_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/40_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/60_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/80_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/100_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/120_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/140_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/160_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/180_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/200_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/250_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/300_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/350_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/400_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/450_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/500_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/550_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/600_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/650_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/700_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/750_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/800_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/850_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/900_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/950_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/1000_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/1100_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/1200_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/1300_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/1400_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/1500_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/1600_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/1700_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/1800_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/1900_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/2000_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/2100_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/2200_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/2300_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/2400_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/2500_USA-road-d.COL.gr",
"ECMA_Instances_2018-2019/20_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/40_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/60_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/80_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/100_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/120_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/140_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/160_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/180_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/200_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/250_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/300_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/350_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/400_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/450_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/500_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/550_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/600_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/650_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/700_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/750_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/800_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/850_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/900_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/950_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/1000_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/1100_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/1200_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/1300_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/1400_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/1500_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/1600_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/1700_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/1800_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/1900_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/2000_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/2100_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/2200_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/2300_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/2400_USA-road-d.NY.gr",
"ECMA_Instances_2018-2019/2500_USA-road-d.NY.gr"]
