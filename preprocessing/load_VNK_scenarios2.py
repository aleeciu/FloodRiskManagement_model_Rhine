from funs_geo import GDAL2Numpy
from funs_generate_network import get_network
import numpy as np
import pandas as pd
import os

' Open VNK flooding scenarios and save it to excel file'

G, branches, upstream_node, dike_nodes, bif_nodes = get_network(Musk_params=True,
                                                                fragcurves=True,
                                                                Losses=False)
damages = pd.read_excel('./data/damages/VNK_damages.xls', index_col=0,
                        usecols=["Locatie", "Schade 2011 (M EUR)", "Unnamed: 3", "Unnamed: 4", "Unnamed: 5",
                                 "Slachtoffers 2011", "Unnamed: 10", "Unnamed: 11", "Unnamed: 12"])

damages.columns = ['{}_Schade'.format(TP) for TP in ['-1', '0', '+1', '+2']
                   ] + ['{}_Slachtoffers'.format(TP) for TP in ['-1', '0', '+1', '+2']]

Npath = r"N:\My Documents\FloodScenarios\\"
data = pd.read_excel(Npath + "Overview.xlsx")  # use the original

#source = './Data/dijkring49/03_Doesburg_noord_49f11_tp/dmax.asc'

data_dic = {name: {'Area': [], 'Volume': [], 'losses': [], 'h': []}
            for name in dike_nodes}


for d in dike_nodes:
    naam = G.node[d]['VNKNAAM']
    print(d, naam)

    for i in np.where(data['Doorbraaklocatie'] == naam)[0]:
        atch = data.loc[i, 'Pad waterdiepte bestand'].split(
            'referentieSchade/')[1]

        source = Npath + atch
        if os.path.isfile(source):
            print(os.path.isfile(source))
            (wd_map, geotransform, projection, nodata) = GDAL2Numpy(source)
            pixel_size = geotransform[1]

            wd_map = np.where(wd_map == nodata, np.nan, wd_map)

            inundated_cells = np.where(wd_map == 0, 0, 1)

            Area = np.nansum(inundated_cells) * pixel_size**2
            Volume = np.nansum(wd_map * pixel_size**2)
            losses = damages.loc[naam] * 1e6

            data_dic[d]['h'].append(Volume / float(Area))
            data_dic[d]['Area'].append(Area)
            data_dic[d]['Volume'].append(Volume)
            data_dic[d]['losses'].append(losses)


# for dike in dikelist:
#        for key in data_dic[dike].keys():
#                data_dic[dike][key] = np.unique(data_dic[dike][key])
#
#        name = 'D:/ciullo/Documents/GitHub/systemrisk/src/preprocessing/losses_lookups/{dikename}_lossestableVNK.xlsx'.format(dikename = dike)
#        pd.DataFrame(data_dic[dike]).to_excel(name)

# xmin=geotransform[0]
# ymax=geotransform[3]
# rows,cols=np.shape(dem)
# ymin=ymax-pixel_size*rows

#import matplotlib.pyplot as plt
# plt.imshow(dem)
# plt.colorbar()
