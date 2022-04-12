

import numpy as np
import proj_a as p

el_grs84 = p.Transformacje(model = "wgs84")
plik = "wsp_inp.txt"

tablica = np.genfromtxt(plik, delimiter = ",", skip_header = 4)
rows,cols = np.shape(tablica)

flh = np.zeros((rows,cols))
xy2000 = np.zeros((rows,3))
xy92 = np.zeros((rows,2))
neu = np.zeros((rows,cols))


tablica_ze_wsp = np.zeros((rows,8))

for i in range(rows):
    flh[i] = el_grs84.xyz2flh(tablica[i,0], tablica[i,1], tablica[i,2])
    xy2000[i] = el_grs84.u2000(flh[i,0], flh[i,1], 0.999923, 21)
    xy92[i] = el_grs84.u92(flh[i,0], flh[i,1], 0.9993)
  
    tablica_ze_wsp[i,0:3] = flh[i]
    tablica_ze_wsp[i,3:6] = xy2000[i]

    
np.savetxt("wsp_FLH.txt", tablica_ze_wsp, delimiter=',')
