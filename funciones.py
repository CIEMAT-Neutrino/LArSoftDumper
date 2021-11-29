### Dependences
import uproot
import numpy as np
import matplotlib.pyplot as plt

def print_picture(asroot,h,w,op=False):

    TPC_picture_root=asroot
    TPC_picture_nump =TPC_picture_root.to_numpy()
    TPC_picture_matr =np.reshape(TPC_picture_nump,(h,w));
    if op:
        TPC_picture_matr[TPC_picture_matr==0]=-100;
        
    plt.imshow(TPC_picture_matr,cmap = 'jet', interpolation='none',aspect=w/h);
    # plt.colorbar()

