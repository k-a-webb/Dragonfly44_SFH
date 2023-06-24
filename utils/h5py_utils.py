import h5py
import numpy as np

# functions for saving data to h5py files
def add_data( hf, name, data, string=False, overwrite=False ):

    if string:
        hf = add_data_string( hf, name, data, overwrite )
    else:
        if name not in hf.keys():
            hf[name] = data
        else:
            if overwrite:
                del hf[name]
                hf[name] = data

    return hf

def add_data_string( hf, name, data, overwrite=False ):
    if (name not in hf.keys()) or overwrite:

        if name in hf.keys():
            del hf[name]

        if np.ndim( data )==1:
            size = len(data)
            sdat = hf.create_dataset(name, (size,), dtype=h5py.special_dtype(vlen=str))
            for i in range(size):
                sdat[i] = data[i]

        elif np.ndim( data )>1:
            shape = np.shape(data)
            sdat = hf.create_dataset(name, shape, dtype=h5py.special_dtype(vlen=str))
            for i in range(shape[0]):
                for j in range(shape[1]):
                    sdat[i,j] = data[i,j]

    return hf

def add_group( hf, name ):
    if name not in hf.keys():
        newgroup = hf.create_group(name)
    else:
        newgroup = hf[name]
    return newgroup
