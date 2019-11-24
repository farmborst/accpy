# -*- coding: utf-8 -*-
''' accpy.dataio.hdf5
author:
    Felix Kramer
'''
from h5py import File as h5pyFile
from h5py._hl.group import Group
from h5py._hl.dataset import Dataset
from time import strftime
from numpy import ndarray


class struct(object):
    def __init__(self, **entries):
        self.__dict__.update(entries)


def h5save(filename, datadict, timestamp=True):
    ''' save dataset to hdf5 format (for load see print(h5load.__doc__))
    input:
        - desired (path/)filename as string
        - dictionary of data
    return:
        - saves data to "(path/)timestamp_filename.hdf5"
        - complete (path/)filename is returned
    usage-example:
                datadict = {'dataset1' : {'x': array(...), 'y': array(...)},
                            'dataset2' : {'x': array(...), 'y': array(...), 'yerr': array(...)},
                            'parameter1' : 1.337,
                            'list1' : [1, 2, 'c']}
                h5save(filename, True. datadict)
    '''
    def dict2h5(datadict, h5id):
        for key, val in datadict.items():
            
            if isinstance(key, bytes):
                key = key.decode().replace('/', '|')
            else:
                key = key.replace('/', '|')
            
            if isinstance(val, (list, tuple, str, bytes, int, float, ndarray)):
                h5id.create_dataset(key, data=val)
            elif isinstance(val, (dict)):
                hdf5_subid = h5id.create_group(key)
                dict2h5(val, hdf5_subid)
            else:
                raise Exception('Data of type {} is not yet supported, sorry for that!'.format(type(val)))
        return

    if timestamp:
        path = '/'.join(filename.split('/')[:-1] + [''])
        filename = strftime('%Y%m%d%H%M%S') + '_' + filename.split('/')[-1]
        filename = path + filename
    if filename[-5:] != '.hdf5':
        filename += '.hdf5'
    hdf5_fid = h5pyFile(filename, 'w')
    dict2h5(datadict, hdf5_fid)
    hdf5_fid.close()
    return filename


def h5load(filename):
    ''' h5load(filename, verbose)
    input:
        - filename (as string) of h5save savedfile
        - desired verbosity
    return:
        - dictionary of saved data
    ALTERNATIVE:
        if the dataset is too large for memory it is also possible to work with it on disk:
        >>> import h5py
        >>> data = h5py.File(filename, 'r')
    '''
    def h52dict(h5id, datadict):
        for key, val in h5id.items():
            if isinstance(val, (Dataset)):
                datadict[key] = h5id[key][()]
            elif isinstance(val, (Group)):
                datadict[key] = {}
                h52dict(h5id[key], datadict[key])
            else:
                raise Exception('Data of type {} is not yet supported, sorry for that!'.format(type(val)))
        return

    if filename[-5:] != '.hdf5':
        filename += '.hdf5'

    data = {}
    hdf5_fid = h5pyFile(filename, 'r')
    h52dict(hdf5_fid, data)
    hdf5_fid.close()
    return data


def confsave(filename, listofvars, listofvals):
    # working with two lists as dictionarys do not accept numpy arrays
    hdf5_fid = h5pyFile(filename, 'w')
    hdf5_fid.create_dataset('listofvars', data=listofvars)
    for var, val in zip(listofvars, listofvals):
        hdf5_fid.create_dataset(var, data=val)
    hdf5_fid.close()


def confload(filename):
    fid = h5pyFile(filename, 'r')
    listofvars = list(fid['listofvars'][()])
    listofvals = []
    for var in listofvars:
        listofvals.append(fid[var][()])
    fid.close()
    return listofvars, listofvals
