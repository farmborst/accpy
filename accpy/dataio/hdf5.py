# -*- coding: utf-8 -*-
''' accpy.dataio.h
author:
    Felix Kramer
'''
from h5py import File as h5pyFile
from time import strftime


class struct(object):
    def __init__(self, **entries):
        self.__dict__.update(entries)


def h5save(filename, verbose=False, timestamp=True, **namesandvariables):
    ''' save dataset to hdf5 format
    input:
        - desired filename as string
        - names = values of variables to be saved
    return:
        - saves data to "timestamp_filename.hdf5 in working directory"
        - complete filename is returned
    usage:
        1.  recommended
                datadict = {'a' : 2,
                            'b' : 'foo',
                            'c' : 1.337,
                            'd' : [1, 2, 'c']}
                h5save(filename, True. **datadict)
        2.  alternative 
                a=2, b='foo', c=1.337, d=[1, 2, 'c']
                h5save(filename, True. a=a, b=b, c=c, d=d)
                accepted datatypes:
                    - int   -> numpy.int64
                    - str   -> str
                    - float -> numpy.float64
                    - list  -> numpy.ndarray of:
                                - np.string__   if >0 string
                                - np.float64    if >0 float
                                - np.int64      if only ints
        
    '''
    if timestamp:
        filename = strftime('%Y%m%d%H%M%S') + '_' + filename + '.hdf5'
    else:
        filename += '.hdf5'
    hdf5_fid = h5pyFile(filename, 'w')              
    if verbose:
        print('\n==========================================================')
        print('Beginning to save to %s ...' % filename)
        print('\n----------------------------------------------------------')
    i = 0
    for key, value in namesandvariables.iteritems():
        i += 1
        if verbose:
            print('{1:0>3} Saving values in {0:} ... '.format(key, i))
        hdf5_fid.create_dataset(key.encode('utf8').replace('/', '|'), data=value)
    if verbose:
        print('\n----------------------------------------------------------')
        print('... finished saving to %s !' % filename)
        print('\n==========================================================')
    hdf5_fid.close()
    return filename


def h5load(filename, verbose):
    ''' h5load(filename, verbose)
    input:
        - filename (as string)
        - names of variables and values to be loaded (as string)
    return:
        - dictionary of saved data
    notice:
        use with files saved with accpy.dataio.save
    '''
    if filename[-5:] != '.hdf5':
        filename += '.hdf5'
    fid = h5pyFile(filename, 'r')
    data = {}
    if verbose:
        print('\n==========================================================')
        print('Beginning to load from %s ...' % filename)
        print('\n----------------------------------------------------------')
    i = 0
    for key in fid:
        i += 1
        try:
            data[key] = fid[key].value
        except:
            subdata = {}
            for subkey in fid[key]:
                subdata[subkey] = fid[key][subkey].value
            data[key] = subdata
        if verbose:
            print('{2:0>3} Loading values from {0:} {1:} ... '.format(key, type(data[key]), i))
    fid.close()
    if verbose:
        print('\n----------------------------------------------------------')
        print('... finished loading from %s !' % filename)
        print('\n==========================================================')
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
    listofvars = list(fid['listofvars'].value)
    listofvals = []
    for var in listofvars:
        listofvals.append(fid[var].value)
    fid.close()
    return listofvars, listofvals
