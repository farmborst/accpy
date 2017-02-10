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


def h5save(filename, showinfo, **namesandvariables):
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
                save(filename, True. **datadict)
        2.  alternative 
                a=2, b='foo', c=1.337, d=[1, 2, 'c']
                save(filename, True. a=a, b=b, c=c, d=d)
                accepted datatypes:
                    - int   -> numpy.int64
                    - str   -> str
                    - float -> numpy.float64
                    - list  -> numpy.ndarray of:
                                - np.string__   if >0 string
                                - np.float64    if >0 float
                                - np.int64      if only ints
        
    '''
    timstamp = strftime('%Y%m%d%H%M%S')
    filename = ''.join([timstamp, '_', filename, '.hdf5'])
    hdf5_fid = h5pyFile(filename, 'w')              
    if showinfo:
        print('\n==========================================================')
        print('Beginning to save to %s ...' % filename)
        print('\n----------------------------------------------------------')
    for key, value in namesandvariables.iteritems():
            if showinfo:
                print('Saving values in %s ... ' % key)
            hdf5_fid.create_dataset(key, data=value)
    if showinfo:
        print('\n----------------------------------------------------------')
        print('... finished saving to %s !' % filename)
        print('\n==========================================================')
    hdf5_fid.close()
    return filename


def h5load(filename, showinfo):
    ''' load(filename,variables(e.g. 'a', 'b', ...))
    input:
        - filename (as string)
        - names of variables and values to be loaded (as string)
    return:
        - dictionary of saved data
    notice:
        use with files saved with accpy.dataio.save
    '''
    if filename[-5:] != '.hdf5':
        filename = ''.join([filename, '.hdf5'])
    fid = h5pyFile(filename, 'r')
    data = {}
    if showinfo:
        print('\n==========================================================')
        print('Beginning to load from %s ...' % filename)
        print('\n----------------------------------------------------------')
        for key in fid:
            data[key] = fid[key].value
            print('Loading values from {0:} {1:} ... '.format(key, type(data[key])))
        print('\n----------------------------------------------------------')
        print('... finished loading from %s !' % filename)
        print('\n==========================================================')
    else:
        for key in fid:
            data[key] = fid[key].value
    fid.close()
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
