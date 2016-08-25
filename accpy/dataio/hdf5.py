# -*- coding: utf-8 -*-
''' accpy.dataio.h
author:
    Felix Kramer
'''
from h5py import File as h5pyFile
from time import strftime


def save(filename, showinfo, **namesandvariables):
    ''' save(filename,variables(e.g. a=a, b=b, ...))
    input:
        - desired filename as string
        - names = values of variables to be saved
    return:
        - saves data to "datetime_filename.hdf5 in working directory"
    notice:
        accepted datatypes:
            -
    '''
    timstamp = strftime('%Y%m%d%H%M%S')
    filename = ''.join([timstamp, '_', filename, '.hdf5'])
    hdf5_fid = h5pyFile(filename, 'w')
    if showinfo:
        print('\n==========================================================')
        print('Beginning to save to %s ...' % filename)
        print('\n----------------------------------------------------------')
        for key, value in namesandvariables.iteritems():
            print('Saving values in %s ... ' % key)
            hdf5_fid.create_dataset(key, data=value)
        print('\n----------------------------------------------------------')
        print('... finished saving to %s !' % filename)
        print('\n==========================================================')
    else:
        for key, value in namesandvariables.iteritems():
            hdf5_fid.create_dataset(key, data=value)
    hdf5_fid.close()
    return filename


def load(filename, showinfo, *varnames):
    ''' load(filename,variables(e.g. 'a', 'b', ...))
    input:
        - desired filename (as string)
        - names of variables and values to be loaded (as string)
    return:
        - wanted variables are loaded from filename
    notice:
        use with files saved with mypy.save
        -
    '''
    if filename[-5:] != '.hdf5':
        filename = ''.join([filename, '.hdf5'])
    fid = h5pyFile(filename, 'r')
    data = []
    if showinfo:
        print('\n==========================================================')
        print('Beginning to load from %s ...' % filename)
        print('\n----------------------------------------------------------')
        for arg in varnames:
            print('Loading values from %s ...' % arg)
            data.append(fid[arg].value)
        print('\n----------------------------------------------------------')
        print('... finished loading from %s !' % filename)
        print('\n==========================================================')
    else:
        for arg in varnames:
            data.append(fid[arg].value)
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
