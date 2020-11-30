#!/usr/bin/env python
"""
Compile candidate list from results_dir
"""
__author__ = "Sidney Mau"

import glob
import yaml

import numpy as np
import csv

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=False,default='config.yaml',
                        help='Config file [.yaml].')
    parser.add_argument('--outfile',type=str,required=False,
                        help='Output file [.csv].')
    args = parser.parse_args()

    with open(args.config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)

    basis_1 = cfg['catalog']['basis_1']
    basis_2 = cfg['catalog']['basis_2']
    candidate_list = cfg['output']['candidate_list']

    if args.outfile:
        outfile = args.outfile
    else:
        outfile = candidate_list

    # Parse results from results_dir into a list of values
    results = []
    for file in glob.glob('{}/*.txt'.format(cfg['output']['results_dir'])):
        with open(file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                results.append([float(val) for val in row])
        csvfile.close()
    
    data = np.asarray(results)
    #data = data[np.unique([data[i][-1] for i in range(len(data))], return_index=True)[1]]
    
    ## Create fits columns
    #c0 = fits.Column(name='SIG',          format='E', array=data[:,0])
    #c1 = fits.Column(name=basis_1,        format='E', array=data[:,1])
    #c2 = fits.Column(name=basis_2,        format='E', array=data[:,2])
    #c3 = fits.Column(name='MODULUS',      format='E', array=data[:,3])
    #c4 = fits.Column(name='R',            format='E', array=data[:,4])
    #c5 = fits.Column(name='N_OBS',        format='E', array=data[:,5])
    #c6 = fits.Column(name='N_OBS_HALF',   format='E', array=data[:,6])
    #c7 = fits.Column(name='N_MODEL',      format='E', array=data[:,7])
    #c8 = fits.Column(name='MC_SOURCE_ID', format='E', array=data[:,8])
    #
    ## Write fits output
    #t = fits.BinTableHDU.from_columns([c0, c1, c2, c3, c4, c5, c6, c7, c8])
    #t.writeto(candidate_list, overwrite=True)

    f = open(outfile, 'ab')
    np.savetxt(f, data, delimiter=',', header='SIG,{},{},MODULUS,R,N_OBS,N_OBS_HALF,N_MODEL,MC_SOURCE_ID'.format(basis_1,basis_2))
    f.close()
