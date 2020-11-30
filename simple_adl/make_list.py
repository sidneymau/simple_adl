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
                        help='Config file [.yaml]')
    parser.add_argument('--outfile',type=str,required=False,
                        help='Output file [.csv]')
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
    for file in glob.glob('{}/*.csv'.format(cfg['output']['results_dir'])):
        with open(file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                results.append([float(val) for val in row])
        csvfile.close()
    
    data = np.asarray(results)
    #data = data[np.unique([data[i][-1] for i in range(len(data))], return_index=True)[1]]
    
    f = open(outfile, 'w')
    np.savetxt(f, data, delimiter=',', header='SIG,{},{},MODULUS,R,N_OBS,N_OBS_HALF,N_MODEL,MC_SOURCE_ID'.format(basis_1,basis_2))
    f.close()
