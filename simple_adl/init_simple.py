#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Sidney Mau"

# Python libraries
import os
import yaml

#------------------------------------------------------------------------------

config = {'survey': {
                     'name'    : 'des',
                     'fracdet' : None
                    },
          'band_1' : 'g',
          'band_2' : 'r',
          'catalog' : {
                       'profile'   : None,
                       'nside'     : 32,
                       'mag_max'   : 23.5,
                       'basis_1'   : 'ra',
                       'basis_2'   : 'dec',
                       'mag'       : '{}mag',
                       'mag_dered' : '{}mag_dered',
                       'mag_err'   : '{}err'
                      },
          'isochrone': { # We may consider leaving this fixed for simplicity
                        'name'        : 'Bressan2012',
                        'survey'      : 'des',
                        'age'         : 12,
                        'metallicity' : 0.0001
                       },
          'output' : {
                      'results_dir'    : 'results_dir',
                      'log_dir'        : 'log_dir',
                      'save_dir'       : 'plot_dir',
                      'candidate_list' : 'candidate_list.csv'
                     },
          'batch' : {'max_jobs' : 20}
         }

def init_dirs(cfg):
    results_dir = os.path.join(os.getcwd(), cfg['output']['results_dir'])
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    log_dir = os.path.join(os.getcwd(), cfg['output']['log_dir'])
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)

    save_dir = os.path.join(os.getcwd(), cfg['output']['save_dir'])
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

#------------------------------------------------------------------------------

if __name__ == '__main__':
    cfgfile = 'config.yaml'
    with open(os.path.join(os.getcwd(), cfgfile), 'w') as outfile:
        yaml.dump(config, outfile, default_flow_style=False)

    init_dirs(config)
