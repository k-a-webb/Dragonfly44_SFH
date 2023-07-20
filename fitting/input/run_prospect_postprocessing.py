
import time
import sys
import os
import argparse

from Dragonfly44_SFH.utils.prospect_postprocessing import save_posterior_draws_models, save_bestfit_model, save_mfrac
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--result_file', type=str, default='', help='Name of Prospector output file')
    parser.add_argument('--overwrite', action="store_true", help='Whether to overwrite values')
    parser.add_argument('--n_seds', type=int, default=500, help='Number of SED models to produce')

    args = parser.parse_args()

    result_file = args.result_file
    overwrite = args.overwrite
    n_seds = args.n_seds

    start_time = time.time()

    print('Result file: ', result_file)
    print('  Getting SED models for {} draws from the posterior probability distributions:'.format(n_seds))
    save_posterior_draws_models( result_file, n_seds=n_seds, overwrite=overwrite )

    #save_bestfit_model( result_file, 'draws',  overwrite=overwrite )

    print('  Calculating mass loss for all draws from the posterior probability distributions:')
    save_mfrac( result_file )

    print("  Duration: {:.3f} min".format( (time.time()-start_time)/60. ))
    print()
