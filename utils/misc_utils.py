import numpy as np

__all__ = [ ]

def flatten( result, showpars, q=[0.025,0.995], sample_weight=None, **extras):
    from .transforms import chain_to_param

    if 'weights' in result.keys(): w = np.copy( result['weights'] )
    else: w = None

    for i,par in enumerate(showpars):

        x = chain_to_param( result['chain'], result['theta_index'], par )
        x = np.squeeze( x )

        qs = weighted_quantile( x, quantiles=q, sample_weight=w )

        if i==0:
            flatchain = x
            brange = qs
        else:
            flatchain = np.vstack([ flatchain, x ])
            brange = np.vstack([ brange, qs ])

    return flatchain.T, w, brange


def weighted_quantile(values, quantiles, sample_weight=None,
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)


def get_residualerror(ys, ys_true):
    resid = np.sqrt( (ys - ys_true)**2 )
    std = np.std(resid, ddof=1)
    return std

def get_data_x3_from_data( table, cols ):
    if type(table) is dict:
        x,em_x,ep_x = [ table[key] for key in cols]
        ex = np.array([[em_x],[ep_x]])
        x = np.array(x)
    else:
        x,em_x,ep_x = table[cols].values.T
        ex = np.array([em_x,ep_x])
    return x, ex


def add_error( arr, error_to_add ):
    new_arr = []
    for xx in arr:
        if np.isnan(xx): new_arr.append( error_to_add )
        else:            new_arr.append( xx + error_to_add )
    return np.array(new_arr)
