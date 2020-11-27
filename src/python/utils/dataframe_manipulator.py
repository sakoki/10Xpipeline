import multiprocessing as mp
import numpy as np
import pandas as pd
import functools

# Functions for commonly used DataFrame manipulations

def combine_columns(TCR_seq, cols, spacer='_', colname=None):
    """Combine multiple columns entries to create a new metalabel
    Entries are automatically converted to string datatype
    The new column is automatically named metalabel unless specified
    
    Parameters
    ----------
    TCR_seq : DataFrame
        TCR sequences with additional metadata
    cols : list of str
        Name of columns to combine
    spacer : str, default '_'
        String used to join different labels
        
    Returns
    -------
    DataFrame
        DataFrame with metalabel added as new column
    """

    TCR_seq = TCR_seq.assign(metalabel = TCR_seq.apply(
        lambda x: f'{spacer}'.join((str(i) for i in x[x.index.isin(cols)] if isinstance(i, str))),
        axis=1)
    )
    if colname:
        TCR_seq.set_axis([*TCR_seq.columns[:-1], colname], axis=1, inplace=True)
    return TCR_seq

def parallelize_dataframe_operation(df, func, n_cores=2, *args, **kwargs):
    """Speed up DataFrame calculations via parallelization
    Breaks the DataFrame into n parts to distribute among n cores
    Runs the function on split dataframes, then obtains the results
    and merges them together.
    
    Note
    ----
    Processes do not share memory and cannot modfiy the same memory concurrently.
    To use parallelization, processes must be independent of others.
    
    Parameters
    ----------
    df : DataFrame
    func : function
    *args
        additional parameters to pass to the parallelized function
    **kwargs
        additional parameters to pass to the parallelized function
        
    Returns
    -------
    DataFrame
    """

    print(f'Parallelizing {func.__name__} on {n_cores} cores...')
    df_split = np.array_split(df, n_cores)  # Break dataframe into n parts for n cores
    pool = mp.Pool(n_cores)  # Spawn n core processes
    func = functools.partial(func, *args, **kwargs)  # Set constant values to function parameters
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df
