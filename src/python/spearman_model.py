import pandas as pd
from scipy.stats import spearmanr
from utils.decorator import timer

def spearman_correlation(query, target, alpha=0.05, no_neg=True):
    """Calculate the query against the target via spearman correlation

    Non-parametric: correlation of ranks

    Correlation coefficient (rho) is in the range between -1 and 1 and the p-value

    Parameters
    ----------
    query : Series
    target : Series
    alpha : float, default 0.05
        Significance threshold for spearmanr
    no_neg : bool, default True
        If True, return correlation values from 0 to 1 (no negatives)

    Returns
    -------
    float
        spearman correlation value
    """

    rho, pval = spearmanr(query, target)
    if no_neg:
        return rho if pval < alpha and rho > 0 else 0
    return rho if pval < alpha else 0

def spearman_matches(query, training_set, training_label, top_N=10, **kwargs):
    """Assigns missing lables using labels from best matched cells found by spearman correlation

    Parameters
    ----------
    query_set : DataFrame
        The object containing entries for which labels will be predicted
    training_set : DataFrame
        The object contraining training data to calculate spearman correlation values
    training_label : Series
        The object containing labels for training set
    top_N : int, default 10
        Specify the number of best matches to return
    **kwargs
        Additional parameters for `spearman_correlation` function

    Returns
    -------
    DataFrame
        spearman correlation values and corresponding labels

    """

    alpha = kwargs.get('alpha', 0.05)
    no_neg = kwargs.get('no_neg', True)

    matches = training_set.apply(
        lambda x: spearman_correlation(query=query, target=x, alpha=alpha, no_neg=no_neg),
        axis=1).sort_values(ascending=False)[:top_N]
    matches.name = 'correlation'

    return matches.to_frame().join(training_label, how='inner')

def score_by_mean(matches, labels, score):
    """Evaluate matches by the average score of some value metric

    Parameters
    ----------
    matches : DataFrame
        Candidate matches with score and label
    labels : str
        Name of column containing labels to group DataFrame by
    score : str
        Name of column containing score values

    Returns
    -------
    str
        label with highest average score
    """

    average = matches.groupby(labels)[score].mean().sort_values(ascending=False)
    if average.index.nunique() > 1 and average.nunique() == 1:
        return 'Tie'
    return average.index[0]

def score_by_count(matches, labels):
    """Evaluates matches by the frequency of a predicted label

    Parameters
    ----------
    matches : DataFrame
        Candidate matches with score and label
    labels : str
        Name of column containing labels to group DataFrame by

    Returns
    -------
    str
        label with highest frequency
    """

    count = matches[labels].value_counts(ascending=False)
    if count.index.nunique() > 1 and count.nunique() == 1:
        return 'Tie'
    return count.index[0]

@timer
def get_spearman_predictions(X_test, X_train, Y_train, cluster=None, score='mean', rm_columns=[], **kwargs):
    """Search for spearman correlation matches

    Parameters
    ----------
    X_test : DataFrame
        Test dataset
    X_train : DataFrame
        Training dataset
    Y_train : Series
        Training dataset labels
    cluster : str, optional
        If provided, limit the rows in X_train by running `groupby` on this column.
        Only use values in X_train that share the same value as X_test for this column.
    score : str, default 'mean'
        Score matches by either average correlation value or count of labels
    rm_columns : list of str, optional
        List of columns to ignore
    **kwargs
        Additional parameters for `spearman_matches` function

    Returns
    -------
    list
        Predicted labels
    """

    top_N = kwargs.get('top_N', 10)
    prediction = []

    for i in range(X_test.shape[0]):
        if cluster:
            cluster_id = X_test.iloc[i][cluster]
            subset = X_train[X_train[cluster] == cluster_id]
        else:
            subset = X_train
        matches = spearman_matches(query=X_test.iloc[i].drop(rm_columns),
                                   training_set=subset.drop(rm_columns, axis=1),
                                   training_label=Y_train, top_N=top_N)
        if score == 'mean':
            prediction.append(score_by_mean(matches, labels='TCR_type', score='correlation'))
        if score == 'count':
            prediction.append(score_by_count(matches, labels='TCR_type'))

    prediction = pd.Series(prediction)
    prediction.index = X_test.index
    return prediction