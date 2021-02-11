import re
import pandas as pd
import numpy as np


def repeat_finder(s):
    #Taken from https://stackoverflow.com/questions/9079797/detect-repetitions-in-string
    r = re.compile(r"(.+?)\1+")
    for match in r.finditer(s, overlapped=True):
        yield (match.group(1), len(match.group(0))/len(match.group(1)))


def parse_input_file(input_file):
    location_list = []
    open_file = open(input_file, 'r')
    for line in open_file:
        chr, start, stop = line.split('\t')
        location_list.append([str(chr), int(start), int(stop)])

    return location_list


def parse_raw_data(repeat_df, sample_name):
    # check to make sure input is a dataframe
    assert type(repeat_df) is pd.core.frame.DataFrame, "Invalid input, must be a pandas dataframe object"

    # make sure we have an expected column

    assert repeat_df.columns[0] == 'Repeat_Length', "Repeat dataframe is missing an expected column, Repeat_Length"

    # z-score transform the dataframe on the sample itself, not on a collection of samples
    # so there is no chance for overfitting.
    mu = repeat_df.mean(axis=0)
    std = repeat_df.std(axis=0)
    df_data = (repeat_df.iloc[:, 1:] - mu) / std  # ignore repeat_length column
    df_data.fillna(0, inplace=True)

    # turn this into a single sample and add the proper feature names
    raveled_data = df_data.T.values.ravel()
    print(df_data)
    columns = np.repeat(repeat_df['Repeat_Length'] + '_', df_data.shape[1]) + np.repeat(df_data.columns,
                                                                                        df_data.shape[0])
    single_sample_df = pd.DataFrame(raveled_data.reshape(1, -1), columns=columns, index=[sample_name])

    return single_sample_df