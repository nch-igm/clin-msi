
"""
        How to use
------------------------------

Python version = 3.7

python plotShapTrace.py zscore_normalized.csv sample_1121-ec

zscore_normalized.csv = path to the design matrix generated by the script parseRaw
sample_1121-ec = the name of the sample in the design matrix

The location of the training data file should be fixed.
It can be found in the deployment_scripts folder

Should take ~ 2 minutes to run

"""





import matplotlib.ticker as ticker
import xgboost as xgb
import pandas as pd
import shap
import pickle
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import numpy as np
import argparse
import sys


# function to obtain the shap and probabilities
def evaluate_test_sample(test_sample, models_list):
    # change the feature names to match jeff's model
    if type(test_sample) is pd.core.frame.Series:
        test_sample = test_sample.to_frame(0).T
    # run each model on a single test_sample
    test_scores = []
    shap_values = []
    for xgb_m in models_list:
        # change the columns to match
        test_sample = test_sample[xgb_m.get_booster().feature_names]
        # compute probability
        prob = xgb_m.predict_proba(test_sample)
        test_scores.append(prob)
        # explain decision via SHAP
        explainer = shap.TreeExplainer(xgb_m)
        s_exp = explainer.shap_values(test_sample)
        shap_values.append(s_exp)
    # transform test_scores by their appropriate axis
    # i.e. if MSS then convert probability to be scaled below 0.5
    arry_pred = np.asarray(test_scores)
    transformed_prb = []
    for i in arry_pred:
        for col in i:
            which_axis_max = col.argmax()
            if which_axis_max == 0:
                transformed_prb.append(1 - col[0])
            else:
                transformed_prb.append(col[1])
    # compute final average probability and other variables
    final_probability = np.mean(transformed_prb)
    standard_dev = np.std(transformed_prb)
    avg_shap_values = np.asarray(shap_values).mean(axis=0).reshape(-1)
    std_shap_values = np.asarray(shap_values).std(axis=0).reshape(-1)
    # output for plotting function
    return final_probability, standard_dev, avg_shap_values, std_shap_values

def fmt(x, pos):
    b = '{:.2f}'.format(x)
    return b





def plot_shap_trace(training_data, test_sample, final_probability, standard_dev, avg_shap_values, std_shap_values, sample_name, output_dir):
    # name the figure
    plt_name = os.path.join(output_dir, 'shape_trace_{}.pdf'.format(sample_name))
    # choose a std
    # shp_thrs = std_shap_values.mean()
    shp_thrs = 1
    # incrementor
    c = 0
    # setup figure
    fig, ax = plt.subplots(7,3)
    ax = ax.flatten()
    # yvalues for dataframe
    y = training_data['y']
    # bound the color coding for the shap values
    _v_min = avg_shap_values.min()
    _v_max = avg_shap_values.max()
    # _v_min = std_shap_values.min()
    # _v_max = std_shap_values.max()
    # store variation for each class
    var_mss = []
    var_msi = []
    # keep track of position
    points_x = []
    points_y_mss = []
    points_y_msi = []
    points_y_test = []
    # custom traversal
    column_walker_0 = np.arange(0,630,30)
    column_walker_1 = np.arange(30,660,30)
    # loop over each feature in the training set and test data
    for col in training_data.columns[training_data.columns != 'y']:
        # house the feature in a new dataframe
        tmp_df = pd.DataFrame({'score': training_data[col],'y':y})
        MSI_mean = 0
        MSS_mean = 0
        MSI_std = 0
        MSS_std = 0
        ymss = np.where(training_data['y'].values == 'MSS')[0]
        ymsi= np.where(training_data['y'].values == 'MSI')[0]
        pos = int(col.split("_")[0])
        marker = "-".join(col.split("_")[1:])
        for g,r in tmp_df.groupby('y'):
            if g == 'MSI':
                MSI_mean = r['score'].mean()
                MSI_std = r['score'].std()
            else:
                MSS_mean = r['score'].mean()
                MSS_std = r['score'].std()
        points_x.append(pos)
        points_y_mss.append(MSS_mean)
        points_y_msi.append(MSI_mean)
        points_y_test.append(test_sample[col])
        var_mss.append(MSS_std)
        var_msi.append(MSI_std)
        if pos == 29:
            # plt the mean of the training data
            try:
                ax[c].plot(points_x,points_y_msi, '-', color='tab:orange',alpha=0.2, label='MSI')
                ax[c].plot(points_x,points_y_mss, '-', color='tab:purple', alpha=0.2, label='MSS')
            except Exception as ex:
                 print('Plotting training data mean failed {}').format(ex)
            # plt the variation of the training data
            try:
                ax[c].fill_between(points_x, np.asarray(points_y_mss)-var_mss, np.asarray(points_y_mss)+var_mss, alpha=0.5,color='tab:purple')
                ax[c].fill_between(points_x, np.asarray(points_y_msi)-var_msi, np.asarray(points_y_msi)+var_msi, alpha=0.5,color='tab:orange')
            except Exception as ex:
                 print('Plotting training data variation failed {}').format(ex)
            # collect the shap data
            try:
                shap_vector = avg_shap_values[column_walker_0[c]:column_walker_1[c]]
                shap_vector_std = std_shap_values[column_walker_0[c]:column_walker_1[c]]
                points_y_test = np.asarray(points_y_test)
                points_x = np.asarray(points_x)
                shap_pos = np.where((shap_vector > 0) & (shap_vector_std <= shp_thrs))[0]
                shap_neg = np.where(shap_vector < 0 & (shap_vector_std <= shp_thrs))[0]
            except Exception as ex:
                print('SHAP vector might be missing information'.format(ex))
            # plot the test data
            try:
                ax[c].plot(points_x, points_y_test, '-',color='k', linewidth=1, snap=False, label='sample')
            except Exception as ex:
                print('Test data trace failed to plot {}'.format(ex))
            # plot shap values
            try:
                ax[c].scatter(points_x,points_y_test, zorder=10,c=shap_vector,cmap=cm.PuOr_r,vmin=_v_min,vmax=_v_max, alpha=1, edgecolors='k',linewidths=0.1, s=5000*np.abs(shap_vector))
                # ax.scatter(points_x[shap_neg],points_y_test[shap_neg], zorder=10, c=shap_vector[shap_neg],cmap=cm.Purples_r,vmin=_v_min,vmax=_v_max, alpha=1, edgecolors='k',linewidths=0.1, s=5000*np.abs(shap_vector[shap_neg]))
            except Exception as ex:
                print('Something is wrong with SHAP values {}'.format(ex))
            # total shap score
            shape_score = avg_shap_values[column_walker_0[c]:column_walker_1[c]].sum()
            # reset for the next marker
            points_x = []
            points_y_mss = []
            points_y_msi = []
            points_y_test = []
            var_mss = []
            var_msi = []
            ax[c].xaxis.set_ticks([])
            ax[c].yaxis.set_ticks([])
            ax[c].set_title('{} Impact: {:.4f}'.format(marker, float(shape_score)) ,fontsize=9, y=1)
            c += 1
    handles, labels = ax[-1].get_legend_handles_labels()
    # add the text information
    norm = mpl.colors.Normalize(vmin=_v_min, vmax=_v_max)
    radii = np.abs(np.linspace(_v_min,_v_max,10))
    color_map = cm.get_cmap('PuOr_r')
    colors = color_map([norm(x) for x in np.linspace(_v_min,_v_max,10)])
    rc = zip(radii,colors)
    circles = []
    spacing = 0
    for x,y in rc:
        circles.append(Circle((0.75+spacing, 0.95), x, facecolor=y,edgecolor='k', linewidth=0.5, alpha=1,transform=fig.transFigure, figure=fig, zorder=1000))
        spacing+=0.01

    # add the circles
    fig.patches.extend(circles)
    fig.legend(handles, labels, loc='upper center')
    fig = mpl.pyplot.gcf()
    fig.set_size_inches(8.5, 11)

    if not sample_name and final_probability >= 0.5:
        plt.suptitle('Score: {:.3f} (MSI)'.format(sample_name, final_probability), x=0.25)
    elif not sample_name and final_probability < 0.5:
        plt.suptitle('Score: {:.3f} (MSS)'.format(sample_name, final_probability), x=0.25)
    elif sample_name and final_probability >= 0.5:
        plt.suptitle('{}\n Score: {:.3f} (MSI)'.format(sample_name, final_probability), x=0.25)
    elif sample_name and final_probability < 0.5:
        plt.suptitle('{}\n Score: {:.3f} (MSS)'.format(sample_name, final_probability), x=0.25)
    plt.text(2.65,2, 'MSS', fontsize=8, transform=ax[0].transAxes)
    plt.text(3.1,2, 'MSI', fontsize=8, transform=ax[0].transAxes)
    plt.text(2.74,1.5, 'Marker Impact', fontsize=8,transform=ax[0].transAxes)
    plt.savefig(plt_name, dpi=500, bbox_inches='tight')
    print('path to figure: {}'.format(plt_name))


def main():
    parser = argparse.ArgumentParser(prog='plot_shap_trace',description='Plot trace of MSI data with SHAP values for each feature.')
    parser.add_argument('path', metavar='/path/to/zscore/normalized/msi/data/csv', type=str, nargs=1, help='path to the test data should be in csv format')
    parser.add_argument('path_to_models', type=str, nargs=1, help='path to model directory')
    parser.add_argument('path_to_training_data', type=str, nargs=1, help='path to training data')
    parser.add_argument('sample_name', type=str, nargs=1, help='name of sample', default='test_sample')
    parser.add_argument('output_dir', type=str, nargs=1, help='path to output dir')
    args = parser.parse_args()
    training_data = pd.read_csv(args.path_to_training_data[0]).set_index('Unnamed: 0')
    test_data = pd.read_csv(args.path[0]).set_index('Unnamed: 0')
    test_sample = test_data.loc[args.sample_name[0]]
    output_dir = args.output_dir[0]

    # prepare model
    path_to_models = args.path_to_models[0]
    model_suffix = 'xgb'

    # generate the file paths for each model
    files = os.listdir(path_to_models)
    globbed_paths = [os.path.join(path_to_models, f) for f in files if f.startswith(model_suffix)]

    # load the model and store in a list
    models_list = []
    for glob in globbed_paths:
        m = pickle.load(open(glob, 'rb'))
        models_list.append(m)
    # create plots
    final_probability, standard_dev, avg_shap_values, std_shap_values = evaluate_test_sample(test_sample, models_list)
    plot_shap_trace(training_data, test_sample, final_probability, standard_dev, avg_shap_values, std_shap_values, args.sample_name[0], output_dir)

def test():
    path_to_training_data = '/Users/axr102/Desktop/Workspace/machine_learning_projects/msi/data/exp4_train_data/marker_matrix_z-normalized.csv'
    path_to_test_data = '/Users/axr102/Desktop/Workspace/machine_learning_projects/msi/data/exp5/test_marker_matrix_z-normalized.csv'
    test_data = pd.read_csv(path_to_test_data).set_index('Unnamed: 0')
    training_data = pd.read_csv(path_to_training_data).set_index('Unnamed: 0')
    test_sample = test_data.loc['Exp5B_C2_1185_tumor_tissue']
    prepare_model()
    final_probability, standard_dev, avg_shap_values, std_shap_values = evaluate_test_sample(test_sample)
    plot_shap_trace(training_data, test_sample, final_probability, standard_dev, avg_shap_values, std_shap_values, 'Exp5B_C2_1185_tumor_tissue')

if __name__ == '__main__':
    main()
