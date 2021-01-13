import os, sys, glob
from configparser import SafeConfigParser

# def read_config(config_file):
#   config = SafeConfigParser()
#   config.read(config_file)
#   return config


def get_read_counts(bam_file_path, sample_name, optec_dir):
    read_count_script = "{}/msi_read_count.py".format(__location__)
    read_count_file = os.path.join(optec_dir, 'msi_count_table_{}.xlsx'.format(sample_name))

    cmd = ['python', read_count_script, bam_file_path, sample_name, optec_dir]
    cmd = ' '.join(cmd)

    return (cmd, read_count_file)


# def normalize_read_counts(read_count_file, sample_name, optec_dir):
#     normalize_read_count_script = "{}/msiPlotter/deployment_scripts/parseRaw.py".format(__location__)
#     normalized_read_count_file = os.path.join(optec_dir, 'zscore_normalized_{}.csv'.format(sample_name))
#
#     cmd = ['python', normalize_read_count_script, read_count_file, sample_name]
#     cmd = ' '.join(cmd)
#
#     return (cmd, normalized_read_count_file)

def plot_msi(normalized_read_count_file, sample_name, optec_dir):
    plot_msi_script = "{}/msiPlotter/deployment_scripts/plotShapTrace.py".format(__location__)
    training_data = '{}/msiPlotter/deployment_scripts/training_data.csv'
    model_files = '{}/shapPlotter/mod_071019_frozen'

    cmd = ['python', plot_msi_script, normalized_read_count_file, model_files, training_data, sample_name, optec_dir]
    cmd = ' '.join(cmd)

    return cmd


def plot_shap(normalized_read_count_file, sample_name, optec_dir):
    plot_shap_script = "{}/shapPlotter/run_MSI_model.sh".format(__location__)

    cmd = ['bash', plot_shap_script, normalized_read_count_file, sample_name, optec_dir, __location__]
    cmd = ' '.join(cmd)

    return cmd

# def plot_qc(mos_depth_path, optec_dir, sample_name):
#     plot_qc_script = "{}/QC/PanelQuality.R".format(__location__)
#     optec_bed = "{}/QC/OPTEC_Genes.bed".format(__location__)
#     sample_mosdepth_path = os.path.join(mos_depth_path, sample_name) + '.Recalibrated'
#
#     cmd = ['Rscript', plot_qc_script, sample_mosdepth_path, optec_bed, '2,3,4', optec_dir]
#     cmd = ' '.join(cmd)
#
#     return cmd

def create_final_report(sample_name, optec_dir):
    final_report_script = "{}/QC/render_report_MSI.Rscript".format(__location__)

    cmd = ['Rscript', final_report_script, optec_dir, sample_name]
    cmd = ' '.join(cmd)

    return cmd



def main(bam_file_path, mos_depth_path, optec_dir):
    sample_name = os.path.basename(bam_file_path).replace('.Recalibrated.bam', '')


    cmd, read_count_file = get_read_counts(bam_file_path, sample_name, optec_dir)
    print (read_count_file)
    os.system(cmd)

    # cmd, normalized_read_count_file = normalize_read_counts(read_count_file, sample_name, optec_dir)
    # print (normalized_read_count_file)
    # os.system(cmd)

    cmd = plot_msi(normalized_read_count_file, sample_name, optec_dir)
    os.system(cmd)

    cmd = plot_shap(normalized_read_count_file, sample_name, optec_dir)
    os.system(cmd)

    # cmd = plot_qc(mos_depth_path, optec_dir, sample_name)
    # os.system(cmd)
    #
    # cmd = create_final_report(sample_name, optec_dir)
    # os.system(cmd)


if __name__ == '__main__':
    __location__ = os.path.abspath(os.path.dirname(__file__))
    #read config
    config_file = sys.argv[1]
    config = read_config(config_file)
    #set directories
    output_dir = config.get("ANALYSIS_PARAMETERS", "OUTPUT_DIR")
    bam_dir = os.path.join(output_dir, 'BAMs')
    qc_dir = os.path.join(output_dir, 'QC')
    mosdepth_dir = os.path.join(qc_dir, 'MosDepth')
    optec_dir = os.path.join(output_dir, 'OPTEC_RESULTS')
    if not os.path.exists(optec_dir):
        os.makedirs(optec_dir)
    # copy necessary files to optec directory
    # os.system('cp {}/QC/MSI-Report.Rmd {}'.format(__location__, optec_dir))
    # os.system('cp {}/QC/NCHLogo-Top.png {}'.format(__location__, optec_dir))
    #create report for each bam file
    for file in glob.glob(os.path.join(bam_dir, '*.Recalibrated.bam')):
        main(file, mosdepth_dir, optec_dir)
        os.system('rm -r {}'.format(os.path.join(optec_dir, 'bam_reads')))