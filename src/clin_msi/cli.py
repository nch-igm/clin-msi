import os

from churchill_design_patterns import EntrypointWorkflow

from .pipelines import train
from .pipelines import predict

class TrainWorkflow(EntrypointWorkflow):
        
    def add_custom_args(self, parser):
        parser.add_argument('--regions-list', type=str, required=True, help="tab separated text file with MSI regions")
        parser.add_argument('--bam-list', type=str, required=True, help="tab separated text file with BAM paths and MSI status")
        parser.add_argument('--reference', type=str, required=True, help="path to reference(.fa or .fasta) file with index (.fai) in the same directory")
        parser.add_argument('--output-dir', type=str, required=True, help="")
        parser.add_argument('--allow-mismatch', action='store_true', help="allows a single base mismatch within the repeat region")
        parser.add_argument('--normalization-scheme', type=str, required=True, choices=['z','std_u'], help="")

    def validate_args(self, regions_list: str, bam_list: str, reference: str, output_dir: str, **kwargs):
        if not os.path.exists(regions_list) and not vcf.startswith("s3://"):
            raise ValueError("regions_list must be a local file or an s3 path.")
        if not os.path.exists(bam_list) and not vcf.startswith("s3://"):
            raise ValueError("bam_list must be a local file or an s3 path.")
        if not os.path.exists(reference) and not vcf.startswith("s3://"):
            raise ValueError("reference must be a local file or an s3 path.")
        if not os.path.exists(output_dir) and not vcf.startswith("s3://"):
            raise ValueError("output_dir must be a local file or an s3 path.")
    
    def workflow(self, wrk_dir: str, regions_list: str, bam_list: str, reference: str, output_dir: str,
                 allow_mismatch: str, normalization_scheme: str) -> None:
        """Download the bam, vcf, and ref fasta and then run the workflow."""
        bam_path = self.flex_input(bam, out_dir=wrk_dir)
        self.flex_input(bam[0:-1] + 'i', out_dir=wrk_dir)
        if bed is not None:
            bed_path = self.flex_input(bed, out_dir=wrk_dir)
        else:
            bed_path = None
        vcf_path = self.flex_input(vcf, out_dir=wrk_dir)
        ref_path = self.flex_input(reference, out_dir=wrk_dir)
        self.flex_input(reference + '.fai', out_dir=wrk_dir)
        # Output name
        if output is not None:
            filtered_vcf = os.path.join(wrk_dir, os.path.basename(output))
            final_dir = os.path.dirname(output)
        elif output_dir is not None:
            filtered_vcf = os.path.join(wrk_dir, "FFPE_filtered_" + os.path.basename(vcf_path))
            final_dir = output_dir
        # run workflow
        parallel_ffpe_filter_workflow.main(
            bam_path=bam_path,
            vcf_path=vcf_path,
            ref_path=ref_path,
            output=filtered_vcf,
            filters=filters,
            bed=bed_path,
            immediate=immediate,
            ignore_sample_names=ignore_sample_names,
            processes=processes)
        # upload to s3 TODO
        self.flex_output(filtered_vcf, out_dir=final_dir)
        #if bed_path is not None:
        #    self.flex_output(filtered_vcf + '.tmb.txt', out_dir=final_dir)
    

class PredictWorkflow(EntrypointWorkflow):
    
    def add_custom_args(self, parser):
        parser.add_argument('--input-file', type=str, required=True, help="path to the tab separated input file with MSI regions")
        parser.add_argument('--bam', type=str, required=True, help="path to BAM file")
        parser.add_argument('--reference', type=str, required=True, help="path to reference(.fa or .fasta) file with index (.fai) in the same directory")
        parser.add_argument('--sample-name', type=str, required=True)
        parser.add_argument('--output-dir', type=str, required=True, help="")
        parser.add_argument('--model-dir', type=str, required=True, help="path to directory containing .pkl files")
        parser.add_argument('--allow-mismatch', action='store_true', help="allows a single base mismatch within the repeat region")
        parser.add_argument('--normalization-scheme', type=str, required=True, choices=['z', 'std_u'], help="")

    def validate_args(self, vcf: str, **kwargs):
        if not os.path.exists(vcf) and not vcf.startswith("s3://"):
            raise ValueError("vcf must be a local file or an s3 path.")
    
    def workflow(self, wrk_dir: str, vcf: str, bam: str, output: str, output_dir: str,
                 bed: str, reference: str, filters: str, immediate: str,
                 ignore_sample_names: bool, processes: int) -> None:
        """Download the bam, vcf, and ref fasta and then run the workflow."""