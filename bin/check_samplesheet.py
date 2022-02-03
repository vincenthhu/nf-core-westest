#!/usr/bin/env python

# TODO nf-core: Update the script to check the samplesheet
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/westest samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


# TODO nf-core: Update the check_samplesheet function
def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    For a normal/tumor pair:

    In this example (example_pair_fastq.tsv), there are 3 read groups for the normal sample and 2 for the tumor sample.

    subject_id,gender,phenotype,paternal_id,maternal_id,sample,lane,matched_sample,fastq_1,fastq_2
    SUBJECT_ID,XX,0,caseyupper,caseyupperlamb,SAMPLE_ID1,1,SAMPLE_ID2,/samples/normal1_R1.fastq.gz,/samples/normal1_R2.fastq.gz
    SUBJECT_ID,XX,0,caseyupper,caseyupperlamb,SAMPLE_ID1,2,SAMPLE_ID2,/samples/normal2_R1.fastq.gz,/samples/normal2_R2.fastq.gz
    SUBJECT_ID,XX,0,caseyupper,caseyupperlamb,SAMPLE_ID1,3,SAMPLE_ID2,/samples/normal3_R1.fastq.gz,/samples/normal3_R2.fastq.gz
    SUBJECT_ID,XX,1,caseyupper,caseyupperlamb,SAMPLE_ID2,1,SAMPLE_ID1,/samples/tumor1_R1.fastq.gz,/samples/tumor1_R2.fastq.gz
    SUBJECT_ID,XX,1,caseyupper,caseyupperlamb,SAMPLE_ID2,2,SAMPLE_ID1,/samples/tumor2_R1.fastq.gz,/samples/tumor2_R2.fastq.gz

    For an example see:
    https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 10
        # TODO nf-core: Update the column names for the input samplesheet
        HEADER = ["subject_id", "gender", "phenotype", "paternal_id", "maternal_id", "sample", "lane", "matched_sample", "fastq_1", "fastq_2"] # 9 standard columns
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            subject_id, gender, phenotype, paternal_id, maternal_id, sample, lane, matched_sample, fastq_1, fastq_2 = lspl[: len(HEADER)]
            subject_id = subject_id.replace(" ", "_")
            sample = sample.replace(" ", "_")
            matched_sample = matched_sample.replace(" ", "_")
            if not subject_id:
                print_error("subject_id entry has not been specified!", "Line", line)
            if not sample:
                print_error("subject_id entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [subject_id, gender, phenotype, paternal_id, maternal_id, single_end, fastq_1, fastq_2]
            if sample and lane and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = [matched_sample, subject_id, gender, phenotype, paternal_id, maternal_id, "0", fastq_1, fastq_2]
            elif sample and lane and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = [matched_sample, subject_id, gender, phenotype, paternal_id, maternal_id, "1", fastq_1, fastq_2]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = { sample: [ single_end, fastq_1, fastq_2 ] }
            if ",".join( sample + lane ) not in sample_mapping_dict:
                sample_mapping_dict[",".join( sample + lane )] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[",".join( sample + lane )]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[",".join( sample + lane )].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "lane", "matched_sample", "subject_id", "gender", "phenotype", "paternal_id", "maternal_id", "single_end", "fastq_1", "fastq_2"]) + "\n")
            for sample_lane in sorted(sample_mapping_dict.keys()):

                ## Check that multiple runs of the same sample (lane) are of the same datatype, i.e. all are PE or all are SE, here is the sample_info[5]
                if not all(x[5] == sample_mapping_dict[sample_lane][0][5] for x in sample_mapping_dict[sample_lane]):
                    print_error("Multiple runs of a sample must be of the same datatype!", "Sample_lane: {}".format(sample_lane))

                for idx, val in enumerate(sample_mapping_dict[sample_lane]):
                    fout.write(",".join(sample_lane + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
