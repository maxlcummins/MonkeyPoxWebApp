#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define parameters
params.input_dir = "fastq"  // Default input directory
params.output_dir = "results"  // Default output directory
params.pattern = "*_{R1,R2}.fastq.gz"  // Default pattern for paired-end FASTQ files
params.s3_paths = []  // Default empty list for S3 paths

// Print workflow header
log.info """
==============================================
FASTQ Processing Workflow
==============================================
Input directory: ${params.input_dir}
Output directory: ${params.output_dir}
Pattern: ${params.pattern}
==============================================
"""

// Function to generate a random letter based on sample ID
def getRandomLetter(sampleId) {
    def random = new Random(sampleId.hashCode())
    def letters = ['A', 'B', 'C']
    return letters[random.nextInt(letters.size())]
}

// Process to download files from S3
process downloadFiles {
    output:
    path "${params.input_dir}/*_{R1,R2}.fastq.gz"

    script:
    """
    mkdir -p ${params.input_dir}
    for s3path in "${params.s3_paths.join(' ')}"; do
        filename=\$(basename \$s3path)
        aws s3 cp \$s3path ${params.input_dir}/\$filename
    done
    """
}

// Define the process to create CSV files
process generateCSV {
    publishDir params.output_dir, mode: 'copy'

    input:
    tuple val(sampleId), path(fastq_files)

    output:
    path "${sampleId}.csv"

    script:
    def letter = getRandomLetter(sampleId)
    """
    echo "${sampleId},${letter}" > ${sampleId}.csv
    """
}

// Define the process to combine all CSV files into a single file
process combineCSV {
    publishDir params.output_dir, mode: 'copy'

    input:
    path csv_files

    output:
    path "combined_results.csv"

    script:
    """
    cat $csv_files > combined_results.csv
    """
}
// Define the workflow
workflow {
    // Download files from S3 and store the output
    downloaded_files = downloadFiles()

    // Create a channel from paired-end FASTQ files using the downloaded files
    fastq_pairs = Channel.fromFilePairs(downloaded_files.collect())
        .ifEmpty { error "Cannot find any files matching the pattern: ${params.pattern}" }
        .ifEmpty { error "Cannot find any files matching the pattern: ${params.pattern}" }

    // Generate individual CSV files
    csv_files = generateCSV(fastq_pairs)

    // Combine all CSV files into one
    combineCSV(csv_files.collect())
}