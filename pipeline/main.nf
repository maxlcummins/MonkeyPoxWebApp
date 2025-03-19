#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define parameters
params.input_dir = "fastq"  // Default input directory
params.output_dir = "results"  // Default output directory
params.pattern = "*_{R1,R2}.fastq.gz"  // Default pattern for paired-end FASTQ files

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
    // Create a channel from paired-end FASTQ files
    fastq_pairs = Channel.fromFilePairs("${params.input_dir}/${params.pattern}")
        .ifEmpty { error "Cannot find any files matching the pattern: ${params.pattern}" }
    
    // Generate individual CSV files
    csv_files = generateCSV(fastq_pairs)
    
    // Combine all CSV files into one
    combineCSV(csv_files.collect())
}