nextflow.enable.dsl = 2

process PREPARE_DATA {
    publishDir "${params.results_dir}/phylo", mode: 'copy', pattern: "metadata.tsv"
    
    input:
    path fhir_files
    path references

    output:
    path "unaligned_sequences.fasta", emit: unaligned
    path "metadata.tsv",              emit: metadata

    script:
    """
    python3 $baseDir/scripts/fhir_phylo.py \\
        --inputs ${fhir_files} \\
        --references ${references}
    """
}

process RUN_MAFFT {
    input:
    path unaligned

    output:
    path "aligned_sequences.fasta", emit: alignment

    script:
    """
    mafft --auto --thread ${task.cpus} ${unaligned} > aligned_sequences.fasta
    """
}

process RUN_IQTREE {
    publishDir "${params.results_dir}/phylo", mode: 'copy'
    
    input:
    path alignment

    output:
    path "phylo_tree.nwk", emit: tree
    path "iqtree.log",     emit: log

    script:
    """
    iqtree2 -s ${alignment} -m GTR+I+G -nt AUTO -bb 1000 -pre denv_analysis
    
    mv denv_analysis.treefile phylo_tree.nwk
    mv denv_analysis.log iqtree.log
    """
}

process CALC_MATRIX {
    publishDir "${params.results_dir}/phylo", mode: 'copy'

    input:
    path alignment

    output:
    path "distance_matrix.tsv", emit: matrix

    script:
    """
    python3 $baseDir/scripts/calc_matrix.py --input ${alignment} --output distance_matrix.tsv
    """
}

workflow PHYLO_ANALYSIS {
    take:
    fhir_files
    references

    main:
    PREPARE_DATA(fhir_files.collect(), references)
    RUN_MAFFT(PREPARE_DATA.out.unaligned)
    RUN_IQTREE(RUN_MAFFT.out.alignment)
    CALC_MATRIX(RUN_MAFFT.out.alignment)

    emit:
    matrix   = CALC_MATRIX.out.matrix
    tree     = RUN_IQTREE.out.tree
    metadata = PREPARE_DATA.out.metadata
}