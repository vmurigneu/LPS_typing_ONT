#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
        Pasteurella multocida LPS analysis pipeline
========================================================================================
 #### Documentation
 #https://github.com/julianzaugg/LPS_typing
 #### Authors
 Valentine Murigneux <v.murigneux@uq.edu.au>
 Julian Zaugg <j.zaugg@uq.edu.au>
========================================================================================
*/

def helpMessage() {
        log.info"""
        =========================================
        Pasteurella multocida LPS analysis pipeline v${workflow.manifest.version}
        =========================================
        Usage:
        Typing workflow
        nextflow main.nf --samplesheet /path/to/samples.csv --outdir /path/to/outdir/ --slurm_account account

        Required arguments:
                --samplesheet                           Path to the samplesheet file
                --outdir                                Path to the output directory to be created
                --slurm_account                         Name of the Bunya account (default='a_uqds')

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

process nanocomp {
        cpus "${params.threads}"
        label "cpu"
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*log'
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*txt'
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*html'
        input:
                val(sampleID_list)
                path(fastq_files)
        output:
                path("2_ONT_nanocomp_report.html"), emit: nanocomp_html
                path("2_ONT_nanostats.txt"), emit: nanocomp_stats
                path("2_ONT_nanocomp.log")
        when:
        !params.skip_nanocomp
        script:
        """
        echo ${fastq_files} > fastq_files.txt
        echo ${sampleID_list} > sampleID_list.txt
        sed  "s/\\[//" sampleID_list.txt | sed "s/\\]//" | sed "s/\\,//g" > sample_list
        sampleID_list_names=\$(cat "sample_list")
        NanoComp -o \$PWD --fastq ${fastq_files} -t ${params.threads} -n \${sampleID_list_names}
        mv NanoStats.txt 2_ONT_nanostats.txt
        mv NanoComp-report.html 2_ONT_nanocomp_report.html
        cp .command.log 2_ONT_nanocomp.log
        """
}

process flye {
        cpus "${params.flye_threads}"
        tag "${sample}"
        label "high_memory" 
        label "cpu"
        publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: "*version.txt", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: "${sample}_assembly.fasta"
        publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: "*info.txt"
        input:
                tuple val(sample), path(fastq)
        output:
                tuple val(sample), path(fastq), path("${sample}_assembly.fasta"), emit: assembly_fasta
                tuple val(sample), path("${sample}_assembly.fasta"), emit: assembly_only
                tuple val(sample), path("*assembly_info.txt"), path("assembly_graph.gfa"),path("assembly_graph.gv"), emit: assembly_graph
                path("*assembly_info.txt"), emit: assembly_info 
                path("flye.log")
                path("flye_version.txt")
        when:
        !params.skip_assembly
        script:
        """
        flye --nano-hq ${fastq} --threads ${params.flye_threads} --out-dir \$PWD ${params.flye_args} --genome-size ${params.genome_size}
        if [ ! -f "assembly.fasta" ]; then
                touch assembly.fasta assembly_info.txt assembly_graph.gfa assembly_graph.gv
        fi
        mv assembly_info.txt ${sample}_assembly_info.txt
        mv assembly.fasta ${sample}_assembly.fasta
        flye -v 2> flye_version.txt
        cp .command.log flye.log
        """  
}

process medaka {
        cpus "${params.medaka_threads}"
        tag "${sample}"
        label "medaka"
        label "cpu"
        publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: "${sample}_flye_polished.fasta"
        publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/3_assembly",  mode: 'copy', pattern: "*_version.txt" 
        input:
                tuple val(sample), path(fastq), path(draft)
        output:
                tuple val(sample), path ("${sample}_flye_polished.fasta"), emit: polished_medaka
        path("medaka.log")
        path("medaka_version.txt")
        when:
        !params.skip_polishing  
        script:
        """
        set +eu
        medaka_consensus -i ${fastq} -d ${draft} -o \$PWD -t ${params.medaka_threads} -m ${params.medaka_model}
        rm consensus_probs.hdf calls_to_draft.bam calls_to_draft.bam.bai
        if [ -f "consensus.fasta" ]; then
                mv consensus.fasta ${sample}_flye_polished.fasta
        else
                touch ${sample}_flye_polished.fasta
        fi
        cp .command.log medaka.log
        medaka --version > medaka_version.txt
        """
}

process quast {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/4_quast",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/4_quast",  mode: 'copy', pattern: '*tsv'
        input:
                tuple val(sample), path(assembly)
        output:
                path("*report.tsv"), emit: quast_results
                path("quast.log")
        when:
        !params.skip_quast
        script:
        """
        quast.py ${assembly} --threads ${params.threads} -o \$PWD
        sed "s/_flye_polished//;s/_assembly//" report.tsv > ${sample}_report.tsv
        rm transposed_report.tsv report.tsv
        cp .command.log quast.log
        """
}

process summary_quast {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
        input:
                path(quast_files)
        output:
                path("4_ONT_quast_report.tsv"), emit: quast_summary
        when:
        !params.skip_quast
        script:
        """
        for file in `ls *report.tsv`; do cut -f2 \$file > \$file.tmp.txt; cut -f1 \$file > rownames.txt; done
        paste rownames.txt *tmp.txt > 4_ONT_quast_report.tsv
        """
}

process summary_flye {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
        input:
                path(flye_info_files)
        output:
                path("3_ONT_flye_stats.tsv"), emit: flye_summary
        when:
        !params.skip_assembly
        script:
        """
        echo -e "sample\tassembly_coverage\tnb_contigs\tassembly_size" > 3_ONT_flye_stats.tsv
        for file in `ls *info.txt`; do
                fileName=\$(basename \$file)
                sample=\${fileName%%_assembly_info.txt}
                grep -v length \$file > tmp
                total_length=`awk '{total_length+=\$2} END {print total_length}' tmp`
                total_cov=`awk '{total_cov+=\$2*\$3} END {print total_cov}' tmp`
                mean_cov=\$((\$total_cov/\$total_length))
                nb_contigs=`grep contig \$file | wc -l`
                echo -e \$sample\\\t\$mean_cov\\\t\$nb_contigs\\\t\$total_length  >> 3_ONT_flye_stats.tsv
        done
        """
}

process download_checkm_db {
        publishDir "$params.outdir/databases/checkm_database",  mode: 'copy'
        output:
                path("checkm_data_2015_01_16"), emit: checkm_db_folder
        when:
        !params.skip_download_checkm_db
        script:
        """
        wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
        mkdir checkm_data_2015_01_16
        tar -xvzf checkm_data_2015_01_16.tar.gz -C checkm_data_2015_01_16
        """
}

process checkm {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        label "high_memory"
        publishDir "$params.outdir/$sample/5_checkm",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/5_checkm",  mode: 'copy', pattern: '*tsv'
        input:
                tuple val(sample), path(assembly), path(checkm_db)
        output:
                path("*checkm_lineage_wf_results.tsv"),  emit: checkm_results
                path("checkm.log")
        when:
        !params.skip_checkm
        script:
        """
        export CHECKM_DATA_PATH=${checkm_db}
        checkm data setRoot ${checkm_db}
        checkm lineage_wf --reduced_tree `dirname ${assembly}` \$PWD --threads ${params.threads} --pplacer_threads ${params.threads} --tab_table -f checkm_lineage_wf_results.tsv -x fasta
        mv checkm_lineage_wf_results.tsv ${sample}_checkm_lineage_wf_results.tsv
        cp .command.log checkm.log
        """
}

process summary_checkm {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
        input:
                path(checkm_files)
        output:
                path("5_ONT_checkm_lineage_wf_results.tsv"), emit: checkm_summary
        when:
        !params.skip_checkm
        script:
        """
        echo -e  sampleID\\\tMarker_lineage\\\tNbGenomes\\\tNbMarkers\\\tNbMarkerSets\\\t0\\\t1\\\t2\\\t3\\\t4\\\t5+\\\tCompleteness\\\tContamination\\\tStrain_heterogeneity > header_checkm
        for file in `ls *checkm_lineage_wf_results.tsv`; do fileName=\$(basename \$file); sample=\${fileName%%_checkm_lineage_wf_results.tsv}; grep -v Bin \$file | \
        sed "s/^flye_polished/\${sample}/;s/_flye_polished//;s/^assembly/\${sample}/;s/_assembly//" >> 5_checkm_lineage_wf_results.tsv.tmp; done
        cat header_checkm 5_checkm_lineage_wf_results.tsv.tmp > 5_ONT_checkm_lineage_wf_results.tsv
        """
}

process sylph_download_db {
        cpus 1
        label "cpu"
        publishDir "$params.outdir/databases/sylph", mode: 'copy', pattern: "*.syldb"
        input:
                val(db)
        output:
                path("*.syldb"), emit: sylph_db
        when:
        !params.skip_download_sylph_db
        script:
        """
        echo "${db}"
        for attempt in 1 2 3 4 5; do
            wget -c -T 60 "${db}" && break
            [ \$attempt -lt 5 ] && sleep 30
        done
        """
}

process sylph {
        cpus "${params.sylph_threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/6_sylph",  mode: 'copy', pattern: "*.tsv", saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(fastq), path(db_files)
        output:
                tuple val(sample), path("*sylph_profile.tsv"), emit: sylph_profile
        when:
        !params.skip_sylph
        script:
        """
        sylph \
        profile ${db_files.join(' ')} \
        -t ${params.sylph_threads} ${fastq} \
        --output-file sylph_profile.tsv
        """
}

process sylph_tax_download_metadata {
        cpus 1
        label "cpu"
        publishDir "$params.outdir/databases/sylph", mode: 'copy', pattern: "*.gz"
        input:
                val(metadata_file)
        output:
                path("*.gz"), emit: sylph_tax_metadata
        when:
        !params.skip_download_sylph_db
        script:
        """
        for attempt in 1 2 3 4 5; do
            wget -c -T 60 "${metadata_file}" -O \$PWD/\$(basename $metadata_file) && break
            [ \$attempt -lt 5 ] && sleep 30
        done
        #wget "https://zenodo.org/records/15314244/files/gtdb_r226_metadata.tsv.gz" -O \$PWD/gtdb_r226_metadata.tsv.gz
        #wget "https://zenodo.org/records/14320496/files/fungi_refseq_2024-07-25_metadata.tsv.gz" -O \$PWD/fungi_refseq_2024-07-25_metadata.tsv.gz
        #sylph-tax download --download-to \$PWD
        """
}

process sylph_tax {
        cpus "${params.sylph_threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/6_sylph",  mode: 'copy', pattern: "*.tsv", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/6_sylph",  mode: 'copy', pattern: "*.sylphmpa", saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(sylph_profile), path(metadata_files)
        output:
                tuple val(sample), path("merged_taxonomic_abundance.tsv"), path("merged_sequence_abundance.tsv"), emit: sylph_tax
        when:
        !params.skip_sylph
        script:
        """
        sylph-tax taxprof \
        ${sylph_profile} \
        -o \$PWD/ \
        -t ${metadata_files.join(' ')}

        # Merge taxonomy outputs
        sylph-tax merge \$PWD/*.sylphmpa \
        --column relative_abundance \
        -o \$PWD/merged_taxonomic_abundance.tsv

        sylph-tax merge \$PWD/*.sylphmpa \
        --column sequence_abundance \
        -o \$PWD/merged_sequence_abundance.tsv
        """
}

process sylph_summary_per_sample {
        //publishDir "$params.outdir/$sample/6_sylph",  mode: 'copy', pattern: "*.tsv"
        input:
               tuple val(sample), path(taxonomic_abundances), path(sequence_abundances)
        output:
               tuple val(sample), path("${sample}_sylph_summary.tsv")
        script:
        """
        taxonomic_abundance_top_species=\$(grep "s__" ${taxonomic_abundances} | grep -v "t__" | sort -t \$'\t' -gr -k 2 | head -n 1 | sed "s/.*s__//g")
        sequence_abundance_top_species=\$(grep "s__" ${sequence_abundances} | grep -v "t__" | sort -t \$'\t' -gr -k 2 | head -n 1 | sed "s/.*s__//g")

        taxonomic_abundance_pasteurella_multocida=\$(grep "s__Pasteurella multocida" ${taxonomic_abundances} | grep -v "t__" | sort -t \$'\t' -gr -k 2 | head -n 1 | sed "s/.*s__//g" | awk -F "\t" '{print \$2}')
        sequence_abundance_pasteurella_multocida=\$(grep "s__Pasteurella multocida" ${sequence_abundances} | grep -v "t__" | sort -t \$'\t' -gr -k 2 | head -n 1 | sed "s/.*s__//g" | awk -F "\t" '{print \$2}')
        echo -e "${sample}\t\$taxonomic_abundance_top_species\t\$sequence_abundance_top_species\tPasteurella_multocida\t\$taxonomic_abundance_pasteurella_multocida\t\$sequence_abundance_pasteurella_multocida" > ${sample}_sylph_summary.tsv
        """
}

process summary_sylph {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
        input:
                path(sylph_summary_files)
        output:
                path("6_ONT_sylph_summary.tsv"), emit: sylph_summary
        when:
        !params.skip_sylph
        script:
        """
        echo -e "sample\ttop_species_by_taxonomic_abundance\ttaxonomic_abundance_for_top_species\ttop_species_by_sequence_abundance\tsequence_abundance_for_top_species\tPasteurella_multocida\ttaxonomic_abundance_for_pasteurella_multocida\tsequence_abundance_for_pasteurella_multocida" > 6_ONT_sylph_summary.tsv
        for file in ${sylph_summary_files.join(' ')}; do
            cat \$file >> 6_ONT_sylph_summary.tsv
        done
        """
}

process kaptive3 {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: '*fna'
        publishDir "$params.outdir/$sample/7_kaptive_v3",  mode: 'copy', pattern: '*tsv'
        input:
                tuple val(sample), path(assembly)
        output:
                tuple val(sample), path("*kaptive_results.tsv"), emit: kaptive_results
                path("*kaptive_results.tsv"),  emit: kaptive_tsv
                path("*fna")
                path("kaptive_v3.log")
        when:
        !params.skip_kaptive3
        script:
        """
        kaptive assembly ${params.kaptive_db_9lps} ${assembly} -f \$PWD/${sample}_kaptive_results.fna -o ${sample}_kaptive_results.tsv
        if [[ -f ${sample}_kaptive_results.fna ]]; then
                sed -i "s/_flye_polished//;s/_assembly//" \$PWD/${sample}_kaptive_results.fna
                sed -i "s/_flye_polished//;s/_assembly//" \$PWD/${sample}_kaptive_results.tsv
        else
                touch ${sample}_kaptive_results.fna
        fi
        cp .command.log kaptive_v3.log
        """
}

process summary_kaptive {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
        input:
                path(kaptive_files)
        output:
                path("7_ONT_kaptive_results.tsv"), emit: kaptive_summary
        when:
        !params.skip_kaptive3
        script:
        """
        echo -e sampleID\\\tBest match locus\\\tBest match type\\\tMatch confidence\\\tProblems\\\tIdentity\\\tCoverage\\\tLength discrepancy\\\tExpected genes in locus\\\tExpected genes in locus, details\\\tMissing expected genes\\\tOther genes in locus\\\tOther genes in locus, details\\\tExpected genes outside locus\\\tExpected genes outside locus, details\\\tOther genes outside locus\\\tOther genes outside locus, details\\\tTruncated genes, details\\\tExtra genes, details >  header_kaptive3
        for file in `ls *_kaptive_results.tsv`; do
        fileName=\$(basename \$file)
        sample=\${fileName%%_kaptive_results.tsv}
        grep -v Assembly \$file | sed "s/^flye_polished/\${sample}/;s/^assembly/\${sample}/;s/_flye_polished//;s/_assembly//"  >> 7_kaptive_results.tsv.tmp; done
        cat header_kaptive3 7_kaptive_results.tsv.tmp > 7_ONT_kaptive_results.tsv
        """
}

process minimap {
        cpus "${params.minimap_threads}"
        tag "${sample}"
        label "cpu"
        label "high_memory"
        publishDir "$params.outdir/$sample/8_clair3",  mode: 'copy', pattern: '*txt', saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/8_clair3",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/8_clair3",  mode: 'copy', pattern: 'minimap2_mapped.ba*', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(fastq), path(kaptive_report)
        output:
                tuple val(sample), path("minimap2_mapped.bam"), path("minimap2_mapped.bam.bai"), path(kaptive_report), emit: minimap_results
                path("minimap.log")
                path("minimap2_flagstat.txt")
        when:
        !params.skip_clair3
        script:
        """

        if [ \$(wc -l < "${kaptive_report}") -le 1 ]; then
                echo "Warning: kaptive_report has only a header for sample ${sample}" > minimap.log
                touch minimap2_mapped.bam minimap2_mapped.bam.bai minimap2_flagstat.txt
        else
                locus=\$(tail -1 "${kaptive_report}" | cut -f3)
                ref_fasta=\$(grep \${locus:0:2} "${params.reference_LPS_directory}/reference_LPS.txt" | cut -f3)
                ref_fasta="${params.reference_LPS_directory}/\$ref_fasta"

                minimap2 -t ${params.minimap_threads} -ax map-ont -k19 -w 19 -U50,500 -g10k \$ref_fasta ${fastq} > minimap2.sam
                samtools sort -o minimap2.bam -@ ${params.minimap_threads} minimap2.sam
                samtools index minimap2.bam
                samtools flagstat minimap2.bam > minimap2_flagstat.txt
                samtools view -b -F 4 minimap2.bam > minimap2_mapped.bam
                samtools index minimap2_mapped.bam
                cp .command.log minimap.log
        fi
        """
}

process clair3 {
        cpus "${params.clair3_threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/8_clair3",  mode: 'copy', pattern: '*vcf', saveAs: { filename -> "${sample}_$filename"}
        publishDir "$params.outdir/$sample/8_clair3",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(bam), path(bai), path(kaptive_report)
        output:
                tuple val(sample), path(bam), path(bai), path ("clair3.vcf"), path(kaptive_report), emit: clair3_results
                path("clair3.log")
        when:
        !params.skip_clair3
        script:
        """
        if [ \$(wc -l < "${kaptive_report}") -le 1 ]; then
                echo "Warning: kaptive_report has only a header for sample ${sample}" > clair3.log
                touch clair3.vcf
        else
                locus=\$(tail -1 "${kaptive_report}" | cut -f3)
                ref_gb=\$(grep \${locus:0:2} "${params.reference_LPS_directory}/reference_LPS.txt" | cut -f2)
                ref_fasta=\$(grep \${locus:0:2} "${params.reference_LPS_directory}/reference_LPS.txt" | cut -f3)
                ref_gb="${params.reference_LPS_directory}/\$ref_gb"
                ref_fasta="${params.reference_LPS_directory}/\$ref_fasta"
                # TODO run clair3 in full alignment (default) and --pileup_only modes, extract high impact mutations missing (or low qual) from full alignment mode and merge in.
                run_clair3.sh --bam_fn=${bam} --ref_fn=\$ref_fasta --threads=${params.clair3_threads} --platform="ont" --model_path=${params.clair3_model} --sample_name=${sample} --output=\$PWD ${params.clair3_args} --no_phasing_for_fa --include_all_ctgs --enable_long_indel
                gunzip -c merge_output.vcf.gz > merge_output.vcf
                #gunzip -c pileup.vcf.gz > merge_output.vcf
                mv merge_output.vcf clair3.vcf
                cp .command.log clair3.log
        fi
        """
}

process snpeff {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/8_clair3",  mode: 'copy', pattern: '*vcf'
        publishDir "$params.outdir/$sample/8_clair3",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(bam), path(bai), path(vcf), path(kaptive_report)
        output:
                tuple val(sample), path ("*clair3.snpeff.vcf"), emit: snpeff_results
                path("snpeff.log")
        when:
        !params.skip_clair3 || !params.skip_snpeff
        script:
        """
        if [ \$(wc -l < "${kaptive_report}") -le 1 ]; then
                echo "Warning: kaptive_report has only a header for sample ${sample}" > snpeff.log
                touch ${sample}_clair3.snpeff.vcf
        else
                locus=\$(tail -1 "${kaptive_report}" | cut -f3)
                ref_gb=\$(grep "\${locus:0:2}" "${params.reference_LPS_directory}/reference_LPS.txt" | cut -f2)
                ref_gb="${params.reference_LPS_directory}/\$ref_gb"
                # mkdir -p LPS_snpeffdb
                mkdir -p snpeff_output/LPS_snpeffdb
                mkdir -p data/LPS_snpeffdb
                cp \${ref_gb} snpeff_output/LPS_snpeffdb/genes.gbk
                # Note: The placeholder genome name used here "LPS_snpeffdb" will be different to that in the GB file. 
                # There will cause a warning "WARNING_REF_DOES_NOT_MATCH_GENOME" to appear in the snpEff output.
                # This can be ignored as the genome will be correct.
                snpEff build -v -configOption 'LPS_snpeffdb'.genome='LPS_snpeffdb' -configOption 'LPS_snpeffdb'.codonTable='Bacterial_and_Plant_Plastid' -genbank -dataDir \$PWD/snpeff_output LPS_snpeffdb
                mv snpeff_output/'LPS_snpeffdb'/*.bin data/'LPS_snpeffdb'
                config_path=\$(find /usr/local/share/ -name "snpEff.config")
                cp \$config_path snpEff.config
                echo 'LPS_snpeffdb.genome : LPS_snpeffdb' >> snpEff.config
                echo 'LPS_snpeffdb.codonTable : Bacterial_and_Plant_Plastid' >> snpEff.config
                
                # Ensure chromosome name in VCF matches that in the SnpEff DB
                current_chromosome_name=\$(grep -v '^#' "${vcf}" | head -n 1 | awk '{print \$1}')
                snpEff dump LPS_snpeffdb 2>/dev/null | grep "Bacterial_and_Plant_Plastid" | awk -F "'" '{print \$2}' > new_name.txt
                new_chromosome_name=\$(cat new_name.txt)
                sed "s/\$current_chromosome_name/\$new_chromosome_name/g" ${vcf} > clair3_name_modified.vcf
                snpEff eff -i vcf -o vcf -c snpEff.config -lof -nodownload -no-downstream -no-intron -no-upstream -no-utr -no-intergenic -v -configOption 'LPS_snpeffdb'.genome='LPS_snpeffdb' -configOption 'LPS_snpeffdb'.codonTable='Bacterial_and_Plant_Plastid' -stats snpeff.html LPS_snpeffdb clair3_name_modified.vcf > clair3.snpeff.vcf

                # Change chromosome name back
                sed -i "s/\$new_chromosome_name/\$current_chromosome_name/g" clair3.snpeff.vcf

                # Remove any non-VCF header lines (like JVM warnings)
                sed -i -n '/^##fileformat/,\$p' clair3.snpeff.vcf

                mv clair3.snpeff.vcf ${sample}_clair3.snpeff.vcf
                cp .command.log snpeff.log
        fi
        """
}

process snpsift {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/8_clair3",  mode: 'copy', pattern: '*log', saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/8_clair3",  mode: 'copy', pattern: '*vcf'
        input:
                tuple val(sample), path(vcf)
        output:
                tuple path(vcf), path ("*clair3.snpeff.high_impact.vcf"), emit: snpsift_results
                path("snpsift.log")
        when:
        !params.skip_clair3 || !params.skip_snpeff
        script:
        """
        SnpSift filter "( EFF[*].IMPACT = 'HIGH' ) && (FILTER = 'PASS')" -f ${vcf} > clair3.snpeff.high_impact.vcf
        
        # Remove any non-VCF header lines (like JVM warnings)
        sed -i -n '/^##fileformat/,\$p' clair3.snpeff.high_impact.vcf
        
        mv clair3.snpeff.high_impact.vcf ${sample}_clair3.snpeff.high_impact.vcf
        cp .command.log snpsift.log
        """
}

process petg_blast {
        cpus "${params.petg_threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/13_petG",  mode: 'copy', pattern: '*_petG_*'
        publishDir "$params.outdir/$sample/13_petG",  mode: 'copy', pattern: 'petg_blast.log', saveAs: { filename -> "${sample}_$filename" }
        input:
                tuple val(sample), path(assembly), path(petg_reference)
        output:
                path("${sample}_petG_hits.fasta"), emit: petg_hits
                path("${sample}_petG_summary.tsv"), emit: petg_summary
                path("${sample}_petG_blast.tsv")
                path("${sample}_petG_blast.filtered.tsv")
                path("petg_blast.log")
        when:
        !params.skip_petg
        script:
        """
        makeblastdb -in ${assembly} -dbtype nucl -parse_seqids -out assembly_db
        blastn \\
                -query ${petg_reference} \\
                -db assembly_db \\
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \\
                -num_threads ${params.petg_threads} \\
                -out ${sample}_petG_blast.tsv

        awk -F '\\t' -v OFS='\\t' -v min_len="${params.petg_min_length}" -v min_ident="${params.petg_min_identity}" '
        {
                start = \$9
                end = \$10
                strand = "plus"
                suffix = "+"
                if (start > end) {
                        start = \$10
                        end = \$9
                        strand = "minus"
                        suffix = "-"
                }
                span = end - start + 1
                if (\$3 >= min_ident && span > min_len) {
                        print \$0, start, end, strand, suffix
                }
        }
        ' ${sample}_petG_blast.tsv > ${sample}_petG_blast.filtered.tsv

        touch ${sample}_petG_hits.fasta
        echo -e "SAMPLE\\tPETG_PRESENT" > ${sample}_petG_summary.tsv
        if [[ -s ${sample}_petG_blast.filtered.tsv ]]; then
                while IFS=\$'\\t' read qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore hit_start hit_end strand suffix; do
                        blastdbcmd -db assembly_db -entry "\$sseqid" -range "\${hit_start}-\${hit_end}" -strand "\$strand" | \\
                        awk -v header=">${sample}|\${sseqid}:\${hit_start}-\${hit_end}\${suffix}" 'BEGIN { print header } /^>/ { next } { print }' >> ${sample}_petG_hits.fasta
                done < ${sample}_petG_blast.filtered.tsv
                echo -e "${sample}\\tyes" >> ${sample}_petG_summary.tsv
        else
                echo -e "${sample}\\t" >> ${sample}_petG_summary.tsv
        fi

        cp .command.log petg_blast.log
        """
}

process empty_mlst_report_input {
        output:
                path("empty_mlst_report_input.txt"), emit: empty_mlst
        when:
        params.skip_mlst
        script:
        """
        touch empty_mlst_report_input.txt
        """
}

process empty_petg_report_input {
        output:
                path("empty_petg_report_input.txt"), emit: empty_petg
        when:
        params.skip_petg
        script:
        """
        touch empty_petg_report_input.txt
        """
}

process report {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*vcf'
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
        input:
                path(clair3_files)
                path(kaptive_summary)
                path(petg_summaries)
                path(mlst_files)
        output:
                tuple path("8_ONT_clair3_snpeff.vcf"), path("8_ONT_clair3_snpeff_high_impact.vcf"), path("10_ONT_subtype_report.tsv"), emit: subtype_report
        when:
        !params.skip_clair3 && !params.skip_kaptive3
        script:
        """
        echo -e  SAMPLEID\\\tCHROM\\\tPOS\\\tID\\\tREF\\\tALT\\\tQUAL\\\tFILTER\\\tINFO\\\tFORMAT\\\tSAMPLE > header_clair3
        for file in `ls *_clair3.snpeff.high_impact.vcf`; do fileName=\$(basename \$file); sample=\${fileName%%_clair3.snpeff.high_impact.vcf}; grep -v "^#" \$file | sed s/^/\${sample}\\\t/  >> 8_clair3_snpeff_high_impact.vcf.tmp; done
        cat header_clair3 8_clair3_snpeff_high_impact.vcf.tmp > 8_ONT_clair3_snpeff_high_impact.vcf
        for file in `ls *_clair3.snpeff.vcf`; do fileName=\$(basename \$file); sample=\${fileName%%_clair3.snpeff.vcf}; grep -v "^#" \$file | sed s/^/\${sample}\\\t/  >> 8_clair3_snpeff.vcf.tmp; done
        cat header_clair3 8_clair3_snpeff.vcf.tmp > 8_ONT_clair3_snpeff.vcf

        echo -e "SAMPLE\\tPETG_PRESENT" > petg_lookup.tsv
        find . -maxdepth 1 -name '*_petG_summary.tsv' | sort | while read petg_file; do
                tail -n +2 "\$petg_file" >> petg_lookup.tsv
        done

        echo -e "SAMPLE\\tMLST" > mlst_lookup.tsv
        find . -maxdepth 1 -name '*_mlst.csv' | sort | while read mlst_file; do
                fileName=\$(basename "\$mlst_file")
                sample=\${fileName%%_mlst.csv}
                mlst_st=\$(awk -F',' 'NR == 1 {print \$3; exit}' "\$mlst_file")
                echo -e "\${sample}\\t\${mlst_st}" >> mlst_lookup.tsv
        done

        awk -F'\t' -v OFS='\t' '
        NR > 1 {
                if (\$4 == "Typeable") {
                        split(\$2, a, "-")
                        gsub("LPS", "L", a[1])
                        print \$1, a[1]
                } else {
                        print \$1, "untypeable"
                }
        }
        ' "${kaptive_summary}" > kaptive_type_lookup.tsv

        subtype_db="${params.reference_LPS_directory}/LPS_subtype_database_v2.txt"
        phenotype_lookup="${params.reference_LPS_directory}/phenotype_lookup.tsv"
        if [[ -f "\$phenotype_lookup" ]]; then
                phenotype_lookup_input="\$phenotype_lookup"
        else
                touch phenotype_lookup_empty.tsv
                phenotype_lookup_input="phenotype_lookup_empty.tsv"
        fi

        awk -F '\\t' -v OFS='\\t' -v subtype_db="\$subtype_db" -v phenotype_lookup="\$phenotype_lookup_input" -v petg_lookup="petg_lookup.tsv" -v mlst_lookup="mlst_lookup.tsv" -v kaptive_type_lookup="kaptive_type_lookup.tsv" '
        function set_header_fields(    i) {
                for (i = 1; i <= NF; i++) {
                        header[\$i] = i
                }
                db_type_col = ("TYPE" in header) ? header["TYPE"] : 1
                db_subtype_col = ("SUBTYPE" in header) ? header["SUBTYPE"] : 2
                db_isolate_col = ("ISOLATE" in header) ? header["ISOLATE"] : 3
                db_chrom_col = ("CHROM" in header) ? header["CHROM"] : 4
                db_pos_col = ("POS" in header) ? header["POS"] : 5
                db_vartype_col = ("VARTYPE" in header) ? header["VARTYPE"] : 6
                db_ref_col = ("REF" in header) ? header["REF"] : 7
                db_alt_col = ("ALT" in header) ? header["ALT"] : 8
                db_gene_col = ("GENE" in header) ? header["GENE"] : 9
                db_pheno_default_col = ("PHENOTYPE_DEFAULT" in header) ? header["PHENOTYPE_DEFAULT"] : 0
                db_pheno_multi_col = ("PHENOTYPE_MULTIPLE_SUBTYPES" in header) ? header["PHENOTYPE_MULTIPLE_SUBTYPES"] : 0
                db_note_col = ("NOTE" in header) ? header["NOTE"] : 10
        }
        function field_value(col) {
                return (col > 0 && col <= NF) ? \$col : ""
        }
        function clean_phenotype(value) {
                return (value == "NA") ? "" : value
        }
        function phenotype_from_rule(type, rule, parts) {
                if (rule == "" || rule == "NA") {
                        return ""
                }
                split(rule, parts, "_")
                return (parts[2] == "") ? rule : type "_" parts[2]
        }
        function choose_phenotype(sample, type, default_phenotype, multi_phenotypes, rules, i, parts) {
                default_phenotype = clean_phenotype(default_phenotype)
                if (multi_phenotypes != "" && multi_phenotypes != "NA") {
                        rule_count = split(multi_phenotypes, rules, ";")
                        for (i = 1; i <= rule_count; i++) {
                                split(rules[i], parts, "_")
                                if (parts[1] != "" && ((sample SUBSEP parts[1]) in sample_subtype)) {
                                        return phenotype_from_rule(type, rules[i])
                                }
                        }
                }
                return default_phenotype
        }
        FILENAME == subtype_db {
                if (FNR == 1) {
                        set_header_fields()
                        next
                }
                key = field_value(db_chrom_col) OFS field_value(db_pos_col) OFS field_value(db_ref_col) OFS field_value(db_alt_col)
                db_count[key]++
                idx = key SUBSEP db_count[key]
                db_type[idx] = field_value(db_type_col)
                db_subtype[idx] = field_value(db_subtype_col)
                db_isolate[idx] = field_value(db_isolate_col)
                db_chrom[idx] = field_value(db_chrom_col)
                db_pos[idx] = field_value(db_pos_col)
                db_vartype[idx] = field_value(db_vartype_col)
                db_ref[idx] = field_value(db_ref_col)
                db_alt[idx] = field_value(db_alt_col)
                db_gene[idx] = field_value(db_gene_col)
                db_pheno_default[idx] = field_value(db_pheno_default_col)
                db_pheno_multi[idx] = field_value(db_pheno_multi_col)
                db_note[idx] = field_value(db_note_col)
                next
        }
        FILENAME == phenotype_lookup {
                if (FNR > 1 && \$1 != "") {
                        phenotype_description[\$1] = \$2
                }
                next
        }
        FILENAME == petg_lookup {
                if (FNR > 1 && \$1 != "") {
                        petg_present[\$1] = \$2
                }
                next
        }
        FILENAME == mlst_lookup {
                if (FNR > 1 && \$1 != "") {
                        mlst_st[\$1] = \$2
                }
                next
        }
        FILENAME == kaptive_type_lookup {
                if (\$1 != "") {
                        kaptive_type[\$1] = \$2
                }
                next
        }
        FNR > 1 {
                sample = \$1
                key = \$2 OFS \$3 OFS \$5 OFS \$6
                for (i = 1; i <= db_count[key]; i++) {
                        idx = key SUBSEP i
                        row_count++
                        row_sample[row_count] = sample
                        row_idx[row_count] = idx
                        sample_subtype[sample SUBSEP db_subtype[idx]] = 1
                }
        }
        END {
                print "SAMPLE", "MLST", "TYPE", "SUBTYPE", "VARTYPE", "ISOLATE_DATABASE", "CHROM", "POS", "REF", "ALT", "GENE", "PREDICTED_PHENOTYPE", "PREDICTED_PHENOTYPE_DESCRIPTION", "PETG_PRESENT", "NOTE"
                for (i = 1; i <= row_count; i++) {
                        sample = row_sample[i]
                        idx = row_idx[i]
                        phenotype = choose_phenotype(sample, db_type[idx], db_pheno_default[idx], db_pheno_multi[idx])
                        description = (phenotype in phenotype_description) ? phenotype_description[phenotype] : ""
                        report_type = (sample in kaptive_type) ? kaptive_type[sample] : db_type[idx]
                        print sample, mlst_st[sample], report_type, db_subtype[idx], db_vartype[idx], db_isolate[idx], db_chrom[idx], db_pos[idx], db_ref[idx], db_alt[idx], db_gene[idx], phenotype, description, petg_present[sample], db_note[idx]
                }
        }
        ' "\$subtype_db" "\$phenotype_lookup_input" petg_lookup.tsv mlst_lookup.tsv kaptive_type_lookup.tsv 8_ONT_clair3_snpeff.vcf > 10_ONT_subtype_report.tsv.tmp
        cp kaptive_type_lookup.tsv kaptive_tmp
        awk -F'\t' 'NR > 1 {print \$1}' 10_ONT_subtype_report.tsv.tmp | sort | uniq > list_samples_clair_exclude
        awk -F'\t' 'FILENAME == ARGV[1] {exclude[\$1] = 1; next} !(\$1 in exclude)' list_samples_clair_exclude kaptive_tmp > kaptive_to_keep
        awk -F'\t' -v OFS='\t' '
        FILENAME == "mlst_lookup.tsv" {
                if (FNR > 1 && \$1 != "") {
                        mlst_st[\$1] = \$2
                }
                next
        }
        FILENAME == "petg_lookup.tsv" {
                if (FNR > 1 && \$1 != "") {
                        petg_present[\$1] = \$2
                }
                next
        }
        {
                print \$1, mlst_st[\$1], \$2, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "", "", petg_present[\$1], ""
        }
        ' mlst_lookup.tsv petg_lookup.tsv kaptive_to_keep > kaptive_to_keep.tsv
        cat 10_ONT_subtype_report.tsv.tmp kaptive_to_keep.tsv > 10_ONT_subtype_report.tsv
        """
}

process html_report {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*.html'
        input:
                path(report_inputs)
                val(pipeline_version)
                val(skipped_steps)
                val(param_note)
        output:
                path("LPS_typing_report.html"), emit: html_report
        when:
        !params.skip_html_report && !params.skip_clair3 && !params.skip_kaptive3
        script:
        """
        generate_lps_report.py \\
                --report-dir . \\
                --lps-db-dir ${params.reference_LPS_directory} \\
                --out LPS_typing_report.html \\
                --pipeline-version "${pipeline_version}" \\
                --skipped "${skipped_steps}" \\
                --params "${param_note}"
        """
}

process prepare_mlst_db {
        container 'docker://quay.io/biocontainers/mlst:2.33.1--hdfd78af_0'
        publishDir "${projectDir}/databases/", mode: 'copy'
        input:
                path(mlst_input_dir)
        output:
                path("mlst_db_prepared"), emit: mlst_db_prepared
        script:
        """
        mkdir -p mlst_db_prepared/blast
        # Copy scheme files, renaming PubMLST .fas allele downloads to .tfa
        for scheme_dir in ${mlst_input_dir}/*/; do
            scheme=\$(basename "\$scheme_dir")
            mkdir -p mlst_db_prepared/\$scheme
            for src in "\$scheme_dir"*; do
                fname=\$(basename "\$src")
                if [[ "\$fname" == *.fas ]]; then
                    cp "\$src" "mlst_db_prepared/\$scheme/\${fname%.fas}.tfa"
                else
                    cp "\$src" "mlst_db_prepared/\$scheme/\$fname"
                fi
            done
        done
        # Build BLAST database with scheme-prefixed sequence headers (as mlst expects)
        for scheme_dir in mlst_db_prepared/*/; do
            scheme=\$(basename "\$scheme_dir")
            for tfa in "\$scheme_dir"*.tfa; do
                [ -f "\$tfa" ] || continue
                sed "s/^>/>\$scheme./" "\$tfa"
            done
        done > mlst_db_prepared/blast/mlst.fa
        makeblastdb -hash_index -in mlst_db_prepared/blast/mlst.fa -dbtype nucl -title "PubMLST" -parse_seqids
        """
}

process mlst {
        cpus "${params.threads}"
        tag "${sample}"
        label "cpu"
        publishDir "$params.outdir/$sample/9_mlst",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/9_mlst",  mode: 'copy', pattern: '*_mlst.csv'
        input:
                tuple val(sample), path(assembly)
                path mlst_db
        output:
                path("*_mlst.csv"),  emit: mlst_results
                path("mlst.log")
        when:
        !params.skip_mlst
        script:
        def custom_db_args = mlst_db ? "--datadir ${mlst_db} --blastdb ${mlst_db}/blast/mlst.fa" : ""
        """
        mlst --scheme ${params.mlst_scheme} ${assembly} --quiet --csv --threads ${params.threads} ${custom_db_args} > mlst.csv
        sed  s/_flye_polished.fasta// mlst.csv > ${sample}_mlst.csv
        cp .command.log mlst.log
        """
}

process summary_mlst {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*csv'
        input:
                path(mlst_files)
        output:
                path("9_ONT_mlst.csv"), emit: mlst_summary
        when:
        !params.skip_mlst
        script:
        """
        for file in `ls *_mlst.csv`; do fileName=\$(basename \$file); sample=\${fileName%%_mlst.csv}; cat \$file >> 9_ONT_mlst.csv; done
        """
}

process download_bakta_db {
    publishDir "$params.outdir/databases/bakta_database",  mode: 'copy'
    output:
        path("bakta_db"), emit: bakta_db
    when:
    !params.skip_download_bakta_db
    script:
    """
    bakta_db download --output bakta_db --type full
    """
}

process bakta {
        cpus "${params.bakta_threads}"
        tag "${sample}"
        publishDir "$params.outdir/$sample/11_bakta",  mode: 'copy', pattern: "*.log"
        publishDir "$params.outdir/$sample/11_bakta",  mode: 'copy', pattern: '*bakta*'
        input:
                tuple val(sample), path(assembly), path(bakta_db)
        output:
                path("*bakta*")
                path("bakta.log")
        when:
        !params.skip_bakta
        script:
        """
        bakta --db ${bakta_db} --threads ${params.bakta_threads} --prefix ${sample}_bakta --output \$PWD/ ${params.bakta_args} ${assembly}
        cp .command.log bakta.log
        """
}

process download_amrfinder_db {
        publishDir "$params.outdir/databases/amrfinderplus_database",  mode: 'copy'
        output:
                path("amrfinderplus_db"), emit: amrfinder_db
        when:
        !params.skip_download_amrfinder_db
        script:
        """
        amrfinder_update -d amrfinderplus_db
        """
}

process amrfinder {
        tag "${sample}"
        publishDir "$params.outdir/$sample/12_amrfinder",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }
        publishDir "$params.outdir/$sample/12_amrfinder",  mode: 'copy', pattern: '*tsv'
        input:
                tuple val(sample), path(assembly), path(amrfinder_db)
        output:
                path("*.tsv"), emit: amrfinder_results
                path("amrfinder.log")
        when:
        !params.skip_amrfinder
        script:
        """
        amrfinder -n ${assembly} -d ${amrfinder_db} -o \$PWD/${sample}_amrfinder.tsv --name ${sample} --threads ${params.threads} --plus ${params.amrfinder_args}
        cp .command.log amrfinder.log
        """
}       

process summary_amrfinder {
        publishDir "$params.outdir/10_report",  mode: 'copy', pattern: '*tsv'
        input:
                path(amrfinder_files)
        output:
                path("12_ONT_amrfinder.tsv"), emit: amrfinder_summary
        when:
        !params.skip_amrfinder
        script:
        """
        echo -e Name\\\tProtein id\\\tContig id\\\tStart\\\tStop\\\tStrand\\\tElement symbol\\\tElement name\\\tScope\\\tType\\\tSubtype\\\tClass\\\tSubclass\\\tMethod\\\tTarget length\\\tReference sequence length\\\t% Coverage of reference\\\t% Identity to reference\\\tAlignment length\\\tClosest reference accession\\\tClosest reference name\\\tHMM accession\\\tHMM description > header_amrfinder
        for file in ${amrfinder_files}; do 
                tail -n +2 "\$file" >> 12_amrfinder.tsv.tmp
        done
        cat header_amrfinder 12_amrfinder.tsv.tmp > 12_ONT_amrfinder.tsv
        """
}

workflow {
        Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
        .splitCsv(header:true, sep:',')
        .map { row -> tuple(row.sample_id, file(row.long_fastq, checkIfExists: true)) }
        .set { ch_samplesheet_ONT }
        Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
        .splitCsv(header:true, sep:',')
        .map { row -> file(row.long_fastq, checkIfExists: true) }
        .set { ch_samplesheet_fastq }
        if (!params.skip_nanocomp) {
                Channel.fromPath( "${params.samplesheet}", checkIfExists:true )
                .splitCsv(header:true, sep:',')
                .map { row -> row.sample_id }
                .collect()
                .set { ch_samplesheet_sampleID }
                ch_samplesheet_sampleID.toList().view()
                nanocomp(ch_samplesheet_sampleID,ch_samplesheet_fastq.collect())
        }
        if (!params.skip_assembly) {
                flye(ch_samplesheet_ONT)
                if (!params.skip_polishing) {
                        medaka(flye.out.assembly_fasta)
                        ch_assembly_for_typing = medaka.out.polished_medaka
                } else if (params.skip_polishing) {
                        ch_assembly_for_typing = flye.out.assembly_only
                }
                summary_flye(flye.out.assembly_info.collect())
        }
        if (!params.skip_quast) {
                if (!params.skip_polishing) {
                        quast(medaka.out.polished_medaka)
                } else if (params.skip_polishing) {
                        quast(flye.out.assembly_only)
                }
                summary_quast(quast.out.quast_results.collect())
        }
        if (!params.skip_checkm) {
                if (!params.skip_download_checkm_db) {
                        download_checkm_db()
                        if (!params.skip_polishing) {
                                checkm(medaka.out.polished_medaka.combine(download_checkm_db.out.checkm_db_folder))
                        }  else if (params.skip_polishing) {
                                checkm(flye.out.assembly_only.combine(download_checkm_db.out.checkm_db_folder))
                        }
                } else{
                        if (!params.skip_polishing) {
                                checkm(medaka.out.polished_medaka.combine(Channel.fromPath( "${params.checkm_db}")))
                        }  else if (params.skip_polishing) {
                                checkm(flye.out.assembly_only.combine(Channel.fromPath( "${params.checkm_db}")))
                        } 
                }
                
                summary_checkm(checkm.out.checkm_results.collect())
        }
        if (!params.skip_sylph) {
                if (!params.skip_download_sylph_db) {
                         // Download the Sylph reference databases
                        ch_sylph_db=Channel.of("${params.sylph_db_gtdb_file}", "${params.sylph_db_fungal_file}")
                        ch_downloaded_dbs = sylph_download_db(ch_sylph_db).sylph_db
                        ch_db_list = ch_downloaded_dbs.collect()
                        ch_db_list.view { "Downloaded DBs: $it" }

                        // Download the Sylph-taxa metadata
                        ch_sylph_metadata=Channel.of("${params.sylph_tax_gtdb_metadata}", "${params.sylph_tax_fungal_metadata}")
                        sylph_tax_download_metadata(ch_sylph_metadata).collect().set{sylph_tax_metadata}

                        // Run sylph
                        ch_db_list.map { dbs -> tuple([dbs]) }.set { ch_db_tuple }

                        ch_samplesheet_ONT
                        .combine(ch_db_tuple)
                        .map { sample_id, fastq, dbs -> tuple(sample_id, fastq, dbs) }
                        .set { ch_sylph_input }

                        sylph(ch_sylph_input)

                        // Run sylph-tax
                        sylph_tax_metadata.map { dbs -> tuple([dbs]) }.set { syph_tax_metadata_tuple }
                        sylph.out.sylph_profile.combine(syph_tax_metadata_tuple).set{ch_sylph_tax_input}
                        sylph_tax(ch_sylph_tax_input)

                } else if (params.skip_download_sylph_db) {
                        ch_sylph_db=Channel.fromPath( "${params.sylph_db}").collect()
                        ch_sylph_metadata=Channel.fromPath("${params.sylph_metadata}").collect()

                        // Run sylph
                        ch_sylph_db.map { dbs -> tuple([dbs]) }.set { ch_db_tuple }

                        ch_samplesheet_ONT
                        .combine(ch_db_tuple)
                        .map { sample_id, fastq, dbs -> tuple(sample_id, fastq, dbs) }
                        .set { ch_sylph_input }

                        sylph(ch_sylph_input)

                        // Run sylph-tax
                        ch_sylph_metadata.map { dbs -> tuple([dbs]) }.set { syph_tax_metadata_tuple }
                        sylph.out.sylph_profile.combine(syph_tax_metadata_tuple).set{ch_sylph_tax_input}
                        sylph_tax(ch_sylph_tax_input)
                }
                
                sylph_tax.out.sylph_tax
                .map { sample, tax, seq -> tuple(sample, tax, seq) }
                .set{sylph_tax_results}

                sylph_summary_per_sample(sylph_tax_results).map{sample, summary_file -> summary_file }.collect().set{all_sylph_summaries}
                summary_sylph(all_sylph_summaries)
        }
        if (params.mlst_datadir) {
                prepare_mlst_db(Channel.value(file(params.mlst_datadir)))
                ch_mlst_db = prepare_mlst_db.out.mlst_db_prepared
        } else {
                ch_mlst_db = Channel.value([])
        }
        if (!params.skip_mlst) {
                mlst(ch_assembly_for_typing, ch_mlst_db)
                summary_mlst(mlst.out.mlst_results.collect())
                mlst_report_ch = mlst.out.mlst_results.collect()
        } else {
                empty_mlst_report_input()
                mlst_report_ch = empty_mlst_report_input.out.empty_mlst
        }
        if (!params.skip_kaptive3) {
                if (!params.skip_petg) {
                        ch_petg_reference = Channel.fromPath("${params.reference_LPS_directory}/petG_X73_NZ_CM001580.fasta", checkIfExists: true)
                        petg_blast(ch_assembly_for_typing.combine(ch_petg_reference))
                        petg_report_ch = petg_blast.out.petg_summary.collect()
                } else {
                        empty_petg_report_input()
                        petg_report_ch = empty_petg_report_input.out.empty_petg
                }
                kaptive3(ch_assembly_for_typing)
                summary_kaptive(kaptive3.out.kaptive_tsv.collect())
                if (!params.skip_clair3) {
                        minimap(ch_samplesheet_ONT.join(kaptive3.out.kaptive_results))
                        clair3(minimap.out.minimap_results)
                        if (!params.skip_snpeff) {
                                snpeff(clair3.out.clair3_results)
                                snpsift(snpeff.out.snpeff_results)
                                clair_tab_ch=snpsift.out.snpsift_results.collect()
                                kaptive_summary_ch=summary_kaptive.out.kaptive_summary
                                report(clair_tab_ch,kaptive_summary_ch,petg_report_ch,mlst_report_ch)
                        }
                }
        }       
        if (!params.skip_bakta) {
                if (!params.skip_download_bakta_db) {
                        download_bakta_db()
                        if (!params.skip_polishing) {
                                bakta(medaka.out.polished_medaka.combine(download_bakta_db.out.bakta_db))
                        } else if (params.skip_polishing) {
                                bakta(flye.out.assembly_only.combine(download_bakta_db.out.bakta_db))
                        }
                } else {
                        ch_bakta_db=Channel.fromPath( "${params.bakta_db}" ).collect()
                        if (!params.skip_polishing) {
                                bakta(medaka.out.polished_medaka.combine(ch_bakta_db))
                        } else if (params.skip_polishing) {
                                bakta(flye.out.assembly_only.combine(ch_bakta_db))
                        }

                }
        }
        if (!params.skip_amrfinder) {
                if (!params.skip_download_amrfinder_db) {
                        download_amrfinder_db()
                        if (!params.skip_polishing) {
                                amrfinder(medaka.out.polished_medaka.combine(download_amrfinder_db.out.amrfinder_db))
                        } else if (params.skip_polishing) {
                                amrfinder(flye.out.assembly_only.combine(download_amrfinder_db.out.amrfinder_db))
                        }
                } else{
                        ch_amrfinder_db=Channel.fromPath( "${params.amrfinder_db}" ).collect()
                        if (!params.skip_polishing) {
                                amrfinder(medaka.out.polished_medaka.combine(ch_amrfinder_db))
                        } else if (params.skip_polishing) {
                                amrfinder(flye.out.assembly_only.combine(ch_amrfinder_db))
                        } 
                }
                summary_amrfinder(amrfinder.out.amrfinder_results.collect())
        }
        // Build the single combined HTML report from the aggregated 10_report outputs.
        // The LPS database directory is read directly via params.reference_LPS_directory,
        // so only the aggregated TSV/HTML files need to be staged here. This is placed at
        // the end of the workflow because it consumes the late summary_* outputs
        // (amrfinder etc.). It is guarded on the same condition that gates the `report`
        // process so report.out.subtype_report is only referenced when it exists.
        if (!params.skip_kaptive3 && !params.skip_clair3 && !params.skip_snpeff) {
                html_inputs_ch = report.out.subtype_report.flatten()
                        .mix(summary_kaptive.out.kaptive_summary)
                if (!params.skip_nanocomp)  html_inputs_ch = html_inputs_ch.mix(nanocomp.out.nanocomp_html)
                if (!params.skip_assembly)  html_inputs_ch = html_inputs_ch.mix(summary_flye.out.flye_summary)
                if (!params.skip_quast)     html_inputs_ch = html_inputs_ch.mix(summary_quast.out.quast_summary)
                if (!params.skip_checkm)    html_inputs_ch = html_inputs_ch.mix(summary_checkm.out.checkm_summary)
                if (!params.skip_sylph)     html_inputs_ch = html_inputs_ch.mix(summary_sylph.out.sylph_summary)
                if (!params.skip_mlst)      html_inputs_ch = html_inputs_ch.mix(summary_mlst.out.mlst_summary)
                if (!params.skip_amrfinder) html_inputs_ch = html_inputs_ch.mix(summary_amrfinder.out.amrfinder_summary)

                skipped_list = []
                if (params.skip_nanocomp)  skipped_list << 'nanocomp'
                if (params.skip_assembly)  skipped_list << 'flye'
                if (params.skip_quast)     skipped_list << 'quast'
                if (params.skip_checkm)    skipped_list << 'checkm'
                if (params.skip_sylph)     skipped_list << 'sylph'
                if (params.skip_mlst)      skipped_list << 'mlst'
                if (params.skip_petg)      skipped_list << 'petG'
                if (params.skip_amrfinder) skipped_list << 'amrfinder'
                if (params.skip_bakta)     skipped_list << 'bakta'
                param_note_str = "MLST scheme: ${params.mlst_scheme}; petG reported present when a hit spans > ${params.petg_min_length} bp at >= ${params.petg_min_identity}% identity"
                html_report(html_inputs_ch.collect(), workflow.manifest.version, skipped_list.join(', '), param_note_str)
        }
}
