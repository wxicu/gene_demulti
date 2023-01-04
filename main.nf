#!/usr/bin/env nextflow
nextflow.enable.dsl=2
'''
BAM Index should in the format BAMFILE.bam.bai

'''
include { data_preprocess } from './modules/samtools'
include { filter_variant } from './modules/bcftools'
include { variant_cellSNP } from './modules/cellsnp'
include { variant_freebayes } from './modules/freebayes'
include { demultiplex_demuxlet } from './modules/demuxlet'
include { demultiplex_freemuxlet } from './modules/freemuxlet'
include { demultiplex_scSplit } from './modules/scsplit'
include { demultiplex_souporcell } from './modules/souporcell'
include { demultiplex_vireo } from './modules/vireo'

process summary{
    publishDir "$params.outdir/compare", mode: 'copy'
    input:
        val demuxlet_result
        val freemuxlet_result
        val vireo_result
        val souporcell_result
        val scsplit_result
    output:
        path '*.csv'

    script:
        def demuxlet_files = ""
        def freemuxlet_files = ""
        def vireo_files = ""
        def souporcell_files = ""
        def scsplit_files = ""
        
        if (demuxlet_result != "no_result"){
            demuxlet_files = "--demuxlet "
            for(r : demuxlet_result) {
                demuxlet_files  = demuxlet_files + r + ":"
            }
        }
        if (freemuxlet_result != "no_result"){
            freemuxlet_files = "--freemuxlet "
            for(r : freemuxlet_result) {
                freemuxlet_files = freemuxlet_files + r + ":"
            }
        }
        if (vireo_result != "no_result"){
            vireo_files = "--vireo "
            for(r : vireo_result) {
                vireo_files = vireo_files + r + ":"
            }
        }
        if (souporcell_result != "no_result"){
            souporcell_files = "--souporcell "
            for(r : souporcell_result) {
                souporcell_files = souporcell_files + r + ":"
            }
        }
        if (scsplit_result != "no_result"){
            scsplit_files = "--scsplit "
            for(r : scsplit_result) {
                scsplit_files = scsplit_files + r + ":"
            }
        }
        
        """
        summary_gene.R $demuxlet_files $vireo_files $souporcell_files $scsplit_files $freemuxlet_files
        """
}


workflow{
    input_bam = params.cellranger == 'True'? Channel.value(params.cellranger_dir).map{ return it + "/outs/possorted_genome_bam.bam"} : Channel.fromPath(params.bam)
    input_bai = params.cellranger == 'True'? Channel.value(params.cellranger_dir).map{ return it + "/outs/possorted_genome_bam.bam.bai"}: Channel.fromPath(params.bam).map{ return it + ".bai"}

    if ((params.demuxlet == "True" & params.demuxlet_preprocess != 'False')| \
       (params.freemuxlet == "True" & params.freemuxlet_preprocess != 'False')| \
       (params.scSplit == "True" & params.scSplit_preprocess != 'False') | \
       (params.souporcell == "True" & params.souporcell_preprocess != 'False')){
        data_preprocess(input_bam)
        qc_bam = data_preprocess.out.map{ return it + "/sorted.bam"}
        qc_bam_bai = data_preprocess.out.map{ return it + "/sorted.bam.bai"}
    }

    if (params.vireo == "True" &  params.vireo_variant != 'False'){
        if (params.vireo_variant == 'cellSNP' | params.vireo_variant == 'Both'){
            if(params.vireo_preprocess != 'False'){
                variant_cellSNP(qc_bam)
            }
            else{
                variant_cellSNP(input_bam)
            }
            cellsnp_vcf = variant_cellSNP.out.map{ return it + "/cellSNP.cells.vcf"}
        }
        else{
            if (params.scSplit == "True" & params.scSplit_variant != 'False'){
                if (params.scSplit_variant == 'cellSNP' | params.scSplit_variant == 'Both'){
                    if(params.scSplit_preprocess != 'False'){
                        variant_cellSNP(qc_bam)
                    }
                    else{
                        variant_cellSNP(input_bam)
                    }
                    cellsnp_vcf = variant_cellSNP.out.map{ return it + "/cellSNP.cells.vcf"}
                }
            }
        }
       
        if (params.vireo_variant == "freebayes" | params.vireo_variant == 'Both'){
            freebayes_region = Channel.from(1..22, "X","Y").flatten()
            if (params.region != "False"){
                freebayes_region = Channel.value(params.region)
            }
            if(params.vireo_preprocess != 'False'){
                variant_freebayes(qc_bam, qc_bam_bai, freebayes_region)
            }
            else{
                variant_freebayes(input_bam, input_bai, freebayes_region)
            }
            filter_variant(variant_freebayes.out, "True","True")
            freebayes_vcf = filter_variant.out.map{ return it + "/filtered_sorted_total_chroms.vcf"}

        }
        else{
            if (params.scSplit == "True" & params.scSplit_variant != 'False'){
                if (params.scSplit_variant == 'freebayes' | params.scSplit_variant == 'Both'){
                    freebayes_region = Channel.from(1..22, "X","Y").flatten()
                    if (params.region != "False"){
                        freebayes_region = Channel.value(params.region)
                    }
                    if(params.scSplit_preprocess != 'False'){
                        variant_freebayes(qc_bam, qc_bam_bai, freebayes_region)
                    }
                    else{
                        variant_freebayes(input_bam, input_bai, freebayes_region)
                    }
                    filter_variant(variant_freebayes.out, "True","True")
                    freebayes_vcf = filter_variant.out.map{ return it + "/filtered_sorted_total_chroms.vcf"}
                }
            }
        }
    }
    else{
        if (params.scSplit == "True" & params.scSplit_variant != 'False'){
            if (params.scSplit_variant == 'cellSNP' | params.scSplit_variant == 'Both'){
                if(params.scSplit_preprocess != 'False'){
                    variant_cellSNP(qc_bam)
                }
                else{
                    variant_cellSNP(input_bam)
                }
                cellsnp_vcf = variant_cellSNP.out.map{ return it + "/cellSNP.cells.vcf"}
            }
            if (params.scSplit_variant == "freebayes" | params.scSplit_variant == 'Both'){
                freebayes_region = Channel.from(1..22, "X","Y","MT").flatten()
                if (params.region != "False"){
                    freebayes_region = Channel.value(params.region)
                }
                if(params.scSplit_preprocess != 'False'){
                    variant_freebayes(qc_bam, qc_bam_bai, freebayes_region)
                }
                else{
                    variant_freebayes(input_bam, input_bai, freebayes_region)
                }
                filter_variant(variant_freebayes.out, "True","True")
                freebayes_vcf = filter_variant.out.map{ return it + "/filtered_sorted_total_chroms.vcf"}
            }
        }
    }
  
    
    if (params.demuxlet == "True"){
        bam = params.demuxlet_preprocess == 'True'? qc_bam: (params.demuxlet_preprocess == 'False'? input_bam : qc_bam.mix(input_bam))
        demultiplex_demuxlet(bam)
        demuxlet_out = demultiplex_demuxlet.out
    }
    else{
        demuxlet_out = channel.value("no_result")
    }
    
    
    if (params.freemuxlet == "True"){
        bam = params.freemuxlet_preprocess == 'True'? qc_bam: (params.freemuxlet_preprocess == 'False'? input_bam : qc_bam.mix(input_bam))
        demultiplex_freemuxlet(bam)
        freemuxlet_out = demultiplex_freemuxlet.out
    }
    else{
        freemuxlet_out = channel.value("no_result")
    }

    
    if (params.vireo == "True"){
        vcf = params.vireo_variant == 'False'? Channel.fromPath(params.celldata): \
            (params.vireo_variant == 'freebayes'? freebayes_vcf: \
            (params.vireo_variant == 'cellSNP'? cellsnp_vcf : \
            (params.vireo_variant == 'Both'? cellsnp_vcf.mix(freebayes_vcf) : \
            (params.vireo_variant == 'noCellSNP'? freebayes_vcf.mix(Channel.fromPath(params.celldata)) : \
            (params.vireo_variant == 'noFreebayes'? cellsnp_vcf.mix(Channel.fromPath(params.celldata)) : \
            cellsnp_vcf.mix(freebayes_vcf, Channel.fromPath(params.celldata)))))))
        demultiplex_vireo(vcf)
        vireo_out = demultiplex_vireo.out
    }
    else{
        vireo_out = channel.value("no_result")
    }

   
    if (params.scSplit == "True"){
        bam = params.scSplit_preprocess == 'True'? qc_bam: (params.scSplit_preprocess == 'False'? input_bam : qc_bam.mix(input_bam))
        bai = params.scSplit_preprocess == 'True'? qc_bam_bai: (params.scSplit_preprocess == 'False'? input_bai : qc_bam_bai.mix(input_bai))
        vcf = params.scSplit_variant == 'False'? Channel.fromPath(params.vcf_mixed): \
            (params.scSplit_variant == 'freebayes'? freebayes_vcf: \
            (params.scSplit_variant == 'cellSNP'? cellsnp_vcf : \
            (params.scSplit_variant == 'Both'? cellsnp_vcf.mix(freebayes_vcf) : \
            (params.scSplit_variant == 'noCellSNP'? freebayes_vcf.mix(Channel.fromPath(params.vcf_mixed)) : \
            (params.scSplit_variant == 'noFreebayes'? cellsnp_vcf.mix(Channel.fromPath(params.vcf_mixed)) : \
            cellsnp_vcf.mix(freebayes_vcf, Channel.fromPath(params.vcf_mixed)))))))
        demultiplex_scSplit(bam,vcf,bai)
        scSplit_out = demultiplex_scSplit.out
    }
    else{
        scSplit_out = channel.value("no_result")
    }

    
    if (params.souporcell == "True"){
        bam = params.souporcell_preprocess == 'True'? qc_bam: (params.souporcell_preprocess == 'False'? input_bam : qc_bam.mix(input_bam))
        demultiplex_souporcell(bam)
        souporcell_out = demultiplex_souporcell.out
    }
    else{
        souporcell_out = channel.value("no_result")
    }

    summary(demuxlet_out, freemuxlet_out, vireo_out, scSplit_out, souporcell_out)
}

