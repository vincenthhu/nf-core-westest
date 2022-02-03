#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/westest
========================================================================================
    Github : https://github.com/nf-core/westest
    Website: https://nf-co.re/westest
    Slack  : https://nfcore.slack.com/channels/westest
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta            = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai        = WoekflowMain.getGenomeAttribute(params, 'fai')
params.dict             = WorkflowMain.getGenomeAttribute(params, 'dict')
params.chromsize        = WorkflowMain.getGenomeAttribute(params, 'chromsize')
params.bwa              = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.dbsnp            = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.g1000snp         = WorkflowMain.getGenomeAttribute(params, 'g1000snp')
params.hapmap           = WorkflowMain.getGenomeAttribute(params, 'hapmap')
params.mills            = WorkflowMain.getGenomeAttribute(params, 'mills')
params.omni             = WorkflowMain.getGenomeAttribute(params, 'omni')
params.axiom            = WorkflowMain.getGenomeAttribute(params, 'axiom')
params.knownindels      = WorkflowMain.getGenomeAttribute(params, 'knownindels')
params.cosmic           = WorkflowMain.getGenomeAttribute(params, 'cosmic')
params.gnomad           = WorkflowMain.getGenomeAttribute(params, 'gnomad')

params.target_bed       = WorkflowMain.getGenomeAttribute(params, 'target_bed')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { WESTEST } from './workflows/westest'

//
// WORKFLOW: Run main nf-core/westest analysis pipeline
//
workflow NFCORE_WESTEST {
    WESTEST ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_WESTEST ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
