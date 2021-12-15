#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/nfcoretest
========================================================================================
    Github : https://github.com/nf-core/nfcoretest
    Website: https://nf-co.re/nfcoretest
    Slack  : https://nfcore.slack.com/channels/nfcoretest
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

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

include { NFCORETEST } from './workflows/nfcoretest'

//
// WORKFLOW: Run main nf-core/nfcoretest analysis pipeline
//
workflow NFCORE_NFCORETEST {
    NFCORETEST ()
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
    NFCORE_NFCORETEST ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
