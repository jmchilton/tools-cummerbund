#!/usr/bin/env python

# Wrapper supports Cuffdiff versions v1.3.0-v2.0

import optparse, os, shutil, subprocess, sys, tempfile

def group_callback( option, op_str, value, parser ):
    groups = []
    flist = []
    for arg in parser.rargs:
        arg = arg.strip()
        if arg[0] is "-":
            break
        elif arg[0] is ",":
            groups.append(flist)
            flist = []
        else:
            flist.append(arg)
    groups.append(flist)

    setattr(parser.values, option.dest, groups)
    
def label_callback( option, op_str, value, parser ):
    labels = []
    for arg in parser.rargs:
        arg = arg.strip()
        if arg[0] is "-":
            break
        else:
            labels.append(arg)

    setattr(parser.values, option.dest, labels)

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()
    
# Copied from sam_to_bam.py:
def check_seq_file( dbkey, cached_seqs_pointer_file ):
    seq_path = ''
    for line in open( cached_seqs_pointer_file ):
        line = line.rstrip( '\r\n' )
        if line and not line.startswith( '#' ) and line.startswith( 'index' ):
            fields = line.split( '\t' )
            if len( fields ) < 3:
                continue
            if fields[1] == dbkey:
                seq_path = fields[2].strip()
                break
    return seq_path

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    
    # Cuffdiff options.
    parser.add_option( '-s', '--inner-dist-std-dev', dest='inner_dist_std_dev', help='The standard deviation for the distribution on inner distances between mate pairs. The default is 20bp.' )
    parser.add_option( '-p', '--num-threads', dest='num_threads', help='Use this many threads to align reads. The default is 1.' )
    parser.add_option( '-m', '--inner-mean-dist', dest='inner_mean_dist', help='This is the expected (mean) inner distance between mate pairs. \
                                                                                For, example, for paired end runs with fragments selected at 300bp, \
                                                                                where each end is 50bp, you should set -r to be 200. The default is 45bp.')
    parser.add_option( '-c', '--min-alignment-count', dest='min_alignment_count', help='The minimum number of alignments in a locus for needed to conduct significance testing on changes in that locus observed between samples. If no testing is performed, changes in the locus are deemed not signficant, and the locus\' observed changes don\'t contribute to correction for multiple testing. The default is 1,000 fragment alignments (up to 2,000 paired reads).' )
    parser.add_option( '--FDR', dest='FDR', help='The allowed false discovery rate. The default is 0.05.' )
    parser.add_option( '-u', '--multi-read-correct', dest='multi_read_correct', action="store_true", help='Tells Cufflinks to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome')

    # Advanced Options:	
    parser.add_option( '--num-importance-samples', dest='num_importance_samples', help='Sets the number of importance samples generated for each locus during abundance estimation. Default: 1000' )
    parser.add_option( '--max-mle-iterations', dest='max_mle_iterations', help='Sets the number of iterations allowed during maximum likelihood estimation of abundances. Default: 5000' )
    
    # Wrapper / Galaxy options.
    parser.add_option( '-f', '--files', dest='groups', action="callback", callback=group_callback, help="Groups to be processed, groups are separated by spaces, replicates in a group comma separated. group1_rep1,group1_rep2 group2_rep1,group2_rep2, ..., groupN_rep1, groupN_rep2" )
    parser.add_option( '-A', '--inputA', dest='inputA', help='A transcript GTF file produced by cufflinks, cuffcompare, or other source.')
    parser.add_option( '-1', '--input1', dest='input1', help='File of RNA-Seq read alignments in the SAM format. SAM is a standard short read alignment, that allows aligners to attach custom tags to individual alignments, and Cufflinks requires that the alignments you supply have some of these tags. Please see Input formats for more details.' )
    parser.add_option( '-2', '--input2', dest='input2', help='File of RNA-Seq read alignments in the SAM format. SAM is a standard short read alignment, that allows aligners to attach custom tags to individual alignments, and Cufflinks requires that the alignments you supply have some of these tags. Please see Input formats for more details.' )

    # Label options
    parser.add_option('-L', '--labels', dest='labels', action="callback", callback=label_callback, help="Labels for the groups the replicates are in.")
    
	# Normalization options.
    parser.add_option( "-N", "--quartile-normalization", dest="do_normalization", action="store_true" )

    # Bias correction options.
    parser.add_option( '-b', dest='do_bias_correction', action="store_true", help='Providing Cufflinks with a multifasta file via this option instructs it to run our new bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates.')
    parser.add_option( '', '--dbkey', dest='dbkey', help='The build of the reference dataset' )
    parser.add_option( '', '--index_dir', dest='index_dir', help='GALAXY_DATA_INDEX_DIR' )
    parser.add_option( '', '--ref_file', dest='ref_file', help='The reference dataset from the history' )

    # Outputs.
    parser.add_option( "--isoforms_fpkm_tracking_output", dest="isoforms_fpkm_tracking_output" )
    parser.add_option( "--genes_fpkm_tracking_output", dest="genes_fpkm_tracking_output" )
    parser.add_option( "--cds_fpkm_tracking_output", dest="cds_fpkm_tracking_output" )
    parser.add_option( "--tss_groups_fpkm_tracking_output", dest="tss_groups_fpkm_tracking_output" )
    parser.add_option( "--isoforms_read_group_tracking_output", dest="isoforms_read_group_tracking_output", default=None)
    parser.add_option( "--genes_read_group_tracking_output", dest="genes_read_group_tracking_output", default=None )
    parser.add_option( "--cds_read_group_tracking_output", dest="cds_read_group_tracking_output", default=None )
    parser.add_option( "--tss_groups_read_group_tracking_output", dest="tss_groups_read_group_tracking_output", default=None )
    parser.add_option( "--isoforms_exp_output", dest="isoforms_exp_output" )
    parser.add_option( "--genes_exp_output", dest="genes_exp_output" )
    parser.add_option( "--tss_groups_exp_output", dest="tss_groups_exp_output" )
    parser.add_option( "--cds_exp_fpkm_tracking_output", dest="cds_exp_fpkm_tracking_output" )
    parser.add_option( "--cds_diff_output", dest="cds_diff_output" )
    parser.add_option( "--isoforms_count_tracking_output", dest="isoforms_count_tracking_output" )
    parser.add_option( "--genes_count_tracking_output", dest="genes_count_tracking_output" )
    parser.add_option( "--cds_count_tracking_output", dest="cds_count_tracking_output" )
    parser.add_option( "--tss_groups_count_tracking_output", dest="tss_groups_count_tracking_output" )
    parser.add_option( "--splicing_diff_output", dest="splicing_diff_output" )
    parser.add_option( "--cds_diff_output", dest="cds_diff_output" )
    parser.add_option( "--promoters_diff_output", dest="promoters_diff_output" )
    parser.add_option( "--run_info_output", dest="run_info_output" )
    parser.add_option( "--read_groups_info_output", dest="read_groups_info_output" )
    parser.add_option( "--cuffdatadir", dest="cuffdatadir", default=None)
    parser.add_option( "--cummeRbund_db_output", dest="cummeRbund_db_output", default=None)
    
    (options, args) = parser.parse_args()
    
    # output version # of tool
    try:
        tmp = tempfile.NamedTemporaryFile().name
        tmp_stdout = open( tmp, 'wb' )
        proc = subprocess.Popen( args='cuffdiff --no-update-check 2>&1', shell=True, stdout=tmp_stdout )
        tmp_stdout.close()
        returncode = proc.wait()
        stdout = None
        for line in open( tmp_stdout.name, 'rb' ):
            if line.lower().find( 'cuffdiff v' ) >= 0:
                stdout = line.strip()
                break
        if stdout:
            sys.stdout.write( '%s\n' % stdout )
        else:
            raise Exception
    except:
        sys.stdout.write( 'Could not determine Cuffdiff version\n' )

    # Make temp directory for output.
    tmp_output_dir = tempfile.mkdtemp()
    cuffdatadir = options.cuffdatadir if options.cuffdatadir else tmp_output_dir
    if not os.path.exists( cuffdatadir ):
        os.makedirs( cuffdatadir )
    
    
    # If doing bias correction, set/link to sequence file.
    if options.do_bias_correction:
        if options.ref_file != 'None':
            # Sequence data from history.
            # Create symbolic link to ref_file so that index will be created in working directory.
            seq_path = os.path.join( cuffdatadir, "ref.fa" )
            os.symlink( options.ref_file, seq_path  )
        else:
            # Sequence data from loc file.
            cached_seqs_pointer_file = os.path.join( options.index_dir, 'sam_fa_indices.loc' )
            if not os.path.exists( cached_seqs_pointer_file ):
                stop_err( 'The required file (%s) does not exist.' % cached_seqs_pointer_file )
            # If found for the dbkey, seq_path will look something like /galaxy/data/equCab2/sam_index/equCab2.fa,
            # and the equCab2.fa file will contain fasta sequences.
            seq_path = check_seq_file( options.dbkey, cached_seqs_pointer_file )
            if seq_path == '':
                stop_err( 'No sequence data found for dbkey %s, so bias correction cannot be used.' % options.dbkey  )            
    
    # Build command.
    
    # Base; always use quiet mode to avoid problems with storing log output.
    cmd = "cuffdiff --no-update-check -q"
    
    # Add options.
    if options.inner_dist_std_dev:
        cmd += ( " -s %i" % int ( options.inner_dist_std_dev ) )
    if options.num_threads:
        cmd += ( " -p %i" % int ( options.num_threads ) )
    if options.inner_mean_dist:
        cmd += ( " -m %i" % int ( options.inner_mean_dist ) )
    if options.min_alignment_count:
        cmd += ( " -c %i" % int ( options.min_alignment_count ) )
    if options.FDR:
        cmd += ( " --FDR %f" % float( options.FDR ) )
    if options.multi_read_correct:
        cmd += ( " -u" )
    if options.num_importance_samples:
        cmd += ( " --num-importance-samples %i" % int ( options.num_importance_samples ) )
    if options.max_mle_iterations:
        cmd += ( " --max-mle-iterations %i" % int ( options.max_mle_iterations ) )
    if options.do_normalization:
        cmd += ( " -N" )
    if options.do_bias_correction:
        cmd += ( " -b %s" % seq_path )
            
    # Add inputs.
    # For replicate analysis: group1_rep1,group1_rep2 groupN_rep1,groupN_rep2
    if options.groups:
        cmd += " --labels "
        for label in options.labels:
            cmd += label + ","
        cmd = cmd[:-1]

        cmd += " " + options.inputA + " "

        for group in options.groups:
            for filename in group:
                cmd += filename + ","
            cmd = cmd[:-1] + " "
    else: 
        cmd += " " + options.inputA + " " + options.input1 + " " + options.input2
        
    # Debugging.
    print cmd

    # Run command.
    try:
        tmp_name = tempfile.NamedTemporaryFile( dir=tmp_output_dir ).name
        tmp_stderr = open( tmp_name, 'wb' )
        proc = subprocess.Popen( args=cmd, shell=True, cwd=cuffdatadir, stderr=tmp_stderr.fileno() )
        returncode = proc.wait()
        tmp_stderr.close()
        
        # Get stderr, allowing for case where it's very large.
        tmp_stderr = open( tmp_name, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stderr.close()
        
        # Error checking.
        if returncode != 0:
            raise Exception, stderr
            
        # check that there are results in the output file
        if len( open( os.path.join( cuffdatadir, "isoforms.fpkm_tracking" ), 'rb' ).read().strip() ) == 0:
            raise Exception, 'The main output file is empty, there may be an error with your input file or settings.'
    except Exception, e:
        stop_err( 'Error running cuffdiff. ' + str( e ) )

        
    # Copy output files from tmp directory to specified files.
    try:
        try:
            if options.isoforms_fpkm_tracking_output and os.path.exists(os.path.join( cuffdatadir, "isoforms.fpkm_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "isoforms.fpkm_tracking" ), options.isoforms_fpkm_tracking_output )
            if options.genes_fpkm_tracking_output and os.path.exists(os.path.join( cuffdatadir, "genes.fpkm_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "genes.fpkm_tracking" ), options.genes_fpkm_tracking_output )
            if options.cds_fpkm_tracking_output and os.path.exists(os.path.join( cuffdatadir, "cds.fpkm_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "cds.fpkm_tracking" ), options.cds_fpkm_tracking_output )
            if options.tss_groups_fpkm_tracking_output and os.path.exists(os.path.join( cuffdatadir, "tss_groups.fpkm_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "tss_groups.fpkm_tracking" ), options.tss_groups_fpkm_tracking_output )

            if options.isoforms_read_group_tracking_output and os.path.exists(os.path.join( cuffdatadir, "isoforms.read_group_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "isoforms.read_group_tracking" ), options.isoforms_read_group_tracking_output )
            if options.genes_read_group_tracking_output and os.path.exists(os.path.join( cuffdatadir, "genes.read_group_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "genes.read_group_tracking" ), options.genes_read_group_tracking_output )
            if options.cds_read_group_tracking_output and os.path.exists(os.path.join( cuffdatadir, "cds.read_group_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "cds.read_group_tracking" ), options.cds_read_group_tracking_output )
            if options.tss_groups_read_group_tracking_output and os.path.exists(os.path.join( cuffdatadir, "tss_groups.read_group_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "tss_groups.read_group_tracking" ), options.tss_groups_read_group_tracking_output )

            if options.isoforms_exp_output and os.path.exists(os.path.join( cuffdatadir, "isoform_exp.diff" )):
                shutil.copyfile( os.path.join( cuffdatadir, "isoform_exp.diff" ), options.isoforms_exp_output )
            if options.genes_exp_output and os.path.exists(os.path.join( cuffdatadir, "gene_exp.diff" )):
                shutil.copyfile( os.path.join( cuffdatadir, "gene_exp.diff" ), options.genes_exp_output )
            if options.cds_exp_fpkm_tracking_output and os.path.exists(os.path.join( cuffdatadir, "cds_exp.diff" )):
                shutil.copyfile( os.path.join( cuffdatadir, "cds_exp.diff" ), options.cds_exp_fpkm_tracking_output )
            if options.tss_groups_exp_output and os.path.exists(os.path.join( cuffdatadir, "tss_group_exp.diff" )):
                shutil.copyfile( os.path.join( cuffdatadir, "tss_group_exp.diff" ), options.tss_groups_exp_output )

            if options.isoforms_count_tracking_output and os.path.exists(os.path.join( cuffdatadir, "isoforms.count_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "isoforms.count_tracking" ), options.isoforms_count_tracking_output )
            if options.genes_count_tracking_output and os.path.exists(os.path.join( cuffdatadir, "genes.count_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "genes.count_tracking" ), options.genes_count_tracking_output )
            if options.cds_count_tracking_output and os.path.exists(os.path.join( cuffdatadir, "cds.count_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "cds.count_tracking" ), options.cds_count_tracking_output )
            if options.tss_groups_count_tracking_output and os.path.exists(os.path.join( cuffdatadir, "tss_groups.count_tracking" )):
                shutil.copyfile( os.path.join( cuffdatadir, "tss_groups.count_tracking" ), options.tss_groups_count_tracking_output )

            if options.cds_diff_output and os.path.exists(os.path.join( cuffdatadir, "cds.diff" )):
                shutil.copyfile( os.path.join( cuffdatadir, "cds.diff" ), options.cds_diff_output )

            if options.splicing_diff_output and os.path.exists(os.path.join( cuffdatadir, "splicing.diff" )):
                shutil.copyfile( os.path.join( cuffdatadir, "splicing.diff" ), options.splicing_diff_output )
            if options.promoters_diff_output and os.path.exists(os.path.join( cuffdatadir, "promoters.diff" )):
                shutil.copyfile( os.path.join( cuffdatadir, "promoters.diff" ), options.promoters_diff_output )    

            if options.run_info_output and os.path.exists(os.path.join( cuffdatadir, "run.info" )):
                shutil.copyfile( os.path.join( cuffdatadir, "run.info" ), options.run_info_output )    
            if options.read_groups_info_output and os.path.exists(os.path.join( cuffdatadir, "read_groups.info" )):
                shutil.copyfile( os.path.join( cuffdatadir, "read_groups.info" ), options.read_groups_info_output )    

        except Exception, e:
            stop_err( 'Error in cuffdiff:\n' + str( e ) ) 
        if options.cummeRbund_db_output:
            try:
                dbFile = 'cuffData.db'
                rscript = tempfile.NamedTemporaryFile( dir=tmp_output_dir,suffix='.r' ).name
                rscript_fh = open( rscript, 'wb' )
                rscript_fh.write('library(cummeRbund)\n')
                if options.inputA and options.ref_file:
                    rscript_fh.write('cuff<-readCufflinks(dir = "%s", dbFile = "%s", gtfFile = "%s", genome = "%s", rebuild = T)\n' % (cuffdatadir,dbFile,options.inputA,options.ref_file))
                else:
                    rscript_fh.write('cuff<-readCufflinks(dir = "%s", dbFile = "%s", rebuild = T)\n' % (cuffdatadir,dbFile))
                rscript_fh.close()
                cmd = ( "Rscript --vanilla %s" % rscript )
                tmp_name = tempfile.NamedTemporaryFile( dir=tmp_output_dir ).name
                tmp_stderr = open( tmp_name, 'wb' )
                proc = subprocess.Popen( args=cmd, shell=True, stderr=tmp_stderr.fileno() )
                #proc = subprocess.Popen( args=cmd, shell=True)
                returncode = proc.wait()
                tmp_stderr.close()
                if os.path.exists(os.path.join( cuffdatadir, dbFile )):
                    shutil.copyfile( os.path.join( cuffdatadir, dbFile ), options.cummeRbund_db_output )    
                    shutil.rmtree(os.path.join( cuffdatadir, dbFile ))
            except Exception, e:
                stop_err( 'Error generating cummeRbund cuffData.db:\n' + str( e ) ) 
    finally:
        # Clean up temp dirs
        if os.path.exists( tmp_output_dir ):
            shutil.rmtree( tmp_output_dir )

if __name__=="__main__": __main__()
