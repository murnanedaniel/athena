#! /usr/bin/env python

# Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration

# Pythonized version of MadGraph steering executables
#    written by Zach Marshall <zach.marshall@cern.ch>
#    updates for aMC@NLO by Josh McFayden <mcfayden@cern.ch>
#  Attempts to remove path-dependence of MadGraph

import os,sys,time,subprocess,shutil,glob,re,difflib,stat
from AthenaCommon import Logging
mglog = Logging.logging.getLogger('MadGraphUtils')

# Magic name of gridpack directory
MADGRAPH_GRIDPACK_LOCATION='madevent'
# Name for the run (since we only have 1, just needs consistency)
MADGRAPH_RUN_NAME='run_01'
# PDF setting (global setting)
MADGRAPH_PDFSETTING=None
from MadGraphUtilsHelpers import *
import MadGraphSystematicsUtils

def setup_path_protection():
    # Addition for models directory
    if 'PYTHONPATH' in os.environ:
        if not 'Generators/madgraph/models' in os.environ['PYTHONPATH']:
            os.environ['PYTHONPATH'] += ':/cvmfs/atlas.cern.ch/repo/sw/Generators/madgraph/models/latest'
    # Make sure that gfortran doesn't write to somewhere it shouldn't
    if 'GFORTRAN_TMPDIR' in os.environ:
        return
    if 'TMPDIR' in os.environ:
        os.environ['GFORTRAN_TMPDIR']=os.environ['TMPDIR']
        return
    if 'TMP' in os.environ:
        os.environ['GFORTRAN_TMPDIR']=os.environ['TMP']
        return


def config_only_check():
    try:
        from __main__ import opts
        if opts.config_only:
            mglog.info('Athena running on config only mode: not executing MadGraph')
            return True
    except:
        pass
    return False


def new_process(process='generate p p > t t~'):
    """ Generate a new process in madgraph.
    Pass a process string.
    Return the name of the process directory.
    """
    if config_only_check(): return

    # Don't run if generating events from gridpack
    if is_gen_from_gridpack():
        return MADGRAPH_GRIDPACK_LOCATION

    # Actually just sent the process card contents - let's make a card
    card_loc='proc_card_mg5.dat'
    mglog.info('Writing process card to '+card_loc)
    a_card = open( card_loc , 'w' )
    a_card.write( process )
    a_card.close()

    madpath=os.environ['MADPATH']
    # Just in case
    setup_path_protection()

    # Check if we have a special output directory
    process_dir = ''
    for l in process.split('\n'):
        # Look for an output line
        if 'output' not in l.split('#')[0].split(): continue
        # Check how many things before the options start
        tmplist = l.split('#')[0].split(' -')[0]
        # if two things, second is the directory
        if len(tmplist.split())==2: process_dir = tmplist.split()[1]
        # if three things, third is the directory (second is the format)
        elif len(tmplist.split())==3: process_dir = tmplist.split()[2]
        # See if we got a directory
        if ''!=process_dir:
            mglog.info('Saw that you asked for a special output directory: '+str(process_dir))
        break

    mglog.info('Started process generation at '+str(time.asctime()))

    generate = subprocess.Popen([madpath+'/bin/mg5_aMC',card_loc],stdin=subprocess.PIPE)
    generate.communicate()

    mglog.info('Finished process generation at '+str(time.asctime()))

    # at this point process_dir is for sure defined - it's equal to '' in the worst case
    if process_dir == '': # no user-defined value, need to find the directory created by MadGraph5
        for adir in sorted(glob.glob( os.getcwd()+'/*PROC*' ),reverse=True):
            if os.access('%s/SubProcesses/subproc.mg'%adir,os.R_OK):
                if process_dir=='':
                    process_dir=adir
                else:
                    mglog.warning('Additional possible process directory, '+adir+' found. Had '+process_dir)
                    mglog.warning('Likely this is because you did not run from a clean directory, and this may cause errors later.')
    else: # user-defined directory
        if not os.access('%s/SubProcesses/subproc.mg'%process_dir,os.R_OK):
            raise RuntimeError('No diagrams for this process in user-define dir='+str(process_dir))
    if process_dir=='':
        raise RuntimeError('No diagrams for this process from list: '+str(sorted(glob.glob(os.getcwd()+'/*PROC*'),reverse=True)))

    # Special catch related to path setting and using afs
    needed_options = ['ninja','collier','fastjet','lhapdf','syscalc_path']
    in_config = open(os.environ['MADPATH']+'/input/mg5_configuration.txt','r')
    option_paths = {}
    for l in in_config.readlines():
        for o in needed_options:
            if o+' =' in l.split('#')[0] and 'MCGenerators' in l.split('#')[0]:
                old_path = l.split('#')[0].split('=')[1].strip().split('MCGenerators')[1]
                old_path = old_path[ old_path.find('/') : ]
                if o =='lhapdf' and 'LHAPATH' in os.environ:
                    # Patch for LHAPDF version
                    version = os.environ['LHAPATH'].split('lhapdf/')[1].split('/')[0]
                    old_version = old_path.split('lhapdf/')[1].split('/')[0]
                    old_path = old_path.replace(old_version,version)
                if o=='ninja':
                    # Patch for stupid naming problem
                    old_path.replace('gosam_contrib','gosam-contrib')
                option_paths[o] = os.environ['MADPATH'].split('madgraph5amc')[0]+old_path
            # Check to see if the option has been commented out
            if o+' =' in l and o+' =' not in l.split('#')[0]:
                mglog.info('Option '+o+' appears commented out in the config file')

    in_config.close()
    for o in needed_options:
        if not o in option_paths: mglog.warning('Path for option '+o+' not found in original config')

    mglog.info('Modifying config paths to avoid use of afs:')
    mglog.info(option_paths)

    # Get the original configuration information
    if not is_NLO_run(process_dir=process_dir):
        config_card=process_dir+'/Cards/me5_configuration.txt'
    else:
        config_card=process_dir+'/Cards/amcatnlo_configuration.txt'
    with file(config_card, 'r') as original:
        data = original.readlines()
    with file(config_card, 'w') as modified:
        for l in data:
            written = False
            for o in needed_options:
                if o+' =' in l.split('#')[0] and o in option_paths:
                    modified.write( o +' = '+option_paths[o]+'\n' )
                    written = True
                    break
            if not written:
                modified.write(l)
    # Done modifying paths

    return process_dir


def get_default_runcard(process_dir=MADGRAPH_GRIDPACK_LOCATION):
    """ Copy the default runcard from one of several locations
    to a local file with name run_card.tmp.dat"""
    output_name = 'run_card.tmp.dat'
    if config_only_check():
        mglog.info('Athena running on config only mode: grabbing run card the old way, as there will be no proc dir')
        mglog.info('Fetching default LO run_card.dat')
        if os.access(os.environ['MADPATH']+'/Template/LO/Cards/run_card.dat',os.R_OK):
            shutil.copy(os.environ['MADPATH']+'/Template/LO/Cards/run_card.dat',output_name)
            return 'run_card.dat'
        elif os.access(os.environ['MADPATH']+'/Template/Cards/run_card.dat',os.R_OK):
            shutil.copy(os.environ['MADPATH']+'/Template/Cards/run_card.dat',output_name)
            return output_name
        else:
            raise RuntimeError('Cannot find default LO run_card.dat!')

    # Get the run card from the installation
    run_card=process_dir+'/Cards/run_card.dat'
    if os.access(run_card,os.R_OK):
        mglog.info('Copying default run_card.dat from '+str(run_card))
        shutil.copy(run_card,output_name)
        return output_name
    else:
        run_card=process_dir+'/Cards/run_card_default.dat'
        mglog.info('Fetching default run_card.dat from '+str(run_card))
        if os.access(run_card,os.R_OK):
            shutil.copy(run_card,output_name)
            return output_name
        else:
            raise RuntimeError('Cannot find default run_card.dat or run_card_default.dat! I was looking here: %s'%run_card)


def generate(process_dir='PROC_mssm_0',grid_pack=False,gridpack_compile=False,cluster_type=None,cluster_queue=None,cluster_nb_retry=None,cluster_temp_path=None,extlhapath=None,required_accuracy=0.01,runArgs=None,reweight_card=None,bias_module=None):
    if config_only_check(): return

    # Just in case
    setup_path_protection()

    # Set consistent mode and number of jobs
    mode = 0
    njobs = 1
    if 'ATHENA_PROC_NUMBER' in os.environ:
        njobs = int(os.environ['ATHENA_PROC_NUMBER'])
        mglog.info('Lucky you - you are running on a full node queue.  Will re-configure for '+str(njobs)+' jobs.')
        mode = 2
    if cluster_type is not None:
        mode = 1

    if is_gen_from_gridpack():
        mglog.info('Running event generation from gridpack (using smarter mode from generate() function)')
        generate_from_gridpack(runArgs=runArgs,extlhapath=extlhapath,gridpack_compile=gridpack_compile,reweight_card=reweight_card)
        return

    # Now get a variety of info out of the runArgs
    beamEnergy,random_seed = get_runArgs_info(runArgs)

    # If we need to get the cards...
    if reweight_card is not None and not os.access(reweight_card,os.R_OK):
        raise RuntimeError('Could not find reweight card '+str(reweight_card))

    # Check if process is NLO or LO
    isNLO=is_NLO_run(process_dir=process_dir)

    # use f2py2 if f2py not available
    if reweight_card is not None:
        from distutils.spawn import find_executable
        if find_executable('f2py2') is not None:
            mglog.info('found f2py2, will update configuration')
            if isNLO:
                config_card=process_dir+'/Cards/amcatnlo_configuration.txt'
            else:
                config_card=process_dir+'/Cards/me5_configuration.txt'
            shutil.move(config_card,config_card+'.old')
            oldcard = open(config_card+'.old','r')
            newcard = open(config_card,'w')
            for line in oldcard:
                if 'f2py_compiler' in line:
                    newcard.write(' f2py_compiler = f2py2\n')
                else:
                    newcard.write(line)
            oldcard.close()
            newcard.close()
        elif find_executable('f2py') is not None:
            mglog.info('Found f2py, will use it for reweighting')
        else:
            raise RuntimeError('Could not find f2py or f2py2, needed for reweighting')

    if grid_pack:
        #Running in gridpack mode
        mglog.info('Started generating gridpack at '+str(time.asctime()))
        mglog.warning(' >>>>>> THIS KIND OF JOB SHOULD ONLY BE RUN LOCALLY - NOT IN GRID JOBS <<<<<<')

        if not isNLO:
            modify_run_card(process_dir=process_dir,settings={'gridpack':'true'})
        else:
            my_settings = {'nevents':'1000','req_acc':str(required_accuracy)}
            modify_run_card(process_dir=process_dir,settings=my_settings)
    else:
        #Running in on-the-fly mode
        mglog.info('Started generating at '+str(time.asctime()))

    mglog.info('Run '+MADGRAPH_RUN_NAME+' will be performed in mode '+str(mode)+' with '+str(njobs)+' jobs in parallel.')

    if reweight_card:
        mglog.info('Running reweighting module. Moving card (%s) into place.'%(reweight_card))
        shutil.copyfile(reweight_card,process_dir+'/Cards/reweight_card.dat')
        check_reweight_card(process_dir+'/Cards/reweight_card.dat')

    # Ensure that things are set up normally
    if reweight_card is not None and not os.access(reweight_card,os.R_OK):
        raise RuntimeError('No reweight card found at '+reweight_card)

    if not os.access(process_dir,os.R_OK):
        raise RuntimeError('No process directory found at '+process_dir)
    if not os.access(process_dir+'/bin/generate_events',os.R_OK):
        raise RuntimeError('No generate_events module found in '+process_dir)

    allow_links = True
    if cluster_type is not None:
        if 'condor' in cluster_type.lower():
            mglog.warning('Condor clusters do not allow links.  Will do more copying rather than linking')
            allow_links = False

    (LHAPATH,origLHAPATH,origLHAPDF_DATA_PATH) = setupLHAPDF(isNLO, process_dir=process_dir, extlhapath=extlhapath, allow_links=allow_links)

    mglog.info('For your information, the libraries available are (should include LHAPDF):')
    mglog.info( sorted( os.listdir( process_dir+'/lib/' ) ) )

    setupFastjet(isNLO, process_dir=process_dir)
    if bias_module!=None:
        setup_bias_module(bias_module,run_card,process_dir)

    mglog.info('Now I will hack the make files a bit.  Apologies, but there seems to be no good way around this.')
    shutil.copyfile(process_dir+'/Source/make_opts',process_dir+'/Source/make_opts_old')
    old_opts = open(process_dir+'/Source/make_opts_old','r')
    new_opts = open(process_dir+'/Source/make_opts','w')
    for aline in old_opts:
        if 'FC=g' in aline:
            mglog.info('Configuring the fancy gfortran compiler instead of g77 / f77')
            new_opts.write('  FC=gfortran\n')
        elif 'FFLAGS+= -ffixed-line-length-132' in aline and 'i686' in os.environ['CMTCONFIG']:
            mglog.info('Setting you up for a 32-bit compilation')
            new_opts.write('FFLAGS+= -ffixed-line-length-132 -m32\n')
        else:
            new_opts.write(aline)
    old_opts.close()
    new_opts.close()
    mglog.info('Make file hacking complete.')

    print_cards_from_dir(process_dir=process_dir)

    currdir=os.getcwd()
    os.chdir(process_dir)

    # Check the run card
    run_card_consistency_check(isNLO=isNLO)

    if mode!=0 and not isNLO:

        if mode==1:
            mglog.info('Setting up parallel running system settings')

            setNCores(process_dir=os.getcwd(), Ncores=njobs)

            config_card='Cards/me5_configuration.txt'
            oldcard = open(config_card,'r')
            newcard = open(config_card+'.tmp','w')

            for line in oldcard:
                if cluster_type!=None and 'cluster_type = ' in line:
                    newcard.write('cluster_type = %s \n'%(cluster_type))
                    mglog.info('Setting cluster type = %s in %s'%(cluster_type,config_card))
                elif cluster_queue!=None and 'cluster_queue = ' in line:
                    newcard.write('cluster_queue = %s \n'%(cluster_queue))
                    mglog.info('Setting cluster queue = %s in %s'%(cluster_queue,config_card))
                else:
                    newcard.write(line)
            oldcard.close()
            newcard.close()
            shutil.move(config_card+'.tmp',config_card)

            mglog.info('New me5_configuration.txt card:')
            configCard = subprocess.Popen(['cat',config_card])
            configCard.wait()

            if cluster_type=='pbs':
                mglog.info('Modifying bin/internal/cluster.py for PBS cluster running')
                os.system("sed -i \"s:text += prog:text += './'+prog:g\" bin/internal/cluster.py")

        generate = subprocess.Popen(['bin/generate_events',str(mode),str(njobs),MADGRAPH_RUN_NAME],stdin=subprocess.PIPE)
        generate.communicate()

    elif not isNLO:

        setNCores(process_dir=os.getcwd(), Ncores=njobs)
        mglog.info('Running serial generation.  This will take a bit more time than parallel generation.')
        generate = subprocess.Popen(['bin/generate_events','0',MADGRAPH_RUN_NAME],stdin=subprocess.PIPE)
        generate.communicate()

    elif isNLO:

        ### Editing config card
        config_card='Cards/amcatnlo_configuration.txt'
        oldcard = open(config_card,'r')
        newcard = open(config_card+'.tmp','w')

        # Make sure params only set once
        run_mode_set=False
        auto_html_set=False
        cltype_set=False
        clqueue_set=False
        nbcore_set=False
        tmppath_set=False
        cluster_nb_retry_set=False

        for line in oldcard:
            if 'run_mode =' in line:
                if not run_mode_set:
                    mglog.info('Setting run_mode = %i'%(mode))
                    newcard.write('run_mode = %i \n'%(mode))
                    run_mode_set=True
            elif 'automatic_html_opening =' in line:
                if not auto_html_set:
                    mglog.info('Setting automatic_html_opening = %s'%('False'))
                    newcard.write('automatic_html_opening = %s \n'%('False'))
                    auto_html_set=True
            elif 'cluster_type = ' in line and mode == 1 and cluster_type:
                if not cltype_set:
                    mglog.info('Setting cluster type = %s in %s'%(cluster_type,config_card))
                    newcard.write('cluster_type = %s \n'%(cluster_type))
                    cltype_set=True
            elif 'cluster_queue = ' in line and mode == 1 and cluster_queue:
                if not clqueue_set:
                    mglog.info('Setting cluster queue = %s in %s'%(cluster_queue,config_card))
                    newcard.write('cluster_queue = %s \n'%(cluster_queue))
                    clqueue_set=True
            elif 'cluster_nb_retry = ' in line and mode == 1 and cluster_nb_retry:
                if not cluster_nb_retry_set:
                    mglog.info('Setting cluster no. retries = %i in %s'%(cluster_nb_retry,config_card))
                    newcard.write('cluster_nb_retry = %i \n'%(cluster_nb_retry))
                    cluster_nb_retry_set=True
            elif 'cluster_temp_path = ' in line and mode == 1 and cluster_temp_path:
                if not tmppath_set:
                    mglog.info('Setting cluster temp path = %s in %s'%(cluster_temp_path,config_card))
                    newcard.write('cluster_temp_path = %s \n'%(cluster_temp_path))
                    tmppath_set=True
            elif 'nb_core = ' in line:
                if not nbcore_set:
                    mglog.info('Setting number of cores = %i in %s'%(njobs,config_card))
                    newcard.write('nb_core = %i \n'%(njobs))
                    nbcore_set=True
            else:
                newcard.write(line)
        oldcard.close()
        newcard.close()
        shutil.move(config_card+'.tmp',config_card)

        mglog.info( "amcatnlo_configuration.txt: "+config_card )
        configCard = subprocess.Popen(['cat',config_card])
        configCard.wait()

        mglog.info('Removing Cards/shower_card.dat to ensure we get parton level events only')
        remove_shower = subprocess.Popen(['rm','Cards/shower_card.dat'])
        remove_shower.wait()

        mygenerate = subprocess.Popen(['bin/generate_events','--name='+MADGRAPH_RUN_NAME],stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
        mygenerate.communicate()

    # Get back to where we came from
    os.chdir(currdir)

    if grid_pack:
        # Name dictacted by https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/PmgMcSoftware
        gridpack_name='mc_'+str(int(beamEnergy*2/1000))+'TeV.'+get_physics_short()+'.GRID.tar.gz'
        mglog.info('Tidying up gridpack (%s)...'%gridpack_name)

        if not isNLO:
            ### LO RUN ###
            if madspin_card:
                shutil.copy((process_dir+'/'+MADGRAPH_RUN_NAME+'_decayed_1_gridpack.tar.gz'),gridpack_name)
            else:
                shutil.copy((process_dir+'/'+MADGRAPH_RUN_NAME+'_gridpack.tar.gz'),gridpack_name)

            if gridpack_compile:
                mkdir = subprocess.Popen(['mkdir','tmp%i/'%os.getpid()])
                mkdir.wait()
                os.chdir('tmp%i/'%os.getpid())
                mglog.info('untar gridpack')
                untar = subprocess.Popen(['tar','xvzf',('../'+gridpack_name)])
                untar.wait()
                mglog.info('compile and clean up')
                os.chdir('madevent/')
                compile = subprocess.Popen(['./bin/compile'])
                compile.wait()
                clean = subprocess.Popen(['./bin/clean4grid'])
                clean.wait()
                os.chdir('../')
                mglog.info('remove old tarball')
                remove_old = subprocess.Popen(['rm',('../'+gridpack_name)])
                remove_old.wait()
                mglog.info('Package up new tarball')
                tar = subprocess.Popen(['tar','cvzf','../'+gridpack_name,'--exclude=lib/PDFsets','.'])
                tar.wait()

                os.chdir('../')
                mglog.info('Remove temporary directory')
                remove_tmp = subprocess.Popen(['rm','-fr','tmp%i/'%os.getpid()])
                remove_tmp.wait()
                mglog.info('Tidying up complete!')

        else:

            ### NLO RUN ###
            mglog.info('Package up process_dir')
            os.rename(process_dir,MADGRAPH_GRIDPACK_LOCATION)
            tar = subprocess.Popen(['tar','czf',gridpack_name,MADGRAPH_GRIDPACK_LOCATION,'--exclude=lib/PDFsets'])
            tar.wait()
            os.rename(MADGRAPH_GRIDPACK_LOCATION,process_dir)

        raise RuntimeError('Gridpack sucessfully created, exiting the transform. IGNORE ERRORS if running gridpack generation!')

    resetLHAPDF(origLHAPATH=origLHAPATH,origLHAPDF_DATA_PATH=origLHAPDF_DATA_PATH)

    mglog.info('Finished at '+str(time.asctime()))
    return 0


def generate_from_gridpack(runArgs=None,reweight_card=None,extlhapath=None, gridpack_compile=None):

    # Get of info out of the runArgs
    beamEnergy,random_seed = get_runArgs_info(runArgs)

    # Just in case
    setup_path_protection()

    isNLO=is_NLO_run(process_dir=MADGRAPH_GRIDPACK_LOCATION)
    LHAPATH=os.environ['LHAPATH'].split(':')[0]

    (LHAPATH,origLHAPATH,origLHAPDF_DATA_PATH) = setupLHAPDF(isNLO, process_dir=MADGRAPH_GRIDPACK_LOCATION, extlhapath=extlhapath)

    setupFastjet(isNLO, process_dir=MADGRAPH_GRIDPACK_LOCATION)

    # Ensure that we only do madspin at the end
    if os.access(MADGRAPH_GRIDPACK_LOCATION+'/Cards/madspin_card.dat',os.R_OK):
        os.rename(MADGRAPH_GRIDPACK_LOCATION+'/Cards/madspin_card.dat',MADGRAPH_GRIDPACK_LOCATION+'/Cards/backup_madspin_card.dat')
        do_madspin=True
    else:
        do_madspin=False

    if reweight_card is not None:
        if os.path.exists(reweight_card):
            shutil.copy( reweight_card , MADGRAPH_GRIDPACK_LOCATION+'/Cards/reweight_card.dat' )
            mglog.info( 'Moved reweight card into place: '+str(reweight_card) )
        else:
            mglog.info( 'Did not find reweight card '+str(reweight_card)+', using the one provided by gridpack' )
        check_reweight_card(MADGRAPH_GRIDPACK_LOCATION+'/Cards/reweight_card.dat')

    # Modify run card, then print
    modify_run_card(process_dir=MADGRAPH_GRIDPACK_LOCATION,settings={'iseed':str(random_seed),'python_seed':str(random_seed)})
    print_cards_from_dir(process_dir=MADGRAPH_GRIDPACK_LOCATION)

    mglog.info('Generating events from gridpack')

    # Ensure that things are set up normally
    if not os.path.exists(MADGRAPH_GRIDPACK_LOCATION):
        raise RuntimeError('Gridpack directory not found at '+MADGRAPH_GRIDPACK_LOCATION)

    nevents = getDictFromCard(MADGRAPH_GRIDPACK_LOCATION+'/Cards/run_card.dat')['nevents']
    mglog.info('>>>> FOUND GRIDPACK <<<<  <- This will be used for generation')
    mglog.info('Generation of '+str(int(nevents))+' events will be performed using the supplied gridpack with random seed '+str(random_seed))
    mglog.info('Started generating events at '+str(time.asctime()))

    #Remove addmasses if it's there
    if os.access(MADGRAPH_GRIDPACK_LOCATION+'/bin/internal/addmasses.py',os.R_OK):
        os.remove(MADGRAPH_GRIDPACK_LOCATION+'/bin/internal/addmasses.py')

    currdir=os.getcwd()

    # Make sure we've set the number of processes appropriately
    setNCores(process_dir=MADGRAPH_GRIDPACK_LOCATION)

    if not isNLO:
        ### LO RUN ###
        if not os.access(MADGRAPH_GRIDPACK_LOCATION+'/bin/run.sh',os.R_OK):
            mglog.error('/bin/run.sh not found at '+MADGRAPH_GRIDPACK_LOCATION)
            raise RuntimeError('Could not find run.sh executable')
        else:
            mglog.info('Found '+MADGRAPH_GRIDPACK_LOCATION+'/bin/run.sh, starting generation.')
        # hack script to add reweighting and systematics, if required
        hack_gridpack_script(reweight_card)

        mglog.info('For your information, ls of '+currdir+':')
        mglog.info( sorted( os.listdir( currdir ) ) )
        mglog.info('For your information, ls of '+MADGRAPH_GRIDPACK_LOCATION+':')
        mglog.info( sorted( os.listdir( MADGRAPH_GRIDPACK_LOCATION ) ) )
        run_card_consistency_check(isNLO=isNLO,process_dir=MADGRAPH_GRIDPACK_LOCATION)
        generate = subprocess.Popen([MADGRAPH_GRIDPACK_LOCATION+'/bin/run.sh',str(int(nevents)),str(int(random_seed))],stdin=subprocess.PIPE)
        generate.communicate()

    else:
        ### NLO RUN ###
        if not os.access(MADGRAPH_GRIDPACK_LOCATION+'/bin/generate_events',os.R_OK):
            mglog.error('bin/generate_events not found at '+MADGRAPH_GRIDPACK_LOCATION)
            raise RuntimeError('Could not find generate_events executable')
        else:
            mglog.info('Found '+MADGRAPH_GRIDPACK_LOCATION+'/bin/generate_events, starting generation.')

        ### Editing config card
        config_card=MADGRAPH_GRIDPACK_LOCATION+'/Cards/amcatnlo_configuration.txt'

        # Make sure params only set once
        cltype_set=False
        clqueue_set=False

        mglog.info( "amcatnlo_configuration.txt: "+config_card )
        configCard = subprocess.Popen(['cat',config_card])
        configCard.wait()

        mglog.info('For your information, ls of '+currdir+':')
        mglog.info( sorted( os.listdir( currdir ) ) )
        mglog.info('For your information, ls of '+MADGRAPH_GRIDPACK_LOCATION+'/Events/:')
        mglog.info( sorted( os.listdir( MADGRAPH_GRIDPACK_LOCATION+'/Events/' ) ) )

        if os.access(MADGRAPH_GRIDPACK_LOCATION+'/Events/'+MADGRAPH_RUN_NAME, os.F_OK):
            mglog.info('Removing %s/Events/%s directory from gridpack generation.'%(MADGRAPH_GRIDPACK_LOCATION,MADGRAPH_RUN_NAME))
            shutil.rmtree(MADGRAPH_GRIDPACK_LOCATION+'/Events/'+MADGRAPH_RUN_NAME)

        # Delete events generated when setting up MadSpin during gridpack generation
        if os.access(MADGRAPH_GRIDPACK_LOCATION+'/Events/'+MADGRAPH_RUN_NAME+'_decayed_1', os.F_OK):
            mglog.info('Removing %s/Events/%s_decayed_1 directory from gridpack generation.'%(MADGRAPH_GRIDPACK_LOCATION,MADGRAPH_RUN_NAME))
            shutil.rmtree(MADGRAPH_GRIDPACK_LOCATION+'/Events/'+MADGRAPH_RUN_NAME+'_decayed_1')

        mglog.info('For your information, ls of '+MADGRAPH_GRIDPACK_LOCATION+'/Events/:')
        mglog.info( sorted( os.listdir( MADGRAPH_GRIDPACK_LOCATION+'/Events/' ) ) )

        run_card_consistency_check(isNLO=isNLO,process_dir=MADGRAPH_GRIDPACK_LOCATION)
        if not gridpack_compile:
            mglog.info('Copying make_opts from Template')
            shutil.copy(os.environ['MADPATH']+'/Template/LO/Source/make_opts',MADGRAPH_GRIDPACK_LOCATION+'/Source/')

            generate = subprocess.Popen([MADGRAPH_GRIDPACK_LOCATION+'/bin/generate_events','--parton','--nocompile','--only_generation','-f','--name=%s'%MADGRAPH_RUN_NAME],stdin=subprocess.PIPE)
            generate.communicate()
        else:
            mglog.info('Allowing recompilation of gridpack')
            if os.path.islink(MADGRAPH_GRIDPACK_LOCATION+'/lib/libLHAPDF.a'):
                mglog.info('Unlinking '+MADGRAPH_GRIDPACK_LOCATION+'/lib/libLHAPDF.a')
                os.unlink(MADGRAPH_GRIDPACK_LOCATION+'/lib/libLHAPDF.a')

            generate = subprocess.Popen([MADGRAPH_GRIDPACK_LOCATION+'/bin/generate_events','--parton','--only_generation','-f','--name=%s'%MADGRAPH_RUN_NAME],stdin=subprocess.PIPE)
            generate.communicate()

    # See if MG5 did the job for us already
    if not os.access('events.lhe.gz',os.R_OK):
        mglog.info('Copying generated events to %s.'%currdir)
        if not os.path.exists(MADGRAPH_GRIDPACK_LOCATION+'Events/GridRun_%i/'%random_seed):
            shutil.copy(MADGRAPH_GRIDPACK_LOCATION+'/Events/'+MADGRAPH_RUN_NAME+'/events.lhe.gz','events.lhe.gz')
    else:
        mglog.info('Events were already in place')

    mglog.info('For your information, ls of '+currdir+':')
    mglog.info( sorted( os.listdir( currdir ) ) )

    mglog.info('Moving generated events to be in correct format for arrange_output().')
    mglog.info('Unzipping generated events.')
    unzip = subprocess.Popen(['gunzip','-f','events.lhe.gz'])
    unzip.wait()

    mglog.info('Moving file over to '+MADGRAPH_GRIDPACK_LOCATION+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe')
    mkdir = subprocess.Popen(['mkdir','-p',(MADGRAPH_GRIDPACK_LOCATION+'/Events/'+MADGRAPH_RUN_NAME)])
    mkdir.wait()
    shutil.move('events.lhe',MADGRAPH_GRIDPACK_LOCATION+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe')

    mglog.info('Re-zipping into dataset name '+MADGRAPH_GRIDPACK_LOCATION+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe.gz')
    rezip = subprocess.Popen(['gzip',MADGRAPH_GRIDPACK_LOCATION+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe'])
    rezip.wait()

    os.chdir(currdir)

    # Now consider MadSpin:
    if do_madspin:
        # Move card back
        os.rename(MADGRAPH_GRIDPACK_LOCATION+'/Cards/backup_madspin_card.dat',MADGRAPH_GRIDPACK_LOCATION+'/Cards/madspin_card.dat')
        mglog.info('Decaying with MadSpin.')
        add_madspin(process_dir=MADGRAPH_GRIDPACK_LOCATION)

    mglog.info('Finished at '+str(time.asctime()))

    resetLHAPDF(origLHAPATH=origLHAPATH,origLHAPDF_DATA_PATH=origLHAPDF_DATA_PATH)

    return 0


def setupFastjet(isNLO, process_dir=None):

    mglog.info('Path to fastjet install dir:%s'%os.environ['FASTJETPATH'])
    fastjetconfig = os.environ['FASTJETPATH']+'/bin/fastjet-config'

    mglog.info('fastjet-config --version:      %s'%str(subprocess.Popen([fastjetconfig, '--version'],stdout = subprocess.PIPE).stdout.read().strip()))
    mglog.info('fastjet-config --prefix:       %s'%str(subprocess.Popen([fastjetconfig, '--prefix'],stdout = subprocess.PIPE).stdout.read().strip()))

    if not isNLO:
        config_card=process_dir+'/Cards/me5_configuration.txt'
    else:
        config_card=process_dir+'/Cards/amcatnlo_configuration.txt'

    oldcard = open(config_card,'r')
    newcard = open(config_card+'.tmp','w')

    for line in oldcard:
        if 'fastjet = ' in line:
            newcard.write('fastjet = %s \n'%(fastjetconfig))
            mglog.info('Setting fastjet = %s in %s'%(fastjetconfig,config_card))
        else:
            newcard.write(line)
    oldcard.close()
    newcard.close()
    shutil.move(config_card+'.tmp',config_card)

    return

def get_LHAPDF_DATA_PATH():
    return get_LHAPDF_PATHS()[1]

def get_LHAPDF_PATHS():
    LHADATAPATH=None
    LHAPATH=None
    for p in os.environ['LHAPATH'].split(':')+os.environ['LHAPDF_DATA_PATH'].split(':'):
        if os.path.exists(p+"/../../lib/") and LHAPATH==None:
            LHAPATH=p
    for p in os.environ['LHAPDF_DATA_PATH'].split(':')+os.environ['LHAPATH'].split(':'):
        if os.path.exists(p) and LHADATAPATH==None and p!=LHAPATH:
            LHADATAPATH=p
    if LHADATAPATH==None:
        LHADATAPATH=LHAPATH
    if LHAPATH==None:
        mglog.error('Could not find path to LHAPDF installation')
    return LHAPATH,LHADATAPATH

# function to get lhapdf id and name from either id or name
def get_lhapdf_id_and_name(pdf):
    pdfname=''
    pdfid=-999
    LHADATAPATH=get_LHAPDF_DATA_PATH()
    pdflist = open(LHADATAPATH+'/pdfsets.index','r')
    if isinstance(pdf,int) or pdf.isdigit():
        pdf=int(pdf)
        pdfid=pdf
        for line in pdflist:
            splitline=line.split()
            if int(splitline[0]) == pdfid:
                pdfname=splitline[1]
                break
    else:
        pdfname=pdf
        for line in pdflist:
            splitline=line.split()
            if splitline[1] == pdfname:
                pdfid=int(splitline[0])
                break
    pdflist.close()

    if pdfname=='':
        err='Couldn\'t find PDF name associated to ID %i in %s.'%(pdfid,LHADATAPATH+'/pdfsets.index')
        mglog.error(err)
        raise RuntimeError(err)
    if pdfid<0:
        err='Couldn\'t find PDF ID associated to name %s in %s.'%(pdfname,LHADATAPATH+'/pdfsets.index')
        mglog.error(err)
        raise RuntimeError(err)

    return pdfid,pdfname

def setupLHAPDF(isNLO, process_dir=None, extlhapath=None, allow_links=True):

    origLHAPATH=os.environ['LHAPATH']
    origLHAPDF_DATA_PATH=os.environ['LHAPDF_DATA_PATH']

    LHAPATH,LHADATAPATH=get_LHAPDF_PATHS()

    pdfname=''
    pdfid=-999

    ### Reading LHAPDF ID from run card
    run_card=process_dir+'/Cards/run_card.dat'
    mydict=getDictFromCard(run_card)

    if mydict["pdlabel"].replace("'","") == 'lhapdf':
        #Make local LHAPDF dir
        mglog.info('creating local LHAPDF dir: MGC_LHAPDF/')
        if os.path.islink('MGC_LHAPDF/'):
            os.unlink('MGC_LHAPDF/')
        elif os.path.isdir('MGC_LHAPDF/'):
            shutil.rmtree('MGC_LHAPDF/')

        newMGCLHA='MGC_LHAPDF/'

        mkdir = subprocess.Popen(['mkdir','-p',newMGCLHA])
        mkdir.wait()

        pdfs_used=[ int(x) for x in mydict['lhaid'].replace(' ',',').split(',') ]
        # included systematics pdfs here
        if 'sys_pdf' in mydict:
            sys_pdf=mydict['sys_pdf'].replace('&&',' ').split()
            for s in sys_pdf:
                if s.isdigit():
                    idx=int(s)
                    if idx>1000: # the sys_pdf syntax is such that small numbers are used to specify the subpdf index
                        pdfs_used.append(idx)
                else:
                    pdfs_used.append(s)
        if 'systematics_arguments' in mydict:
            systematics_arguments=MadGraphSystematicsUtils.parse_systematics_arguments(mydict['systematics_arguments'])
            if 'pdf' in systematics_arguments:
                sys_pdf=systematics_arguments['pdf'].replace(',',' ').replace('@',' ').split()
                for s in sys_pdf:
                    if s.isdigit():
                        idx=int(s)
                        if idx>1000: # the sys_pdf syntax is such that small numbers are used to specify the subpdf index
                            pdfs_used.append(idx)
                    else:
                        pdfs_used.append(s)
        for pdf in pdfs_used:
            if isinstance(pdf,str) and (pdf.lower()=='errorset' or pdf.lower()=='central'):
                continue
            # new function to get both lhapdf id and name
            pdfid,pdfname=get_lhapdf_id_and_name(pdf)
            mglog.info("Found LHAPDF ID=%i, name=%s!"%(pdfid,pdfname))

            if not os.path.exists(newMGCLHA+pdfname) and not os.path.lexists(newMGCLHA+pdfname):
                if not os.path.exists(LHADATAPATH+'/'+pdfname):
                    mglog.warning('PDF not installed at '+LHADATAPATH+'/'+pdfname)
                if allow_links:
                    mglog.info('linking '+LHADATAPATH+'/'+pdfname+' --> '+newMGCLHA+pdfname)
                    os.symlink(LHADATAPATH+'/'+pdfname,newMGCLHA+pdfname)
                else:
                    mglog.info('copying '+LHADATAPATH+'/'+pdfname+' --> '+newMGCLHA+pdfname)
                    shutil.copytree(LHADATAPATH+'/'+pdfname,newMGCLHA+pdfname)

        if allow_links:
            mglog.info('linking '+LHADATAPATH+'/pdfsets.index --> '+newMGCLHA+'pdfsets.index')
            os.symlink(LHADATAPATH+'/pdfsets.index',newMGCLHA+'pdfsets.index')

            atlasLHADATAPATH=LHADATAPATH.replace('sft.cern.ch/lcg/external/lhapdfsets/current','atlas.cern.ch/repo/sw/Generators/lhapdfsets/current')
            mglog.info('linking '+atlasLHADATAPATH+'/lhapdf.conf --> '+newMGCLHA+'lhapdf.conf')
            os.symlink(atlasLHADATAPATH+'/lhapdf.conf',newMGCLHA+'lhapdf.conf')
        else:
            mglog.info('copying '+LHADATAPATH+'/pdfsets.index --> '+newMGCLHA+'pdfsets.index')
            shutil.copy2(LHADATAPATH+'/pdfsets.index',newMGCLHA+'pdfsets.index')

            atlasLHADATAPATH=LHADATAPATH.replace('sft.cern.ch/lcg/external/lhapdfsets/current','atlas.cern.ch/repo/sw/Generators/lhapdfsets/current')
            mglog.info('copying '+atlasLHADATAPATH+'/lhapdf.conf -->'+newMGCLHA+'lhapdf.conf')
            shutil.copy2(atlasLHADATAPATH+'/lhapdf.conf',newMGCLHA+'lhapdf.conf')


        LHADATAPATH=os.getcwd()+'/MGC_LHAPDF'

    else:
        mglog.info('Not using LHAPDF')
        return (LHAPATH,origLHAPATH,origLHAPDF_DATA_PATH)


    if isNLO:
        os.environ['LHAPDF_DATA_PATH']=LHADATAPATH

    mglog.info('Path to LHAPDF install dir:%s'%LHAPATH)
    mglog.info('Path to LHAPDF data dir: %s'%LHADATAPATH)
    if not os.path.isdir(LHADATAPATH):
        mglog.error('LHAPDF data dir: %s is not accesible'%LHADATAPATH)
        return 1
    if not os.path.isdir(LHAPATH):
        mglog.error('LHAPDF path dir: %s is not accesible'%LHAPATH)
        return 1

    # Dealing with LHAPDF
    if extlhapath:
        lhapdfconfig=extlhapath
        if not os.access(lhapdfconfig,os.X_OK):
            mglog.error('Failed to find valid external lhapdf-config at %s'%lhapdfconfig)
            return 1
        LHADATAPATH=subprocess.Popen([lhapdfconfig, '--datadir'],stdout = subprocess.PIPE).stdout.read().strip()
        mglog.info('Changing LHAPDF_DATA_PATH to %s'%LHADATAPATH)
        os.environ['LHAPDF_DATA_PATH']=LHADATAPATH
    else:
        getlhaconfig = subprocess.Popen(['get_files','-data','lhapdf-config'])
        getlhaconfig.wait()
        #Get custom lhapdf-config
        if not os.access(os.getcwd()+'/lhapdf-config',os.X_OK):
            mglog.error('Failed to get lhapdf-config from MadGraphControl')
            return 1
        lhapdfconfig = os.getcwd()+'/lhapdf-config'

    mglog.info('lhapdf-config --version:      %s'%str(subprocess.Popen([lhapdfconfig, '--version'],stdout = subprocess.PIPE).stdout.read().strip()))
    mglog.info('lhapdf-config --prefix:       %s'%str(subprocess.Popen([lhapdfconfig, '--prefix'],stdout = subprocess.PIPE).stdout.read().strip()))
    mglog.info('lhapdf-config --libdir:       %s'%str(subprocess.Popen([lhapdfconfig, '--libdir'],stdout = subprocess.PIPE).stdout.read().strip()))
    mglog.info('lhapdf-config --datadir:      %s'%str(subprocess.Popen([lhapdfconfig, '--datadir'],stdout = subprocess.PIPE).stdout.read().strip()))
    mglog.info('lhapdf-config --pdfsets-path: %s'%str(subprocess.Popen([lhapdfconfig, '--pdfsets-path'],stdout = subprocess.PIPE).stdout.read().strip()))

    if not isNLO:
        config_card=process_dir+'/Cards/me5_configuration.txt'
    else:
        config_card=process_dir+'/Cards/amcatnlo_configuration.txt'

    oldcard = open(config_card,'r')
    newcard = open(config_card+'.tmp','w')

    for line in oldcard:
        if 'lhapdf = ' in line:
            newcard.write('lhapdf = %s \n'%(lhapdfconfig))
            mglog.info('Setting lhapdf = %s in %s'%(lhapdfconfig,config_card))
        else:
            newcard.write(line)
    oldcard.close()
    newcard.close()
    shutil.move(config_card+'.tmp',config_card)

    mglog.info('Creating links for LHAPDF')
    if os.path.islink(process_dir+'/lib/PDFsets'):
        os.unlink(process_dir+'/lib/PDFsets')
    elif os.path.isdir(process_dir+'/lib/PDFsets'):
        shutil.rmtree(process_dir+'/lib/PDFsets')
    if allow_links:
        os.symlink(LHADATAPATH,process_dir+'/lib/PDFsets')
    else:
        shutil.copytree(LHADATAPATH,process_dir+'/lib/PDFsets')
    mglog.info('Available PDFs are:')
    mglog.info( sorted( [ x for x in os.listdir(process_dir+'/lib/PDFsets') if not ".tar.gz" in x ] ) )

    return (LHAPATH,origLHAPATH,origLHAPDF_DATA_PATH)

# Function to set the number of cores and the running mode in the run card
def setNCores(process_dir=None, Ncores=None):
    if process_dir is None:
        mglog.warning('Cannot setNCores because no process dir was provided')
        return
    my_Ncores = Ncores
    if Ncores is None and 'ATHENA_PROC_NUMBER' in os.environ:
        my_Ncores = int(os.environ['ATHENA_PROC_NUMBER'])
    if my_Ncores is None:
        mglog.info('Setting up for serial run')
        my_Ncores = 1

    if not is_NLO_run(process_dir=process_dir):
        config_card=process_dir+'/Cards/me5_configuration.txt'
    else:
        config_card=process_dir+'/Cards/amcatnlo_configuration.txt'

    import fileinput
    # core configuration with the one we need
    oldcard = open(config_card,'r')
    newcard = open(config_card+'.tmp','w')
    for line in oldcard.readlines():
        if 'nb_core = ' in line:
            mglog.info('Setting number of cores = %i in %s'%(my_Ncores,config_card))
            newcard.write('nb_core = %i \n'%(my_Ncores))
        elif 'run_mode = ' in line:
            mglog.info('Setting run mode = %i in %s'%(0 if my_Ncores==1 else 2,config_card))
            newcard.write('run_mode = %i \n'%(0 if my_Ncores==1 else 2))
        elif 'automatic_html_opening =' in line:
            mglog.info('Setting automatic_html_opening = %s'%('False'))
            newcard.write('automatic_html_opening = %s \n'%('False'))
        else:
            newcard.write(line)
    oldcard.close()
    newcard.close()
    shutil.move(config_card+'.tmp',config_card)


def resetLHAPDF(origLHAPATH='',origLHAPDF_DATA_PATH=''):
    mglog.info('Restoring original LHAPDF env variables:')
    os.environ['LHAPATH']=origLHAPATH
    os.environ['LHAPDF_DATA_PATH']=origLHAPDF_DATA_PATH
    mglog.info('LHAPATH=%s'%os.environ['LHAPATH'])
    mglog.info('LHAPDF_DATA_PATH=%s'%os.environ['LHAPDF_DATA_PATH'])


def get_mg5_executable():
    madpath=os.environ['MADPATH']
    if not os.access(madpath+'/bin/mg5_aMC',os.R_OK):
        raise RuntimeError('mg5_aMC executable not found in '+madpath)
    return madpath+'/bin/mg5_aMC'


def add_lifetimes(process_dir,threshold=None):
    """ Add lifetimes to the generated LHE file.  Should be
    called after generate_events is called.
    """
    if config_only_check(): return

    me_exec=get_mg5_executable()

    if len(glob.glob(process_dir+'/Events/*'))<1:
        mglog.error('Process dir %s does not contain events?'%process_dir)
    run = glob.glob(process_dir+'/Events/*')[0].split('/')[-1]

    # Note : This slightly clunky implementation is needed for the time being
    # See : https://answers.launchpad.net/mg5amcnlo/+question/267904

    tof_c = open('time_of_flight_exec_card','w')
    tof_c.write('launch '+process_dir+''' -i
add_time_of_flight '''+run+((' --threshold='+str(threshold)) if threshold is not None else ''))
    tof_c.close()

    mglog.info('Started adding time of flight info '+str(time.asctime()))

    generate = subprocess.Popen([me_exec,'time_of_flight_exec_card'],stdin=subprocess.PIPE)
    generate.communicate()

    mglog.info('Finished adding time of flight information at '+str(time.asctime()))

    # Re-zip the file if needed
    lhe_gz = glob.glob(process_dir+'/Events/*/*lhe.gz')[0]
    if not os.access(lhe_gz,os.R_OK):
        mglog.info('LHE file needs to be zipped')
        lhe = glob.glob(process_dir+'/Events/*/*lhe.gz')[0]
        rezip = subprocess.Popen(['gzip',lhe])
        mglog.info('Zipped')
    else:
        mglog.info('LHE file zipped by MadGraph automatically. Nothing to do')

    return True


def add_madspin(madspin_card=None,process_dir=MADGRAPH_GRIDPACK_LOCATION):
    """ Run madspin on the generated LHE file.  Should be
    run when you have inputGeneratorFile set.
    Only requires a simplified process with the same model that you are
    interested in (needed to set up a process directory for MG5_aMC)
    """
    if config_only_check(): return

    me_exec=get_mg5_executable()

    if madspin_card is not None:
        shutil.copyfile(madspin_card,process_dir+'/Cards/madspin_card.dat')

    if len(glob.glob(process_dir+'/Events/*'))<1:
        mglog.error('Process dir %s does not contain events?'%process_dir)
    run = glob.glob(process_dir+'/Events/*')[0].split('/')[-1]

    # Note : This slightly clunky implementation is needed for the time being
    # See : https://answers.launchpad.net/mg5amcnlo/+question/267904

    ms_c = open('madspin_exec_card','w')
    ms_c.write('launch '+process_dir+''' -i
decay_events '''+run)
    ms_c.close()

    mglog.info('Started running madspin at '+str(time.asctime()))

    generate = subprocess.Popen([me_exec,'madspin_exec_card'],stdin=subprocess.PIPE)
    generate.communicate()

    mglog.info('Finished running madspin at '+str(time.asctime()))

    # Re-zip the file if needed
    lhe_gz = glob.glob(process_dir+'/Events/*/*lhe.gz')[0]
    if not os.access(lhe_gz,os.R_OK):
        mglog.info('LHE file needs to be zipped')
        lhe = glob.glob(process_dir+'/Events/*/*lhe.gz')[0]
        rezip = subprocess.Popen(['gzip',lhe])
        mglog.info('Zipped')
    else:
        mglog.info('LHE file zipped by MadGraph automatically. Nothing to do')


def arrange_output(process_dir='PROC_mssm_0',lhe_version=None,saveProcDir=False,runArgs=None,madspin_card=None,fixEventWeightsForBridgeMode=False):
    if config_only_check(): return

    # NLO is not *really* the question here, we need to know if we should look for weighted or
    #  unweighted events in the output directory.  MadSpin (above) only seems to give weighted
    #  results for now?
    #isNLO=is_NLO_run(process_dir=process_dir)
    hasUnweighted = os.access(process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe.gz',os.R_OK)

    hasRunMadSpin=False
    madspinDirs=sorted(glob.glob(process_dir+'/Events/'+MADGRAPH_RUN_NAME+'_decayed_*/'))
    if len(madspinDirs): hasRunMadSpin=True
    if hasRunMadSpin and not hasUnweighted:
        # check again:
        hasUnweighted = os.access(madspinDirs[-1]+'/unweighted_events.lhe.gz',os.R_OK)

    if madspin_card or hasRunMadSpin:
        if len(madspinDirs):
            if hasUnweighted:
                # so this is a bit of a mess now...
                # if madspin is run from an NLO grid pack the correct lhe events are at both
                #      madevent/Events/run_01/unweighted_events.lhe.gz
                # and  madevent/Events/run_01_decayed_1/events.lhe.gz
                # so there are unweighted events but not in the madspinDir...
                if os.path.exists(madspinDirs[-1]+'/unweighted_events.lhe.gz'):
                    shutil.move(madspinDirs[-1]+'/unweighted_events.lhe.gz',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe.gz')
                    mglog.info('Moving MadSpin events from %s to %s.'%(madspinDirs[-1]+'/unweighted_events.lhe.gz',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe.gz'))
                elif os.path.exists(madspinDirs[-1]+'/events.lhe.gz'):
                    shutil.move(madspinDirs[-1]+'/events.lhe.gz',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe.gz')
                    mglog.info('Moving MadSpin events from %s to %s.'%(madspinDirs[-1]+'/events.lhe.gz',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe.gz'))
                else:
                    raise RuntimeError('MadSpin was run but can\'t find files :(')

            else:
                shutil.move(madspinDirs[-1]+'/events.lhe.gz',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/events.lhe.gz')
                mglog.info('Moving MadSpin events from %s to %s.'%(madspinDirs[-1]+'/events.lhe.gz',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/events.lhe.gz'))

        else:
            mglog.error('MadSpin was run but can\'t find output folder %s.'%(process_dir+'/Events/'+MADGRAPH_RUN_NAME+'_decayed_1/'))
            raise RuntimeError('MadSpin was run but can\'t find output folder %s.'%(process_dir+'/Events/'+MADGRAPH_RUN_NAME+'_decayed_1/'))

        if fixEventWeightsForBridgeMode:
            mglog.info("Fixing event weights after MadSpin... initial checks.")

            # get the cross section from the undecayed LHE file
            spinmodenone=False
            MGnumevents=-1
            MGintweight=-1

            if hasUnweighted:
                eventsfilename="unweighted_events"
            else:
                eventsfilename="events"
            unzip = subprocess.Popen(['gunzip','-f',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/%s.lhe.gz' % eventsfilename])
            unzip.wait()

            for line in open(process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/%s.lhe'%eventsfilename):
                if "Number of Events" in line:
                    sline=line.split()
                    MGnumevents=int(sline[-1])
                elif "Integrated weight (pb)" in line:
                    sline=line.split()
                    MGintweight=float(sline[-1])
                elif "set spinmode none" in line:
                    spinmodenone=True
                elif "</header>" in line:
                    break

            if spinmodenone and MGnumevents>0 and MGintweight>0:
                mglog.info("Fixing event weights after MadSpin... modifying LHE file.")
                newlhe=open(process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/%s_fixXS.lhe'%eventsfilename,'w')
                initlinecount=0
                eventlinecount=0
                inInit=False
                inEvent=False

                # new default for MG 2.6.1+ (https://its.cern.ch/jira/browse/AGENE-1725)
                # but verified from LHE below.
                event_norm_setting="average" 

                for line in open(process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/%s.lhe'%eventsfilename):

                    newline=line
                    if "<init>" in line:                         
                        inInit=True
                        initlinecount=0
                    elif "</init>" in line:
                        inInit=False
                    elif inInit and initlinecount==0:
                        initlinecount=1
                        # check event_norm setting in LHE file, deteremines how Pythia interprets event weights
                        sline=line.split()
                        if abs(int(sline[-2])) == 3:
                            event_norm_setting="sum"
                        elif abs(int(sline[-2])) == 4:
                            event_norm_setting="average"
                    elif inInit and initlinecount==1:
                        sline=line.split()
                        # update the global XS info
                        relunc=float(sline[1])/float(sline[0])
                        sline[0]=str(MGintweight)                
                        sline[1]=str(float(sline[0])*relunc)     
                        if event_norm_setting=="sum":
                            sline[2]=str(MGintweight/MGnumevents)
                        elif event_norm_setting=="average":
                            sline[2]=str(MGintweight)            
                        newline=' '.join(sline)
                        newline+="\n"
                        initlinecount+=1
                    elif inInit and initlinecount>1:
                        initlinecount+=1
                    elif "<event>" in line:                      
                        inEvent=True
                        eventlinecount=0
                    elif "</event>" in line:
                        inEvent=False
                    elif inEvent and eventlinecount==0:
                        sline=line.split()
                        # next change the per-event weights
                        if event_norm_setting=="sum":
                            sline[2]=str(MGintweight/MGnumevents)
                        elif event_norm_setting=="average":
                            sline[2]=str(MGintweight)            
                        newline=' '.join(sline)
                        newline+="\n"
                        eventlinecount+=1
                    newlhe.write(newline)
                newlhe.close()

                mglog.info("Fixing event weights after MadSpin... cleaning up.")
                shutil.copyfile(process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/%s.lhe' % eventsfilename,
                                process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/%s_badXS.lhe' % eventsfilename)

                shutil.move(process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/%s_fixXS.lhe' % eventsfilename,
                            process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/%s.lhe' % eventsfilename)

                rezip = subprocess.Popen(['gzip',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/%s.lhe' % eventsfilename])
                rezip.wait()

                rezip = subprocess.Popen(['gzip',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/%s_badXS.lhe' % eventsfilename])
                rezip.wait()

    # Clean up in case a link or file was already there
    if os.path.exists(os.getcwd()+'/events.lhe'): os.remove(os.getcwd()+'/events.lhe')

    mglog.info('Unzipping generated events.')
    if hasUnweighted:
        unzip = subprocess.Popen(['gunzip','-f',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe.gz'])
        unzip.wait()
    else:
        unzip = subprocess.Popen(['gunzip','-f',process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/events.lhe.gz'])
        unzip.wait()

    mglog.info('Putting a copy in place for the transform.')
    if hasUnweighted:
        orig_input = process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/unweighted_events.lhe'
        mod_output = open(os.getcwd()+'/events.lhe','w')
    else:
        orig_input = process_dir+'/Events/'+MADGRAPH_RUN_NAME+'/events.lhe'
        mod_output = open(os.getcwd()+'/events.lhe','w')

    #Removing empty lines in LHE
    nEmpty=0
    with open(orig_input,'r') as fileobject:
        for line in fileobject:
            if line.strip():
                mod_output.write(line)
            else:
                nEmpty=nEmpty+1
    mod_output.close()

    mglog.info('Removed %i empty lines from LHEF.'%nEmpty)

    if lhe_version:
        mod_output2 = open(os.getcwd()+'/events.lhe','r')
        test=mod_output2.readline()
        if 'version="' in test:
            mglog.info('Applying LHE version hack')
            final_file = open(os.getcwd()+'/events.lhe.copy','w')
            final_file.write('<LesHouchesEvents version="%i.0">\n'%lhe_version)
            shutil.copyfileobj(mod_output2, final_file)
            final_file.close()
            shutil.copy(os.getcwd()+'/events.lhe.copy',os.getcwd()+'/events.lhe')
        mod_output2.close()

    # Actually move over the dataset - this first part is horrible...
    outputTXTFile = None
    if runArgs is None:
        raise RuntimeError('Must provide runArgs to arrange_output')

    if hasattr(runArgs,'outputTXTFile'):
        outputDS = runArgs.outputTXTFile
    else:
        outputDS = 'tmp_LHE_events'

    mglog.info('Moving file over to '+outputDS.split('.tar.gz')[0]+'.events')

    shutil.move(os.getcwd()+'/events.lhe',outputDS.split('.tar.gz')[0]+'.events')

    mglog.info('Re-zipping into dataset name '+outputDS)
    rezip = subprocess.Popen(['tar','cvzf',outputDS,outputDS.split('.tar.gz')[0]+'.events'])
    rezip.wait()

    if not saveProcDir:
        mglog.info('Blasting away the process directory')
        shutil.rmtree(process_dir,ignore_errors=True)

        if os.path.isdir('MGC_LHAPDF/'):
            shutil.rmtree('MGC_LHAPDF/',ignore_errors=True)

    # shortening the outputDS in the case of an output TXT file
    if outputTXTFile is not None:
        outputDS = outputDS.split('.TXT')[0]
    # Do some fixing up for them
    if runArgs is not None:
        mglog.debug('Setting inputGenerator file to '+outputDS)
        runArgs.inputGeneratorFile=outputDS

    mglog.info('All done with output arranging!')
    return outputDS


def setup_bias_module(bias_module,run_card,process_dir):
    if isinstance(bias_module,tuple):
        mglog.info('Using bias module '+bias_module[0])
        the_run_card = open(run_card,'r')
        for line in the_run_card:
            if 'bias_module' in line and not bias_module[0] in line:
                mglog.error('You need to add the bias module '+bias_module[0]+' to the run card to actually run it')
                return 1
        the_run_card.close()
        if len(bias_module)!=3:
            mglog.error('Please give a 3-tuple of strings containing bias module name, bias module, and makefile. Alternatively, give path to bias module tarball.')
            return 1
        bias_module_newpath=process_dir+'/Source/BIAS/'+bias_module[0]
        os.makedirs(bias_module_newpath)
        bias_module_file=open(bias_module_newpath+'/'+bias_module[0]+'.f','w')
        bias_module_file.write(bias_module[1])
        bias_module_file.close()
        bias_module_make_file=open(bias_module_newpath+'/Makefile','w')
        bias_module_make_file.write(bias_module[2])
        bias_module_make_file.close()
    else:
        mglog.info('Using bias module '+bias_module)
        bias_module_name=bias_module.split('/')[-1].replace('.gz','')
        bias_module_name=bias_module_name.replace('.tar','')
        the_run_card = open(run_card,'r')
        for line in the_run_card:
            if 'bias_module' in line and not bias_module_name in line:
                mglog.error('You need to add the bias module '+bias_module_name+' to the run card to actually run it')
                return 1
        the_run_card.close()

        if os.path.exists(bias_module+'.tar.gz'):
            bias_module_path=bias_module+'.tar.gz'
        elif os.path.exists(bias_module+'.gz'):
            bias_module_path=bias_module+'.gz'
        elif os.path.exists(bias_module):
            bias_module_path=bias_module
        else:
            mglog.error('Did not find bias module '+bias_module+' , this path should point to folder or tarball.  Alternatively give a tuple of strings containing module name, module, and makefile')
            return 1
        bias_module_newpath=process_dir+'/Source/BIAS/'+bias_module_path.split('/')[-1]
        mglog.info('Copying bias module into place: '+bias_module_newpath)
        shutil.copy(bias_module_path,bias_module_newpath)
        mglog.info('Unpacking bias module')
        if bias_module_newpath.endswith('.tar.gz'):
            untar = subprocess.Popen(['tar','xvzf',bias_module_newpath,'--directory='+process_dir+'/Source/BIAS/'])
            untar.wait()
        elif bias_module_path.endswith('.gz'):
            gunzip = subprocess.Popen(['gunzip',bias_module_newpath])
            gunzip.wait()


def check_reweight_card(reweight_card):
    shutil.move(reweight_card,reweight_card+'.old')
    oldcard = open(reweight_card+'.old','r')
    newcard = open(reweight_card,'w')
    changed = False
    nrwgt = 1
    for line in oldcard:
        if not line.strip().startswith('launch') :
            newcard.write(line)
        else:
            # if both rwgt_info and name are given, fine
            if '--rwgt_info' in line and '--rwgt_name' in line:
                newcard.write(line)
            # if only one of the two is defined, set the other to the same value
            elif '--rwgt_info' in line:
                m=re.match('launch\s*--rwgt_info\s*=\s*([\s\S]+)',line.strip())
                if m==None or len(m.groups())!=1:
                    raise RuntimeError('Unexpected format of reweight card')
                else:
                    newcard.write(line.strip()+' --rwgt_name='+m.group(1).strip()+'\n')
                    changed=True
            elif '--rwgt_name' in line:
                m=re.match('launch\s*--rwgt_name\s*=\s*([\s\S]+)',line.strip())
                if m==None or len(m.groups())!=1:
                    raise RuntimeError('Unexpected format of reweight card')
                else:
                    newcard.write(line.strip()+' --rwgt_info='+m.group(1).strip()+'\n')
                    changed=True
            # if none is defined: complain
            else:
                mglog.warning('Every reweighting launch needs a --rwgt_name (see https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/Reweight), please update your reweight_card accordingly. Added a dummy name for now.')
                newcard.write(line.strip()+' --rwgt_name=dummy_rwgt_name_{0}  --rwgt_info=dummy_rwgt_info_{0}\n'.format(nrwgt))
                nrwgt+=1

    if changed:
        mglog.info('Updated reweight_card')
    newcard.close()
    oldcard.close()


def helpful_definitions():
    return """
# Define multiparticle labels
define p = g u c d s u~ c~ d~ s~
define j = g u c d s u~ c~ d~ s~
define pb = g u c d s b u~ c~ d~ s~ b~
define jb = g u c d s b u~ c~ d~ s~ b~
define l+ = e+ mu+
define l- = e- mu-
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
define fu = u c e+ mu+ ta+
define fu~ = u~ c~ e- mu- ta-
define fd = d s ve~ vm~ vt~
define fd~ = d~ s~ ve vm vt
define susystrong = go ul ur dl dr cl cr sl sr t1 t2 b1 b2 ul~ ur~ dl~ dr~ cl~ cr~ sl~ sr~ t1~ t2~ b1~ b2~
define susyweak = el- el+ er- er+ mul- mul+ mur- mur+ ta1- ta1+ ta2- ta2+ n1 n2 n3 n4 x1- x1+ x2- x2+ sve sve~ svm svm~ svt svt~
define susylq = ul ur dl dr cl cr sl sr
define susylq~ = ul~ ur~ dl~ dr~ cl~ cr~ sl~ sr~
define susysq = ul ur dl dr cl cr sl sr t1 t2 b1 b2
define susysq~ = ul~ ur~ dl~ dr~ cl~ cr~ sl~ sr~ t1~ t2~ b1~ b2~
define susysl = el- el+ er- er+ mul- mul+ mur- mur+ ta1- ta1+ ta2- ta2+
define susyv = sve svm svt
define susyv~ = sve~ svm~ svt~
"""

def strong_process_dict( njets = 1 , gentype = 'GG' ):
    header = "import model MSSM_SLHA2\n"
    header += helpful_definitions()
    header += """
# Specify process(es) to run

"""
    footer = """
# Output processes to MadEvent directory
output -f
"""

    jetloop = [ '' ]
    if 0<njets: jetloop += [ 'j' ]
    if 1<njets: jetloop += [ 'j j' ]
    if 2<njets: jetloop += [ 'j j j' ]
    if 3<njets: jetloop += [ 'j j j j' ]
    if 4<njets: jetloop += [ 'j j j j j' ]
    if 5<njets: jetloop += [ 'j j j j j j' ]
    if 6<njets: jetloop += [ 'j j j j j j j' ]
    if 7<njets: jetloop += [ 'j j j j j j j j' ]
    if 8<njets:
        mglog.error('That is just ridiculous.  Do not generate more than eight jets off of your main process in MadGraph.  The world will implode.')
        return {}

    starter = 'generate'
    squarks = ['ul','dl','cl','sl','ur','dr','cr','sr']
    ss_string = ''
    sg_string = ''
    ssb_string = ''
    gg_string = ''
    n1n1_string = ''
    n2n3_string = ''
    c1n2_string = ''
    cc_string = ''
    bb_string = ''
    tt_string = ''
    Scharm_string = ''
    first = True
    for jets in jetloop:
        bb_string += '%s p p > b1 b1~ %s $ go ul ur dl dr cl cr sl sr t1 t2 b2 ul~ ur~ dl~ dr~ cl~ cr~ sl~ sr~ t1~ t2~ b2~ \n'%(starter,jets)
        tt_string += '%s p p > t1 t1~ %s $ go ul ur dl dr cl cr sl sr t2 b1 b2 ul~ ur~ dl~ dr~ cl~ cr~ sl~ sr~ b1~ t2~ b2~ \n'%(starter,jets)
        n1n1_string += '%s p p > n1 n1 %s $ susystrong \n'%(starter,jets)
        n2n3_string += '%s p p > n2 n3 %s $ susystrong \n'%(starter,jets)
        gg_string += '%s p p > go go %s $ susysq \n'%(starter,jets)
        for sign1 in ['+','-']:
            c1n2_string += '%s p p > x1%s n2 %s $ susystrong \n'%(starter,sign1,jets)
            for sign2 in ['+','-']:
                cc_string += '%s p p > x1%s x1%s %s $ susystrong \n'%(starter,sign1,sign2,jets)
                starter = 'add process'
        if first: starter = 'generate'
        for s in ['cl','cr']:
            for t in ['cl','cr']:
                Scharm_string += '%s p p > %s %s~ %s $ go \n'%(starter,s,t,jets)
                starter = 'add process'
        if first: starter = 'generate'
        for s in squarks:
            sg_string += '%s p p > %s go %s \n'%(starter,s,jets)
            sg_string += '%s p p > %s~ go %s \n'%(starter,s,jets)
            for t in squarks:
                if 'phot' in gentype:
                    ss_string += '%s p p > %s %s~ a %s $ go \n'%(starter,s,t,jets)
                    starter = 'add process'
                else:
                    ss_string += '%s p p > %s %s~ %s $ go \n'%(starter,s,t,jets)
                    starter = 'add process'
                ssb_string += '%s p p > %s %s~ %s \n'%(starter,s,t,jets)
        first = False

    processes = {
        'SS' : ss_string ,
        'SSphot' : ss_string,
        'GG' : gg_string ,
        'SG' : sg_string ,
        'GS' : sg_string ,
        'CC' : cc_string ,
        'BB' : bb_string ,
        'Scharm' : Scharm_string,
        'TT' : tt_string ,
        'N2N3' : n2n3_string ,
        'N1N1' : n1n1_string ,
        'C1N2' : c1n2_string ,
        'SSB' : ssb_string ,
        'ALL' : ss_string+'\n'+gg_string+'\n'+ssb_string+'\n'+sg_string
      }

    if not gentype in processes.keys():
        mglog.error('No idea how to deal with the simplified models for'+gentype+'.  Sorry!')
        return ''
    if processes[gentype] == '':
        mglog.error('No processes found for the set up you requested: '+str(gentype))
        return ''
    mglog.info('Will run MadGraph over the process:')
    mglog.info(str(processes[gentype]))

    return header+processes[gentype]+footer


def get_variations( gentype , masses , syst_mod , xqcut = None ):
    if xqcut is None:
        my_mass=2000 # default
        if 'Scharm'==gentype:
            my_mass = min([masses[x] for x in ['1000004','2000004'] if x in masses ])
        elif 'N2N3'==gentype:
            my_mass = min([masses[x] for x in ['1000023','1000025'] if x in masses ])
        elif 'N1N1'==gentype:
            my_mass = min([masses[x] for x in ['1000022'] if x in masses ])
        elif 'C1N2'==gentype:
            my_mass = min([masses[x] for x in ['1000023','1000024'] if x in masses ])
        elif 'Stau'==gentype:
            my_mass = min([masses[x] for x in ['1000015','2000015'] if x in masses ])
        elif 'StauStau'==gentype:
            my_mass = min([masses[x] for x in ['1000015','2000015'] if x in masses ])
        elif 'SlepSlep'==gentype:
            my_mass = min([masses[x] for x in ['1000011','1000013','1000015','2000011','2000013','2000015'] if x in masses ])
        elif 'T2' in gentype:
            my_mass = min([masses[x] for x in ['2000006'] if x in masses ])
        else:
            if 'G' in gentype or 'ALL' in gentype:
                my_mass = masses['1000021']
            if 'S' in gentype or 'ALL' in gentype:
                my_mass = min([masses[x] for x in ['1000001','1000002','1000003','1000004','2000001','2000002','2000003','2000004'] if x in masses])
            if 'T' in gentype:
                my_mass = masses['1000006']
            if 'B' in gentype:
                my_mass = masses['1000005']
            if 'C' in gentype:
                my_mass = masses['1000024']
            if 'D' in gentype:
                my_mass = masses['2000001']
        xqcut = min(my_mass*0.25,500)
        if syst_mod is not None and 'qup' in syst_mod.lower(): xqcut = xqcut*2.
        elif syst_mod is not None and 'qdown' in syst_mod.lower(): xqcut = xqcut*0.5
    mglog.info('For matching, will use xqcut of '+str(xqcut))

    alpsfact = 1.0
    scalefact = 1.0
    if syst_mod is not None and 'alpsfactup' in syst_mod.lower():
        alpsfact = 2.0
    elif syst_mod is not None and 'alpsfactdown' in syst_mod.lower():
        alpsfact = 0.5

    if syst_mod is not None and 'scalefactup' in syst_mod.lower(): scalefact = 2.0
    elif syst_mod is not None and 'scalefactdown' in syst_mod.lower(): scalefact = 0.5

    return abs(xqcut) , alpsfact , scalefact


def SUSY_process(process=''):
    # Generate the new process!
    if 'import model' in process:
        mglog.info('Assuming that you have specified the model in your process string already')
        full_proc = ''
        for l in process.split('\n'):
            if 'import model' in l:
                full_proc += l+'\n'
                break
        full_proc+=helpful_definitions()
        for l in process.split('\n'):
            if 'import model' not in l:
                full_proc += l+'\n'
        full_proc+="""
# Output processes to MadEvent directory
output -f
"""
    else:
        full_proc = "import model MSSM_SLHA2\n"+helpful_definitions()+"""
# Specify process(es) to run

"""+process+"""
# Output processes to MadEvent directory
output -f
"""
    return full_proc


def SUSY_Generation(runArgs = None, process=None, gentype='SS', njets=1, decaytype='direct',\
                    nevents=10000, syst_mod=None, xqcut=None, keepOutput=False, param_card=None, writeGridpack=False,\
                    madspin_card=None, extras={}, params={}, fixEventWeightsForBridgeMode=False):

    xqcut , alpsfact , scalefact = get_variations( gentype , params['MASS'] , syst_mod , xqcut=xqcut )

    process_dir = MADGRAPH_GRIDPACK_LOCATION
    if not is_gen_from_gridpack():
        if process is not None: full_proc = SUSY_process(process)
        else: full_proc = strong_process_dict(njets,gentype)
        process_dir = new_process(full_proc)
    mglog.info('Using process directory '+str(process_dir))

    # Grab the param card and move the new masses into place
    modify_param_card(param_card_input=param_card,process_dir=process_dir,params=params)

    # Set up the extras dictionary
    settings = {'drjj':0.0,'lhe_version':'3.0','cut_decays':'F','pdflabel':"'cteq6l1'",'lhaid':10042,'ktdurham':xqcut,'nevents':nevents}
    settings.update(extras)
    settings.update({'scalefact':scalefact,'alpsfact':alpsfact})

    # Generate events!
    if is_gen_from_gridpack():
        modify_run_card(process_dir=process_dir,runArgs=runArgs,settings=settings)
        generate_from_gridpack(runArgs=runArgs)
    else:
        # Grab the run card and move it into place
        modify_run_card(process_dir=process_dir,runArgs=runArgs,settings=settings)
        generate(process_dir=process_dir,grid_pack=writeGridpack)

    # Move output files into the appropriate place, with the appropriate name
    the_spot = arrange_output(process_dir=process_dir,saveProcDir=keepOutput,runArgs=runArgs,fixEventWeightsForBridgeMode=fixEventWeightsForBridgeMode)

    mglog.info('All done generating events!!')
    return [xqcut,the_spot]


def update_lhe_file(lhe_file_old,param_card_old=None,lhe_file_new=None,masses={},delete_old_lhe=True):
    """Build a new LHE file from an old one and an updated param card.
    The masses of some particles can be changed via the masses dictionary.  No particles that appear in the events
    may have their masses changed.
    If the param card is provided, the decay block in the LHE file will be replaced with the one in the param card.
    By default, the old LHE file is removed.
    If None is provided as a new LHE file name, the new file will replace the old one."""
    # If we want to just use a temp file, then put in a little temp holder
    lhe_file_new_tmp = lhe_file_new if lhe_file_new is not None else lhe_file_old+'.tmp'
    # Make sure the LHE file is there
    if not os.access(lhe_file_old,os.R_OK):
        mglog.error('Could not access old LHE file at '+str(lhe_file_old)+'. Please check the file location.')
        return -1
    # Grab the old param card
    paramcard = subprocess.Popen(['get_files','-data',param_card_old])
    paramcard.wait()
    if not os.access(param_card_old,os.R_OK):
        mglog.info('Could not get param card '+param_card_old)
    # Don't overwrite old param cards
    if os.access(lhe_file_new_tmp,os.R_OK):
        mglog.error('Old file at'+str(lhe_file_new_tmp)+' in the current directory. Dont want to clobber it. Please move it first.')
        return -1

    oldparam = open(param_card_old,'r')
    newlhe = open(lhe_file_new_tmp,'w')
    blockName = None
    decayEdit = False
    eventRead = False
    particles_in_events = []
    # Decay block ends with </slha>

    with open(lhe_file_old,'r') as fileobject:
        for line in fileobject:
            if decayEdit and not '</slha>' in line: continue
            if decayEdit and '</slha>' in line:
                decayEdit = False
            if line.strip().upper().startswith('BLOCK') or line.strip().upper().startswith('DECAY')\
                        and len(line.strip().split()) > 1:
                pos = 0 if line.strip().startswith('DECAY') else 1
                blockName = line.strip().upper().split()[pos]

            akey = None
            if blockName != 'DECAY' and len(line.strip().split()) > 0:
                akey = line.strip().split()[0]
            elif blockName == 'DECAY' and len(line.strip().split()) > 1:
                akey = line.strip().split()[1]

            # Replace the masses with those in the dictionary
            if akey != None and blockName == 'MASS'  and akey in masses:
                newlhe.write('   %s    %s  # \n'%(akey,str(masses[akey])))
                mglog.info('   %s    %s  #'%(akey,str(masses[akey])))
                decayEdit = False
                continue

            # Replace the entire decay section of the LHE file with the one from the param card
            if blockName == 'DECAY' and param_card_old is not None:
                # We are now reading the decay blocks!  Take them from the param card
                oldparam = open(param_card_old,'r')
                newDecays = False
                for old_line in oldparam.readlines():
                    newBlockName = None
                    if old_line.strip().upper().startswith('DECAY') and len(old_line.strip().split()) > 1:
                        newBlockName = line.strip().upper().split()[pos]
                    if newDecays:
                        newlhe.write(old_line)
                    elif newBlockName == 'DECAY':
                        newDecays = True
                        newlhe.write(old_line)
                oldparam.close()
                # Done adding the decays
                decayEdit = True
                blockName = None
                continue

            # Keep a record of the particles that are in the events
            if not eventRead and '<event>' in line: eventRead = True
            if eventRead:
                if len(line.split())==11:
                    aparticle = line.split()[0]
                    if not aparticle in particles_in_events: particles_in_events += [aparticle]

            # Otherwise write the line again
            newlhe.write(line)

    # Check that none of the particles that we were setting the masses of appear in the LHE events
    for akey in masses:
        if akey in particles_in_events:
            mglog.error('Attempted to change mass of a particle that was in an LHE event!  This is not allowed!')
            return -1

    # Close up and return
    newlhe.close()

    # Move the new file to the old file location
    if lhe_file_new is None:
        os.remove(lhe_file_old)
        shutil.move(lhe_file_new_tmp,lhe_file_old)
        lhe_file_new_tmp = lhe_file_old
    # Delete the old file if requested
    elif delete_old_lhe:
        os.remove(lhe_file_old)

    return lhe_file_new_tmp


def modify_param_card(param_card_input=None,param_card_backup=None,process_dir=MADGRAPH_GRIDPACK_LOCATION,params={}):
    """Build a new param_card.dat from an existing one.
    Params should be a dictionary of dictionaries. The first key is the block name, and the second in the param name.
    Keys can include MASS (for masses) and DECAY X (for decays of particle X)"""
    # Grab the old param card and move it into place

    # Check for the default run card location
    if param_card_input is None:
        param_card_input=process_dir+'/Cards/param_card.dat'
    elif param_card_input is not None and not os.access(param_card_input,os.R_OK):
        paramcard = subprocess.Popen(['get_files','-data',param_card_input])
        paramcard.wait()
        if not os.access(param_card_input,os.R_OK):
            raise RuntimeError('Could not get param card '+param_card_input)
        mglog.info('Using input param card at '+param_card_input)

    #ensure all blocknames and paramnames are upper case
    for blockName in params:
       params[blockName.upper()] = params.pop(blockName)
       for paramName in params[blockName.upper()]:
          params[blockName.upper()][paramName.upper()] = params[blockName.upper()].pop(paramName)

    if param_card_backup is not None:
        mglog.info('Keeping backup of original param card at '+param_card_backup)
        param_card_old = param_card_backup
    else:
        param_card_old = param_card_input+'.old_to_be_deleted'
    if os.path.isfile(param_card_old): os.unlink(param_card_old) # delete old backup
    os.rename(param_card_input, param_card_old) # change name of original card

    oldcard = open(param_card_old,'r')
    newcard = open(process_dir+'/Cards/param_card.dat','w')
    decayEdit = False #only becomes true in a DECAY block when specifying the BR
    blockName = ""
    doneParams = {} #tracks which params have been done
    for line in oldcard:
        if line.strip().upper().startswith('BLOCK') or line.strip().upper().startswith('DECAY')\
                    and len(line.strip().split()) > 1:
            if decayEdit and blockName == 'DECAY': decayEdit = False # Start a new DECAY block
            pos = 0 if line.strip().startswith('DECAY') else 1
            blockName = line.strip().upper().split()[pos]
        if decayEdit: continue #skipping these lines because we are in an edit of the DECAY BR

        akey = None
        if blockName != 'DECAY' and len(line.strip().split()) > 0:
            akey = line.strip().split()[0]
        elif blockName == 'DECAY' and len(line.strip().split()) > 1:
            akey = line.strip().split()[1]
        if akey==None:
           newcard.write(line)
           continue

        #check if we have params for this block
        if not params.has_key(blockName):
           newcard.write(line)
           continue
        blockParams = params[blockName]

        # look for a string key, which would follow a #
        stringkey = None
        if '#' in line: #ignores comment lines
           stringkey = line.strip()[line.strip().find('#')+1:].strip()
           if len(stringkey.split()) > 0: stringkey = stringkey.split()[0].upper()

        if not akey in blockParams and not (stringkey != None and stringkey in blockParams):
           newcard.write(line)
           continue

        if akey in blockParams and (stringkey != None and stringkey in blockParams):
           raise RuntimeError('Conflicting use of numeric and string keys %s and %s' % (akey,stringkey))
        theParam = blockParams.get(akey,blockParams[stringkey])
        if not blockName in doneParams: doneParams[blockName] = {}
        if akey in blockParams: doneParams[blockName][akey]=True
        elif stringkey != None and stringkey in blockParams: doneParams[blockName][stringkey]=True

        #do special case of DECAY block
        if blockName=="DECAY":
           if theParam.splitlines()[0].split()[0]=="DECAY":
               #specifying the full decay block
               for newline in theParam.splitlines():
                    newcard.write(newline+'\n')
                    mglog.info(newline)
               decayEdit = True
           else: #just updating the total width
              newcard.write('DECAY   %s    %s  # %s\n'%(akey,str(theParam),line.strip()[line.strip().find('#')+1:] if line.strip().find('#')>0 else ""))
              mglog.info('DECAY   %s    %s  # %s\n'%(akey,str(theParam),line.strip()[line.strip().find('#')+1:] if line.strip().find('#')>0 else ""))
        # second special case of QNUMBERS
        elif blockName=='QNUMBERS':
           #specifying the full QNUMBERS block
           for newline in theParam.splitlines():
                newcard.write(newline+'\n')
                mglog.info(newline)
           decayEdit = True
        else: #just updating the parameter
           newcard.write('   %s    %s  # %s\n'%(akey,str(theParam),line.strip()[line.strip().find('#')+1:] if line.strip().find('#')>0 else ""))
           mglog.info('   %s    %s  # %s\n'%(akey,str(theParam),line.strip()[line.strip().find('#')+1:] if line.strip().find('#')>0 else ""))
        # Done editing the line!

    #check that all specified parameters have been updated (helps to catch typos)
    for blockName in params:
       if not blockName in doneParams:
          raise RuntimeError('Did not find any of the parameters for block %s in param_card' % blockName)
       for paramName in params[blockName]:
          if not paramName in doneParams[blockName]:
            raise RuntimeError('Was not able to replace parameter %s in param_card' % paramName)

    # Close up and return
    oldcard.close()
    newcard.close()
    return param_card_new


def modify_run_card(run_card_input=None,run_card_backup=None,process_dir=MADGRAPH_GRIDPACK_LOCATION,runArgs=None,settings={}):
    """Build a new run_card.dat from an existing one.
    This function can get a fresh runcard from DATAPATH or start from the process directory.
    Settings is a dictionary of keys (no spaces needed) and values to replace.
    """
    # Check for the default run card location
    if run_card_input is None:
        run_card_input=get_default_runcard(process_dir)
    elif run_card_input is not None and not os.access(run_card_input,os.R_OK):
        runcard = subprocess.Popen(['get_files','-data',run_card_input])
        runcard.wait()
        if not os.access(run_card_input,os.R_OK):
            raise RuntimeError('Could not get run card '+run_card_input)

    # guess NLO
    isNLO=is_NLO_run(process_dir=process_dir)
    # add gobal PDF and scale uncertainty config to extras, except PDF or weights for syscal config are explictly set
    MadGraphSystematicsUtils.setup_pdf_and_systematic_weights(MADGRAPH_PDFSETTING,settings,isNLO)

    # Get some info out of the runArgs
    if runArgs is not None:
        beamEnergy,rand_seed = get_runArgs_info(runArgs)
        if not 'iseed' in settings: settings['iseed']=rand_seed
        if not isNLO and not 'python_seed' in settings: settings['python_seed']=rand_seed
        if 'beamEnergy' in settings:
            mglog.warning('Do not set beamEnergy in MG settings. The variables are ebeam1 and ebeam2. Will use your setting of '+str(settings['beamEnergy']))
            beamEnergy=settings['beamEnergy']
            settings.pop('beamEnergy')
        if not 'ebeam1' in settings: settings['ebeam1']=beamEnergy
        if not 'ebeam2' in settings: settings['ebeam2']=beamEnergy

    mglog.info('Modifying run card located at '+run_card_input)
    if run_card_backup is not None:
        mglog.info('Keeping backup of original run card at '+run_card_backup)
        run_card_old = run_card_backup
    else:
        run_card_old = run_card_input+'.old_to_be_deleted'
    mglog.debug('Modifying runcard settings: '+str(settings))
    if os.path.isfile(run_card_old): os.unlink(run_card_old) # delete old backup
    os.rename(run_card_input, run_card_old) # change name of original card

    oldCard = open(run_card_old, 'r')
    newCard = open(process_dir+'/Cards/run_card.dat', 'w')
    used_settings = []
    for line in iter(oldCard):
        if not line.strip().startswith('#'): # line commented out
            command = line.split('!', 1)[0]
            comment = line.split('!', 1)[1] if '!' in line else ''
            if '=' in command:
                setting = command.split('=')[-1] #.strip()
                stripped_setting = setting.strip()
                oldValue = '='.join(command.split('=')[:-1])
                if stripped_setting in settings:
                    # if setting set to 'None' it will be removed from run_card
                    if settings[stripped_setting]==None:
                        line=''
                        mglog.info('Removing '+stripped_setting+'.')
                        used_settings += [ stripped_setting ]
                    else:
                        line = oldValue.replace(oldValue.strip(), str(settings[stripped_setting]))+'='+setting
                        if comment != '': line += '  !' + comment
                        mglog.info('Setting '+stripped_setting+' = '+str(settings[stripped_setting])+'.')
                        used_settings += [ stripped_setting ]
        newCard.write(line)

    # Clean up unused options
    for asetting in settings:
        if asetting in used_settings: continue
        if settings[asetting]==None: continue
        mglog.warning('Option '+asetting+' was not in the default run_card.  Adding by hand a setting to '+str(settings[asetting]) )
        newCard.write( ' '+str(settings[asetting])+'   = '+str(asetting)+'\n')
    # close files
    oldCard.close()
    newCard.close()
    mglog.info('Finished modification of run card.')
    if run_card_backup is None: os.unlink(run_card_old)


def print_cards_from_dir(process_dir=MADGRAPH_GRIDPACK_LOCATION):
    card_dir=process_dir+'/Cards/'
    print_cards(proc_card=card_dir+'proc_card_mg5.dat',run_card=card_dir+'run_card.dat',param_card=card_dir+'param_card.dat',\
                madspin_card=card_dir+'madspin_card.dat',reweight_card=card_dir+'reweight_card.dat',warn_on_missing=False)


def print_cards(proc_card='proc_card_mg5.dat',run_card=None,param_card=None,madspin_card=None,reweight_card=None,warn_on_missing=True):
    if os.access(proc_card,os.R_OK):
        mglog.info("proc_card:")
        procCard = subprocess.Popen(['cat',proc_card])
        procCard.wait()
    elif warn_on_missing:
        mglog.warning('No proc_card: '+proc_card+' found')

    if run_card is not None and os.access(run_card,os.R_OK):
        mglog.info("run_card:")
        runCard = subprocess.Popen(['cat',run_card])
        runCard.wait()
    elif run_card is not None and warn_on_missing:
        mglog.warning('No run_card: '+run_card+' found')
    else:
        mglog.info('Default run card in use')

    if param_card is not None and os.access(param_card,os.R_OK):
        mglog.info("param_card:")
        paramCard = subprocess.Popen(['cat',param_card])
        paramCard.wait()
    elif param_card is not None and warn_on_missing:
        mglog.warning('No param_card: '+param_card+' found')
    else:
        mglog.info('Default param card in use')

    if madspin_card is not None and os.access(madspin_card,os.R_OK):
        mglog.info("madspin_card:")
        madspinCard = subprocess.Popen(['cat',madspin_card])
        madspinCard.wait()
    elif madspin_card is not None and warn_on_missing:
        mglog.warning('No madspin_card: '+madspin_card+' found')
    else:
        mglog.info('No madspin card in use')

    if reweight_card is not None and os.access(reweight_card,os.R_OK):
        mglog.info("reweight_card:")
        madspinCard = subprocess.Popen(['cat',reweight_card])
        madspinCard.wait()
    elif reweight_card is not None and warn_on_missing:
        mglog.warning('No reweight_card: '+reweight_card+' found')
    else:
        mglog.info('No reweight card in use')


def is_gen_from_gridpack():
    """ Simple function for checking if there is a grid pack.
    Relies on the specific location of the unpacked gridpack (madevent)
    which is here set as a global variable. The gridpack is untarred by
    the transform (Gen_tf.py) and no sign is sent to the job itself
    that there is a gridpack in use except the file's existence"""
    if os.access(MADGRAPH_GRIDPACK_LOCATION,os.R_OK):
        mglog.info('Located input grid pack area')
        return True
    return False


def is_NLO_run(process_dir=MADGRAPH_GRIDPACK_LOCATION):
    isNLO=False

    lo_config_card=process_dir+'/Cards/me5_configuration.txt'
    nlo_config_card=process_dir+'/Cards/amcatnlo_configuration.txt'

    if os.access(lo_config_card,os.R_OK) and not os.access(nlo_config_card,os.R_OK):
        isNLO=False
    elif os.access(nlo_config_card,os.R_OK) and not os.access(lo_config_card,os.R_OK):
        isNLO=True
    else:
        mglog.error("Neither configuration card found in "+process_dir+". Unable to determine LO or NLO process!")
        RuntimeError('Unable to locate configuration card')

    return isNLO


def run_card_consistency_check(isNLO=False,process_dir='.'):
    cardpath=process_dir+'/Cards/run_card.dat'
    mydict=getDictFromCard(cardpath)

    for k,v in mydict.iteritems():
        mglog.info( '"%s" = %s'%(k,v) )

    # We should always use event_norm = average [AGENE-1725] otherwise Pythia cross sections are wrong
    if not checkSetting('event_norm','average',mydict):
        modify_run_card(process_dir=process_dir,settings={'event_norm':'average'})
    # Only needed for 2.5.0 to 2.6.0 to battle problem with inconsistent event weights [AGENE-1542]
    # Will likely cause Pythia to calculate the wrong cross section [AGENE-1725]

    if not isNLO:
        #Check CKKW-L setting
        if float(mydict['ktdurham']) > 0 and int(mydict['ickkw']) != 0:
            log='Bad combination of settings for CKKW-L merging! ktdurham=%s and ickkw=%s.'%(mydict['ktdurham'],mydict['ickkw'])
            mglog.error(log)
            raise RuntimeError(log)

        # Check if user is trying to use deprecated syscalc arguments with the other systematics script
        if not 'systematics_program' in mydict or mydict['systematics_program']=='systematics':
            syscalc_settings=['sys_pdf', 'sys_scalefact', 'sys_alpsfact', 'sys_matchscale']
            found_syscalc_setting=False
            for s in syscalc_settings:
                if s in mydict:
                    mglog.warning('Using syscalc setting %s with new systematics script. Systematics script is default from 2.6.2 and steered differently (https://cp3.irmp.ucl.ac.be/projects/madgraph/wiki/Systematics#Systematicspythonmodule)'%(s))
                    found_syscalc_setting=True
            if found_syscalc_setting:
                syst_arguments=MadGraphSystematicsUtils.convertSysCalcArguments(mydict)
                mglog.info('Converted syscalc arguments to systematics arguments: '+syst_arguments)
                syst_settings_update={'systematics_arguments':syst_arguments}
                for s in syscalc_settings:
                    syst_settings_update[s]=None
                modify_run_card(process_dir=process_dir,settings=syst_settings_update)


    # usually the pdf and systematics should be set during modify_run_card
    # but check again in case the user did not call the function or provides a different card here
    mglog.info('Checking PDF and systematics settings')
    if not MadGraphSystematicsUtils.base_fragment_setup_check(MADGRAPH_PDFSETTING,mydict,isNLO):
        # still need to set pdf and systematics
        syst_settings=MadGraphSystematicsUtils.get_pdf_and_systematic_settings(MADGRAPH_PDFSETTING,isNLO)
        modify_run_card(process_dir=process_dir,settings=syst_settings)

    mydict_new=getDictFromCard(cardpath)
    if 'systematics_arguments' in mydict_new:
        systematics_arguments=MadGraphSystematicsUtils.parse_systematics_arguments(mydict_new['systematics_arguments'])
        if not 'weight_info' in systematics_arguments:
            mglog.info('Enforcing systematic weight name convention')
            systematics_arguments['weight_info']=MadGraphSystematicsUtils.SYSTEMATICS_WEIGHT_INFO
            modify_run_card(process_dir=process_dir,settings={'systematics_arguments':MadGraphSystematicsUtils.write_systematics_arguments(systematics_arguments)})

    if not isNLO:
        if not 'python_seed' in mydict:
            mglog.warning('No python seed set in run_card -- adding one with same value as iseed')
            modify_run_card(process_dir=process_dir,settings={'python_seed':mydict['iseed']})

    mglog.info('Finished checking run card - All OK!')


def hack_gridpack_script(reweight_card):
    need_to_add_rwgt=reweight_card!=None

    run_card_dict=getDictFromCard(MADGRAPH_GRIDPACK_LOCATION+'/Cards/run_card.dat',lowercase=True)

    systematics_program=None
    if settingIsTrue(run_card_dict['use_syst']):
        systematics_program='systematics'
        if checkSettingExists('systematics_program',run_card_dict):
            if checkSetting('systematics_program','systematics',run_card_dict):
                systematics_program='systematics'
            if checkSetting('systematics_program','syscalc',run_card_dict):
                systematics_program='syscalc'
            if checkSetting('systematics_program','none',run_card_dict):
                systematics_program=None
    need_to_add_syst=systematics_program!=None

    systematics_arguments=''
    if checkSettingExists('systematics_arguments',run_card_dict):
        sys_dict=MadGraphSystematicsUtils.parse_systematics_arguments(run_card_dict['systematics_arguments'])
        for s in sys_dict:
            systematics_arguments+=' --'+s+'='+sys_dict[s]

    # add systematics calculation and reweighting to run.sh
    runscript=MADGRAPH_GRIDPACK_LOCATION+'/bin/run.sh'
    oldscript = open(runscript,'r')
    newscript = open(runscript+'.tmp','w')
    # in older MG versions the gridpack is run with the command below
    gridrun_line_old='./bin/gridrun $num_events $seed'
    syst_line_old=''
    reweight_line_old='./bin/madevent reweight GridRun_${seed} -f\n'
    # in new versions it is run like this
    gridrun_line_new='${DIR}/bin/gridrun $num_events $seed $gran'
    syst_line_new=''
    reweight_line_new='${DIR}/bin/madevent reweight GridRun_${seed} -f\n'

    for line in oldscript:
        if (need_to_add_rwgt or need_to_add_syst) and gridrun_line_old in line:
            newscript.write(line)
            # run systematics
            if need_to_add_syst:
                newscript.write(syst_line_old)
                need_to_add_syst=False
            # reweight
            if need_to_add_rwgt:
                newscript.write(reweight_line_old)
                need_to_add_rwgt=False
        elif (need_to_add_rwgt or need_to_add_syst) and gridrun_line_new in line:
            newscript.write(line)
            # run systematics
            if need_to_add_syst:
                newscript.write(syst_line_new)
                need_to_add_syst=False
            # reweight
            if need_to_add_rwgt:
                newscript.write(reweight_line_new)
                need_to_add_rwgt=False

        else:
            newscript.write(line)
    oldscript.close()
    newscript.close()
    mglog.info('created '+runscript+'.tmp')

    if reweight_card and need_to_add_rwgt:
        raise RuntimeError('Could not add reweighting to gridpack script: '+runscript+' maybe line to generate events changed')
    shutil.move(runscript+'.tmp',runscript)
    st = os.stat(runscript)
    os.chmod(runscript, st.st_mode | stat.S_IEXEC)
