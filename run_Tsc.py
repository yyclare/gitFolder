#!/usr/bin/env python

import argparse
import os
import shutil
import sys
import subprocess
import time
import platform

from Bio import SeqIO
from config import Config
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def main():
    """
        Master script of running shotgun metagenomics pipeline.
    """

    # parse arguments
    parser = argparse.ArgumentParser(
        description='''A pipeline to analyze shotgun metagenomics NGS Illumina data.''',
        epilog='''Please read the documentation for more details.''')
    parser.add_argument('-c', '--config', nargs=1, type=file, required=True,
                        help='configuration file for the pipeline', metavar='shotgun_config.cfg',
                        dest='config')
    parser.add_argument('-q', '--fastqc', action='store_true',
                        help='run fastqc to check the sequence quality', dest='q')
    parser.add_argument('-t', '--trimmomatic', action='store', nargs='*',
                        help='run trimmomatic to clean the sequences', choices=["n", "r", "f"], dest='t')
    parser.add_argument('-r', '--randomlysample', action='store_true',
                        help='run randomlysample to rarefy the reads number to the lowest among the samples', dest='r')
    parser.add_argument('-a', '--assembly', action='store', nargs='*',
                        help='run assembly tools to assemble the reads to contigs', choices=["n", "f", "r", "rf", "soap",
                                                                                             "idba", "ray"], dest='a')
    parser.add_argument('-p', '--ppsp', action='store', nargs='*', help='run ppsp to perform taxonomic binning',
                        choices=["m", "n", "f"], dest='p')
    parser.add_argument('-x', '--taxator', action='store_true', help='run taxator-tk to perform taxaonomic binning',
                        dest='x')
    parser.add_argument('-f', '--annotation', action='store',  help='run functional annotation',
                        choices=[None, "h"], dest='f', nargs='*')


    args = parser.parse_args()
    

    # read configuration file
    config = Config(args.config[0], 'SMP')


    # pipeline directory
    pipelineDir = os.path.normpath(config.get('pipelineDir'))
    if not os.path.isdir(pipelineDir):
        print("Pipeline directory doesn't exist:", pipelineDir)
        return

    # create the following dirs in the pipeline directory if they do not exist
    workingDir = os.path.join(pipelineDir, 'working')
    outputDir = os.path.join(pipelineDir, 'output')
    profileDir = os.path.join(pipelineDir, 'profile')
    fastqcDir = os.path.join(outputDir, 'fastqc')
    trimmomaticDir = os.path.join(outputDir, 'trimmomatic')
    trimmomaticProfileDir = os.path.join(profileDir, 'trimmomatic')
    randomlysampleDir = os.path.join(outputDir, 'randomlysample')
    #randomlysampleProfileDir = os.path.join(profileDir, 'randomlysample')
    randomlysampleWorkingDir = os.path.join(workingDir, 'randomlysample')
    filterreadsDir = os.path.join(outputDir, 'filterreads')
    filterreadsWorkingDir = os.path.join(workingDir, 'filterreads')
    forAssemblyWorkingDir = os.path.join(workingDir, 'forAssembly')
    soapdenovoDir = os.path.join(outputDir, 'soapdenovo')
    soapdenovoWorkingDir = os.path.join(forAssemblyWorkingDir, 'soapdenovo')
    idbaDir = os.path.join(outputDir,'IDBA-UD')
    idbaWorkingDir = os.path.join(forAssemblyWorkingDir,'IDBA-UD')
    rayDir = os.path.join(outputDir,'Ray')
    rayWorkingDir = os.path.join(forAssemblyWorkingDir, 'Ray')
    assemblyDir = os.path.join(outputDir, 'assembly')
    ppspDir = os.path.join(outputDir, 'PPSP')
    #taxatorDir = os.path.join(outputDir, 'taxator')
    annotationDir = os.path.join(outputDir, 'annotation')
    orfDir = os.path.join(annotationDir, 'ORFs')
    hmmDir = os.path.join(annotationDir, 'HMM')



    dirArray = [workingDir, outputDir, fastqcDir, trimmomaticDir, profileDir, trimmomaticProfileDir,
                randomlysampleDir, randomlysampleWorkingDir, filterreadsDir, filterreadsWorkingDir,
                forAssemblyWorkingDir, soapdenovoDir, idbaDir, rayDir, assemblyDir, ppspDir, annotationDir, orfDir, hmmDir]


    for dirPath in dirArray:
        if not os.path.isdir(dirPath):
            try:
                os.mkdir(dirPath)
            except OSError:
                print ("Can't create the directory", dirPath)
                return

    # set up log file
    logDir = os.path.join(outputDir, 'log')
    if not os.path.exists(logDir):
        try:
            os.mkdir(logDir)
        except OSError:
            print("Can't create log directory:", logDir)
            return

    # read multiple input raw fastqc reads
    #def readInputForward():
    #    count = 1
    #    list1 = []
        #list2 = []
        #for i in range(10000):
        #    inputFastqFileForward = config.get('inputFastqFileForward'+str(count))
            #inputFastqFileReverse = config.get('inputFastqFileReverse'+str(count))
            #if inputFastqFileForward is not None:
            #    count += 1
            #    list1.append(inputFastqFileForward)
                #list2.append(inputFastqFileReverse)
        #listf = list1
        #listr = list2
        #if listf != '':
        #    return listf

    #Judge if you need to split the original File, only work for fastq File
    #inputFileJudge = config.get('pairedOrNot')
    #if inputFileJudge == 'True':
    #    SplitPairedFile(os.path.normpath(config.get('inputPiredFile')), os.path.dirname)



    count = 1
    list1 = []
    list2 = []
    for i in range(10000):
        inputFastqFileForward = config.get('inputFastqFileForward'+str(count))
        inputFastqFileReverse = config.get('inputFastqFileReverse'+str(count))
        if inputFastqFileForward is None:
            break
        else:
            list1.append(inputFastqFileForward)
            list2.append(inputFastqFileReverse)
            count += 1
    listf = list1
    listr = list2

    # run fastqc
    if args.q is not None:
        fastqcInstallDir = os.path.normpath(config.get('fastqcInstallDir'))

        for i, j in zip(listf, listr):
            if j is not None:
                cmd = str(os.path.join(fastqcInstallDir, 'fastqc') + ' -o ' + fastqcDir + ' -f fastq ' +
                      os.path.normpath(i) + ' ' + os.path.normpath(j))
            else: #for single reads
                cmd = str(os.path.join(fastqcInstallDir, 'fastqc') + ' -o ' + fastqcDir + ' -f fastq ' +
                      os.path.normpath(i))
            proc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=fastqcDir)
            proc.wait() #check proc.wait



    # run trimmomatic-0.32 to trim the adapters and clean the sequences according to quality score and the seq length.
    # The program first cut the adapters then cleans the seqs with quality scores and in the end discard the seqs with
    # length shorter than the minimum length.
    # Yao's ubuntu's version, needs to change to server version.
    if args.t is not None:
        trimmomaticInstallDir = os.path.normpath(config.get('trimmomaticInstallDir'))
        #lista=[]
        #listb=[]
        for i, j in zip(listf, listr):
            cmd = str('java -jar' + ' ' + os.path.join(trimmomaticInstallDir, 'trimmomatic-0.32.jar') + ' '
                      + config.get('seqType') + ' -phred' + config.get('phred') + ' -trimlog ' +
                      os.path.join(logDir, os.path.basename(i).split('.')[0]+'trimmomatic.log')
                      + ' ' + os.path.normpath(i) + ' ' + os.path.normpath(j) + ' ' + os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]
                      + '_paired.fq') + ' ' + os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0] + '_unpaired.fq') + ' ' + os.path.join(trimmomaticDir, os.path.basename(j).split('.')[0]
                      + '_paired.fq') + ' ' + os.path.join(trimmomaticDir, os.path.basename(j).split('.')[0] + '_unpaired.fq') + ' ' +'ILLUMINACLIP:'
                      + trimmomaticProfileDir + '/' + config.get('adapter') + ':' + config.get('seedMismatches') + ':'
                      + config.get('palindromeClipThreshold') + ':' + config.get('simpleClipThreshold') + ' '
                      + 'LEADING:' + config.get('leading') + ' ' + 'TRAILING:' + config.get('trailing') + ' '
                      + 'SLIDINGWINDOW:' + config.get('windowSize') + ':' + config.get('requiredQuality') + ' '
                      + 'MINLEN:' + config.get('minLength'))
            #print cmd
            proc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=trimmomaticDir)
            proc.wait()

            # unify the sequences ids pattern (there are old and new versions of seq ids, after this step, the id pattern will be "A /1" or "A 1:N:0:GA..." )
            seqidSuffix((os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_paired.fq')),
                                 (os.path.join(trimmomaticDir, os.path.basename(j).split('.')[0]+'_paired.fq')),
                                 (os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_temp_paired.fq')),
                                 (os.path.join(trimmomaticDir, os.path.basename(j).split('.')[0]+'_temp_paired.fq')))

            # combine separated forward and reverse files to one single fastaq file (/1,/2,/1,/2).
            combinePairedEndFileFastq((os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_paired.fq')),
                                 (os.path.join(trimmomaticDir, os.path.basename(j).split('.')[0]+'_paired.fq')),
                                 (os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_combined.fq')))

            if 'n' in args.t:
                fastqTofasta((os.path.join(trimmomaticDir, (os.path.basename(i).split('.')[0] + '_paired.fq'))), 1,
                             (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyTrim_forward.fa')))
                fastqTofasta((os.path.join(trimmomaticDir, (os.path.basename(j).split('.')[0] + '_paired.fq'))), 2,
                             (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyTrim_reverse.fa')))
                removeN((os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyTrim_forward.fa')),
                        (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyTrim_forward_noN.fa')), config.get('removeNlen'))
                removeN((os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyTrim_reverse.fa')),
                        (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyTrim_reverse_noN.fa')), config.get('removeNlen'))
                getConsensusSeqs((os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyTrim_forward_noN.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyTrim_reverse_noN.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyTrim_forward_consensus.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyTrim_reverse_consensus.fa')))

            #print args.t
            #First run RandomlySample then run FilterReads-P(MPI-version available) -
            # yao's ubuntu's version, needs to change to server version.
            elif 'r' in args.t and ('f' in args.t):
                # combine separated forward and reverse files to one single fastaq file (/1,/2,/1,/2).
                combinePairedEndFileFastq((os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_paired.fq')),
                                 (os.path.join(trimmomaticDir, os.path.basename(j).split('.')[0]+'_paired.fq')),
                                 (os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_combined.fq')))

                randomlysampleInstallDir = os.path.normpath(config.get('randomlysampleInstallDir'))
                filterreadsInstallDir = os.path.normpath(config.get('filterreadsInstallDir'))
                cmd2 = str(os.path.join(randomlysampleInstallDir, 'RandomlySample') + ' --ignore-quality=1 --by-pair=1 --num-samples=' +
                           config.get('randomlySampleNum-samples') + ' --format-output=1' + ' --fastq-base-quality=33 ' + ' --temp-dir=' +
                           os.path.join(randomlysampleWorkingDir, 'randomlysampletmp1') + ' --keep-temp-dir=' +
                           os.path.join(randomlysampleWorkingDir, 'randomlysampletmp2') +
                           ' --input-file ' + os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_combined.fq') +
                           ' --output-file ' + os.path.join(randomlysampleDir, os.path.basename(i).split('.')[0]+'_combined_random.fa'))
                proc2 = subprocess.Popen(cmd2, shell=True, bufsize=-1, cwd=randomlysampleDir)
                proc2.wait()
                print "RandomlySample return code:", proc2.returncode
                if proc2.returncode != 0:
                    raise Exception("Command returned with non-zero %s status: %s" % (proc2.returncode, cmd2))
                print "Filtering reads based on kmer counting is going to start."
                cmd3 = str('export LD_LIBRARY_PATH=/home/yao/Softwares/openmpi-1.6.5/lib:$LD_LIBRARY_PATH;mpirun -np ' + config.get('numberOfProcess_mpi') + ' ' + os.path.join(filterreadsInstallDir, 'FilterReads-P')
                           + ' --histogram-file ' + os.path.join(filterreadsWorkingDir, os.path.basename(i).split('.')[0]+'_histogram')
                           + ' --size-history-file ' + os.path.join(filterreadsWorkingDir, os.path.basename(i).split('.')[0]+'_sizeHistory')
                           + ' --ignore-quality=1 --format-output=1 --keep-read-comment=1' + ' --kmer-size=' + config.get('kmerSize_filterReads')
                           + ' --min-depth=' + config.get('minDepth_filterReads') + ' --separate-output=0 --min-passing-in-pair=2'
                           + ' --min-read-length=' + config.get('minReadLen_filterReads') + ' --max-kmer-output-depth=-1 '
                           + ' --skip-artifact-filter=1 --dedup-single=1 --dedup-edit-distance=0' + ' --input-file '
                           + os.path.join(randomlysampleDir, os.path.basename(i).split('.')[0]+'_combined_random.fa')
                           + ' --output-file ' + os.path.join(filterreadsDir, os.path.basename(i).split('.')[0]+'_combined_random_filter.fa')
                           + ' --log-file ' + os.path.join(logDir, os.path.basename(i).split('.')[0]+'_filterReads_log'))
                proc3 = subprocess.Popen(cmd3, shell=True, bufsize=-1, cwd=filterreadsDir)
                proc3.wait()
                print "Filterreads return code:", proc3.returncode
                if proc3.returncode != 0:
                    raise Exception("Command returned with non-zero %s status: %s" % (proc3.returncode, cmd3))
                #
                splitPairedReads((os.path.join(filterreadsDir, os.path.basename(i).split('.')[0]+'_combined_random_filter.fa')),
                                 (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0]+'_random_filter_forward.fa')),
                                 (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0]+'_random_filter_reverse.fa')))
                removeN((os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_random_filter_forward.fa')),
                        (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_random_filter_forward_noN.fa')), config.get('removeNlen'))
                removeN((os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_random_reverse.fa')),
                        (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_random_filter_reverse_noN.fa')), config.get('removeNlen'))
                getConsensusSeqs((os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_random_filter_forward_noN.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_random_filter_reverse_noN.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_random_filter_forward_consensus.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_random_filter_reverse_consensus.fa')))



            # Only run RandomlySample (without FilterReads-P)
            elif 'r' in args.t:
                # combine separated forward and reverse files to one single fastaq file (/1,/2,/1,/2).
                combinePairedEndFileFastq((os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_paired.fq')),
                                 (os.path.join(trimmomaticDir, os.path.basename(j).split('.')[0]+'_paired.fq')),
                                 (os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_combined.fq')))

                randomlysampleInstallDir = os.path.normpath(config.get('randomlysampleInstallDir'))
                cmd2 = str(os.path.join(randomlysampleInstallDir, 'RandomlySample') + ' --ignore-quality=1 --by-pair=1 --num-samples=' +
                           config.get('randomlySampleNum-samples') + ' --format-output=1' + ' --fastq-base-quality=33 ' + ' --temp-dir=' +
                           os.path.join(randomlysampleWorkingDir, 'randomlysampletmp1') + ' --keep-temp-dir=' +
                           os.path.join(randomlysampleWorkingDir, 'randomlysampletmp2') +
                           ' --input-file ' + os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_combined.fq') +
                           ' --output-file ' + os.path.join(randomlysampleDir, os.path.basename(i).split('.')[0]+'_combined_onlyRandom.fa'))
                proc2 = subprocess.Popen(cmd2, shell=True, bufsize=-1, cwd=randomlysampleDir)
                proc2.wait()
                print "RandomlySample return code:", proc2.returncode
                if proc2.returncode != 0:
                    raise Exception("Command returned with non-zero %s status: %s" % (proc2.returncode, cmd2))
                print "You only rarefied the reads and you did not want to filter the reads, right?"

                splitPairedReads((os.path.join(randomlysampleDir, os.path.basename(i).split('.')[0]+'_combined_onlyRandom.fa')),
                                 (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0]+'_onlyRandom_forward.fa')),
                                 (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0]+'_OnlyRandom_reverse.fa')))
                removeN((os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyRandom_forward.fa')),
                        (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyRandom_forward_noN.fa')), config.get('removeNlen'))
                removeN((os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyRandom_reverse.fa')),
                        (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyRandom_reverse_noN.fa')), config.get('removeNlen'))
                getConsensusSeqs((os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyRandom_forward_noN.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyRandom_reverse_noN.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyRandom_forward_consensus.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyRandom_reverse_consensus.fa')))

            #Only run FilterReads-P(MPI-version available)--yao's ubuntu's version, needs to change to server version.
            elif 'f' in args.t:
                # combine separated forward and reverse files to one single fastaq file (/1,/2,/1,/2).
                combinePairedEndFileFastq((os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_paired.fq')),
                                          (os.path.join(trimmomaticDir, os.path.basename(j).split('.')[0]+'_paired.fq')),
                                          (os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_combined.fq')))

                filterreadsInstallDir = os.path.normpath(config.get('filterreadsInstallDir'))
                #cmd3 = str('export LD_LIBRARY_PATH=/home/yao/Softwares/openmpi-1.6.5/lib:$LD_LIBRARY_PATH;mpirun -np '
                cmd3 = str('mpirun -np '
                           + config.get('numberOfProcess_mpi') + ' ' + os.path.join(filterreadsInstallDir, 'FilterReads-P')
                           + ' --histogram-file ' + os.path.join(filterreadsWorkingDir, os.path.basename(i).split('.')[0]+'_onlyFiltering_histogram')
                           + ' --size-history-file ' + os.path.join(filterreadsWorkingDir, os.path.basename(i).split('.')[0]+'_onlyFiltering_sizeHistory')
                           + ' --ignore-quality=1 --format-output=1 --keep-read-comment=1' + ' --kmer-size=' + config.get('kmerSize_filterReads')
                           + ' --min-depth=' + config.get('minDepth_filterReads') + ' --separate-output=0 --min-passing-in-pair=2'
                           + ' --min-read-length=' + config.get('minReadLen_filterReads') + ' --max-kmer-output-depth=-1 '
                           + ' --skip-artifact-filter=1 --dedup-single=1 --dedup-edit-distance=0' + ' --input-file '
                           + os.path.join(os.path.join(trimmomaticDir, os.path.basename(i).split('.')[0]+'_combined.fq'))
                           + ' --output-file ' + os.path.join(filterreadsDir, os.path.basename(i).split('.')[0]+'_combined_onlyFilter.fa')
                           + ' --log-file ' + os.path.join(logDir, os.path.basename(i).split('.')[0]+'_onlyFilterReadslog'))
                proc3 = subprocess.Popen(cmd3, shell=True, bufsize=-1, cwd=filterreadsDir)
                proc3.wait()
                print "Filterreads return code:", proc3.returncode
                if proc3.returncode != 0:
                    raise Exception("Command returned with non-zero %s status: %s" % (proc3.returncode, cmd3))
                print "You only filtered the reads, you did not want to rarefy the reads, right?"

                splitPairedReads((os.path.join(filterreadsDir, os.path.basename(i).split('.')[0]+'_combined_onlyFilter.fa')),
                                 (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0]+'_onlyFilter_forward.fa')),
                                 (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0]+'_OnlyFilter_reverse.fa')))
                removeN((os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyFilter_forward.fa')),
                        (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyFilter_forward_noN.fa')), config.get('removeNlen'))
                removeN((os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyFilter_reverse.fa')),
                        (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyFilter_reverse_noN.fa')), config.get('removeNlen'))
                getConsensusSeqs((os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyFilter_forward_noN.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyFilter_reverse_noN.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyFilter_forward_consensus.fa')),
                                  (os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyFilter_reverse_consensus.fa')))

            #lista.append(os.path.basename(i).split('.')[0] + '_paired.fq ')
            #listb.append(os.path.basename(j).split('.')[0] + '_paired.fq ')
        #print lista, listb

                #sys.stdout.write(cmd3)


    fileHolder = FileHolder(outputTrimForward=[], outputTrimReverse=[], outputRandomForward=[], outputRandomReverse=[],
                            outputFilterForward=[], outputFilterReverse=[], outputRandomFilterForward=[],
                            outputRandomFilterReverse=[])

    if args.a:
        if "soap" in args.a:
            def soapAssembly(file1,file2):
                counter = 1
                for i in range(1000):
                    readsNumberPerSample = config.get('sample'+str(counter))

                    if readsNumberPerSample is not None:
                        soapdenovoPerSampleDir = os.path.join(soapdenovoDir, ('sample'+str(counter)))
                        if not os.path.exists(soapdenovoPerSampleDir):
                            try:
                                os.mkdir(soapdenovoPerSampleDir)
                            except OSError:
                                print("Can't create soapdenovoPerSample directory:", soapdenovoPerSampleDir)
                        soapdenovoPerSampleWorkingDir = os.path.join(forAssemblyWorkingDir, ('soap_'+'sample'+str(counter)))
                        if not os.path.exists(soapdenovoPerSampleWorkingDir):
                            try:
                                os.mkdir(soapdenovoPerSampleWorkingDir)
                            except OSError:
                                print("Can't create soapdenovoPerSample directory:", soapdenovoPerSampleWorkingDir)

                        soapConfigFile = open(os.path.join(soapdenovoPerSampleDir, ('sample'+str(counter) + '.config')), 'w')
                        soapConfigFile.write('max_rd_len=' + config.get('sample'+str(counter)+'_' + 'maxReadLen') + '\n'
                                         + '[LIB]' + '\n' + '#average insert size'
                                         + '\n' + 'avg_ins=' + config.get('sample'+str(counter)+'_'+'averageInsertSize')
                                         + '\n' +'#if sequence needs to be reversed'
                                         + '\n' + 'reverse_seq=' + config.get('sample'+str(counter)+'_' + 'reverseSeq')
                                         + '\n' + '#in which part(s) the reads are used'
                                         + '\n' + 'asm_flags=' + config.get('sample'+str(counter)+'_' + 'asmFlag')
                                         + '\n' + '#paired data files' + '\n')

                        for k, l in zip(file1[(int(readsNumberPerSample.split('-')[0])-1)/2:int(readsNumberPerSample.split('-')[-1])/2],
                            (file2[(int(readsNumberPerSample.split('-')[0])-1)/2:int(readsNumberPerSample.split('-')[-1])/2])):

                            shutil.copy(os.path.join(forAssemblyWorkingDir, k), soapdenovoPerSampleWorkingDir)
                            shutil.copy(os.path.join(forAssemblyWorkingDir, l), soapdenovoPerSampleWorkingDir)
                        #print 'p=',k
                        #print 'p=',l
                            soapConfigFile.write('p='+k+'\n'+'p='+l+'\n')
                        soapConfigFile.close()


                        soapdenovoInstallDir = os.path.normpath(config.get('soapdenovoInstallDir'))
                        runSoapdenovo = str('export PATH=/net/programs/Debian-6.0.3-x86_64/soapdenovo-2.04:$PATH'+';' + ' '
                                            +os.path.join(soapdenovoInstallDir, config.get('soapdenovoMethod')) + ' all -s '
                                + os.path.join(os.path.normpath(soapdenovoPerSampleDir), 'sample'+str(counter) + '.config')
                                + ' -K '+config.get('assemblyKmerSize') + ' -R' + ' -p ' + config.get('n_cpu') + ' -o '
                                + 'soapdenovo_sample' + str(counter) + '_' + ' ' + '1>'+'sample'+str(counter)+'_assembly.log' + ' '
                                + '2>'+'sample'+str(counter) + '_assembly.error')
                        #print runSoapdenovo
                        procRunSoapdenovo = subprocess.Popen(runSoapdenovo, shell=True, bufsize=-1, cwd=soapdenovoPerSampleDir)
                        procRunSoapdenovo.wait()
                        print "Soapdenovo return code:", procRunSoapdenovo.returncode
                        if procRunSoapdenovo.returncode != 0:
                            raise Exception("Command returned with non-zero %s status: %s" % (procRunSoapdenovo.returncode, runSoapdenovo))
                        print "Soapdenovo is finished."
                        shutil.copy(os.path.join(soapdenovoPerSampleDir, 'soapdenovo_sample' + str(counter) + '_'+'.contig'),
                                    os.path.join(assemblyDir, 'soapdenovo_sample'+str(counter)+'.contig.fa'))
                        counter += 1

        if "idba" in args.a:
            def idbaAssembly(file1, file2):
                counter = 1
                for i in range(1000):
                    readsNumberPerSample = config.get('sample'+str(counter))

                    if readsNumberPerSample is not None:
                        idbaPerSampleDir = os.path.join(idbaDir, ('sample'+str(counter)))
                        if not os.path.exists(idbaPerSampleDir):
                            try:
                                os.mkdir(idbaPerSampleDir)
                            except OSError:
                                print("Can't create idbaPerSample directory:", idbaPerSampleDir)
                    if readsNumberPerSample is not None:
                        idbaPerSampleWorkingDir = os.path.join(forAssemblyWorkingDir, ('idba_'+'sample'+str(counter)))
                        if not os.path.exists(idbaPerSampleWorkingDir):
                            try:
                                os.mkdir(idbaPerSampleWorkingDir)
                            except OSError:
                                print("Can't create idbaPerSample directory:", idbaPerSampleWorkingDir)

                        outputFile1 = open(os.path.join(idbaPerSampleWorkingDir, 'sample' + str(counter) + '.forwardMerge.fa'),'w')
                        outputFile2 = open(os.path.join(idbaPerSampleWorkingDir, 'sample' + str(counter) + '.reverseMerge.fa'),'w')

                        for k, l in zip((file1[(int(readsNumberPerSample.split('-')[0])-1)/2:int(readsNumberPerSample.split('-')[-1])/2]),
                                       (file2[(int(readsNumberPerSample.split('-')[0])-1)/2:int(readsNumberPerSample.split('-')[-1])/2])):

                            forwardDir = os.path.join(idbaPerSampleWorkingDir, 'forward')
                            reverseDir = os.path.join(idbaPerSampleWorkingDir, 'reverse')
                            if not os.path.exists(forwardDir):
                                try:
                                    os.mkdir(forwardDir)
                                except OSError:
                                    print ("Can't create forward directory:", forwardDir)
                            shutil.copy(os.path.join(forAssemblyWorkingDir, k), forwardDir)
                            inputFile1 = open(os.path.normpath(k), 'r')
                            lines = inputFile1.readlines()
                            #while True:
                            if len(lines) == 0:
                                break
                            for line in lines:
                                outputFile1.write(line)

                            if not os.path.exists(reverseDir):
                                try:
                                    os.mkdir(reverseDir)
                                except OSError:
                                    print ("Can't create reverse directory:", reverseDir)

                            shutil.copy(os.path.join(forAssemblyWorkingDir, l), reverseDir)
                            inputFile2 = open(os.path.normpath(l), 'r')
                            lines = inputFile2.readlines()
                            #while True:
                            if len(lines) == 0:
                                break
                            for line in lines:
                                outputFile2.write(line)

                        outputFile1.close()
                        outputFile2.close()

                        combinePairedEndFileFasta(os.path.join(idbaPerSampleWorkingDir, 'sample' + str(counter) + '.forwardMerge.fa'),
                                                  os.path.join(idbaPerSampleWorkingDir, 'sample' + str(counter) + '.reverseMerge.fa'),
                                                  os.path.join(idbaPerSampleWorkingDir, 'sample' + str(counter) + '.mergePaired.fa'))

                        idbaInstallDir = os.path.normpath(config.get('idbaInstallDir'))

                        runIdba = str('export PATH=/net/programs/Debian-6.0.3-x86_64/idba-1.1.1/bin:$PATH'+';' + os.path.join(idbaInstallDir, 'idba_ud') + ' '
                                      + '-' + config.get('idbaReadsLen') + ' '
                                      + os.path.join(idbaPerSampleWorkingDir, 'sample' + str(counter) + '.mergePaired.fa')
                                      + ' -o ' + idbaPerSampleDir + ' --num_threads ' + config.get('idbaNumThreads'))
                        #print runIdba
                        procRunIdba = subprocess.Popen(runIdba, shell=True, bufsize=-1, cwd=idbaPerSampleDir)
                        procRunIdba.wait()
                        print "IDBA_UD return code:", procRunIdba.returncode
                        if procRunIdba.returncode != 0:
                            raise Exception("Command returned with non-zero %s status: %s" % (procRunIdba.returncode, runIdba))
                        print "sample" + str(counter) + "IDBA_UD is finished."
                        shutil.copy(os.path.join(idbaPerSampleDir, 'contig.fa'), os.path.join(idbaPerSampleDir, 'idba_sample'+str(counter)+'.contig.fa'))
                        shutil.copy(os.path.join(idbaPerSampleDir, 'idba_sample'+str(counter)+'.contig.fa'), assemblyDir)
                        counter += 1

        if "ray" in args.a:
            def rayAssembly(file1,file2):
                counter = 1
                for i in range(1000):
                    readsNumberPerSample = config.get('sample'+str(counter))

                    if readsNumberPerSample is not None:
                        rayPerSampleDir = os.path.join(rayDir, ('sample'+str(counter)))
                        if not os.path.exists(rayPerSampleDir):
                            try:
                                os.mkdir(rayPerSampleDir)
                            except OSError:
                                print("Can't create idbaPerSample directory:", rayPerSampleDir)
                    #if readsNumberPerSample is not None:
                        rayPerSampleWorkingDir = os.path.join(forAssemblyWorkingDir, ('ray_'+'sample'+str(counter)))
                        if not os.path.exists(rayPerSampleWorkingDir):
                            try:
                                os.mkdir(rayPerSampleWorkingDir)
                            except OSError:
                                print("Can't create idbaPerSample directory:", rayPerSampleWorkingDir)

                        outputFile1 = open(os.path.join(rayPerSampleWorkingDir, 'sample' + str(counter) + '.forwardMerge.fa'),'w')
                        outputFile2 = open(os.path.join(rayPerSampleWorkingDir, 'sample' + str(counter) + '.reverseMerge.fa'),'w')

                        for k, l in zip((file1[(int(readsNumberPerSample.split('-')[0])-1)/2:int(readsNumberPerSample.split('-')[-1])/2]),
                                       (file2[(int(readsNumberPerSample.split('-')[0])-1)/2:int(readsNumberPerSample.split('-')[-1])/2])):

                            inputFile1 = open(os.path.normpath(k), 'r')
                            lines = inputFile1.readlines()
                            #while True:
                            if len(lines) == 0:
                                break
                            for line in lines:
                                outputFile1.write(line)

                            inputFile2 = open(os.path.normpath(l), 'r')
                            lines = inputFile2.readlines()
                            #while True:
                            if len(lines) == 0:
                                break
                            for line in lines:
                                outputFile2.write(line)

                        outputFile1.close()
                        outputFile2.close()

                        rayInstallDir = os.path.normpath(config.get('rayInstallDir'))
                        rayPerSampleDirOutput = os.path.join(rayPerSampleDir, 'output')
                        runRay = str('export PATH=/net/programs/Debian-6.0.3-x86_64/Ray-2.3.0:$PATH'+';' + 'mpirun -np '
                                     + config.get('numProcess_mpi_ray') + ' ' + os.path.join(rayInstallDir, 'Ray')
                                     + ' -k ' + config.get('kmerSize_ray') + ' -p '
                                     + os.path.join(rayPerSampleWorkingDir, 'sample' + str(counter) + '.forwardMerge.fa')
                                     + ' ' + os.path.join(rayPerSampleWorkingDir, 'sample' + str(counter) + '.reverseMerge.fa')
                                     + ' -o ' + rayPerSampleDirOutput)
                        #print runRay
                        procRunRay = subprocess.Popen(runRay, shell=True, bufsize=-1, cwd=rayPerSampleDir)
                        procRunRay.wait()
                        print "Ray return code:", procRunRay.returncode
                        if procRunRay.returncode != 0:
                            raise Exception("Command returned with non-zero %s status: %s" % (procRunRay.returncode, runRay))
                        print "Ray is finished."
                        shutil.copy(os.path.join(rayPerSampleDirOutput, 'Contigs.fasta'), os.path.join(rayPerSampleDir, 'ray_sample'+str(counter)+'.contig.fa'))
                        shutil.copy(os.path.join(rayPerSampleDir, 'ray_sample'+str(counter)+'.contig.fa'), assemblyDir)
                        counter += 1

        if "n" in args.a:
            for i, j in zip(listf, listr):
                fileHolder.outputTrimForward.append(os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyTrim_forward_consensus.fa'))
                fileHolder.outputTrimReverse.append(os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyTrim_reverse_consensus.fa'))
            file1 = fileHolder.outputTrimForward
            file2 = fileHolder.outputTrimReverse

            if "soap" in args.a:
                soapAssembly(file1, file2)
            if "idba" in args.a:
                #starttime = time.clock()
                idbaAssembly(file1, file2)
                #endtime = time.clock()
                #print (endtime-starttime)
            if "ray" in args.a:
                rayAssembly(file1, file2)

        if "r" in args.a:
            for i, j in zip(listf,listr):
                fileHolder.outputRandomForward.append(os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyRandom_forward_consensus.fa'))
                fileHolder.outputRandomReverse.append(os.path.join(forAssemblyWorkingDir, os.path.basename(j).split('.')[0] + '_onlyRandom_forward_consensus.fa'))
            file1 = fileHolder.outputRandomForward
            file2 = fileHolder.outputRandomReverse
            if "soap" in args.a:
                soapAssembly(file1, file2)
        #    return fileHolder.outputRandomForward, fileHolder.outputRandomReverse
            if "idba" in args.a:
                idbaAssembly(file1, file2)
            if "ray" in args.a:
                rayAssembly(file1, file2)
        if "f" in args.a:
            for i, j in zip(listf, listr):
                fileHolder.outputFilterForward.append(os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyFilter_forward_consensus.fa'))
                fileHolder.outputFilterForward.append(os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] + '_onlyFilter_reverse_consensus.fa'))
            file1 = fileHolder.outputFilterForward
            file2 = fileHolder.outputFilterReverse
            if "soap" in args.a:
                soapAssembly(file1, file2)
            if "idba" in args.a:
                idbaAssembly(file1, file2)
            if "ray" in args.a:
                rayAssembly(file1, file2)
        #    return fileHolder.outputFilterForward, fileHolder.outputFilterReverse
        if "rf" in args.a:
            for i, j in zip(listf, listr):
                fileHolder.outputRandomFilterForward.append(os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] +'_random_filter_forward_consensus.fa'))
                fileHolder.outputRandomFilterReverse.append(os.path.join(forAssemblyWorkingDir, os.path.basename(i).split('.')[0] +'_random_filter_reverse_consensus.fa'))
            file1 = fileHolder.outputRandomFilterForward
            file2 = fileHolder.outputRandomFilterReverse
            if "soap" in args.a:
                soapAssembly(file1, file2)
            if "idba" in args.a:
                idbaAssembly(file1, file2)
            if "ray" in args.a:
                rayAssembly(file1, file2)
        #    return fileHolder.outputRandomFilterForward, fileHolder.outputRandomFilterReverse

    if args.p:
        counter = 1
        for i in range(1000):
            ppspPerSampleDir = os.path.join(ppspDir, ('sample'+str(counter)))
            if not os.path.exists(ppspPerSampleDir):
                try:
                    os.mkdir(ppspPerSampleDir)
                except OSError:
                    print("Can't create ppspPerSample directory:", ppspPerSampleDir)

            filterReadLen(os.path.join(assemblyDir, config.get('ppspAssemblyMethod') + '_sample' + str(counter)+'.contig.fa'),
                          os.path.join(assemblyDir, config.get('ppspAssemblyMethod') + '_sample' + str(counter)+'.contig.ppspfilter.fa'), 1000)

            ppspConfigFile = open(os.path.join(ppspPerSampleDir, ('sample'+str(counter) + '.cfg')), 'w')
            ppspConfigFile.write('[PhyloPythiaS_Plus]' + '\n' + '# Configuration file of the PhyloPythiaS Plus pipeline.'
                                 + '\n' + 'pipelineDir=' + os.path.normpath(ppspPerSampleDir) + '\n' + '# INPUT FILES' + '\n' +
                                 'inputFastaFile=' + os.path.join(assemblyDir, config.get('ppspAssemblyMethod') + '_sample'
                                 + str(counter)+'.contig.ppspfilter.fa') + '\n' + '# Fasta file with scaffolds (optional)'
                                 + '\n' + 'inputFastaScaffoldsFile=' + noneToEmptyStr(config.get('inputFastaScaffoldsFile')) + '\n'
                                 + '# Scaffold-contig mapping, tab separated file; map: scaffold -> contigs; (optional)'
                                 + ' \n' + 'scaffoldsToContigsMapFile=' + noneToEmptyStr(config.get('scaffoldsToContigsMapFile')) + '\n'
                                 + '# Reference prediction file in PPS (*.out) format (optional)' + '\n'
                                 + 'referencePlacementFileOut=' + noneToEmptyStr(config.get('referencePlacementFileOut')) + '\n'
                                 + '# REFERENCE\n' + 'databaseFile=/net/metagenomics/projects/PPSmg/taxonomy/20121122\n'
                                 + 'refSeq=/net/metagenomics/projects/PPSmg/data/nobackup/NCBI20121122/sequences\n'
                                 + '# Exclude reference sequences according to the reference prediction file at a given '
                                   'rank (optional)\n# allowed ranks are: phylum, class, order, family, genus, species, strain\n'
                                 + 'excludeRefSeqRank=' + noneToEmptyStr(config.get('excludeRefSeqRank')) + '\n'
                                 + 's16Database=/net/metagenomics/projects/PPSmg/database/silva111' + '\n'
                                 + 'mgDatabase=/net/metagenomics/projects/PPSmg/database/mg3' + '\n' + 'excludeRefMgRank='
                                 + noneToEmptyStr(config.get('excludeRefMgRank')) + '\n' + '# TOOLS\n'
                                 + 'ppsInstallDir=/net/metagenomics/projects/PPSmg/tools/PhyloPythiaS/vm/1_3\n'
                                 + 'mothurInstallDir=/net/metagenomics/projects/PPSmg/tools/mothur/mothur_1_23\n'
                                 + 'hmmerBinDir=/net/metagenomics/projects/PPSmg/tools/hmmer-3.0/binaries\n'
                                 + 'rnaHmmInstallDir=/net/metagenomics/projects/PPSmg/tools/rna_hmm3\n'
                                 + '#BASIC SETTINGS\n' + '# e.g.: 1 ~ phylum, 2 ~ class, 3 ~ order, 4 ~ family, 5 ~ genus, 6 ~ species\n'
                                 + 'rankIdCut=' + config.get('rankIdCut') + '\n' + 'maxLeafClades=' + config.get('maxLeafClades')
                                 + '\n' + 'minPercentInLeaf=' + config.get('minPercentInLeaf') + '\n'
                                 + 'minSeqLen=' + config.get('minSeqLen') + '\n' + 'parallelPPSmodels=' + config.get('parallelPPSmodels')
                                 + '\n' + '# ADVANCED SETTINGS\n' + '# It is not recommended to change the advanced settings, ' +
                                 'however it is possible to give you more control.\n' + 'placeContigsFromTheSameScaffold=True\n'
                                 + 'agThreshold=0.3\n' + 'assignedPartThreshold=0.5\n' + 'rankIdAll=0\n' + 'weightStayAll=60.0\n'
                                 + 'rankIdCutMinBp=100\n' + 'overrideWithPPSPlacements=True\n' + 'minBpToModel=100000\n'
                                 + 'minSSDfileSize=100\n' + 'maxSSDfileSize=400000\n' + 'minGenomesWgs=1\n' + 'minBpPerSpecies=2000000\n'
                                 + 'candidatePlTopPercentThreshold=0.1\n' + 'outputFileContigSubPattern=^(.*)\n'
                                 + 'mothurClassifyParamOther=method=bayesian, cutoff=80, iters=300' + '\n' + 'configPPS=' + '\n'
                                 + 'recallMinFracClade=0.001' + '\n' + 'precisionMinFracPred=0.001' + '\n' + 'correctLabelThreshold=0.9')
            ppspConfigFile.close() #Close file is important

            #only run maker gene analysis not the whole ppsp pipeline
            if "m" in args.p:
                runPpsp = str('export PYTHONPATH=/net/metagenomics/projects/PPSmg/scripts/current:$PYTHONPATH'+';' + ' '
                              + 'python /net/metagenomics/projects/PPSmg/scripts/current/algbioi/core/run.py'
                              + ' -c ' + 'sample'+str(counter) + '.cfg' + ' -n -g -o s16 mg -s')
            #run whole ppsp pipeline without inputFastaScaffoldsFile setting in ppsp configuration file
            if "n" in args.p:
                runPpsp = str('export PYTHONPATH=/net/metagenomics/projects/PPSmg/scripts/current:$PYTHONPATH'+';' + ' ' + 'python /net/metagenomics/projects/PPSmg/scripts/current/algbioi/core/run.py'
                              + ' -c ' + os.path.join(ppspPerSampleDir,'sample'+str(counter) + '.cfg') + ' -n -g -o s16 mg -t -a -p c -r -s')
                print runPpsp
            #run whole ppsp pipeline with inputFastaScaffoldsFile setting in ppsp configuration file
            if "f" in args.p:
                runPpsp = str('export PYTHONPATH=/net/metagenomics/projects/PPSmg/scripts/current:$PYTHONPATH'+';' + ' ' + 'python /net/metagenomics/projects/PPSmg/scripts/current/algbioi/core/run.py'
                              + ' -c ' + 'sample'+str(counter) + '.cfg' + ' -n -g -o s16 mg -t -a -p c s v -r -s')

            procRunPpsp = subprocess.Popen(runPpsp, shell=True, bufsize=-1, cwd=ppspPerSampleDir)
            procRunPpsp.wait()
            print "PPSP return code:", procRunPpsp.returncode
            if procRunPpsp.returncode != 0:
                raise Exception("Command returned with non-zero %s status: %s" % (procRunPpsp.returncode, runPpsp))
            else:
                print "sample"+str(counter),"PPSP is finished."

            counter += 1
            if os.path.exists(os.path.join(assemblyDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.contig.fa')):
                pass
            else:
                break

    if args.x:
        counter = 1
        for i in range(1000):
            #taxatorPerSampleDir = os.path.join(taxatorDir, ('sample'+str(counter)))
            # if not os.path.exists(taxatorPerSampleDir):
            #     try:
            #         os.mkdir(taxatorPerSampleDir)
            #     except OSError:
            #         print("Can't create taxatorPerSample directory:", taxatorPerSampleDir)
            renameReads(os.path.join(assemblyDir, config.get('taxatorAssemblyMethod') + '_sample' + str(counter)+'.contig.fa'),
                          os.path.join(assemblyDir, config.get('taxatorAssemblyMethod') + '_sample' + str(counter)+'.contig.rename.fa'))
            taxatorInstallDir = os.path.normpath(config.get('taxatorInstallDir'))
            #shutil.copy(os.path.join(taxatorInstallDir, 'binning-working-fasta-blast.bash'), taxatorPerSampleDir)
            # runTaxator = str(os.path.join(taxatorInstallDir, 'binning-workflow-fasta-blast.bash') + ' '
            runTaxator = str(os.path.join(taxatorInstallDir, 'binning-workflow-fasta-last.sh') + ' '
                             + os.path.join(assemblyDir, config.get('taxatorAssemblyMethod') + '_sample' + str(counter)+'.contig.rename.fa'))
            # print runTaxator
            procRunTaxator = subprocess.Popen(runTaxator, shell=True, bufsize=-1, cwd=taxatorInstallDir)
            procRunTaxator.wait()
            print "Taxator-tk return code:", procRunTaxator.returncode
            if procRunTaxator.returncode != 0:
                raise Exception("Command returned with non-zero %s status: %s" % (procRunTaxator.returncode, runTaxator))
            else:
                print "sample"+str(counter), "Taxator-tk is finished."

            counter += 1
            if os.path.exists(os.path.join(assemblyDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.contig.fa')):
                pass
            else:
                break

    if args.f is not None:
        counter = 1
        for i in range(1000):
            orfPerSampleDir = os.path.join(orfDir, ('sample'+str(counter)))
            hmmPerSampleDir = os.path.join(hmmDir, ('sample'+str(counter)))
            if not os.path.exists(orfPerSampleDir):
                try:
                    os.mkdir(orfPerSampleDir)
                except OSError:
                    print("Can't create orfPerSampleDir:", orfPerSampleDir)
            if not os.path.exists(hmmPerSampleDir):
                try:
                    os.mkdir(hmmPerSampleDir)
                except OSError:
                    print("Can't create hmmPerSampleDir:", hmmPerSampleDir)

            filterReadLen(os.path.join(assemblyDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.contig.fa'),
                          os.path.join(assemblyDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.contig.annofilter.fa'), int(config.get('annoSeqLen')))

            MetaGeneMarkInstallDir = os.path.normpath(config.get('MetaGeneMarkInstallDir'))
            runMetaGeneMark = str(os.path.join(MetaGeneMarkInstallDir, 'gmhmmp') + ' -m '
                                  + os.path.join(MetaGeneMarkInstallDir, 'MetaGeneMark_v1.mod') + ' '
                                  + os.path.join(assemblyDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.contig.annofilter.fa')
                                  + ' -f G -a -d -r -o ' + os.path.join(orfPerSampleDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.ntaa.gff')
                                  + ';' + os.path.join(MetaGeneMarkInstallDir, 'aa_from_gff.pl') + ' < '
                                  + os.path.join(orfPerSampleDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.ntaa.gff')
                                  + ' > ' + os.path.join(orfPerSampleDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.gff.faa')
                                  + ';' + os.path.join(MetaGeneMarkInstallDir, 'nt_from_gff.pl') + ' < '
                                  + os.path.join(orfPerSampleDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.ntaa.gff')
                                  + ' > ' + os.path.join(orfPerSampleDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.gff.fna'))
            procRunMetaGeneMark = subprocess.Popen(runMetaGeneMark, shell=True, bufsize=-1, cwd=orfPerSampleDir)
            procRunMetaGeneMark.wait()
            print "ORFs calling return code:", procRunMetaGeneMark.returncode
            if procRunMetaGeneMark.returncode != 0:
                raise Exception("Command returned with non-zero %s status: %s" % (procRunMetaGeneMark.returncode, runMetaGeneMark))
            else:
                print "sample"+str(counter), "ORFs calling is finished."
            f1 = open(os.path.join(orfPerSampleDir, config.get('assemblyMethod') + '_sample' + str(counter)+'ORF_aa.txt'), 'w')
            f1.write(os.path.join(orfPerSampleDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.gff.faa'))
            f1.close()

            if None in args.f:
                counter += 1
            #
                if os.path.exists(os.path.join(assemblyDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.contig.fa')):
                    pass
                else:
                    break



            if "h" in args.f:
                hmmScriptInstallDir = os.path.normpath(config.get('hmmScriptInstallDir'))
                runAnnotationScriptPath = os.path.normpath("/net/metagenomics/") #at the moment, Aaron's pipeline only support to run under the /net/metagenomics/
                arrayOutDir = os.path.join(hmmPerSampleDir, (config.get('assemblyMethod') + '_sample' + str(counter) + '_array_out'))
                if not os.path.exists(arrayOutDir):
                    try:
                        os.mkdir(arrayOutDir)
                    except OSError:
                        print("Can't create arrayOutDir:", arrayOutDir)
                runHmm = str('importpackage parallel; python ' + os.path.join(hmmScriptInstallDir, 'write_annot_array_script.py') + ' -d '
                             + os.path.join(orfPerSampleDir, config.get('assemblyMethod') + '_sample' + str(counter)+'ORF_aa.txt')[len(runAnnotationScriptPath)+1:]
                             + ' -t ' + os.path.normpath(hmmPerSampleDir)[len(runAnnotationScriptPath)+1:] + ' -f '
                             + os.path.normpath(arrayOutDir)[len(runAnnotationScriptPath)+1:] + ' -s ' + config.get('annotationSource') + ' -h ')
                # print runHmm
                procRunHmm = subprocess.Popen(runHmm, shell=True, bufsize=-1, cwd=runAnnotationScriptPath)
                procRunHmm.wait()
                print "Hmm annotation return code:", procRunHmm.returncode
                if procRunHmm.returncode != 0:
                    raise Exception("Command returned with non-zero %s status : %s" % (procRunHmm.returncode, runHmm))
                else:
                    print "sample"+str(counter), "Hmm annotation is finished."
                #
                counter += 1
            #
            if os.path.exists(os.path.join(assemblyDir, config.get('assemblyMethod') + '_sample' + str(counter)+'.contig.fa')):
                pass
            else:
                break
            #if "s" in args.f:









class FileHolder():
    """
    Store and generate the input and output files names.
    """
    def __init__(self, outputTrimForward=None, outputTrimReverse=None, outputRandomForward=None, outputRandomReverse=None,
                outputFilterForward=None, outputFilterReverse=None, outputRandomFilterForward=None, outputRandomFilterReverse=None):
        self.outputTrimForward = outputTrimForward
        self.outputTrimReverse = outputTrimReverse
        self.outputRandomForward = outputRandomForward
        self.outputRandomReverse = outputRandomReverse
        self.outputFilterForward = outputFilterForward
        self.outputFilterReverse = outputFilterReverse
        self.outputRandomFilterForward = outputRandomFilterForward
        self.outputRandomFilterReverse = outputRandomFilterReverse
        #self.outputTrimPaired = outputTrimPaired


    def getTrimAssemblyNameForward(self, assembler, filetype='fa'):
        return self.outputTrimForward + '_' + assembler + '_' + '.' + filetype
    def getTrimAssemblyNameReverse(self, assembler, filetype='fa'):
        return self.outputTrimReverse + '_' + assembler + '_' + '.' + filetype
    def getRandomAssemblyNameForward(self, assembler, filetype='fa'):
        return self.outputRandomForward + '_' + assembler + '_' + '.' + filetype
    def getRandomAssemblyNameReverse(self, assembler, filetype='fa'):
        return self.outputRandomReverse + '_' + assembler + '_' + '.' + filetype
    def getFilterAssemblyNameForward(self, assembler, filetype='fa'):
        return self.outputFilterForward + '_' + assembler + '_' + '.' + filetype
    def getFilterAssemblyNameReverse(self, assembler, filetype='fa'):
        return self.outputFilterReverse + '_' + assembler + '_' + '.' + filetype
    def getRandomFilterAssemblyNameForward(self,assembler, filetype='fa'):
        return self.outputRandomFilterForward + '_' + assembler + '_' + '.' + filetype
    def getRandomFilterAssemblyNameReverse(self,assembler, filetype='fa'):
        return self.outputRandomFilterReverse + '_' + assembler + '_' + '.' + filetype


def seqidSuffix(inputForward, inputReverse, outputFile1, outputFile2):
    f1 = open(outputFile1, "w")
    f2 = open(outputFile2, "w")

    for name, seq, qual in FastqGeneralIterator(open(inputForward)):
        if name.endswith('/1'):
            id1 = name[:-2]
            f1.write("@%s\n%s\n+\n%s\n" % (id1+' /1', seq, qual))
        elif '1:N:0' in name:
            pass
        else:
            print "Forward sequences ids patterns are not fit for the pipeline."
    for name, seq, qual in FastqGeneralIterator(open(inputReverse)):
        if name.endswith('/2'):
            id2 = name[:-2]
            f2.write("@%s\n%s\n+\n%s\n" % (id2+' /2', seq, qual))
        elif '2:N:0' in name:
            pass
        else:
            print "Reverse sequences ids patterns are not fit for the pipeline."
    f1.close()
    f2.close()

    if os.path.exists(os.path.normpath(outputFile1)):
        os.remove(inputForward)
        os.rename(outputFile1, inputForward)
    if os.path.exists(outputFile2):
        os.remove(inputReverse)
        os.rename(outputFile2, inputReverse)

#def combinePairedEndFileFastq(file1, file2, outputFile):
#    outputFile = open(outputFile, 'w')
#    f1 = SeqIO.parse(os.path.normpath(file1), "fastq")
#    f2 = SeqIO.parse(os.path.normpath(file2), "fastq")
#    reverse_ids = set()
#    paired_ids = set()
#    for record1, record2 in zip(f1, f2):
#        reverse_ids.add(record2.id)
#        name = record1.id
#        if name in reverse_ids:
#            paired_ids.add(name)
#            reverse_ids.remove(name)
            #outputFile.write("@" + (str(paired_ids)+f_suffix) +'\n' +str(record1.seq) + '\n'+ '+' + '\n'  + record1.\n" %(,record1.seq,record1)
        #if record2.id in paired_ids:
        #    SeqIO.write([record1, record2], outputFile, "fastq")
    #outputFile.close()

def combinePairedEndFileFastq(file1, file2, outputFile):
    f1 = open(file1, "r")
    f2 = open(file2, "r")
    f3 = open(outputFile, 'w')

    g1 = SeqIO.parse(f1, "fastq")
    g2 = SeqIO.parse(f2, "fastq")

    try:
        while True:
            SeqIO.write([g1.next(), g2.next()], f3, "fastq")
    except StopIteration:
        pass
    f3.close()

#def combinePairedEndFileFastq(file1, file2, outputFile1):
#    outputFile1 = open(outputFile1, 'w')
#    reverse_ids = set()
#    paired_ids = set()
#    for record2 in SeqIO.parse(os.path.normpath(file2), "fastq"):
#        if record2.id.endswith('/2'):
#            reverse_ids.add(record2.id[:-2])
#        else:
#            reverse_ids.add(record2.id)
#    record2_dict = SeqIO.index(os.path.normpath(file2), "fastq")
#    for record1 in SeqIO.parse(os.path.normpath(file1), "fastq"):
#        if record1.id.endswith('/1'):
#            name = record1.id[:-2]
#        else:
#            name = record1.id
#        if name in reverse_ids:
#            paired_ids.add(name)
#            reverse_ids.remove(name)
                ##print paired_ids
            #SeqIO.write(record1, outputFile1, "fastq")
            #if record2.id.endswith('/2'):
            #    SeqIO.write(record2_dict[name+'/2'], outputFile1,"fastq")
            #else:
            #    SeqIO.write(record2_dict[name], outputFile1, "fastq")



#def combinePairedEndFileFasta(file1, file2, outputFile):
#    outputFile = open(outputFile, 'w')
#    f1 = SeqIO.parse(os.path.normpath(file1), "fasta")
#    f2 = SeqIO.parse(os.path.normpath(file2), "fasta")
#    reverse_ids = set()
#    paired_ids = set()

    #for record1, record2 in zip(f1, f2):
    #    reverse_ids.add(record2.id[:-2])
    #    name = record1.id[:-2]
    #    if name in reverse_ids:
    #        paired_ids.add(name)
    #    if record2.id[:-2] in paired_ids:
    #        outputFile.write('>' + str(record1.id) + '\n' + str(record1.seq) + '\n' + '>' + str(record2.id) + '\n' + str(record2.seq) + '\n')
    #outputFile.close()

def combinePairedEndFileFasta(file1, file2, outputFile):
    f1 = open(file1, "r")
    f2 = open(file2, "r")
    f3 = open(outputFile, 'w')

    g1 = SeqIO.parse(f1, "fasta")
    g2 = SeqIO.parse(f2, "fasta")

    try:
        while True:
            SeqIO.write([g1.next(), g2.next()], f3, "fasta")
    except StopIteration:
        pass
    f3.close()

def splitPairedReads(inPairedFasta, outOddFasta, outEvenFasta):
    outEvenFasta = open(outEvenFasta, 'w')
    outOddFasta = open(outOddFasta, 'w')
    counter = 1
    inPairedFasta = SeqIO.parse(os.path.normpath(inPairedFasta), 'fasta')
    for record in inPairedFasta:
        #entry = '>' + str(record.id) + '\n' + str(record.seq) + '\n'
        if counter % 2 != 0:
            outOddFasta.write('>' + str(record.id) + '/1' + '\n' + str(record.seq) + '\n')
        else:
            outEvenFasta.write('>' + str(record.id) + '/2' + '\n' + str(record.seq) + '\n')
        counter += 1

    outOddFasta.close()
    outEvenFasta.close()

def fastqTofasta(inputFile, direction, outputFile):
    f1 = SeqIO.parse(os.path.normpath(inputFile), "fastq")
    f2 = open(outputFile, 'w')
    for record in f1:
        entry = '>' + str(record.id) + '/' + str(direction) + '\n' + str(record.seq) + '\n'
        f2.write(entry)


def removeN(inPutFile, outputFile, removeNlength):
    outputFile = open(outputFile, 'w')
    f1 = SeqIO.parse(os.path.normpath(inPutFile), "fasta")
    for record in f1:
        seqs = record.seq.split('N')[0]
        if len(seqs) > int(removeNlength):
            outputFile.write('>' + str(record.id) + '\n' + str(seqs) + '\n')

#def getConsensusSeqs(file1, file2, file3, file4):
    #using {} see # content
    #records = {}
    #reverse_ids = set()
    #paired_ids = set()
    #forward =
    #reverse =
    #consensusf = open(file3, 'w')
    #consensusr = open(file4, 'w')
    #for record1, record2 in zip(forward, reverse):
    #for record2 in SeqIO.parse(os.path.normpath(file2), "fasta"):
    #    if record2.id.endswith('/2/2'):
    #        reverse_ids.add(record2.id[:-4])
    #    else:
    #        reverse_ids.add(record2.id[:-2])
        #record1_id = record1.id[:-2]
        #reverse_ids.add(record2.id[:-2])
        #print record1_id
        #records[record1_id] = record1.seq
        #print record2.id[:-2]
        #if record2.id[:-2] in records:
        #if record1.id[:-2] in reverse_ids:
    #for record1 in SeqIO.parse(os.path.normpath(file1), "fasta"):
    #    if record1.id.endswith('/1/1'):
    #        namef = record1.id_new[:-4]
    #    else:
    #        namef = record1.id[:-2]
    #    print namef
    #    if namef in reverse_ids:
    #        paired_ids.add(namef)
    #        reverse_ids.remove(namef)
    #        if namef.endswith('/1/1'):
    #            consensusf.write('>' + str(record1.id[:-2]) + '\n' + str(record1.seq) + '\n')
    #        else:
    #            consensusf.write('>' + str(record1.id) + '\n' + str(record1.seq) + '\n')
    #print paired_ids

    #for record2 in SeqIO.parse(os.path.normpath(file2), "fasta"):
    #    if record2.id.endswith('/2/2'):
    #        namer = record2.id[:-4]
    #    else:
    #        namer = record2.id[:-2]
    #    print namer
        #if namer in paired_ids:
        #    if namer.endswith('/2/2'):
        #        consensusr.write('>' + str(record2.id[:-2]) + '\n' + str(record2.seq) + '\n')
        #    else:
        #        consensusr.write('>' + str(record2.id) + '\n' + str(record2.seq) + '\n')
            #print namer

def getConsensusSeqs(file1, file2, file3, file4):
    #records = {}
    reverse_ids = set()
    paired_ids = set()
    #forward =
    #reverse =
    consensusf = open(file3, 'w')
    consensusr = open(file4, 'w')
    #for record1, record2 in zip(forward, reverse):
    for record2 in SeqIO.parse(os.path.normpath(file2), "fasta"):
        reverse_ids.add(record2.id[:-2])
        #record1_id = record1.id[:-2]
        #reverse_ids.add(record2.id[:-2])
        #print record1_id
        #records[record1_id] = record1.seq
        #print record2.id[:-2]
        #if record2.id[:-2] in records:
        #if record1.id[:-2] in reverse_ids:
    for record1 in SeqIO.parse(os.path.normpath(file1), "fasta"):
        namef = record1.id[:-2]
        if namef in reverse_ids:
            paired_ids.add(namef)
            reverse_ids.remove(namef)
            consensusf.write('>' + str(record1.id) + '\n' + str(record1.seq) + '\n')
    for record2 in SeqIO.parse(os.path.normpath(file2), "fasta"):
        namer = record2.id[:-2]
        if namer in paired_ids:
            consensusr.write('>' + str(record2.id) + '\n' + str(record2.seq) + '\n')


def filterReadLen(inputFile, outputFile, readLen):
    f1 = SeqIO.parse(os.path.normpath(inputFile), "fasta")
    f2 = open(outputFile, 'w')
    for record in f1:
        if len(record.seq) >= readLen:
            # entry = '>' + str(record.id) + '\n' + str(record.seq) + '\n'
            f2.write('>' + str(record.id) + '\n' + str(record.seq) + '\n')
    f2.close()

def renameReads(inputFile, outputFile):
    f1 = SeqIO.parse(os.path.normpath(inputFile), "fasta")
    f2 = open(outputFile, 'w')
    for record in f1:
        f2.write('>' + str(record.id) + '\n' + str(record.seq) + '\n')
    f2.close()

def SplitPairedFile(inputFile,outputFile1,outputFile2):
    """Split paired fastq file to two separate files"""

    inputFile1 = open(inputFile, 'r')
    outputFile1 = open(outputFile1, 'w')
    outputFile2 = open(outputFile2, 'w')

    for record in SeqIO.parse(inputFile1, "fastq"):
        if record.id.endswith('/1'):
            SeqIO.write(record, outputFile1,"fastq")
        elif record.id.endswith('/2'):
            SeqIO.write(record, outputFile2,"fastq")
        elif ('1:N:0') in record.description:
            SeqIO.write(record, outputFile1,"fastq")
        elif ('2:N:0') in record.description:
            SeqIO.write(record, outputFile2,"fastq")

    inputFile1.close()
    outputFile1.close()
    outputFile2.close()

def noneToEmptyStr(x):
    if x is None:
        return ""
    else:
        return x


if __name__ == "__main__":
    main()
