#!/home/hkim/anaconda3/bin/python

import os, sys, pickle, regex
import subprocess as sp
import multiprocessing as mp

sBASE_DIR       = '/extdata1/Jinman/pe_screening'
sFLASH          = '/home/hkim/bin/FLASH-1.2.11-Linux-x86_64/flash'
nLINE_CNT_LIMIT = 500000 ## split files, must be multiples of 4. file size < 1G recommendation:40000, size > 1G recommendation:400000

class cPEData:pass

def load_NGS_files (sInFile):
    InFile       = open(sInFile, 'r')
    list_sOutput = [sReadLine.strip('\n') for sReadLine in InFile if not sReadLine.startswith('#')]
    InFile.close()

    dict_sOutput = {}
    for sFile in list_sOutput:

        sSample = '_'.join(sFile.split('_')[:3])

        if sSample not in dict_sOutput:
            dict_sOutput[sSample] = []
        dict_sOutput[sSample].append(sFile)
    #loop END: sFile

    return dict_sOutput
#def END: load_NGS_files


def load_PE_input (sInFile):

    dict_sOutput = {}
    InFile       = open(sInFile, 'r')
    list_sTest   = []
    for sReadLine in InFile:
        ## File Format ##
        ## Target#  | Barcode | WT_Target | Edit_Target
        ## 181      | TTT.... | CTGCC..   | CTGCC...

        if sReadLine.startswith('Target'): continue ## SKIP HEADER

        list_sColumn = sReadLine.strip('\n').split('\t')

        cPE              = cPEData()
        cPE.nTargetNo    = int(list_sColumn[0])
        cPE.sBarcode     = list_sColumn[1]
        cPE.sWTSeq       = list_sColumn[2].upper()
        cPE.sAltSeq      = list_sColumn[3].upper()
        sKey             = cPE.sBarcode

        list_sTest.append(cPE.sBarcode)

        if sKey not in dict_sOutput:
            dict_sOutput[sKey] = ''
        dict_sOutput[sKey] = cPE
    #loop END:
    '''
    print('list_sTest', len(list(set(list_sTest))))

    list_sDup  = [item for item, count in collections.Counter(list_sTest).items() if count > 1]

    print(len(list_sDup))
    for sDup in list_sDup:
        print(sDup)
    '''
    return dict_sOutput
#def END: load_PE_input


def run_FLASH (sAnalysis, sWorkDir, dict_sFiles, bTestRun):

    sInDir   = '%s/raw'    % sWorkDir
    sLogDir  = '%s/log/Jinman.RunFLASH.%s'  % (sBASE_DIR, sAnalysis)
    os.makedirs(sLogDir, exist_ok=True)

    list_sProjects = []
    for sSample in dict_sFiles:

        sOutDir    = '%s/flash/%s'      % (sWorkDir, sSample)
        os.makedirs(sOutDir, exist_ok=True)

        sLogFile   = '%s/%s.fastq.log'  % (sLogDir, sSample)
        sInFile1   = '%s/%s_1.fastq.gz' % (sInDir, sSample)
        sInFile2   = '%s/%s_2.fastq.gz' % (sInDir, sSample)

        print(sInFile1, sInFile2)
        assert os.path.isfile(sInFile1) and os.path.isfile(sInFile2)

        sScript      = '%s %s %s ' % (sFLASH, sInFile1, sInFile2)
        sScript     += '-M 400 '   # max overlap
        sScript     += '-m 5 '     # min overlap
        sScript     += '-O '       # allow "outies"  Read ends overlap
                                   # Read 1: <-----------
                                   # Read 2:       ------------>

        sScript     += '-d %s '    % sOutDir
        sScript     += '-o %s '    % sSample

        sCmd         = '%s | tee %s ;'                       % (sScript, sLogFile)
        sCmd        += 'rm -rf %s/*hist* %s/*notCombined* ;' % (sOutDir, sOutDir)

        sFileName    = '%s.extendedFrags' % sSample
        list_sProjects.append([sSample, sFileName])

        if bTestRun: print(sCmd)
        else: os.system(sCmd)
    #loop END: sSample

    #VS Check
    if not list_sProjects: sys.exit('Empty List : run_FLASH : list_sProjects= %s' % len(list_sProjects))

    sOutFile = '%s/Project_list_%s.txt' % (sWorkDir, sAnalysis)
    OutFile  = open(sOutFile, 'w')
    for sProject, sFileName in list_sProjects:
        sOut = '%s\t%s\n' % (sProject, sFileName)
        OutFile.write(sOut)
    #loop END: sProject, sFileName
#def END: run_FLASH


def split_fq_file (sWorkDir, sFastqTag, bTestRun):

    sOutDir    = '%s/split'                  % sWorkDir
    os.makedirs(sOutDir,exist_ok=True)

    sInFile    = '%s/%s.fastq'               % (sWorkDir, sFastqTag)
    print(sInFile)
    assert os.path.isfile(sInFile)

    sOutTag    = '%s/%s_fastq'               % (sOutDir, sFastqTag)

    sScript     = 'split --verbose '         # For Logging Purposes
    sScript    += '-l %s '                   % nLINE_CNT_LIMIT
    sScript    += '-a 3 '                    # Number of suffice places e.g. 001.fq

    sScript    += '--numeric-suffixes=1 '    # Start with number 1
    sScript    += '--additional-suffix=.fq ' # Add suffix .fq'
    sScript    += '%s %s_'                   % (sInFile, sOutTag)

    if bTestRun: print(sScript)
    else:
        os.system(sScript)
        list_sFiles = os.listdir(sOutDir)
        sOutFile    = '%s/%s.split.txt'       % (sWorkDir, sFastqTag)
        OutFile     = open(sOutFile, 'w')
        for sFile in list_sFiles:
            sOut = '%s\n' % sFile
            OutFile.write(sOut)
        #loop END: sFile
        OutFile.close()
    #if END:
#def END: split_fq_file


def get_split_list (sWorkDir, sFastqTag):
    sInFile = '%s/%s.split.txt' % (sWorkDir, sFastqTag)
    InFile = open(sInFile, 'r')
    list_sSplits = [sReadLine.strip('\n') for sReadLine in InFile if not sReadLine.startswith('#')]
    InFile.close()
    return list_sSplits
#def END: get_split_list


def get_line_cnt (sWorkDir, sFastqTag):

    sInFile  = '%s/%s.fastq'       % (sWorkDir, sFastqTag)
    sOutFile = '%s/%s.linecnt.txt' % (sWorkDir, sFastqTag)

    if os.path.isfile(sOutFile):
        InFile   = open(sOutFile, 'r')
        nLineCnt = int([sReadLine.strip('\n').split(' ') for sReadLine in InFile][0][0])
        print(nLineCnt)
        return nLineCnt
    else:
        sCmd     = 'wc -l %s > %s'     % (sInFile, sOutFile)
        os.system(sCmd)
        InFile   = open(sOutFile, 'r')
        nLineCnt = int([sReadLine.strip('\n').split(' ') for sReadLine in InFile][0][0])
        return nLineCnt
#def END: get_line_cnt


def load_read_data (dict_sOutFreq, dict_sOutRead, sInFile):

    InFile       = open(sInFile, 'r')
    for sReadLine in InFile:

        list_sColumn      = sReadLine.strip('\n').split('\t')
        sBarcode          = list_sColumn[0]
        list_WT_reads     = [sReadID for sReadID in list_sColumn[1].split(',') if sReadID]
        list_edited_reads = [sReadID for sReadID in list_sColumn[2].split(',') if sReadID]
        list_other_reads  = [sReadID for sReadID in list_sColumn[3].split(',') if sReadID]

        dict_sOutFreq[sBarcode][0] += len(list_WT_reads)
        dict_sOutFreq[sBarcode][1] += len(list_edited_reads)
        dict_sOutFreq[sBarcode][2] += len(list_other_reads)

        dict_sOutRead[sBarcode][0] += list_WT_reads
        dict_sOutRead[sBarcode][1] += list_edited_reads
        dict_sOutRead[sBarcode][2] += list_other_reads
    #loop END: sReadLine
    InFile.close()
#def END: load_freq_data


def run_seqkit_grep (sParameters):
    sReadIDFile = sParameters[0]
    sInFastq    = sParameters[1]
    sOutFile    = sParameters[2]
    sScript     = 'seqkit grep -f %s %s > %s' % (sReadIDFile, sInFastq, sOutFile)
    os.system(sScript)
#def END: run_seqkit_grep


def mp_sort_by_barcode (nCores, sSample, sInDir, sOutDir, sBarcodeFile, list_sSplits):

    list_sParameters = []
    for sSplitFile in list_sSplits:
        # HKK_191230_1.extendedFrags_|fastq_01|.fq
        sSplitTag = '_'.join(sSplitFile.split('.')[1].split('_')[-2:])
        sInFile   = '%s/split/%s' % (sInDir, sSplitFile)
        assert os.path.isfile(sInFile)

        sTempDir = '%s/temp/%s/%s' % (sOutDir, sSample, sSplitTag)
        os.makedirs(sTempDir, exist_ok=True)
        list_sParameters.append([sSplitTag, sInFile, sTempDir, sBarcodeFile])
    #loop END: sSplitFile

    p = mp.Pool(nCores)
    p.map_async(sort_by_barcode, list_sParameters).get()
    p.map_async(determine_output, list_sParameters).get()
#def END: mp_sort_by_barcode


def sort_by_barcode (list_sParameters):

    print('Processing %s' % list_sParameters[1])

    sSplitTag      = list_sParameters[0]
    sInFastq       = list_sParameters[1]
    sTempOut       = list_sParameters[2]
    sBarcodeFile   = list_sParameters[3]
    sRE            = '[T]{7}'
    nBarcode3Cut   = 3
    nBarcode5Ext   = 18
    dict_sBarcodes = load_PE_input(sBarcodeFile)

    dict_sOutput   = {}
    InFile = open(sInFastq, 'r')
    for i, sReadLine in enumerate(InFile):

        if i % 4 == 0: sReadID = sReadLine.replace('\n', '')
        if i % 4 != 1: continue

        sNGSSeq = sReadLine.replace('\n', '').upper()

        for sReIndex in regex.finditer(sRE, sNGSSeq, overlapped=True):
            nIndexStart = sReIndex.start()
            nIndexEnd   = sReIndex.end()
            sBarcode    = sNGSSeq[nIndexStart+nBarcode3Cut:nIndexEnd+nBarcode5Ext]

            if nIndexStart > (len(sNGSSeq) / 2): continue # SKIP barcode in back of read

            ### Skip Non-barcodes ###
            try: cPE = dict_sBarcodes[sBarcode]
            except KeyError: continue
            #########################

            if sBarcode not in dict_sOutput:
                dict_sOutput[sBarcode] = []
            dict_sOutput[sBarcode].append([sReadID, sNGSSeq])

        #loop END: i, sReadLine
    #loop END: cPE
    InFile.close()
    ## Pickle Out ##
    sOutFile = '%s/%s.data' % (sTempOut, sSplitTag)
    OutFile = open(sOutFile, 'wb')
    pickle.dump(dict_sOutput, OutFile)
    OutFile.close()
#def END: check_barcode


def determine_output (list_sParameters):

    print('Processing %s' % list_sParameters[2])

    sSplitTag      = list_sParameters[0]
    sInFastq       = list_sParameters[1]
    sTempOut       = list_sParameters[2]
    sBarcodeFile   = list_sParameters[3]
    dict_cPE       = load_PE_input(sBarcodeFile)

    ## Pickle Load ##
    sInFile        = '%s/%s.data' % (sTempOut, sSplitTag)
    InFile         = open(sInFile, 'rb')
    dict_sBarcodes = pickle.load(InFile)
    InFile.close()
    print('%s dict_sBarcodes =%s' % (sSplitTag, len(dict_sBarcodes)))

    dict_sOutput   = {}
    for sBarcode in dict_sBarcodes:

        cPE      = dict_cPE[sBarcode]

        if sBarcode not in dict_sOutput:
            dict_sOutput[sBarcode] = {'WT': [], 'Alt': [], 'Other': []}

        nWTSize  = len(cPE.sWTSeq)
        nAltSize = len(cPE.sAltSeq)

        for sReadID, sNGSSeq in dict_sBarcodes[sBarcode]:

            nBarcodeS      = sNGSSeq.find(sBarcode)
            nBarcodeE      = nBarcodeS + len(sBarcode)
            sWTSeqCheck    = sNGSSeq[nBarcodeE:nBarcodeE+nWTSize]
            sAltSeqCheck   = sNGSSeq[nBarcodeE:nBarcodeE+nAltSize]

            if sWTSeqCheck == cPE.sWTSeq:
                dict_sOutput[cPE.sBarcode]['WT'].append(sReadID)

            elif sAltSeqCheck == cPE.sAltSeq:
                dict_sOutput[cPE.sBarcode]['Alt'].append(sReadID)

            elif sWTSeqCheck != cPE.sWTSeq and sAltSeqCheck != cPE.sAltSeq:
                dict_sOutput[cPE.sBarcode]['Other'].append(sReadID)
            #if END:
        #loop END: sReadID, sNGSSeq
    #loop END: sBarcode
    list_sKeys = ['WT', 'Alt', 'Other']

    sOutFile   = '%s/%s.reads.txt' % (sTempOut, sSplitTag)
    OutFile    = open(sOutFile, 'w')

    for sBarcode in dict_sOutput:
        sOut = '%s\t%s\n' % (sBarcode, '\t'.join([','.join(dict_sOutput[sBarcode][sType]) for sType in list_sKeys]))
        OutFile.write(sOut)
    # loop END: sBarcode
    OutFile.close()
#def END: determine_output_WTandEdited


def combine_output (sSample, sOutDir, sBarcodeFile, list_sSplits):

    dict_cPE      = load_PE_input(sBarcodeFile)
    sCombineOut   = '%s/combined' % sOutDir
    sOthersOut    = '%s/others'   % sOutDir
    os.makedirs(sCombineOut, exist_ok=True)
    os.makedirs(sOthersOut, exist_ok=True)

    dict_sOutFreq = {sBarcode: [0, 0, 0]    for sBarcode in dict_cPE}
    dict_sOutRead = {sBarcode: [[], [], []] for sBarcode in dict_cPE}
    nTotal        = len(list_sSplits)
    for i, sSplitFile in enumerate(list_sSplits):

        print('%s/%s -- %s' % ((i+1), nTotal, sSplitFile))

        sSplitTag  = '_'.join(sSplitFile.split('.')[1].split('_')[-2:])
        sTempOut   = '%s/temp/%s/%s'    % (sOutDir, sSample, sSplitTag)
        sInFile    = '%s/%s.reads.txt'  % (sTempOut, sSplitTag)
        assert os.path.isfile(sInFile)
        load_read_data (dict_sOutFreq, dict_sOutRead, sInFile)
    #loop END: sSplitFile

    sHeader   = '%s\t%s\t%s\t%s\n' % ('Barcode', 'WT', 'Alt', 'Other')
    sOutFile  = '%s/%s.combinedFreq.txt'    % (sCombineOut, sSample)
    OutFile   = open(sOutFile, 'w')
    OutFile.write(sHeader)
    for sBarcode in dict_cPE:
        nWT, nEdited, nOther = dict_sOutFreq[sBarcode]
        sOut = '%s\t%s\t%s\t%s\n' % (sBarcode, nWT, nEdited, nOther)
        OutFile.write(sOut)
    #loop END: sBarcode
    OutFile.close()
    sOutFile2 = '%s/%s.combinedReadID.data' % (sCombineOut, sSample)
    OutFile2 = open(sOutFile2, 'wb')
    pickle.dump(dict_sOutRead, OutFile2)
    OutFile2.close()

    ## Others Out ##
    sOutDir    = '%s/%s' % (sOthersOut, sSample)
    os.makedirs(sOutDir, exist_ok=True)
    sOutFile3  = '%s/%s.readIDs.txt' % (sOthersOut, sSample)
    OutFile3   = open(sOutFile3, 'w')

    for sBarcode in dict_cPE:
        list_sOther   = dict_sOutRead[sBarcode][2]
        sOutFile      = '%s/%s.readIDs.txt' % (sOutDir, sBarcode)
        OutFile       = open(sOutFile, 'w')
        for sReadID in list_sOther:

            sReadName = sReadID.replace('@','').split(' ')[0]
            sOut      = '%s\n' % sReadName
            OutFile.write(sOut)
        #loop END: sReadID
        OutFile.close()

        sOut = '%s\t%s\n' % (sBarcode, ','.join(list_sOther))
        OutFile3.write(sOut)
    #loop END: sBarcode
    OutFile3.close()
#def END: combine_output


def gzip_fastq (sInDir, sFastqTag, bTestRun):

    sScript = 'pigz --fast -c %s/%s.fastq > %s/%s.fastq.gz' % (sInDir, sFastqTag, sInDir, sFastqTag)
    if bTestRun: print(sScript)
    else:        os.system(sScript)
#def END: gzip_fastq


def extract_reads_by_readID_for_indelsearcher (sSample, sInDir, sFastqTag, sOutDir, nCores, nTop):

    sInFastq       = '%s/%s.fastq.gz'       % (sInDir, sFastqTag)
    sOthersDir     = '%s/others'            % sOutDir
    sInFile        = '%s/%s.readIDs.txt'    % (sOthersDir, sSample)
    InFile         = open(sInFile, 'r')
    dict_sReadData = {}
    for sReadLine in InFile:

        list_sColumn = sReadLine.strip('\n').split('\t')

        sBarcode     = list_sColumn[0]
        list_sReadID = list_sColumn[1].split(',')

        if sBarcode not in dict_sReadData:
            dict_sReadData[sBarcode] = ''
        dict_sReadData[sBarcode] = list_sReadID
    #loop END: sReadLine
    InFile.close()

    list_sBasicStats  = [[sBarcode, len(dict_sReadData[sBarcode])] for sBarcode in dict_sReadData]
    list_sBasicStats  = sorted(list_sBasicStats, key=lambda e:e[1], reverse=True)

    sIndelTemp = '%s/%s_forIndel'       % (sOthersDir, sSample)
    os.makedirs(sIndelTemp, exist_ok=True)

    sOutFile_fulllist = '%s/%s_full.txt' % (sOthersDir, sSample)
    OutFile           = open(sOutFile_fulllist, 'w')
    list_sParameters  = []
    for sBarcode, nReadCnt in list_sBasicStats:

        sOut = '%s\t%s\n' % (sBarcode, nReadCnt)
        OutFile.write(sOut)

        sReadIDFile = '%s/%s/%s.readIDs.txt'   % (sOthersDir, sSample, sBarcode)
        sOutFile    = '%s/%s.seqs.txt'         % (sIndelTemp, sBarcode)
        list_sParameters.append([sReadIDFile, sInFastq, sOutFile])
    #loop END: sBarcode, nReadCnt
    OutFile.close()

    p = mp.Pool(nCores)
    p.map_async(run_seqkit_grep, list_sParameters).get()

    sIndelOut  = '%s/forIndelSearcher/%s-PE'  % (sOutDir, sSample)
    os.makedirs(sIndelOut, exist_ok=True)

    sOutFastq = '%s/%s.others.fastq' % (sIndelOut, sSample)
    ## Partially Cat files into Final File ##
    nBins       = 100
    list_sFiles = [sData[2] for sData in list_sParameters]
    nTotalJobs  = len(list_sFiles)
    list_nIndexKey = [[i + 1, int(nTotalJobs * (i) / nBins), int(nTotalJobs * (i + 1) / nBins)] for i in
                      range(nBins)]
    for nIndex, nFileS, nFileE in list_nIndexKey:
        list_sInFiles = list_sFiles[nFileS:nFileE]
        sOutFile      = '%s/%s.part%s.fastq' % (sIndelOut, sSample, nIndex)
        print(sOutFile)

        sp.call('cat %s >  %s' % (' '.join(list_sInFiles), sOutFile), shell=True)
    #loop END: nIndex, nFileS, nFileE
    sCmd      = 'cat %s/*part*.fastq > %s'    % (sIndelOut, sOutFastq)
    sp.call(sCmd, shell=True)

    print('%s-PE\t%s.others' % (sSample, sSample))
#def END: extract_reads_by_readID_for_indelsearcher


def combined_output_fastq (sSample, sInDir, sFastqTag, sOutDir, nCores, nTop):

    sInFastq       = '%s/%s.fastq.gz'       % (sInDir, sFastqTag)
    sOthersDir     = '%s/others'            % sOutDir
    sInFile        = '%s/%s.readIDs.txt'    % (sOthersDir, sSample)
    InFile         = open(sInFile, 'r')
    dict_sReadData = {}
    for sReadLine in InFile:

        list_sColumn = sReadLine.strip('\n').split('\t')

        sBarcode     = list_sColumn[0]
        list_sReadID = list_sColumn[1].split(',')

        if sBarcode not in dict_sReadData:
            dict_sReadData[sBarcode] = ''
        dict_sReadData[sBarcode] = list_sReadID
    #loop END: sReadLine
    InFile.close()

    list_sBasicStats  = [[sBarcode, len(dict_sReadData[sBarcode])] for sBarcode in dict_sReadData]
    list_sBasicStats  = sorted(list_sBasicStats, key=lambda e:e[1], reverse=True)

    sIndelTemp = '%s/%s_forIndel'       % (sOthersDir, sSample)
    os.makedirs(sIndelTemp, exist_ok=True)

    list_sParameters  = []
    for sBarcode, nReadCnt in list_sBasicStats:

        sOut = '%s\t%s\n' % (sBarcode, nReadCnt)

        sReadIDFile = '%s/%s/%s.readIDs.txt'   % (sOthersDir, sSample, sBarcode)
        sOutFile    = '%s/%s.seqs.txt'         % (sIndelTemp, sBarcode)
        list_sParameters.append([sReadIDFile, sInFastq, sOutFile])
    #loop END: sBarcode, nReadCnt

    sIndelOut  = '%s/forIndelSearcher/%s-PE'  % (sOutDir, sSample)
    os.makedirs(sIndelOut, exist_ok=True)

    sOutFastq = '%s/%s.others.fastq' % (sIndelOut, sSample)
    ## Partially Cat files into Final File ##
    nBins       = 100
    list_sFiles = [sData[2] for sData in list_sParameters]
    nTotalJobs  = len(list_sFiles)
    list_nIndexKey = [[i + 1, int(nTotalJobs * (i) / nBins), int(nTotalJobs * (i + 1) / nBins)] for i in
                      range(nBins)]
    for nIndex, nFileS, nFileE in list_nIndexKey:
        list_sInFiles = list_sFiles[nFileS:nFileE]
        sOutFile      = '%s/%s.part%s.fastq' % (sIndelOut, sSample, nIndex)
        print(sOutFile)

        sp.call('cat %s >  %s' % (' '.join(list_sInFiles), sOutFile), shell=True)
    #loop END: nIndex, nFileS, nFileE
    sCmd      = 'cat %s/*part*.fastq > %s'    % (sIndelOut, sOutFastq)
    sp.call(sCmd, shell=True)

    print('%s-PE\t%s.others' % (sSample, sSample))
#def END: extract_reads_by_readID_for_indelsearcher


def main():

    sAnalysis    = 'HKK_191230'
    sQueue       = 'workq'
    bTestRun     = False
    sTimeLimit   = '250:00:00'
    nCores       = 20

    sInputDir    = '%s/input'     % sBASE_DIR
    sOutputDir   = '%s/output/%s' % (sBASE_DIR, sAnalysis)
    sBarcodeFile = '%s/191216_Barcode_Targets.txt' % sInputDir

    ## Run FLASH ##
    sDataDir     = '%s/%s'       % (sInputDir, sAnalysis)
    sFileList    = '%s/%s.txt'   % (sInputDir, sAnalysis)
    dict_sFiles  = load_NGS_files (sFileList)
    run_FLASH (sAnalysis, sDataDir, dict_sFiles, bTestRun)

    ## Run Analysis ##
    list_sSamples = list(dict_sFiles.keys())
    for sSample in list_sSamples:

        print('Processing %s' % sSample)
        sInDir    = '%s/%s/flash/%s'   % (sInputDir, sAnalysis, sSample)
        sFastqTag = '%s.extendedFrags' % sSample
        split_fq_file(sInDir, sFastqTag, bTestRun)
        nLineCnt  = get_line_cnt (sInDir, sFastqTag)
        list_sSplitFile = get_split_list (sInDir, sFastqTag)
        mp_sort_by_barcode       (nCores, sSample, sInDir, sOutputDir, sBarcodeFile, list_sSplitFile)
        combine_output           (sSample, sOutputDir, sBarcodeFile, list_sSplitFile)

        ## Examining Top "Other" Reads ##
        nTop = 50
        gzip_fastq (sInDir, sFastqTag, bTestRun)
        extract_reads_by_readID_for_indelsearcher (sSample, sInDir, sFastqTag, sOutputDir, nCores, nTop)
        combined_output_fastq                     (sSample, sInDir, sFastqTag, sOutputDir, nCores, nTop)

        #analyze_top_barcodes (sSample, sBarcodeFile, sOutputDir, nTop, nLineCnt)
    #loop END: sSample

#def END: main
