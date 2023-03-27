from Bio import Entrez

Entrez.email = 'kimura.kazuki44@gmail.com'
path_23andme = 'C:/Users/kimur/Downloads/SNP.txt'
path_VEP = 'C:/Users/kimur/Downloads/VEP.txt'

def ToVEP(path_23andme, path_VEP):
    f = open(path_23andme, 'r')
    g = open(path_VEP, 'w')
    line = f.readline().splitlines()[0].split('\t')
    while line:
        if line[0][0] == '#':
            line = f.readline().split('\t')
            continue
        chr = line[1]
        if chr == 'MT':
            continue
        pos = line[2]
        ref = line[3][0]
        if chr == 'Y':
            alt == line[3][0]
        else:
            alt = line[3][1]
        if ref == '-':
            ref = 'A'
            alt = 'A'
        VEPline = chr + ' ' + pos + ' ' + pos + ' ' + ref + '/' + alt + ' + \n'
        g.write(VEPline)
        line = f.readline().split('\t')

def SearchSNPwithVEP(path_23andme, path_SNPs): ### path_SNPs : temporally assume as data table downloaded from NCBI Genome Data Viewer
    f = open(path_23andme, 'r')
    g = open(path_SNPs, 'r')
    line = f.readline().split('\t')
    RSIDs = []
    while line != ['']:
        if line[0][0] == '#':
            line = f.readline().split('\t')
            continue
        RSID = line[0]
        chr = line[1]
        if chr == 'X':
            RSIDs.append(RSID)
        line = f.readline().split('\t')
    rawline = g.readline()
    SNPs = []
    while rawline != '':
        line = rawline.split('\t')
        if line[0][0] == '#':
            rawline = g.readline()
            continue
        RSID = line[0]
        if RSID in RSIDs:
            SNPs.append(RSID)
            print(rawline)
        rawline = g.readline()
    print('SNP count = ', len(SNPs))

def SearchSNP(path_23andme, gene):
    term = gene + '[GENE]'
    handle = Entrez.esearch(db = 'SNP', term = term, retmax = '4000')
    record = Entrez.read(handle)
    f = open(path_23andme, 'r')
    line = f.readline()
    count = 0
    while line:
        rsid = line.split('\t')[0][2:]
        if rsid in record['IdList']:
            print(rsid)
        count += 1
        if count % 10000 == 0:
            print(count, rsid) 
        line = f.readline()


