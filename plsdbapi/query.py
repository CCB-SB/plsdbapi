import requests
import pandas
import logging
import re
from tqdm import tqdm

###########################
### variables
###########################
URL = 'https://ccb-microbe.cs.uni-saarland.de/plsdb2025/api/'

###########################
### logger
###########################

# file handler
fh = logging.FileHandler('plsdbapi.log')
fh.setLevel(logging.INFO)
# console handler
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# formatter
formatter = logging.Formatter(
    fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%y:%m:%d %H:%M:%S'
)
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# logger
logger = logging.getLogger('plsdbapi_logger')
logger.setLevel(logging.INFO)
logger.addHandler(fh)
logger.addHandler(ch)

###########################
### helper
###########################

def test_value(value, pct):
    if value < 0.0:
        return False
    if pct:
        if value > 100:
            return False
    else:
        if value > 1.0:
            return False
    return True

###########################
### API requests
###########################

def download_fasta(ids):
    """
    :param ids: List of NCBI sequence accession ids
    :return pandas dataframe with id and fasta sequence
    """
    import time

    # join ids to string
    fastas = ''
    for id in ids:
        fastas += id + ';'
    PARAMS = {'fastas': fastas}

    # send first request to start preparing fasta file
    response = requests.post(url=URL + 'fasta', data=PARAMS)
    assert response.status_code == 200, "Response status is %s" % response.status_code
    d = response.json()
    if 'job_id' in d:
        job_id = d['job_id']
    else:
        raise Exception("Wrong response")

    logger.info("START with job id: %s" % job_id)

    # use returned job id to request job status
    PARAMS = {'job_id':job_id}
    response = requests.get(url=URL + 'fasta', params=PARAMS)
    assert response.status_code == 200, "Response status is %s" % response.status_code
    
    while response.status_code == 200:

        # if json is returned job is not finished -> request staus again after some time
        if response.headers['Content-Type'] == 'application/json':
            d = response.json()
            if d['label'] == 'running':
                logger.info("preparing fasta file")
                time.sleep(5)
                response = requests.get(url=URL + 'fasta', params=PARAMS, stream=True)
            else: 
                if 'notfound' in d:
                    notfound = ', '.join(d['notfound'])
                    raise Exception("Request failed. No plasmid in PLSDB with the following ID(s): %s" % notfound)
                else:
                    raise Exception("Request failed. \"%s\" was returned." % d['label'])
        # fasta file is returned
        else:
            header = response.headers['content-disposition']
            fname = re.findall("filename=(.+)", header)[0]
            logger.info("finished")

            f = open(fname, 'wb')
            logger.info("starting file download")
            for chunk in tqdm(response.iter_content(chunk_size=1024)):
                if chunk:
                    f.write(chunk)

            logger.info("DONE")
            return True            

    return False     
    
        
def summary(ids, fasta=False):
    """
    This function processes the plasmid IDs and returns their information stored in PLSDB as .tsv and (optional) .fasta file.
    :param ids: List of plasmid NCBI sequence accession ids
    :param fasta: When set to True fasta file is generated
    :return: Pandas dataframe containing plasmid info from PLSDB plus fasta sequence if fasta=True
    """
    if not isinstance(ids, list):
        ids = ids.split()
    
    assert len(ids) > 0, "List of ids is empty."

    found = []
    notfound = []
    fastas = []
    results = []

    # search ids in PLSDB
    logger.info('start searching for plasmids')
    for id in ids:
        PARAMS = {'NUCCORE_ACC': id}
        response = requests.get(url=URL+'summary/', params=PARAMS)
        if response.status_code == 200:
            d = response.json()
            if 'searched' in d:
                # more than one accessions found -  no unique accession for id
                d['label'] = 'notfound'
                notfound.append(d)
            elif 'Metadata_annotations' in d:
                # exactly one matching accession was found for id
                d['searched'] = id
                d['label'] = 'found'
                found.append(d)
                # if fasta=True its sequence will be downloaded
                if fasta:
                    fastas.append(d['searched'])
            results.append(d)
            
        else:
            raise Exception('Http response is not OK. Response status code is %s.' % response.status_code)

    logger.info('search is finished')
    logger.info('%s of %s ids were found' % (str(len(found)), str(len(ids))))

    # download fastas for found accessions and merge dataframes
    if fasta:
        seqs = download_fasta(fastas)
        if not seqs:
            logger.info('fasta download failed')

    return results


def query_plasmid_sequence(search_type, ifile='', iseq='', seqname='sequence',
                           mash_max_v=0.1, mash_max_d=0.1, mash_min_i=0.99,
                           mash_screen_w=False, mash_dist_i=False,
                           blastn_min_i=60, blastn_min_c=90, tblastn_min_c=90):

    """
    Queries the sequence(s) from ifile or iseq according to search_type
    :param search_type: Either 'mash_screen', 'mash_dist', 'blastn' or 'tblastn'
    :param ifile: input file. For mash search the maximal file size is 500MB. For blast sequences between 100-5000bp are allowed with a maximum of 10 sequences.
    :param iseq: input sequence
    :param seqname: name for sequence iseq
    :param mash_max_v: maximum p-value for mash screen and mash dist: between 0 ad 1
    :param mash_max_d: maximum distance for mash dist: between 0 and 1
    :param mash_min_i: minimum identity for mash screen: between 0 and 1
    :param mash_screen_w: winner takes all strategy for mash screen: Hashes found in mutiple queries will be removed except for the query with highest identity. Removes redundancy from output.
    :param mash_dist_i: individually strategy for mash dist: Process each sequence individually rather than processing them together.
    :param blastn_min_i: minimum identity for blastn: between 0 and 100
    :param blastn_min_c: minimum query coverage/HSP for blastn: between 0 and 100
    :param tblastn_min_c: minimum query coverage/hSP for tblastn: between 0 and 100
    :return: Pandas dataframe containing search output and info of plasmids
    """
    import time
    from tempfile import NamedTemporaryFile

    allowed_search_types = ['mash_dist', 'mash_screen', 'blastn', 'tblastn']
    assert (search_type in allowed_search_types), 'Invalid search type.'
    assert test_value(mash_max_v, False), 'Value of mash_max_v is out of range.'
    assert test_value(mash_max_d, False), 'Value of mash_max_d is out of range.'
    assert test_value(mash_min_i, False), 'Value of mash_min_i is out of range.'
    assert test_value(blastn_min_i, True), 'Value of blastn_min_i is out of range.'
    assert test_value(blastn_min_c, True), 'Value of blastn_min_c is out of range.'
    assert test_value(tblastn_min_c, True), 'Value of tblastn_min_c is out of range.'

    # set parameters for search
    PARAMS = {'search_type': search_type,
              'mash_max_v': mash_max_v, 'mash_max_d': mash_max_d, 'mash_min_i': mash_min_i,
              'mash_screen_w': mash_screen_w, 'mash_dist_i': mash_dist_i,
              'blastn_min_i': blastn_min_i, 'blastn_min_c': blastn_min_c, 'tblastn_min_c': tblastn_min_c}

    logger.info('START sequence search')
    if not ifile == '':
        # search with ifile
        FILES = {'fasta_file': open(ifile, 'rb')}
        response = requests.post(url=URL+'sequence/', files=FILES, data=PARAMS)
    elif not iseq == '':
        # search with sequence string
        file = NamedTemporaryFile()
        file.write(str.encode('>' + seqname + '\n'))
        file.write(str.encode(iseq))
        file.read()
        FILES = {'fasta_file': open(file.name, 'rb')}
        response = requests.post(url=URL+'sequence/', files=FILES, data=PARAMS)
        file.close()
    else:
        raise Exception("Either ifile or iseq is required.")

    errors = ['failed', 'error', 'no job id']
    if response.status_code == 200:
        data = response.json()
        while data['label'] != 'finished':
            if data['label'] in errors:
                # error occurred on server side
                raise Exception("Returned response: " + data['label'] + ' Please check again your input file and the input specifications as well as the chosen search_type.')
            if data['label'] == 'invalid post':
                raise Exception("Submission failed. " + data['error'])
            # wait before asking if finished
            logger.info("search job is running")
            time.sleep(5)
            job_id = data['job_id']
            response = requests.get(url=URL+'sequence/', params={'job_id': job_id})
            data = response.json()

        # server side is done
        df = pandas.read_json(data['results'], orient='split')
        logger.info('DONE')
        return df
    else:
        raise Exception('Http response is not OK. Response status code is %s.' % response.status_code)


def filter_nuccore(fasta=False, NUCCORE_Source="", NUCCORE_Topology="", NUCCORE_has_identical="",
                   AMR_genes="", BGC_types=""):
    """
    Filter the PLSDB plasmids
    :param fasta: if True generates fasta file of filtered plasmids
    :param NUCCORE_Source: Search for plasmid with selected source. Available values : RefSeq, INSDC
    :param NUCCORE_Source: Report all plasmids with selected 'topology'. Available values : circular, linear
    :param NUCCORE_has_identical: Report all plasmids with detected Identical sequences. Available values : yes, empty
    :param AMR_genes: Search for plasmid that contains the select AMR gene. Str with comma separated values
    :param BGC_types: Search for plasmid that contains the select BGC type. Str with comma separated values
    :return:
    """

    local = locals()
    PARAMS = dict()
    for var in local:
        if var == 'fasta': continue
        if not local[var]: continue
        PARAMS[var] = local[var]
    
    print(PARAMS)
    logger.info('START search for plasmids')
    response = requests.get(url=URL + 'filter_nuccore', params=PARAMS)

    if response.status_code == 200:
        data = response.json()

        if len(data) == 0:
            logger.info('No plasmid found for filter %s.' % PARAMS)
            return

        if fasta:
            fastas = data['NUCCORE_ACC']
            seqs = download_fasta(fastas)
            if not seqs:
                logger.info("fasta download failed")

        logger.info("DONE")       
        return data
    else:
        raise Exception('Http response is not OK. Response status code is %s.' % response.status_code)

def filter_biosample(fasta=False, BIOSAMPLE_UID="", LOCATION_name="", ECOSYSTEM_tags="",
                   ECOSYSTEM_taxid_name="", ECOSYSTEM_taxid="", DISEASE_ontid_name="",
                   DISEASE_ontid=""):
    """
    Filter the PLSDB plasmids
    :param fasta: if True generates fasta file of filtered plasmids
    :param BIOSAMPLE_UID: Report all plasmids found with the selected BIOSAMPLE_ACC
    :param LOCATION_name: Report all plasmids found at the given location, e.g. Germany
    :param ECOSYSTEM_tags: Report all plasmids with entered isolation_source, e.g. fecal.
    :param ECOSYSTEM_taxid_name: Report all plasmids from host-associated ecosystems using NCBI Taxonomy taxids, e.g Homo sapiens = 9606
    :param DISEASE_ontid_name: Report all plasmids from host-associated ecosystems with associated diseases, e.g Aspiration pneumonia
    :param DISEASE_ontid: Report all plasmids from host-associated ecosystems with associated diseases by Disease/Symptom ID, e.g Aspiration pneumonia = 'DOID:0050152'
    :return:
    """

    local = locals()
    PARAMS = dict()
    for var in local:
        if var == 'fasta': continue
        if not local[var]: continue
        PARAMS[var] = local[var]
    
    print(PARAMS)
    logger.info('START search for plasmids')
    response = requests.get(url=URL + 'filter_biosample', params=PARAMS)

    if response.status_code == 200:
        data = response.json()

        if len(data) == 0:
            logger.info('No plasmid found for filter %s.' % PARAMS)
            return

        if fasta:
            fastas = data['NUCCORE_ACC']
            seqs = download_fasta(fastas)
            if not seqs:
                logger.info("fasta download failed")

        logger.info("DONE")       
        return data
    else:
        raise Exception('Http response is not OK. Response status code is %s.' % response.status_code)

def filter_taxonomy(fasta=False, TAXONOMY_strain="", TAXONOMY_strain_id="", TAXONOMY_species="",
                   TAXONOMY_species_id="", TAXONOMY_genus="", TAXONOMY_genus_id="",
                   TAXONOMY_family="", TAXONOMY_family_id="", TAXONOMY_order="",
                   TAXONOMY_order_id="", TAXONOMY_class="", TAXONOMY_class_id="",
                   TAXONOMY_phylum="", TAXONOMY_phylum_id="", TAXONOMY_superkingdom="",
                   TAXONOMY_superkingdom_id=""
                     ):
    """
    Filter the PLSDB plasmids
    :param fasta: if True generates fasta file of filtered plasmids
    :param TAXONOMY_strain: Report all plasmids according to given host strain (NCBI Taxonomy name), e.g. Escherichia coli B7A
    :param TAXONOMY_strain_id: Report all plasmids according to given host strain ID (NCBI Taxonomy taxid), e.g. Escherichia coli B7A = 340184
    :param TAXONOMY_species: Report all plasmids according to given host species (NCBI Taxonomy name), e.g. Escherichia coli.
    :param TAXONOMY_species_id: Report all plasmids according to given host species ID (NCBI Taxonomy taxid), e.g. Escherichia coli = 562.
    :param TAXONOMY_genus_id: Report all plasmids according to given host genus ID (NCBI Taxonomy taxid), e.g. Escherichia = 561.
    :param TAXONOMY_family: Report all plasmids according to given host family (NCBI Taxonomy name), e.g. Enterobacteriaceae.
    :param TAXONOMY_family_id: Report all plasmids according to given host family ID (NCBI Taxonomy taxid), e.g. Enterobacteriaceae = 543.
    :param TAXONOMY_order: Report all plasmids according to given host order (NCBI Taxonomy name), e.g. Enterobacterales.
    :param TAXONOMY_order_id: Report all plasmids according to given host order ID (NCBI Taxonomy taxid), e.g. Enterobacterales = 91347.
    :param TAXONOMY_class: Report all plasmids according to given host class (NCBI Taxonomy name), e.g. Gammaproteobacteria.
    :param TAXONOMY_class_id: Report all plasmids according to given host class ID (NCBI Taxonomy taxid), e.g. Gammaproteobacteria = 1236.
    :param TAXONOMY_phylum: Report all plasmids according to given host phylum (NCBI Taxonomy name), e.g. Pseudomonadota.
    :param TAXONOMY_phylum_id: Report all plasmids according to given host phylum ID (NCBI Taxonomy taxid), e.g. Pseudomonadota = 1224.
    :param TAXONOMY_superkingdom: Report all plasmids according to given host superkigdom (NCBI Taxonomy name), e.g. Bacteria.
    :param TAXONOMY_superkingdom_id: Report all plasmids according to given host superkingdom ID (NCBI Taxonomy taxid), e.g. Bacteria = 2.
    :return:
    """

    local = locals()
    PARAMS = dict()
    for var in local:
        if var == 'fasta': continue
        if not local[var]: continue
        PARAMS[var] = local[var]
    
    print(PARAMS)
    logger.info('START search for plasmids')
    response = requests.get(url=URL + 'filter_taxonomy', params=PARAMS)

    if response.status_code == 200:
        data = response.json()

        if len(data) == 0:
            logger.info('No plasmid found for filter %s.' % PARAMS)
            return

        logger.info("DONE")       
        return data
    else:
        raise Exception('Http response is not OK. Response status code is %s.' % response.status_code)


    
if __name__ == "__main__":

    import query
    # query.summary(['NZ_CP031107.1'], fasta=True)
    query.filter_nuccore(NUCCORE_Source="RefSeq", NUCCORE_Topology="circular", NUCCORE_has_identical='yes', AMR_genes="espP,toxB,ehxA,katP")
    query.filter_biosample(ECOSYSTEM_tags="fecal", ECOSYSTEM_taxid=9606, DISEASE_ontid_name='Aspiration pneumonia')
    query.filter_taxonomy(TAXONOMY_strain_id=340184)






