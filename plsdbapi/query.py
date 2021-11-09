import requests
import pandas
import logging
import re
from tqdm import tqdm

###########################
### variables
###########################
URL = 'https://www.ccb.uni-saarland.de/plsdb/plasmids/api/'

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
    
        
def query_plasmid_id(ids, fasta=False):
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

    # search ids in PLSDB
    logger.info('start searching for plasmids')
    for id in ids:
        PARAMS = {'plasmid_id': id}
        response = requests.get(url=URL+'plasmid/', params=PARAMS)
        if response.status_code == 200:
            d = response.json()
            if 'found' in d:
                # exactly one matching accession was found for id
                d['found']['searched'] = id
                d['found']['label'] = 'found'
                found.append(d['found'])
                # if fasta=True its sequence will be downloaded
                if fasta:
                    fastas.append(d['found']['searched'])
            elif 'searched' in d:
                # more than one accessions found -  no unique accession for id
                d['label'] = 'notfound'
                notfound.append(d)
            else:
                raise Exception("Unexpected response for plasmid %s: %s" % (id, d))
            
        else:
            raise Exception('Http response is not OK. Response status code is %s.' % response.status_code)

    logger.info('search is finished')
    logger.info('%s of %s ids were found' % (str(len(found)), str(len(ids))))

    # convert arrays to pd dataframe and merge results
    df1 = pandas.DataFrame(data=found)
    df2 = pandas.DataFrame(data=notfound)
    df = df1.append(df2, ignore_index=True, sort=False)

    # download fastas for found accessions and merge dataframes
    if fasta:
        seqs = download_fasta(fastas)
        if not seqs:
            logger.info('fasta download failed')

    # re-order columns 
    cols = list(df)
    cols.insert(0, cols.pop(cols.index('searched')))
    cols.insert(0, cols.pop(cols.index('label')))   
    df = df.loc[:, cols]

    return df


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


def query_plasmid_filter(fasta=False, plasmid='all', plasmid_strategy='contains',
                         topology='all', topology_strategy='is',
                         earliest_date='1982/01/01', latest_date=None,
                         location='all', location_strategy='contains',
                         location_mapped='all', location_mapped_strategy='contains',
                         min_latitude=-90, max_latitude=90,
                         min_longitude=-180, max_longitude=180,
                         isolation_source='all', isolation_source_strategy='contains',
                         host='all', host_strategy='contains',
                         sample_type='all', sample_type_strategy='contains',
                         plasmidfinder='all', plasmidfinder_strategy='contains',
                         pmlst='all', pmlst_strategy='contains',
                         min_length=0, max_length=10000000,
                         min_gc=0, max_gc=100,
                         taxon='all', taxon_strategy='contains',
                         species='all', species_strategy='contains',
                         genus='all', genus_strategy='contains',
                         family='all', family_strategy='contains',
                         order='all', order_strategy='contains',
                         bclass='all', bclass_strategy='contains',
                         phylum='all', phylum_strategy='contains',
                         kingdom='all', kingdom_strategy='contains'
                         ):
    """
    Filter the PLSDB plasmids
    :param fasta: if True generates fasta file of filtered plasmids
    :param plasmid: NCBI sequence accession
    :param plasmid_strategy: 'contains', 'begins', 'ends', 'is' for plasmid
    :param topology: 'circular', 'linear' or 'not-set'
    :param topology_strategy: 'is'
    :param earliest_date: format YYYY/MM/DD
    :param latest_date: format YYYY/MM/DD
    :param location: Country (:city)
    :param location_strategy: 'contains', 'begins', 'ends', 'is' for location
    :param location_mapped: Country (,city)
    :param location_mapped_strategy: 'contains', 'begins', 'ends', 'is' for location_mapped
    :param min_latitude: values between -90 and 90
    :param max_latitude: values between -90 and 90
    :param min_longitude: values between -180 and 180
    :param max_longitude: values between -180 and 180
    :param isolation_source:
    :param isolation_source_strategy: 'contains', 'begins', 'ends', 'is' for isolation_source
    :param host:
    :param host_strategy: 'contains', 'begins', 'ends', 'is' for host
    :param sample_type:
    :param sample_type_strategy: 'contains', 'begins', 'ends', 'is' for sample_type
    :param plasmidfinder:
    :param plasmidfinder_strategy: 'contains', 'begins', 'ends', 'is' for plasmidfinder
    :param pmlst:
    :param pmlst_strategy: 'contains', 'begins', 'ends', 'is' for plasmidfinder
    :param min_length:
    :param max_length:
    :param min_gc: values between 0 and 100
    :param max_gc: values between 0 and 100
    :param taxon:
    :param taxon_strategy: 'contains', 'begins', 'ends', 'is' for for taxon
    :param species:
    :param species_strategy: 'contains', 'begins', 'ends', 'is' for species
    :param genus:
    :param genus_strategy: 'contains', 'begins', 'ends', 'is' for genus
    :param family:
    :param family_strategy: 'contains', 'begins', 'ends', 'is' for family
    :param order:
    :param order_strategy: 'contains', 'begins', 'ends', 'is' for order
    :param bclass:
    :param bclass_strategy: 'contains', 'begins', 'ends', 'is' for class
    :param phylum:
    :param phylum_strategy: 'contains', 'begins', 'ends', 'is' for phylum
    :param kingdom:
    :param kingdom_strategy: 'contains', 'begins', 'ends', 'is' for kingdom
    :return:
    """

    local = locals()
    PARAMS = dict()
    for var in local:
        if var == 'fasta': continue
        PARAMS[var] = local[var]
    
    logger.info('START search for plasmids')
    response = requests.get(url=URL + 'filter/', params=PARAMS)

    if response.status_code == 200:
        data = response.json()
        df = pandas.read_json(data, orient='split')

        if len(data) == 0:
            logger.info('No plasmid found for filter %s.' % PARAMS)
            return

        if fasta:
            fastas = df['ACC_NUCCORE']
            seqs = download_fasta(fastas)
            if not seqs:
                logger.info("fasta download failed")

        logger.info("DONE")       
        return df
    else:
        raise Exception('Http response is not OK. Response status code is %s.' % response.status_code)



    









