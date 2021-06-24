#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 09:50:38 2021

@author: antonio
"""
from Bio import Entrez
import codecs
import unicodedata
import os
import datetime
import re
import argparse
import pickle

_RE_COMBINE_WHITESPACE = re.compile(r"\s+")


def search(query):
    # Get PubMed IDS
    Entrez.email = 'PLACEHOLDER@gmail.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax='2000000',
                            retmode='xml',
                            term=query) 
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    # Get full records
    ids = ','.join(id_list)
    Entrez.email = 'PLACEHOLDER@gmail.com'
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
        
def get_PMIDs(queries, fnames, base_path, logfile, all_time=False):
    sets = []
    for query,fname in zip(queries, fnames):
        complete_query = f'{query} AND (hasabstract)'
        if all_time==False:
             complete_query = complete_query + 'AND (2011/01/01:2021/12/31[dp])'
        r = search(complete_query)
        logfile.write(fname + '\t' + str(len(set(r['IdList']))) + '\n')
        fout = open(os.path.join(base_path, f'pmids/{fname}.txt'), 'w')
        for item in sorted(set(r['IdList'])):
            fout.write(item + '\n')
        fout.close()
        sets.append(r['IdList'])
        
    # Set unions
    setsv2 = set().union(*sets)
    fout = open(os.path.join(base_path, 'pmids/union.txt'), 'w')
    for item in sorted(setsv2):
        fout.write(item + '\n')
    fout.close()
    
    return setsv2


def basic_text_cleaning(s, _RE_COMBINE_WHITESPACE):
    #s = s.replace("±", "+/-") #there are times in which the web version has '±' and other '+/-'. I cannot standardize this
    s = unicodedata.normalize(u'NFC',s.strip()).encode("utf-8").decode("utf-8")
    # if "</" in s:
    #     s = re.sub('<[^<]+>', "", s)
    
    s_final = _RE_COMBINE_WHITESPACE.sub(" ", s).strip()
    return s_final

def parse_abstract(abstract, PMID, logfile):
    
    # Parse structured abstracts
    if len(abstract)>1:
        abstract_str = ''
        for item in abstract:
            if item.attributes == {}:
                if abstract_str == '':
                    abstract_str = str(item)
                else:
                    abstract_str = abstract_str+' '+ str(item)
                continue
                
            if 'Label' not in item.attributes.keys():
                logfile.write(f'WARNING: Skipping {PMID}: Abstract length > 1 but "Label" not in item attributes\n')
                break
            if abstract_str == '':
                if item.attributes['Label'] == 'UNLABELLED':
                    abstract_str = str(item)
                else:
                    abstract_str = item.attributes['Label']+': '+ str(item)
            else:
                if item.attributes['Label'] == 'UNLABELLED':
                    abstract_str = abstract_str+' '+ str(item)
                else:
                    abstract_str = abstract_str+' '+item.attributes['Label']+': '+ str(item)
        
    else:
        abstract_str = abstract[0]
                
    return abstract_str
    
def get_full_records(base_path, setsv2, logfile):
    
    c = 0
    # check = []
    logfile.write('Chunks:\n')
    for item in sorted(list(chunks(list(setsv2), 10000))):
        c=c+1
        fout = codecs.open(os.path.join(base_path,f'{c}_background_set.txt'), 'w', 'utf-8')
        logfile.write(f'{c}\n')
        res = fetch_details(item)
        
        for result in res['PubmedArticle']:
            PMID = str(result['MedlineCitation']['PMID'])
            #logfile.write(f'{PMID}\n') # Added
            title = result['MedlineCitation']['Article']['ArticleTitle'].strip()
            
            # Skip records with no abstract
            if 'Abstract' not in result['MedlineCitation']['Article'].keys():
                logfile.write(f'WARNING: Skipping {PMID}: Abstract key not in record\n')
                continue
            if 'AbstractText' not in result['MedlineCitation']['Article']['Abstract'].keys():
                logfile.write(f'WARNING: Skipping {PMID}: AbstractText key not in record\n')
                continue
            
            abstract = result['MedlineCitation']['Article']['Abstract']['AbstractText']
            
            abstract_str = parse_abstract(abstract, PMID, logfile)
                
            # Basic text cleaning
            abstract_str = basic_text_cleaning(abstract_str, _RE_COMBINE_WHITESPACE)
            title = basic_text_cleaning(title, _RE_COMBINE_WHITESPACE)
            
            """if any(x in abstract_str for x in ['\t', '\n']):
                logfile.write(f'WARNING: Skipping {PMID}: Abstract contains invalid characters\n')
                continue"""
            
            fout.write(f'{PMID}\t{title}\t{abstract_str}\n')
        fout.close()
    # return check


def argparser():
    '''
    DESCRIPTION: Parse command line arguments
    '''
    
    parser = argparse.ArgumentParser(description='process user given parameters')
    parser.add_argument("-i", "--input", required = True, dest = "input", 
                        help = "path to input set")
    parser.add_argument("-o", "--output", required =  True, dest="output", 
                        help = "path to output folder")
    parser.add_argument("-l", "--logfile", required =  True, dest="logfile", 
                help = "logfilepath")
    args = parser.parse_args()

    
    return args.input, args.output, args.logfile

def load_obj(directory, name):
    '''Helper function using pickle to save and load objects'''
    with open(os.path.join(directory,name), "rb") as f:
        print(os.path.join(directory,name))
        return pickle.load(f)


if __name__ == '__main__':

    mode = 'fetch'
    if mode == 'fetch':
        
        setpath, basepathoutput, logfilepath = argparser()
    
        setsv2 = set(map(lambda x: str(x.strip('\n')), open(setpath)))
        logfile = open(logfilepath, 'w')
        ## Get Full records
        logfile.write(f'{datetime.datetime.now()}: Get Full records\n')
        logfile.write(f'Store full records in {basepathoutput}\n')
        get_full_records(basepathoutput, setsv2, logfile)
        logfile.write(f'{datetime.datetime.now()}: Process finished')
        logfile.close()
        
    elif mode=='search_pmids':
        
        '''
        base_path, rules_name, logfilepath = argparser()
        # rules_name = 'criteria10-rules.txt''
        # base_path = '/home/antonio/Documents/Projects/Tasks/biocreative/background-set/my-downloaded-files/criteria10-all-time'
        logfile = open(logfilepath, 'w')
        
        queries = list(map(lambda x: x.strip('\n'), open(os.path.join(base_path, rules_name)).readlines()))
        fnames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                  'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
                  'AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH', 'II', 'JJ', 'KK',
                  'LL', 'MM', 'NN', 'OO', 'PP', 'QQ', 'RR', 'SS', 'TT', 'UU', 'VV',
                  'WW', 'XX', 'YY', 'ZZ', 'AAA', 'BBB', 'CCC', 'DDD', 'EEE', 'FFF',
                  'GGG', 'HHH', 'III', 'JJJ', 'KKK', 'LLL', 'MMM', 'NNN', 'OOO', 
                  'PPP', 'QQQ', 'RRR', 'SSS', 'TTT', 'UUU', 'VVV', 'WWW', 'XXX', 
                  'YYY', 'ZZZ','AAAA', 'BBBB', 'CCCC', 'DDDD', 'EEEE', 'FFFF', 
                  'GGGG', 'HHHH', 'IIII', 'JJJJ', 'KKKK', 'LLLL', 'MMMM', 'NNNN', 
                  'OOOO','PPPP', 'QQQQ', 'RRRR', 'SSSS', 'TTTT', 'UUUU', 'VVVV', 
                  'WWWW', 'XXXX', 'YYYY', 'ZZZZ','AAAAA', 'BBBBB', 'CCCCC', 'DDDDD', 
                  'EEEEE', 'FFFFF', 'GGGGG', 'HHHHH', 'IIIII', 'JJJJJ', 'KKKKK', 
                  'LLLLL', 'MMMMM', 'NNNNN', 'OOOOO','PPPPP', 'QQQQQ', 'RRRRR', 
                  'SSSSS', 'TTTTT', 'UUUUU', 'VVVVV', 'WWWWW', 'XXXXX', 'YYYYY', 
                  'ZZZZZ','AAAAAA', 'BBBBBB', 'CCCCCC', 'DDDDDD', 'EEEEEE', 'FFFFFF', 
                  'GGGGGG', 'HHHHHH', 'IIIIII', 'JJJJJJ', 'KKKKKK', 'LLLLLL', 'MMMMMM', 
                  'NNNNNN', 'OOOOOO','PPPPPP', 'QQQQQQ', 'RRRRRR', 'SSSSSS', 'TTTTTT', 
                  'UUUUUU', 'VVVVVV', 'WWWWWW', 'XXXXXX', 'YYYYYY', 'ZZZZZZ','AAAAAAA', 
                  'BBBBBBB', 'CCCCCCC', 'DDDDDDD', 'EEEEEEE', 'FFFFFFF', 'GGGGGGG', 
                  'HHHHHHH', 'IIIIIII', 'JJJJJJJ', 'KKKKKKK', 'LLLLLLL', 'MMMMMMM', 
                  'NNNNNNN', 'OOOOOOO','PPPPPPP', 'QQQQQQQ', 'RRRRRRR', 'SSSSSSS', 
                  'TTTTTTT', 'UUUUUUU', 'VVVVVVV', 'WWWWWWW', 'XXXXXXX', 'YYYYYYY', 
                  'ZZZZZZZ','AAAAAAAA']
        
        ## Get PubMed IDs
        logfile.write(f'{datetime.datetime.now()}: Get PMIDs\n')
        logfile.write(f'Store PMIDs in {os.path.join(base_path, "pmids")}\n')
        setsv2 = get_PMIDs(queries, fnames, base_path, logfile, all_time=True)
        logfile.write(f'{datetime.datetime.now()}: Total number of unique PMIDs: {str(len(setsv2))}\n')
        logfile.close()'''

