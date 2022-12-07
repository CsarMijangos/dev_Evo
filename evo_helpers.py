import evo_constants as CTS


                ######################## 
                ####   function to ##
                ####   say welcome ##
                ########################

def welcome():
    """ This function presents evomining to the user. 
    """
    print("Welcome to EvoMinig\n")
    input("Press enter to continune ...")
    return

    
    
######################## 
####     function to 
###  show the options 
########################
#OPTION = 0
def option_menu():
    """This function displays the option for execute evomining.
    """
    option = 0
    while option not in [1,2]:
        print("Press 1 to continue with precharged example.")
        print("Press 2 to continue with your own databases.")
        try: 
            option = get_option()
        except:
            print("Invalid option")
            continue
    return option
     

####################### 
####     function to 
###  manage the user option 
########################

def get_option():
    """ This function requires the user to choose
    an option for the user to run evomining. 
    Precharged example or the user data.
    """
    option = int(input("\n Enter option: "))
    if option not in [1,2]:
        raise ValueError("Not a valid option, try again")
    else:
        return option
    
    
######################## 
####     function to 
###  obtain the user databases 
########################

def get_dbs():
    """This function requires the databases paths from 
    the user.
    """
    central_db = input("\n Enter path to the directory that contains the central_DB: ")
    genomes_db = input("\n Enter path to the directory that contains the genomes_DB (.faa and .txt files): ")
    rast_ids = input("\n Enter path to the file rast.ids: ")
    nat_prods_db = input("\n Enter path to the directory that contains the nat_prods_DB: ")
    DBs = { "central" : central_db,
           "genomes" : genomes_db,
           "rast_ids" : rast_ids,
           "nat_prods" : nat_prods_db }
    return DBs

#####################################################################
##############################   Functions to obtain ################
###############################   the headers from     ##############
###############################  the genomes' files  ###############
####################################################################
#####################################################################

import os


def dict_ids_to_names_genomes(file_rast_ids):
    """ This function returns a dictionary with the genome id's as keys
    and the genome's names as values. The input file_rast_ids is the Rast.ids
    file within the genomes' database.
    """
    genomes_ids_names_dict = {}
    with open(file_rast_ids, "r") as file:
        for line in file:
            genomes_ids_names_dict[line.split("\t")[1]] = line.split("\t")[2].rstrip()
            
    return genomes_ids_names_dict


def make_EvoHeader(txt_file_line, dict_genome):
    """ This function obtains the headers in the
    evomining format. The arguments are the header
    in the .txt file in the genomes db and a
    dictionary of ids to names of the genomes."""
    x = txt_file_line.split('\t')[0:8] 
    y = x[1].split('|')[1].split('.')
    #eliminamos la cadena del type: 'peg':
    featureId = [y[i] for i in range(0,len(y)) if i != 2] 
    genomeId = ".".join(featureId[0:2])
    genomeIdPeg = '.'.join(featureId)
    function = '_'.join(x[7].split(" "))
    header = [">gi",genomeIdPeg,genomeId,dict_genome[genomeId].split(" ")[-1], 
             function,"_".join(dict_genome[genomeId].split(" "))]
    evo_header = "|".join(header)+"\n"
    return evo_header


def dict_faaHead_seqLines(faa_file):
    """ This function returns a dictionary
    with the headers in the .faa file of a genome
    as keys and a list with the potition of the
    first and last line of its sequence in the file.
    """
    dict_faaHead_seqlines = dict()
    range_seq_lines = []
    with open(faa_file, 'r') as file:
        current_line = 0
        prev_head = file.readline()
        dict_faaHead_seqlines[prev_head] = [current_line + 1]
        for line in file:
            current_line+=1 
            if ">" in line:
                dict_faaHead_seqlines[prev_head].append(current_line - 1)
                dict_faaHead_seqlines[line] = [current_line+1]
                prev_head = line
               
        dict_faaHead_seqlines[prev_head].append(current_line)
    return dict_faaHead_seqlines


def header_dict_EvoFaa(txt_file, dict_genome):
    """ This function returns a dictionary  with keys
    equal to the headers in the .faa file and values 
    equal to the evomining headers. Its arguments are
    the .txt file and a dictionary with the genomes ids
    as keys and the genomes names as values.
    """
   
    header_dict_evoheader = {}
    with open(txt_file,'r') as file:
        next(file)
        for line in file:
            line = line.strip()
            #Nos aseguramos de solo considerar los "peg":
            if "peg" in line:
                evo_header = make_EvoHeader(line,dict_genome)
                aux = evo_header.split("|")
                aux2 = aux[1].split(".")
                header_faa = ">fig"+ "|" + aux[2] + ".peg." + aux2[2]+"\n"
                header_dict_evoheader[header_faa] = evo_header 
            else: pass
    return header_dict_evoheader



def make_evo_headers_file(faa_file, txt_file, rast_ids_file):
    """ This function obtains the .faa
    file with the headers in the evomining format.
    The arguments are the paths to the .faa and .txt files in 
    the genomes data base, and the path to the rastids.
    """
    # Dictionary of faa headers to lines of its sequence:
    dicc_headers_seqs = dict_faaHead_seqLines(faa_file)
    
   
    # Dictionary of genomes ids to genomes names:
    dict_genomes_names = dict_ids_to_names_genomes(rast_ids_file)
    
   
    # Dictionary of faa headers to evomining fmt headers:
    header_faa_evo = header_dict_EvoFaa(txt_file, dict_genomes_names)
    
    # list to filter the headers with "peg":
    peg_keys = [key for key in list(dicc_headers_seqs.keys()) if "peg" in key]
    
    # Path to store the file with the .faa files with evomining headers:
    evo_header_file = CTS.GENOMES_EvoFmt + "evo_" + faa_file.split("/")[-1]
    
    read_file = open(faa_file, 'r')
    write_file = open(evo_header_file, "w")
    for pos, line in enumerate(read_file):
        if line in peg_keys:
            init = dicc_headers_seqs[line][0]
            end = dicc_headers_seqs[line][1] + 1
            line_list = list(range(init,end))
            evo_header = header_faa_evo[line]
            write_file.write(evo_header)
        if pos in line_list:
            write_file.write(line)
    read_file.close()
    write_file.close()
    return 

def make_all_evo_headers(faa_files_path, rast_ids_file):
    """ This function creates all the .faa files
    with the headers in the evomining's format.
    """
    faa_files_list = [x for x in os.listdir(faa_files_path) if ".faa" in x]
    
    if not faa_files_path.endswith("/"):
        faa_files_path = faa_files_path + "/" 
        
    for file in faa_files_list:
        faa_file = faa_files_path + file
        txt_file = faa_files_path + file.split(".")[0] + ".txt"
        
        make_evo_headers_file(faa_file, txt_file, rast_ids_file)
    return
    
## Make the evo_genomes_db: the data base of genomes with 
## the evomining's headers.
    
def join_evo_headers_files():   
    """This function joins the .faa files that have the evomining's headers
    in a single file named evo_genomes_db.faa.
    """
    files_list = [x for x in os.listdir(CTS.GENOMES_EvoFmt) if ".faa" in x]
    ## Path to the location of the evo_genomes_db:
    evo_genomes_db = CTS.EVO_GENOMES_DB + "evo_genomes_db.fasta"
    
    with open(evo_genomes_db, "w") as new_file:
        for file in files_list:
            name = CTS.GENOMES_EvoFmt + file
            f = open(name, "r")
            #with open(name) as f:
            for line in f:
                new_file.write(line)
            f.close()
                
    return 
    

################################################################################
########################## Delete the files once used ##########################
################################################################################
    
from os import remove
def remove_evo_headers_files():
    """This function removes the .faa files
    used to obtain the evo_genomes_db.faa.
    """
    files_list = [x for x in os.listdir(CTS.GENOMES_EvoFmt) if ".faa" in x]
    for file in files_list:
        name = CTS.GENOMES_EvoFmt + file
        remove(name)
    return
    
    

    
    
####################################################################
####################################################################
#################### Helpers to apply blast ########################
####################################################################
####################################################################

import matplotlib.pyplot as plt


######################### Makeblastdbs ##############################

from Bio.Blast.Applications import NcbimakeblastdbCommandline

def makeblast_db(path_input_db):
    """ This function obtains the blastdb of the db passed as 
    its argument.
    """
    
    aux = [x for x in os.listdir(path_input_db) if ".fasta" in x]
    
    input_db = path_input_db + aux[0]
    
    output_name = (CTS.BLASTDBs_PATH + aux[0].replace(".fasta", "")+
                   "_blastdb/"+aux[0].replace(".fasta", "") + "_blastdb")
    #print(output_name)
    cline_db = NcbimakeblastdbCommandline(cmd = CTS.MAKEBLASTDB_CMD,
                                    dbtype="prot",
                                   input_file= input_db,
                                     out = output_name)
    cline_db()
    return
    


###########################   Blastp    #################################

from Bio.Blast.Applications import NcbiblastpCommandline

def apply_blastp(query_path, blastdb_path):
    """This function applies blastp with query
    equal to the file with path given by query_path,
    and blast database in the directory given by
    blastdb_path.
    """
    ## obtaining query name:
    if query_path.endswith(".fasta"):
        query_file = query_path
    else:
        aux = [x for x in os.listdir(query_path) if ".fasta" in x]
        query_file = query_path + aux[0]
    
    ## obtaining the blastdb's path:
    blastDB = blastdb_path + blastdb_path.split("/")[-2]
    
    
    ## obtaining output's path:
    if "central" in aux[0]:
        output_blastp = (CTS.BLASTp_PATH +
                         "central_to_genomes/central_to_genomes.blast")
    elif "evo_genomes_db" in aux[0]:
        output_blastp = (CTS.BLASTp_PATH +
                         "genomes_to_central/genomes_to_central.blast")
    elif "exp_fam" in aux[0]:
        output_blastp = (CTS.BLASTp_PATH +
                         "exp_fam_to_nat_prods/exp_fam_to_nat_prods.blast")
    #print(blastDB)
    cline_blastp = NcbiblastpCommandline(cmd = CTS.BLASTp_CMD,
                                     query= query_file,
                                     db= blastDB,
                                     evalue=0.0001, 
                                     remote= False, 
                                     ungapped=False,
                                     outfmt = 6,
                                     max_target_seqs = 10000,
                                     out = output_blastp)
    cline_blastp()
    return 



################################################################################
###########################  Helpers to obtain de BBH  #########################
################################################################################


















    