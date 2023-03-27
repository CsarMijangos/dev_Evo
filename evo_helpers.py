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

####-----------------------------------------------------------------------------------------------------------------------#####

#####################################################################
##############################   Functions to obtain ################
###############################   the headers from     ##############
###############################  the genomes' files  ###############
####################################################################
#####################################################################

import os

#Esta función también se utiliza para obtener los bbh's
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
    evomining format. The arguments are the headers
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
    as keys and a list with the position of the
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
    in a single file named evo_genomes_db.fasta.
    """
    files_list = [x for x in os.listdir(CTS.GENOMES_EvoFmt) if ".faa" in x]
    
    
    new_files_list = []
    for file in files_list:
        old_name = CTS.GENOMES_EvoFmt + file 
        F = open(old_name, "r")
        new_file = F.readline().split("|")[2] + ".faa"
        new_files_list.append(new_file)
        new_name = CTS.GENOMES_EvoFmt + new_file
        F.close()
        os.rename(old_name, new_name)
    
    ordered_files_list = sorted(new_files_list, key = lambda t : float(t.split(".faa")[0]))
    ## Path to the location of the evo_genomes_db:
    evo_genomes_db = CTS.EVO_GENOMES_DB + "evo_genomes_db.fasta"
    
    with open(evo_genomes_db, "w") as new_file:
        for file in ordered_files_list:
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
    
    

####-----------------------------------------------------------------------------------------------------------------------#####    
    
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
    
    cline_db = NcbimakeblastdbCommandline(cmd = CTS.MAKEBLASTDB_CMD,
                                    dbtype="prot",
                                   input_file= input_db,
                                     out = output_name)
    cline_db()
    return
    

#########################################################################
###########################   Blastp    #################################
#########################################################################


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
    
    #if "central" in aux[0]:
    if "central.fasta" in query_file:
        output_blastp = (CTS.BLASTp_PATH +
                         "central_to_genomes/central_to_genomes.blast")
    #elif "evo_genomes_db" in aux[0]:
    elif "evo_genomes_db.fasta" in query_file:
        output_blastp = (CTS.BLASTp_PATH +
                         "genomes_to_central/genomes_to_central.blast")
    elif "expanded_fam" in query_file:
        input_file = query_file.split("/")[-1].split(".")[0]
        output_name = input_file.split(".")[0] + "_to_mibig.blast"
        output_blastp = (CTS.BLASTp_PATH +
                         "exp_fam_to_nat_prods/" + output_name)
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



def blast_vs_mibig():
    """ This function applies blastp with queries the copies in each organism, and 
    subject the mibig database.
    """
    expanded_families_files = [x for x in os.listdir(CTS.EXPANDED_FAMS) if x.endswith(".fasta")]
    if expanded_families == []:
        return print("There are no expanded families")
    else:
        mibig_db = CTS.BLASTDBs_PATH + "nat_prods_blastdb/"
        for file in expanded_families_files:
            query = CTS.EXPANDED_FAMS + file
            apply_blastp(query, mibig_db)
        return 
####-----------------------------------------------------------------------------------------------------------------------#####

################################################################################
####################### Make json files ########################################
################################################################################















################################################################################
###################  Helpers to obtain the expanded families  ##################
################################################################################

import numpy as np
import pandas as pd
    


    
    
############ Create auxiliary files with the enzimes and organisms information ###########
def make_central_names_file(central_fasta_file_path):
    """This function creates a file with the names of all
    the elements in central_db.
    """
    output_path = CTS.BBH_aux_files + "central_enzimes.txt"
    with open(central_fasta_file_path,"r") as fasta_file:
        with open(output_path, "w") as output:
            for line in fasta_file:
                if ">" in line:
                    Line = line.split(">")[1]
                    output.write(Line)
    return
    
    
##### Create the ids to orgs names auxiliary file

def ids_to_names_genomes_file(file_rast_ids):
    """ This function returns a file with the genome id's as keys
    and the genome's names as values. The input file_rast_ids is the Rast.ids
    file within the genomes' database.
    """
    
    genomes_ids_names_dict = {}
    output_path = CTS.BBH_aux_files + "ids_to_orgs_names.txt"
    with open(file_rast_ids, "r") as file:
        with open(output_path, "w") as output:
            for line in file:
                Line = line.split("\t")[1] + "-->" + line.split("\t")[2].rstrip()
                output.write(Line+"\n")
            
    return 
    


################################################################################
######################  Helpers to obtain the best hit  ########################
######################    from central to organisms      #########################
################################################################################

def best_hit_central_to_orgs(blast_file):
    """ This function will obtain the best hit from central to organisms. Its output 
    is a file with these hits. From the output we will obtain the expanded families.
    """
    
    # Charge the dataframe with the blastp file:
    df = pd.read_csv(blast_file, sep = '\t', names = CTS.BLAST_COLS)
    
    ### SE NECESITA QUE LOS ARCHIVOS AUXILIARES ESTÉN CREADOS ###
    #Files with the enzimes and orgs info:
    
    # Make a list with the enzimes:
    central_enzimes_file = CTS.BBH_aux_files + 'central_enzimes.txt'       #Agregar ruta a constantes
    central_enzimes = []
    with open(central_enzimes_file, "r") as central_file:
        for line in central_file:
            if line:
                central_enzimes.append(line.rstrip("\n"))
    
    
    # Make a dictionary with the ids and names of organisms:
    ids_to_orgs_names_file = CTS.BBH_aux_files + 'ids_to_orgs_names.txt'   #Agregar ruta a constantes
    genomes = dict()
    with open(ids_to_orgs_names_file, "r") as orgs_file:
        for line in orgs_file:
            if line:
                key = line.split("-->")[0]
                value = line.split("-->")[1]
                genomes[key] = value
    
    
    
    #Define and initialize the output file:
    
    # Path to best hits file:
    
    central_to_orgs_bh = CTS.BBH_aux_files + "central_to_orgs_BH.uniq"  #Agregar ruta a constantes
    
    
    # Finding the best hit:
    best_hits_file = open(central_to_orgs_bh, "w")
    for enzime in central_enzimes:
        for genome_id in genomes.keys():
            aux_list = []
            
            for index in range(len(df)):
                if (aux_list == []):
                    if ((enzime == df["query"][index])&(genome_id in df["subject"][index])):
                        aux_list = list(df.iloc[index])
                    else:
                        continue
                else:
                    if ((enzime == df["query"][index]) & (genome_id in df["subject"][index]) & (df["bitscore"][index] > aux_list[11])):
                        aux_list = list(df.iloc[index])
            if aux_list != []: 
                line = ["{} ".format(item) for item in aux_list]
                best_hits_file.writelines(line)
                best_hits_file.write("\n")
                
    best_hits_file.close()
    return
    
################################################################################
######################  Helpers to obtain the best hit  ########################
######################    from organisms to central      #########################
################################################################################


def best_hit_orgs_to_central(blast_file):
    """ This function will obtain the best hit from organisms to central, i.e., 
    assigns to each position id of a blast hit the corresponding enzime that best 
    hits in that particular place.
    """
    
    #charge dataframe:
    df = pd.read_csv(blast_file, sep = '\t', names = CTS.BLAST_COLS, index_col=False)
    
    #Obtain dictionary with the copy id position to the line of blast
    
    uniq2 = dict()
    for index in df.index:
        actual_line = df.iloc[index]
        copy_id = actual_line["query"].split("|")[1]  #id of the hit position
        if copy_id in uniq2:
            if actual_line["bitscore"] > uniq2[copy_id][11]:
                uniq2[copy_id] = list(actual_line)
        else:
            uniq2[copy_id] = list(actual_line)
    
    
    # Write the file with the best hits:
    
    orgs_to_central_bh = CTS.BBH_aux_files + "orgs_to_central_BH.uniq"  #Agregar ruta a constantes
    best_hits_file = open(orgs_to_central_bh, "w")
     
    for key in uniq2:
        line = ["{} ".format(item) for item in uniq2[key]] 
        best_hits_file.writelines(line)
        best_hits_file.write("\n")
    
    
    best_hits_file.close()
    return




################################################################################
###########################  Helpers to obtain de BBH  #########################
################################################################################


def bbh_detector(central_to_orgs_bhs, orgs_to_central_bhs):
    """ This function obtains the best bidirectional hits between
    the blast of central to organisms and organisms to central respectively.
    """

    # charge dataframes:
    
    df = pd.read_csv(central_to_orgs_bhs, sep = ' ', names = CTS.BLAST_COLS, index_col = False)
    
    df1 = pd.read_csv(orgs_to_central_bhs, sep = ' ', names = CTS.BLAST_COLS, index_col = False)
    
    # Write the file with the bbh:
    
    BBH_file = CTS.BBH_aux_files + "BBH.bbh"  #Agregar ruta a constantes
    bbh_file = open(BBH_file, "w")
     
    for indx in df1.index:
        copy_id = df1.iloc[indx]["query"].split("|")[1]
        enzime = df1.iloc[indx][1]
        for indxx in list(df.index):
            if copy_id not in df.iloc[indxx]["subject"]:
                continue
            elif enzime not in df.iloc[indxx]["query"]:
                continue
            else:
                #bbh[copy_id] = enzime
                bbh_file.write(copy_id + " -> ")
                bbh_file.write(enzime +"\t")
                bbh_file.write(df1.iloc[indx]["query"]+"\n")
           
    
    bbh_file.close()
    return



#######################################################################################################################
############################################ Helpers to obtain the ######################################################
###########################################   expansions and the ######################################################
###########################################    recruited enzimes ######################################################
#######################################################################################################################


#No se utilizó:
def bitscore_filter(dataframe, threshold):
    """ This function returns a dataframe with only those
    entries with a bitscore >= threshold.
    """
    dataframe.query("bitscore >= @threshold", inplace = True)
    dataframe.reset_index(drop=True, inplace=True)
    return dataframe





##*********************************************************************************************************************
########################################### Obtain how many copies are ###########################################
########################################### in each organism            ###########################################
#*********************************************************************************************************************


def copy_count(bistcore_threshold = 100, evalue_threshold = 0.001):
    """This function returns a dictionary with
    the family's number of the enzimes as keys. Each key
    has as a value another dictionary that has as keys the
    organisms names and as values the list of the copies ids into 
    this organism of each enzime in the corresponding family.
    """
    
    
    # charge the dataframe:
    # blast_file is the path to the blast of central to genomes
    blast_file = CTS.BLASTp_PATH + '/central_to_genomes/central_to_genomes.blast'
    df = pd.read_csv(blast_file, sep = '\t', names = CTS.BLAST_COLS, index_col=False)
    
    
    ## Obtain the dictionary of fam to enzimes:
    queries = df["query"].unique()
    fam_keys = set([item.split("|")[1] for item in queries])
    enz_fams = dict()
    for key in fam_keys:
        enz_fams[key] = []
    for qry in queries:
         for key in fam_keys:
            if key == qry.split("|")[1]:
                enz_fams[key].append(qry)
    
    ## Obtain the list of organisms ids:
    orgs_ids_list = list()
    with open(CTS.BBH_aux_files+"ids_to_orgs_names.txt", "r") as file:
        for line in file:
            org_id = line.split("-->")[0]
            orgs_ids_list.append(org_id)
    
    #### Filtering by bitscore and evalue:
    
    DF1 = df[(df["bitscore"]>=bitscore_threshold) & (df["e_value"]<= evalue_threshold)]
    DF1.reset_index(drop=True, inplace=True)
    
    #### Obtain the output dictionary:
    
    # For the i-th family we obtain a dictionary with the organisms as keys and the position id
    # of the copies of each member of that family.
    copies_by_family = dict()
    for fam in fam_keys:
        
        copies_into = dict()
        for org in orgs_ids_list:
            copies_into[org] = []
        for indx in DF1.index:
            actual_sbjt = DF1.iloc[indx]["subject"]
            actual_org = actual_sbjt.split("|")[2]
            actual_copy = actual_sbjt.split("|")[1]
            actual_enzm = DF1.iloc[indx]["query"]
            actual_fam = actual_enzm.split("|")[1]
            for org in orgs_ids_list: #the last condition is for avoid duplicates
                if ((actual_fam == fam) & (actual_org == org) & (actual_copy not in copies_into[org])):
                    copies_into[org].append(actual_copy)
        for org in orgs_ids_list:
            if copies_into[org] == []:
                del copies_into[org]
        copies_by_family[fam] = copies_into
    
    return copies_by_family


# obtaining the thresholds for the number of copies present for considering that
# an expansion event has ocurred. 

###The input is the output of the copy_count function.

def expansion_thresholds(copies_by_fam):
    """ This function returns a dictionary with the families numbers ids
    as keys, and the value of each key is a dictionary with the mean, median,
    mean + std, and mean + 2*std of the number of copies of the members
    of the family into each organism. 
    """
    
    thresholds_by_family = dict()
    thresholds_types = ["mean","median", "mean + std", "mean + 2*sdt"]
    for fam in copies_by_fam:
        counting_array = np.array([len(copies_by_fam[fam][org]) for org in copies_by_fam[fam] if copies_by_fam[fam][org]!=[]])
        mean =  counting_array.mean()
        median = np.median(counting_array)
        std = counting_array.std()
        thresholds_by_family[fam] = {"mean": mean,
                                     "median": median,
                                    "mean+std": mean+std,
                                    "mean+2*std": mean + 2*std}
            
                
    return thresholds_by_family


##################################################################################################################################
########################################### Selecting the expansions ####################################################
##################################################################################################################################
        
def detect_expansions(bitscore_threshold = 100, evalue_threshold=0.001):
    """ This function returns a dictionary with keys equal to the names
    of each organism that has expansions, and a list with the positions
    of the expansion copies as values.
    """
    
    copies_by_family = copy_count(bitscore_threshold, evalue_threshold)
    
    expansion_thrlds = expansion_thresholds(copies_by_family)
    
    #dictionary of org_ids to org_names:
    orgs_ids_names = dict()
    with open(CTS.BBH_aux_files+"ids_to_orgs_names.txt", "r") as file:
        for line in file:
            org_id = line.split("-->")[0]
            org_name = line.split("-->")[1]
            orgs_ids_names[org_id] = org_name
    
    #dictionary with the lengths of the number of copies in each organism by family:
    number_of_copies = dict()
    for fam in copies_by_family:
        leng = dict()
        for org in copies_by_family[fam]:
            leng[org] = len(copies_by_family[fam][org])
        number_of_copies[fam] = leng
    # obtaining the expansions by criteria:
    expansions =dict()
    for fam in copies_by_family:
        copies = {"mean": [],
                  "median": [],
                 "mean+std": [],
                 "mean+2*std": []}
        for org in copies_by_family[fam]:
            if number_of_copies[fam][org] > expansion_thrlds[fam]["mean"]:
                copies["mean"].append(copies_by_family[fam][org])
            elif number_of_copies[fam][org] > expansion_thrlds[fam]["median"]:
                copies["median"].append(copies_by_family[fam][org])
            elif number_of_copies[fam][org] > expansion_thrlds[fam]["mean+std"]:
                copies["mean+std"].append(copies_by_family[fam][org])
            elif number_of_copies[fam][org] > expansion_thrlds[fam]["mean+2*std"]:
                copies["mean+2*std"].append(copies_by_family[fam][org])
            else:
                continue
        expansions[fam] = copies
    return expansions

##################################################################################################
    ############################     Version 2  ##########################################
####################################################################################

import json

def best_hits_blast1(path_blast_1, rast_ids_file):
    """ This function creates json files with dictionaries that contain the info of the 
    families (families.json), the ids to names of organisms (orgs_ids_names.json), 
    a dictionary with the copies of the enzimes of each family into each organism (Dictio_1.json),
    and a dictionary with the best hit of the copies.
    """
    
    # charge dataframe:
    df = pd.read_csv(path_blast_1, sep = "\t", names = CTS.BLAST_COLS, index_col = False )
    
    #Create dictionary of families to enzimes
    uniq_qries = df["query"].unique()
    families = dict()
    for enz in uniq_qries:
        fam = enz.split("|")[1]
        key = "fam_"+fam
        families[key] = []
    
    for enz in uniq_qries:
        fam = enz.split("|")[1]
        key = "fam_"+fam
        families[key].append(enz)
    
    
    fam_file = open(CTS.JSON_FILES_PATH + "families.json", "w")
    json.dump(families,fam_file)
    fam_file.close()

    # Create dictionary of orgs_ids_names:
    orgs_ids_names = dict()
    rast_file = open(rast_ids_file, "r")
    for line in rast_file:
        org_id = line.split("\t")[1]
        org_name = line.split("\t")[2]
        orgs_ids_names[org_id] = org_name
    rast_file.close()
    
    orgs_file = open(CTS.JSON_FILES_PATH + "orgs_ids_names.json", "w")
    json.dump(orgs_ids_names, orgs_file)
    orgs_file.close()
    
    ##### Create Dictionary with the copies by org:
    
    # Keys for the dictio
    keys_dict = []
    for fam in families:
        for org in orgs_ids_names:
            key = fam + "|" + org
            keys_dict.append(key)
            
    # Make dictionary   
    Dictio = dict()
    for key in keys_dict:
        Dictio[key] = dict()
    
    
    for key in keys_dict:
        fam = key.split("|")[0]
        for enz in families[fam]:
            Dictio[key][enz] = []
    
    
    for indx in df.index:
        if ((df.iloc[indx]["bitscore"]>=100) & (df.iloc[indx]["e_value"]<0.001)):
           
            fam = "fam_" + df.iloc[indx]["query"].split("|")[1]
            org = df.iloc[indx]["subject"].split("|")[2]
            key = fam + "|" + org
            enz = df.iloc[indx]["query"]
            enz_num = enz.split("_")[-1].split("|")[0]  # It is important that the enzimes are numered 
            copy = df.iloc[indx]['subject'].split("|")[1]  
            bitsc = df.iloc[indx]["bitscore"]
            element = [copy, bitsc, enz_num]
            Dictio[key][enz].append(element)            
        else:
            pass
    
    # Make file:
    dictio_file = open(CTS.JSON_FILES_PATH + "Dictio_1.json", "w")
    json.dump(Dictio, dictio_file)
    dictio_file.close()
    
    
    ### Create dictionary with the best hits:
    
    Dictio_bh = dict()
    for key in keys_dict:
        Dictio_bh[key] = []
    
        
    for org in Dictio:
        total_copies = []
        for enzs in Dictio[org]:
            total_copies += Dictio[org][enzs]
        Dictio_bh[org] = total_copies
    
    for key in Dictio_bh:
        best_hits = dict()
        for hit in Dictio_bh[key]:    # hit = [copy_id, bitscore, enz_num]
            if best_hits.keys() == []:
                best_hits[hit[0]] = hit  # best_hits ={ copy_id: hit}
            else:
                if hit[0] not in best_hits:
                    best_hits[hit[0]] = hit
                else:
                    if hit[1] > best_hits[hit[0]][1]:
                        best_hits[hit[0]] = hit
        Dictio_bh[key] = sorted(list(best_hits.values()), key = lambda v : int(v[0].split(".")[-1]))
    
    dict_bh_file = open(CTS.JSON_FILES_PATH + "Dictio_expanded_fams.json", "w")
    json.dump(Dictio_bh, dict_bh_file)
    dict_bh_file.close()
    
    return
    
    


    
    
########################################################################################################
################################## Obtaining the expanded families #####################################
########################################################################################################

def expanded_families(bitscore=100, evalue=0.001, criteria="mean"):
    """ This function returns the expanded families in the form of a
    dictionary with the families as keys, and by values has a list with 
    organisms ids that had expansions.
    """
    
    expansions = detect_expansions(bitscore, evalue)
    expanded_families = dict()
    for fam in expansions:
        if expansions[fam][criteria]:
            expanded_family_orgs_ids = []
            for copy_list in expansions[fam][criteria]:
                numbers_id = copy_list[0].split(".")
                org_id = numbers_id[0]+numbers_id[1]
                expanded_family_orgs_ids.append(org_id)
            expanden_families[fam] = expanded_family_orgs_ids 
        else:
             continue 
    return expanden_families
  
  

def obtain_seq_of_expFam(evo_db_path):
    """This function returns files with sequences of the copies in each
    organism. One file by each family of central enzimes that undergoes expansions.
    """
    import json
    
    ## Issue: REvisar si los archivos ya existen!!
    input_file = open(evo_db_path, "r")
    output_path = CTS.EXPANDED_FAMS
    input_file.seek(0,2)
    file_size = input_file.tell()
    input_file.seek(0)
    
    dictio_file = open(CTS.EXP_FAMS_JSON, "r")
    dictio_copies = json.load(dictio_file)
    dictio_file.close()
    
    orgs_file = open(CTS.ORGS_IDS_NMS, "r")
    orgs_dict = json.load(orgs_file)
    orgs_ids = list(orgs_dict.keys())
    orgs_file.close()
    
    fams_file = open(CTS.FAMILIES, "r")
    fams_dict = json.load(fams_file)
    families = list(fams_dict.keys())
    fams_file.close()
    
    sorted_orgs = sorted(orgs_ids, key = float)
    sorted_families = sorted(families, key= lambda f : int(f.split("_")[1]))
    
    
    for fam in sorted_families:
        output_file = open(output_path+"expanded_"+fam+".fasta", "a")
        for org in sorted_orgs:
            key = fam + "|" + org
            
            #return print(key)
            for copies_info in dictio_copies[key]:
                copy_id = copies_info[0]
                
                ptrn = ">gi|" + copy_id
                line = input_file.readline()
                #for line in input_file:
                while line:
                    if ptrn in line:
                        ptrn = ""
                        output_file.write(line)
                        next_line = input_file.readline()
                        while(">gi" not in next_line):
                            output_file.write(next_line)
                            next_line = input_file.readline()
                        break
                    if (input_file.tell() == file_size) & (ptrn != ''):
                        input_file.seek(0)
                        line = input_file.readline()
                        continue
                    line = input_file.readline()
                    
        output_file.close()        
    input_file.close()
    
    return 


########################################################################################################
#####################################  Obtaining the heat map ##########################################
########################################################################################################

#####################################  ************************     #####################################
#####################################  * 
def mean_highlighter(x):
    #style for less than mean
    style_lt = "background-color: white; color: black; font-weight: bold; border: solid; text-align: center; font-size:1.2em;"
    #style for greater than mean
    style_gt = "background-color: #C04000; color: black; font-weight: bold; border: solid; text-align:center;"
    gt_mean = x > x.mean()
    return [style_gt if i else style_lt for i in gt_mean]

def median_highlighter(x):
    style_lt = "background-color: white; color: black; font-weight: bold; border: solid; text-align: center; font-size:1.2em;"
    style_gt = "background-color: #C04000; color: black; font-weight: bold; border: solid; text-align:center;"
    gt_median = x > x.mean()
    return [style_gt if i else style_lt for i in gt_median]

def std_highlighter(x):
    style_lt = "background-color: white; color: black; font-weight: bold; border: solid; text-align: center; font-size:1.2em;"
    style_gt = "background-color: #C04000; color: black; font-weight: bold; border: solid; text-align:center;"
    gt_std = x > x.mean() + x.std()
    return [style_gt if i else style_lt for i in gt_std]

def std2_highlighter(x):
    style_lt = "background-color: white; color: black; font-weight: bold; border: solid; text-align: center; font-size:1.2em;"
    style_gt = "background-color: #C04000; color: black; font-weight: bold; border: solid; text-align:center;"
    gt_std = x > x.mean() + 2*x.std()
    return [style_gt if i else style_lt for i in gt_std]

def mode_highlighter(x):
    import statistics as st
    style_lt = "background-color: white; color: black; font-weight: bold; border: solid; text-align: center; font-size:1.2em;"
    style_gt = "background-color: #C04000; color: black; font-weight: bold; border: solid; text-align:center;"
    gt_mode = x > st.mode(x)
    return [style_gt if i else style_lt for i in gt_mode]


def genomes_highlighter(x):
    style_genome = "background-color: white; color: black; font-weight: bold; border: solid; text-align: center; font-size:1.2em;"
    return [style_genome for i in x]

def apply_highlighter2(DF_heat,function, style_highlighter, headers_style):
    return DF_heat.style.apply(function,subset=["fam_1", "fam_2"]).apply(style_highlighter, subset="orgs_ids").hide(axis="index").set_table_styles([headers_style], overwrite = False)


def make_heat_map(criteria):
    import json
    
    json_file_Expfams = open(CTS.EXP_FAMS_JSON, "r")
    exp_fam_dictio = json.load(json_file_Expfams)
    json_file_Expfams.close()
    
    json_file_org = open(CTS.ORGS_IDS_NMS, "r")
    orgs_ids_dict = json.load(json_file_org)
    orgs_ids_list = list(orgs_ids_dict.keys())
    json_file_org.close()
    
    json_file_fams = open(CTS.FAMILIES, "r")
    fams_dict = json.load(json_file_fams)
    families = list(fams_dict.keys())
    json_file_fams.close()
    
    Dictio_DF = {"orgs_ids": orgs_ids_list}
    for fam in families:
        Dictio_DF[fam] = []
    for org in Dictio_DF["orgs_ids"]:
        for key in exp_fam_dictio:
            if org in key:
                fmly = key.split("|")[0]
                Dictio_DF[fmly].append(len(exp_fam_dictio[key]))
    
    #Dictio_DF
    DF_heat = pd.DataFrame.from_dict(Dictio_DF)
    #for indx in DF_heat.index:
    #    Id = DF_heat["organisms"][indx]
    #    DF_heat.at[indx,"organisms"] = orgs_ids_dict[Id]
    
    
    #DF_heat
    headers_style = {
    'selector': 'th.col_heading',
    'props': "background-color: OldLace; color: black; font-weight: bold; border: solid; text-align: center; font-size:1.2em;"
    }
    if criteria == "mean":
        return apply_highlighter2(DF_heat,mean_highlighter, genomes_highlighter, headers_style)
    elif criteria == "median":
        return apply_highlighter2(DF_heat,median_highlighter, genomes_highlighter, headers_style)
    elif criteria == "mean + std":
        return apply_highlighter2(DF_heat,std_highlighter, genomes_highlighter, headers_style)
    elif criteria == "mean + 2*std":
        return apply_highlighter2(DF_heat,std2_highlighter, genomes_highlighter, headers_style)
    elif criteria == "mode":
        return apply_highlighter2(DF_heat,mode_highlighter, genomes_highlighter, headers_style)
    else: 
        return print(f"{criteria} is not a valid criteria\n")
    
    

    
    
##################################### ##################################### #####################################      
##################################### Obtain the BGC and writing the sequences ##################################### 
#####################################   to the expanded families files       ##################################### 
##################################### ##################################### ##################################### 
    
    
def obtain_bgcs():
    """ This function obtains the bgc that were detected by blast
    from expanded families to MIBiG.
    """
    
    blast_files = [x for x in os.listdir(CTS.BLASTp_PATH + "exp_fam_to_nat_prods/") if x.endswith(".blast")]
    for blst_file in blast_files:
        blast_path = "/home/csar/Proyectos/Posdoc/Proyecto_pos/dev_package/data/data_bases/blastp/exp_fam_to_nat_prods/"

        df = pd.read_csv(blast_path + blst_file,sep = '\t', names = CTS.BLAST_COLS, index_col = False)
        bgc_recruited = list(df[df["bitscore"]>100]["subject"].unique())
        bgc_recruited.sort()
        
        output_file_name = CTS.EXPANDED_FAMS + blst_file.split("_to_")[0]+ ".fasta"
        mibig_db = CTS.MIBiG_NAT_PRODS_DB + "nat_prods.fasta"
        
        mibig_file = open(mibig_db, "r")
        out_file = open(output_file_name, "a")
        for line in mibig_file:
            li = line.strip()
            if ">" in li:
                li = li.split(">")[1]
            if li in bgc_recruited:
                #print("Se imprimio una línea\n")
                out_file.write(line)
                next_line = mibig_file.readline()
                out_file.write(next_line)
                bgc_recruited.remove(li)
            if bgc_recruited == []:
                #print(f"Se imprimieron todas las líneas del {blst_file}\n")
                out_file.close()
                mibig_file.close()
                break
                 
        if bgc_recruited != []:
            out_file.close()
            mibig_file.close()
            print(f"Algo salió mal con {blst_file} \n")
            print(li+"\n")
            print(bgc_recruited)
    return 



###############################################################################################################
##################################### Muscle Alignment ########################################################
#####################################                  #####################################
##########################################################################

def apply_muscle():
    """ This function alignes the sequences of the expanded families + mibig recruits.
    """
    
    import subprocess

    expanded_families_files = [x for x in os.listdir(CTS.EXPANDED_FAMS) if x.endswith(".fasta")]
    
    for exp_fam in expanded_families_files:
        in_file = CTS.EXPANDED_FAMS + exp_fam
        out_file = CTS.MUSCLE_OUTPUT + "MUSCLE_" + exp_fam
        subprocess.run([CTS.MUSCLE_EXE, "-align", in_file, "-output", out_file])
    return 
        

def obtain_tree():
    """This function returns a tree file with the aligned fasta. 
    """ 
    from Bio.Phylo.Applications import _Fasttree
    #import subprocess
     
    aligned_fastas = [x for x in os.listdir(CTS.MUSCLE_OUTPUT) if x.endswith(".fasta")]

    for aligned_fasta in aligned_fastas:
        in_file = CTS.MUSCLE_OUTPUT + aligned_fasta
        out_file = CTS.MUSCLE_OUTPUT + "TREE_" + aligned_fasta
        cmd = _Fasttree.FastTreeCommandline(CTS.FastTree_EXE, input = in_file, out = out_file)
        cmd()
        #subprocess.run([CTS.MUSCLE_EXE, "muscle -maketree -in %s -output %s" %(in_file, out_file)])
        #subprocess.run([CTS.FastTree_EXE, "FastTree out",out_file, in_file])
    return 
############################################################################################################
###################################### VIEW THE TREES ######################################################
############################################################################################################

def tree_view():
    """ This function shows the tree from the newick format file obtained with
    the obtain_tree function.
    """
    from ete3 import PhyloTree, TreeStyle
    
    tree_files = [x for x in os.listdir(CTS.MUSCLE_OUTPUT) if (x.startswith("TREE") and x.endswith(".fasta"))]
    for t in tree_files:
        algt_file = t.split("TREE_")[1]
        tree_obj = PhyloTree(t, alignment = algt_file, alg_format = "fasta")
        ts = TreeStyle()
        t.show(tree_style = ts)
    return 


