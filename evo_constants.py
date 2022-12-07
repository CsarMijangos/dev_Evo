#################################################################                        
               ############# DATA BASES #########
#################################################################
    
# Data bases of the precharged example:

#Enzime family DB:
EXMPL_CENTRAL_DB = "data/data_bases/example/central/central.fasta"   

# Genome DB:
EXMPL_GENOMES_DB = "data/data_bases/example/genomes/los17/"  

# Rast_ids file:
EXMPL_RAST_IDS = "data/data_bases/example/genomes/los17Rast.ids" 

#Natural products DB(default =MIBiG):
MIBiG_NAT_PRODS_DB = "data/data_bases/example/mibig/nat_prods.fasta" 

##Dictionary with the all the paths:
EXMPL_DB = { "central" : EXMPL_CENTRAL_DB,
            "genomes" : EXMPL_GENOMES_DB,
            "rast_ids" : EXMPL_RAST_IDS,
            "nat_prods" : MIBiG_NAT_PRODS_DB}




###################################################################
                ############ Paths to auxiliary files ######
                ############ used by evomining ############
###################################################################


# This directory will contain the auxiliary files 
# to obtain the genomes .faa files with the evomining format:
        
genomes_EvoFmt = "data/aux_data/genomes_EvoFmt/"  





###################################################################
#########################    Paths to the   #######################
#########################  commands in the  #######################
#########################    local blast    #######################
#########################   installation    #######################
###################################################################

## path to the command for makeblastdb:
makeblastdb_cmd = "/home/csar/ncbi-blast-2.13.0+/bin/makeblastdb"

## path to the command for blastp:
blastp_cmd = "/home/csar/ncbi-blast-2.13.0+/bin/blastp"



