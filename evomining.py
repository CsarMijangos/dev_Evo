import evo_constants as CTS
import evo_helpers as hlp

from Bio.Blast.Applications import NcbimakeblastdbCommandline

class Evomining:
    """This class creates an evomining object.
    Its methods permit to access to the evomining
    funtionality.
    """

    def __init__(self, centralDB,
                evo_genomesDB,
                nat_prodsDB):
        self.central = centralDB
        self.genomesDB = evo_genomesDB
        self.nat_prods = nat_prodsDB
        
    def _make_genomes_blast_db(self):
        aux = [x for x in os.listdir(self.genomesDB) if ".fasta" in x]
        
        input_db = self.genomesDB + aux[0]
        output_name = aux[0].replace(".fasta", "") + "_BlastDB"
        cline_db = NcbimakeblastdbCommandline(cmd = CTS.MAKEBLASTDB_CMD,
                                    dbtype="prot",
                                   input_file= input_db,
                                     out = output_name)
        return
        
      
        
    





if __name__ == "__main__":
    hlp.welcome()
    opt = hlp.option_menu()
    if opt == 1:
        hlp.make_all_evo_headers(CTS.EXMPL_DB["genomes"],
                                 CTS.EXMPL_DB["rast_ids"])
        hlp.join_evo_headers_files()
        evo_genomes_path = CTS.GENOMES_EvoFmt
        
        evomining = EvoMining(CTS.EXMPL_DB["central"],
                          evo_genomes_path,
                          CTS.EXMPL_DB["nat_prods"])
        evomining.run()
        
    else: 
        user_dbs = hlp.get_dbs()
        hlp.make_all_evo_headers(user_dbs["genomes"], 
                             user_dbs["rast_ids"])
        hlp.join_evo_headers_files()
        
    
        user_dbs["genomes"] = CTS.GENOMES_EvoFmt
    
        evomining = EvoMining(user_dbs["central"],
                         user_dbs["genomes"],
                         user_dbs["nat_prods"])
        evomining.run()
        
        
    