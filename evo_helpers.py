import constant as cts


######################## 
####     function to 
###  say welcome 
########################

def welcome():
    print("Welcome to EvoMinig\n")
    input("Press enter to continune ...")
    

    
    
######################## 
####     function to 
###  show the options 
########################
#OPTION = 0
def option_menu():
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
    central_db = input("\n Enter path to central_DB: ")
    genomes_db = input("\n Enter path to genomes_DB: ")
    rast_ids = input("\n Enter path to rast_ids: ")
    nat_prods_db = input("\n Enter path to nat_prods_DB: ")
    DBs = { "central" : central_db,
           "genomes" : genomes_db,
           "rast_ids" : rast_ids,
           "nat_prods" : nat_prods_db }
    return DBs

















    