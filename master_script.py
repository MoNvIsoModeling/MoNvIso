"""
MoNvIso, the Modeling eNvirnment for Isoforms.

 Francesco Oliva
"""

# Libraries for the main
import time, sys, os
import library as lib

start_time = time.time()
master_directory = os.getcwd()

if "\\" in master_directory:  # if running on windows
    slash = "\\"


if "/" in master_directory:   # if running on mac/linux
    slash = "/"

AA_dict = dict(ALA="A", ARG="R", ASN="N", ASP="D", ASX="B", CYS="C", GLU="E", GLN="Q", GLX="Z", GLY="G", HIS="H", ILE="I",
          LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V", MSE="M")

reportname = "report.log"
csvname = "report.csv"


mutations_list = "mutations.txt"    # name of the file where gene names and mutations are listed
blocks, protein_list = lib.parse_input(mutations_list)   # process the mutations.txt file

parameters = lib.get_parameters()   # processes the parameters
print(parameters)

lib.make_directories(blocks, master_directory, slash)   # build the directories
lib.get_isoforms_from_db(protein_list, master_directory, str(parameters["DB_LOCATION"]), slash)  # scan the DB for the possible isoforms

genes = []
for gene in protein_list:  # starts the real work
    newpath = master_directory + slash + str(gene) + slash
    os.chdir(newpath)  # change to the folder of the gene in process

    lib.master_isoform(mutations_list, newpath)
    check_output = lib.run_hmm(master_directory, newpath, gene, slash, float(parameters["RESOLUTION"]),
                               float(parameters["SEQID"]), str(parameters["HMMER_HOME"]),
                               str(parameters["COBALT_HOME"]), int(parameters["PDB_TO_USE"]))

    report = open(reportname, "a+")
    csvfile = open(csvname, "a+", newline='')
    if os.stat("usable_isoforms.txt").st_size != 0:
        with open(reportname, "a") as report, open(csvname, "a", newline='') as csvfile:
            winning_iso = lib.score_isoform(AA_dict, gene, mutations_list, report, csvfile, slash)
    else:
        winning_iso = 42
    report.close()
    csvfile.close()
    if winning_iso != 42:
        print(f"The winning isoform for gene {gene} is: {winning_iso}")
        pdbs_and_chains = lib.prepare_modeller(winning_iso, gene)
        pdbs = lib.prepare_alignments(pdbs_and_chains, winning_iso, gene, int(parameters["MODEL_CUTOFF"]))
        lib.write_modeller(pdbs, gene, str(parameters["NUM_OF_MOD_WT"]))
        lib.run_modeller(gene, slash)
        os.chdir("..")
        lib.cleaning(gene, pdbs)
        os.chdir(newpath+winning_iso)
        lib.insert_mutations(AA_dict, winning_iso, str(parameters["NUM_OF_MOD_MUT"]), slash, gene)
        os.chdir("../..")
    os.chdir(master_directory)
    lib.final_cleaning(gene, master_directory)
    genes.append(gene)
lib.problematic_genes(genes)
lib.final_report()

print("Execution time: %s seconds" % ("{:.2f}".format(time.time() - start_time)))
