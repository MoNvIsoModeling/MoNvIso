import shutil
import re
import subprocess
import numpy as np
import sys, os
from os import path
import time
from Bio.PDB import *  # tbd
from Bio import AlignIO, Align, pairwise2, SeqIO  # tbd
from Bio.Blast import NCBIWWW as blastq
from Bio.Blast import NCBIXML as blastparser
import psutil  # tbd
import csv


class Chain():
    # how to use the class and access the info:
    # for pdb in data: # run through the pdbs
    #     print(pdb) # print the pdb name
    #     for chain in range(len(data[pdb])): # run through the chains in a pdb
    #         print(data[pdb][chain].chainid) # print the chain id
    #         print(data[pdb][chain].score/coverage/seqid) # print the score/coverage/seqid of the chain after alignment
    #         for interval in range(len(data[pdb][chain].ranges)): # run through the intervals covered by the chain
    #             print(str(data[pdb][chain].ranges[interval].start)) # print the beginning of the interval
    #             print(str(data[pdb][chain].ranges[interval].end)) # print the end of the interval

    def __init__(self, chainame, score=None, coverage=None, seqid=None):
        self.chainid = chainame
        self.score = score
        self.coverage = coverage
        self.seqid = seqid
        self.ranges = []

    def add_range(self, range):
        self.ranges.append(range)

    def add_score(self, score):
        self.score = score

    def add_coverage(self, coverage):
        self.coverage = coverage

    def add_seqid(self, seqid):
        self.seqid = seqid


class range1():
    def __init__(self, str):
        x = str.split("-")
        self.start = int(x[0])
        self.end = int(x[1])


class Iso_score:
    def __init__(self, mutation_hits, notfound, neigh_list, name, modelmut, cantmodel, final_mut_score, cover_score,
                 totalscore):
        self.hits = mutation_hits  # AA mapped on the isoform
        self.notfound = notfound  # AA NOT mapped on the isoform
        self.neighbours = neigh_list  # First neighbours ( AA +- 2)
        self.name = name  # Isoform name (i.e. isoform1)
        self.modelmut = modelmut  # AA that can be modelled
        self.totalscore = totalscore  # total score of the isoform = ((modelmut/total mut found in all isoforms) + (total covered AA in the isoform/total AA in the isoform))*10
        self.cantmodel = cantmodel  # AA that are mapped but cannot be modelled


def get_parameters():
    keywords = {
        "RESOLUTION": None,
        "SEQID": None,
        "PDB_TO_USE": None,
        "DB_LOCATION": None,
        "HMMER_HOME": None,
        "COBALT_HOME": None,
        "MODEL_CUTOFF": None,
        "NUM_OF_MOD_WT": None,
        "NUM_OF_MOD_MUT": None
    }
    keys = list(keywords)
    filename = "parameters.dat"
    if os.path.isfile(filename) == False:
        print("PARAMETERS FILE NOT FOUND!")
        sys.exit()
    with open(filename, "r") as f:
        lines = f.readlines()

    for key in keys:
        for line in lines:
            if key in line:
                value = line[line.find("=") + 1:].strip()
                keywords[key] = value
    return keywords


def parse_input(mutation_list):
    with open(mutation_list, "r") as f:
        content = f.read()
    blocks = []
    protein_list = []
    for block in content.split("\n\n"):
        blocks.append(block.splitlines())

    block = 0

    while block < len(blocks):
        protein_list.append(blocks[block][0].upper())
        block += 1

    return blocks, protein_list


def make_directories(blocks, master_directory, slash):
    block = 0

    while block < len(blocks):
        gene_name = blocks[block][0].upper()
        if os.path.exists(gene_name):
            shutil.rmtree(gene_name)
        os.mkdir(gene_name)

        path = master_directory + slash + str(gene_name)

        newfile = open(str(path) + slash + "mutations.txt", "w")

        mutation = 1

        while mutation < len(blocks[block]):
            newfile.write(str(blocks[block][mutation]) + "\n")
            mutation += 1

        block += 1


def get_isoforms_from_db(protein_list, master_directory, db_location, slash):
    # the first part deals with the canonical isoform which is located in uniprot_sprot.fasta
    canonical_db = str(db_location) + "uniprot_sprot.fasta"
    isoforms_db = str(db_location) + "uniprot_sprot_varsplic.fasta"
    isoforms = parse_isoform(canonical_db)

    for gene in protein_list:
        name = "GN=" + str(gene.upper()) + " "
        counter = 0
        while counter < len(isoforms):
            if name in str(isoforms[counter][0]) and "Homo" in str(isoforms[counter][0]):
                path = master_directory + slash + str(gene) + slash
                filename = str(path) + "isoform0.fasta"
                files = open(filename, "w")
                files.write(">" + str(isoforms[counter][0]) + "\n")
                for line in isoforms[counter][1:]:
                    files.write(line)
                    if line != isoforms[counter][-1]:
                        files.write("\n")
                files.close()
            counter += 1

    # the second part deals with all the other isoforms which are in the file called uniprot_sprot_varsplic.fasta
    isoforms = parse_isoform(isoforms_db)

    for gene in protein_list:
        name = "GN=" + str(gene.upper())
        numb = len(name)
        counter = 0
        while counter < len(isoforms):
            if isoforms[counter][0][-numb:] == name and "Homo" in isoforms[counter][0]:
                path = master_directory + slash + str(gene) + slash
                words = isoforms[counter][0].split(" ")
                filename = str(path) + "isoform" + str(words[2]) + ".fasta"
                files = open(filename, "w")
                files.write(">" + str(isoforms[counter][0]) + "\n")
                for line in isoforms[counter][1:]:
                    files.write(line)
                    if line != isoforms[counter][-1]:
                        files.write("\n")
                files.close()

            counter += 1


def parse_isoform(infile):
    with open(infile, "r") as f:
        content = f.read()

    return [block.splitlines() for block in content.split(">") if "sp|" in block]


def master_isoform(mutation_list, path):
    parse_input(mutation_list)  # I'm not sure this does anything anymore
    isolist = os.listdir(path=path)
    for i in range(len(isolist)):
        if isolist[i] == "mutations.txt":
            isolist.pop(i)
            break
    isolist.sort()
    miso = open("master_isoform.txt", "w")
    numbers = []
    for i in isolist:
        if str(i[:7]) == "isoform":
            try:
                numbers.append(int(i[7:-6]))
            except:
                numbers.append(i[7:-6])
    try:
        numbers.sort()
    except:
        pass
    for i in numbers:
        miso.write(f"isoform{str(i)}.fasta\n")


def standardize(AA_dict, removes_multiple, mutations_list):  # this function could be improved using functions like "isnumber" or "isalpha" or "isalphanum"
    g = open(mutations_list, "r")
    shadowlist = g.readlines()
    mutlist = []
    keys = list(AA_dict)
    for i in shadowlist:
        appo = re.split('(\d+)', i)
        mutlist.append(appo)
    for i in range(len(mutlist)):
        mutlist[i][0] = re.sub(r'[^a-zA-Z]', "", mutlist[i][0])  # leaves only letters
        mutlist[i][2] = re.sub(r'[^a-zA-Z]', "", mutlist[i][2])
        x = len(mutlist[i][0])
        if x > 1:
            for j in keys:
                if j == mutlist[i][0].upper():
                    mutlist[i][0] = AA_dict[j]
                if j == mutlist[i][2].upper():
                    mutlist[i][2] = AA_dict[j]
    g.close()
    sorted(mutlist, key=lambda x: int(x[1]))

    # removes multiple mutations of the same residue
    shortlist = []
    if removes_multiple == 1:
        shortlist = [mutlist[0]]
        line = 1
        while line < len(mutlist):
            if mutlist[line][1] != mutlist[line - 1][1]:
                shortlist.append(mutlist[line])
            line += 1
    return mutlist, shortlist


def create_aalist(isoin):  # could be updated
    isowork = isoin.readlines()[1:]  # load the sequence (but it's uncomfortable to use)
    for i in range(len(isowork)):
        isowork.append(isowork[i])  # appends the AA as elements of the list (comfortable to use)
    todel = len(isowork) // 2
    del isowork[:todel]

    isowork = "".join(isowork).replace("\n", "")  # removes the "\n" elements
    return list(isowork)


def generate_mut_dict(mutlist):
    return {'%s' % (i[0] + i[1]): 0 for i in mutlist}


def get_neighbours(mutlist, aalist, aa):
    neigh = []
    maxnumber = len(aalist)
    if 2 < int(aa[1]) < int(maxnumber) - 2:
        neigh.append(aalist[int(aa[1]) - 3])
        neigh.append(aalist[int(aa[1]) - 2])
        neigh.append(aalist[int(aa[1])])
        neigh.append(aalist[int(aa[1]) + 1])
    if int(aa[1]) == 1:
        neigh.append(aalist[int(aa[1])])
        neigh.append(aalist[int(aa[1]) + 1])
    if int(aa[1]) == 2:
        neigh.append(aalist[int(aa[1]) - 2])
        neigh.append(aalist[int(aa[1])])
        neigh.append(aalist[int(aa[1]) + 1])
    if int(aa[1]) == int(maxnumber):
        neigh.append(aalist[int(aa[1]) - 3])
        neigh.append(aalist[int(aa[1]) - 2])
    if int(aa[1]) == int(maxnumber) - 1:
        neigh.append(aalist[int(aa[1]) - 3])
        neigh.append(aalist[int(aa[1]) - 2])
        neigh.append(aalist[int(aa[1])])

    return neigh


def get_isonames(isolist, slash):
    isonames = []
    for i in isolist:
        i.replace(".fasta", "")
        isonames.append(i)
    return isonames


def parser_aligned(alignment, slash):
    with open(alignment, "r") as f:
        content = f.read()

    # keeps only the pdb entries and the target sequence

    newblock = []
    for block in content.split(">"):
        block = ">" + block
        if "|1 sp|" in block:
            newblock.append(block.splitlines())
        elif "_clean" in block:
            newblock.append(block.splitlines())

    # writes a file with the target sequence at the bottom (since the rest of the code previously written assumes it at the bottom)
    newfile = "alignment_pdbs.fasta"
    with open(newfile, "w") as f:
        blockcounter = 1
        while blockcounter < len(newblock):
            for line in newblock[blockcounter]:
                f.write(str(line))
                f.write("\n")
            f.write("\n")
            blockcounter += 1
        for line in newblock[0]:
            f.write(str(line))
            f.write("\n")

    # open the alignment_pdbs.fasta file and reads it

    with open(newfile, "r") as f:
        content = f.read()

    newblock = []
    for block in content.split("\n\n"):
        newblock.append(block.splitlines())

    # replace the extra AAs in the PDB with -

    o = open("temp_align.fasta", "w")
    #    newblock.pop()
    for block in newblock:
        if block != newblock[-1]:
            o.write(block[0] + "\n")
            line = 1
            while line < len(block):
                aa = 0
                while aa < len(block[line]):
                    safe = 0
                    if block[line][aa] != "-" and newblock[-1][line][aa] == "-":
                        o.write("-")
                        safe = 1
                    if safe == 0:
                        o.write(block[line][aa])
                    aa += 1
                o.write("\n")
                line += 1
            o.write("\n")
    o.close()

    # Remove every - that appear in both the PDBs and the reference
    with open("temp_align.fasta", "r") as i:
        newcontent = i.read()
    o = open("cleaned_align.fasta", "w")
    newopen = open("aligned_oneline.fasta", "w")
    thirdopen = open("aligned_oneline_withgaps.fasta", "w")
    theblock = []
    for block in newcontent.split("\n\n"):
        theblock.append(block.splitlines())
    theblock.pop()
    for block in theblock:
        o.write(block[0] + "\n")
        newopen.write(block[0] + "\n")
        thirdopen.write(block[0] + "\n")
        line = 1
        while line < len(block):
            aa = 0
            while aa < len(block[line]):
                safe = 0
                thirdopen.write(block[line][aa])
                if block[line][aa] == "-" and newblock[-1][line][aa] == "-":
                    safe = 1
                if safe == 0:
                    o.write(block[line][aa])
                    newopen.write(block[line][aa])
                aa += 1
            line += 1
            if block != newblock[-1]:
                o.write("\n")

        o.write("\n")
        newopen.write("\n\n")
        thirdopen.write("\n\n")

    newline = 1
    o.write(newblock[-1][0] + "\n")
    newopen.write(newblock[-1][0] + "\n")
    thirdopen.write(newblock[-1][0] + "\n")
    while newline < len(newblock[-1]):
        newaa = 0
        while newaa < len(newblock[-1][newline]):
            thirdopen.write(newblock[-1][newline][newaa])
            if newblock[-1][newline][newaa] != "-":
                o.write(newblock[-1][newline][newaa])
                newopen.write(newblock[-1][newline][newaa])
            newaa += 1
        o.write("\n")
        newline += 1
    o.close()


def evaluate_modellability(path, isoform, mut, slash):
    openfile = str(path) + slash + str(isoform) + slash + "covered_aa.txt"
    with open(openfile, "r") as f:
        content = f.read()

    newblock = []
    begins = []
    ends = []
    for block in content.split("\n"):
        newblock.append(block.split("\t"))
    i = 0
    while i < len(newblock):
        begins.append(newblock[i][1])
        ends.append(newblock[i][2])
        i += 1
    i = 0
    ranger = []
    while i < len(begins):
        test = range(int(begins[i]), int(ends[i]))
        for x in test:
            ranger.append(x)
        i += 1

    if int(mut) in ranger:
        modelmut = "yes"
    else:
        modelmut = "no"
    return modelmut


def sort_list(sub_list):  # this sort a list according to the second element (i.e. a list that list the mutations has the res number as second element of the sublist)
    return sorted(sub_list, key=lambda x: int(x[1]))


def coverage_score(path, aa_list, isoform, slash):  # the coverage score is defined as the number of AA covered by the templates divided by the total amount of AA in the isoform and multiplied by 10
    total_aa = len(aa_list)
    namein = str(path) + slash + "isoform" + str(isoform) + slash + "covered_aa.txt"
    with open(namein, "r") as f:
        content = f.read()

    newblock = []
    for block in content.split("\n"):
        newblock.append(block.split("\t"))
    newblock = sort_list(newblock)
    begins = []
    ends = []
    for block in newblock:
        begins.append(block[1])
        ends.append(block[2])

    covered_aa = 0
    counter = 0
    aa = 1
    global_interval = []
    while counter < len(begins):
        interval = np.arange(int(begins[counter]), int(ends[counter]) + 1)
        global_interval.append(interval)
        counter += 1
    global_interval = np.unique(np.hstack(global_interval))
    while aa <= total_aa:
        if aa in global_interval:
            covered_aa += 1
        aa += 1
    score = "{:.1f}".format(10 * covered_aa / total_aa)
    return score


def write_report(isonames, protein, max_mut, mut_score, max_mut_found, report, winning_iso, mutlist):
    residues = [str(res[0]) + str(res[1]) for res in mutlist]
    if str(winning_iso.upper()) == "ISOFORM0":
        report.write(
            f"## {str(protein)} ##"
            + "\n\nWinning isoform: CANONICAL (isoform 0)\n\n"
        )

    else:
        report.write(
            f"## {str(protein)} ##"
            + "\n\nWinning isoform: "
            + str(winning_iso.upper())
            + "\n\n"
        )

    count_res = 1
    report.write("Possible mutating residues: \t")
    for res in residues:
        report.write(str(res))
        count_res += 1
        if res != residues[-1]:
            report.write(", ")
        else:
            report.write("\n\n")

        if count_res == 10:
            report.write("\n\t\t\t\t")
            count_res = 1

    writeit = "no"
    first = 0
    for i in mut_score:
        if mut_score[i] == 0:
            if first == 0:
                report.write("\nResidue(s) ")
                first = 1
            report.write(f'{i}, ')
            writeit = "yes"
    if writeit == "yes":
        report.write(
            (
                (
                    f"is (are) never found. Thus only {str(max_mut_found)}"
                    + " residues are considered for the score calculation instead of "
                )
                + str(max_mut)
                + "\n\n"
            )
        )


    for prot in isonames:
        no_bueno = open(f"missing_mutations_{str(prot.name)}.dat", "w")

        if str(prot.name.upper()) == "ISOFORM0":
            report.write("      CANONICAL ISOFORM\n-----------------------------\n\n")
        else:
            report.write(f"      {str(prot.name.upper())}" + "\n--------------------\n\n")
        report.write("\t\tMutating residues mapped\t\t\t\t" + str(prot.hits) + "/" + str(max_mut_found) + "\n")
        report.write(
            "\t\tMutating residues that can be modelled\t\t\t" + str(prot.modelmut) + "/" + str(prot.hits) + "\n")

        if prot.notfound != []:
            report.write("\n\t\t")
            count_res = 1
            report.write("NOT mapped on the isoform:\n\t\t")
            for nofound in prot.notfound:
                report.write(f'{str(nofound)}, ')
                count_res += 1
                if count_res == 10:
                    report.write("\n\t\t")
                    count_res = 1
                no_bueno.write(str(nofound) + "\n")
        if int(prot.hits) != 0:
            report.write("\n\n\t\tResidues MAPPED on the isoform:\n\t\t")
        count_res = 1
        for res in residues:
            if prot.notfound == []:
                report.write(str(res))
                count_res += 1
                if res != residues[-1]:
                    report.write(", ")


            else:
                nofound = 0
                if str(res) not in str(prot.notfound):
                    report.write(str(res))
                    if res != residues[-1]:
                        report.write(", ")
                    else:
                        report.write("\n")
                    count_res += 1
            if count_res == 10:
                report.write("\n\t\t")
                count_res = 1

        if prot.cantmodel != []:
            report.write("\n")
            counter = 0
            report.write("\t\tNOT found in the crystal structures: ")
            count_res = 1
            for notmodellable in prot.cantmodel:
                report.write(f'{str(notmodellable)}, ')
                count_res += 1
                if count_res == 10:
                    report.write("\n\t\t")
                    count_res = 1
                no_bueno.write(str(notmodellable))
                if counter < int(len(prot.cantmodel)) - 1:
                    no_bueno.write("\n")
            counter += 1
            report.write("\n\n")
            no_bueno.close()

        report.write("\n\tSCORES:\n\t\tMutations: " + str("{:.1f}".format(prot.final_mut_score)) + "\tCoverage: " + str(
            prot.cover_score) + "\tTotal: " + str("{:.1f}".format(prot.totalscore)) + "\n")
        report.write("\n\n")


def write_csv(isonames, protein, mutlist, csvfile):
    writer = csv.writer(csvfile)
    writer.writerow([protein])
    listofiso = [""]
    for prot in isonames:
        if str(prot.name) == "isoform0":
            listofiso.append("canonical")
        else:
            listofiso.append("isoform " + str(prot.name[7:]))

    writer.writerow(listofiso)

    for mut in mutlist:
        mutline = []
        mutation = mut[0] + mut[1]
        mutline.append(mutation)
        for prot in isonames:
            if mutation in prot.notfound:
                mutline.append("")
            else:
                mutline.append("x")
        writer.writerow(mutline)

    writer.writerow("")


def score_isoform(AA_dict, protein, mutations_list, report, csvfile, slash):
    removes_multiple = 1
    mutlist = standardize(AA_dict, removes_multiple, mutations_list)[1]
    f = open("usable_isoforms.txt", "r")
    # f = open("master_isoform.txt", "r")
    iso_list = f.read().strip().splitlines()
    isonames = get_isonames(iso_list, slash)  # name of the isoform
    path = os.getcwd()
    mut_score = generate_mut_dict(mutlist)  # keep track if a mutation is found (not how many times)
    count = 0
    aa_lists = []
    numbers = []
    for i in isonames:
        numbers.append(i[7:-6])
        cantmodel = []
        cantfind = []
        neigh_list = []
        iso_score = 0
        g = open(i, "r")  # it opens the fasta files of the isoform
        aalist = create_aalist(g)
        aa_lists.append(aalist)
        modelmut = 0
        for j in mutlist:
            safe = 0
            try:
                aalist[int(j[1]) - 1]
            except IndexError:
                safe = 1
            if safe == 0 and str(j[0]) == str(aalist[int(j[1]) - 1]):
                iso_score += 1
                mutname = j[0] + j[1]
                mut_score[mutname] = 1
                neigh_list.append(get_neighbours(mutlist, aalist, j))
                mutmodel = evaluate_modellability(path, i[:-6], j[1], slash)
                if mutmodel == "yes":
                    modelmut += 1
                if mutmodel == "no":
                    cantmodel.append(j[0] + j[1])
            else:
                cantfind.append(j[0] + j[1])

        isonames[count] = Iso_score(iso_score, cantfind, neigh_list, i[:-6], modelmut, cantmodel, None, None, None)
        # isonames is now a list that contains the objects with the scores and neighbours list
        count += 1

    # TODO Aggiungere i vari tipi di warnings

    max_mut = len(list(mut_score))

    never_found = 0
    for i in mut_score:
        if mut_score[i] == 0:
            never_found += 1

    max_mut_found = max_mut - never_found

    for i in range(len(isonames)):
        if max_mut_found != 0:
            isonames[i].final_mut_score = 10 * int(isonames[i].modelmut) / max_mut_found
        else:
            isonames[i].final_mut_score = 0
        isonames[i].cover_score = coverage_score(path, aa_lists[i], numbers[i], slash)
        isonames[i].totalscore = float(isonames[i].cover_score) + isonames[i].final_mut_score

    ordered_dict = {}
    count = 0
    for i in range(len(isonames)):
        ordered_dict['%s' % isonames[i].name] = isonames[
            i].totalscore  # this part orders the isoform in decrescent order from the one with the highest (totalscore) score to the lowest to select the winning isoform
    for w in sorted(ordered_dict, key=ordered_dict.get, reverse=True):
        if count == 0:
            winning_iso = w
        count += 1

    # if os.stat("usable_isoform.txt").st_size == 0:
    #     winning_iso = "There are no usable isoform for this gene. You may want a meatbag to look at it"

    write_report(isonames, protein, max_mut, mut_score, max_mut_found, report, winning_iso, mutlist)
    write_csv(isonames, protein, mutlist, csvfile)

    # if os.stat("usable_isoform.txt").st_size == 0:
    #     winning_iso = 42

    return winning_iso


def write_modeller2(mutation, gene, pdbs, number_of_models, slash):
    # local_slave may need an update to local_worker
    pdb_location = os.getcwd() + slash
    directory = pdb_location + "mutation_" + str(mutation) + slash
    o = open(str(directory) + "run_modeller.py", "w")
    o.write(
        "from modeller import *\nfrom modeller.automodel import *\nfrom modeller.scripts import complete_pdb\nfrom modeller.parallel import *\nimport os\n\n")
    o.write("log.verbose()\nj=job()\ncores_number = 8\nthread = 0\n")
    o.write("while thread < cores_number:\n    j.append(local_slave())\n    thread += 1\n")
    o.write("env = environ()\nenv.io.atom_files_directory = ['" + str(pdb_location) + "']\nenv.io.hetatm = True\n")
    o.write("a = automodel(env,\n    alnfile =\"" + str(directory) + "input_modeller.dat\",\n    knowns = (" + str(
        pdbs) + ") ,\n    sequence = \"" + str(
        gene) + "_" + mutation + "\",\n    assess_methods=(assess.DOPE, assess.GA341))\n")
    o.write("a.starting_model= 1\na.ending_model  = " + number_of_models + "\na.use_parallel_job(j)\na.make()\n")
    o.write(
        "ok_models = filter(lambda x: x['failure'] is None, a.outputs)\ntoscore = 'DOPE score'\nok_models = sorted(ok_models, key=lambda k: k[toscore])\nm = ok_models[0]\nmyout = open(\"MYOUT.dat\", \"w\")")
    o.write("\nmyout.write(\"Top model: \" + str(m['name']) + \" (DOPE SCORE: %.3f)\" % (m[toscore]))\n")
    o.write(
        "env.libs.topology.read(file='$(LIB)/top_heav.lib')\nenv.libs.parameters.read(file='$(LIB)/par.lib')\nmdl = complete_pdb(env, m['name'])\ns = selection(mdl)\n")
    o.write("s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=\"" + str(
        gene.upper()) + mutation + "_profile.dat\", normalize_profile=True, smoothing_window=15)")


def insert_mutations(AA_dict, winning_iso, number_of_models, slash, gene):

    with open("final_msa.fasta", "r") as f:
        block = f.readlines()
    isoform = []
    for line in block:
        if ">lcl|2 sp|" in line:
            break
        if ">lcl" not in line:
            isoform.append(line.strip())
    original_fasta = list(''.join(isoform))
    removes_multiple = 0
    mutations_list = "../mutations.txt"
    mutlist = standardize(AA_dict, removes_multiple, mutations_list)[0]
    f = open("../missing_mutations_" + winning_iso + ".dat")
    bad_mutlist = f.readlines()
    cool_mutlist = []
    f.close()

    # step 1: select the mutations that can be modeled
    for mut in mutlist:
        found = "No"
        good = str(mut[0]) + str(mut[1])
        for bad in bad_mutlist:
            if good == bad.strip():
                found = "Yes"
                break
        if found == "No":
            cool_mutlist.append(mut)

    # step 2 apply the mutation
    with open("aligned_oneline_withgaps.fasta", "r") as f:
        content = f.read()
    blocks = []
    for block in content.split("\n\n"):
        blocks.append(block.splitlines())

    with open("input_modeller.dat", "r") as f:
        content = f.read()
    blocks = []
    for block in content.split("\n\n"):
        blocks.append(block.splitlines())

    input_modeller = blocks[-1][2:]

    ref_sequence = []
    for lines in input_modeller:
        for aa in lines:
            ref_sequence.append(aa.strip())
    ref_sequence.pop()


    for j in cool_mutlist:
        fastain = original_fasta.copy()
        mutname = str(j[0]) + str(j[1]) + str(j[2])
        os.mkdir("mutation_" + mutname)
        count_aa = 1
        safe_counter = 0
        while safe_counter < len(fastain):
            if count_aa == int(j[1]):
                fastain[safe_counter] = str(j[2])
                break
            if fastain[safe_counter] != "-":
                count_aa += 1
            safe_counter += 1

        counter = 0
        while counter < len(ref_sequence):
            if ref_sequence[counter] == "/":
                fastain.insert(counter, "/")
            counter += 1
        counter = 0
        while counter < len(ref_sequence):
            if ref_sequence[counter] == "-":
                fastain[counter] = "-"
            counter += 1
        fastain.append("*")

        # copy old input and replace the wt with the mutant
        newinput = open("mutation_" + mutname + slash + "input_modeller.dat", "w")
        with open("input_modeller.dat", "r") as f:
            content = f.read()
            for block in content.split("\n\n"):
                blocks.append(block)

        for block in blocks:
            if ">P1;" + str(gene) in block[0]:
                break
            else:
                for line in block:
                    newinput.write(str(line))
                    newinput.write("\n")
                newinput.write("\n\n")

        newinput.write(">P1;" + gene + "_" + mutname + "\n" + "sequence:" + gene + "_" + mutname + ":xxx:x:xxx:x::::\n")

        for aa in fastain:
            newinput.write(aa)

        newinput.close()

        all_files = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]
        pdbs = []
        for name in all_files:
            if "_clean.pdb" in name:
                pdbs.append(name[:-4])
        write_modeller2(mutname, gene, pdbs, number_of_models, slash)
        directory = os.getcwd() + slash + "mutation_" + mutname
        os.chdir(directory)
        run_modeller(gene, slash)
        os.chdir("..")


def chainids(chains):  # unused
    chain_ID = []
    for i in chains:
        test = str(i)
        index1 = test[test.find('='):]
        index2 = index1[:index1.find('>')]
        index3 = index2[index2.find('Waldo')]
        chain_ID.append(index3)
    return chain_ID


def run_hmm(master_directory, newpath, gene, slash, res_cutoff, seqid_cut, hmmer_home, cobalt_home, max_pdbs):
    with open("master_isoform.txt", "r") as f:
        iso_list = f.read().strip().splitlines()
    isonames = get_isonames(iso_list, slash)
    with open("usable_isoforms.txt", "w") as f:
        for iso in isonames:
            if os.stat(iso).st_size == 0:
                os.remove(iso)
            else:
                outdir = newpath + str(iso[:-6])
                os.mkdir(outdir)  # make the output directories
                os.chdir(outdir)
                myfasta = import_sequence(f"../{str(iso)}")
                start_query(myfasta, gene)
                time.sleep(10)
                check_cobalt = run_cobalt(gene, cobalt_home)
                if not check_cobalt:
                    build_hmm(gene, hmmer_home)
                    check_output = search_templates(gene, max_pdbs)
                    if not check_output:
                        pdbs, chains = get_templates()
                        toremove = []
                        for pdbs_to_remove_positions_counter, (pdb, chain) in enumerate(zip(pdbs, chains)):
                            pdb_exist = download_pdb(pdb)
                            if pdb_exist:
                                check_pdb = manage_pdb(pdb, chain)  # NMR structures (no resolution) will be removed
                            else:
                                check_pdb = False
                            if not check_pdb:
                                toremove.append(pdbs_to_remove_positions_counter)
                        if toremove:
                            for counter, position in enumerate(toremove):
                                pdbs.pop(position-counter)
                                chains.pop(position-counter)
                        for pdb, chain in zip (pdbs, chains):
                             get_fasta(pdb, chain)
                        merge_fastas(pdbs, chains, gene)
                        align_templates_to_msa(gene, cobalt_home)
                        parser_aligned("final_msa.fasta", slash)
                        used_beginnings, used_ends, used_chains, used_pdbs = run_templates_analysis(res_cutoff, seqid_cut,
                                                                            slash)
                        if used_beginnings == 42:
                            check_output = True
                    os.chdir("..")
                    if not check_output:
                        f.write(iso)
                        if iso != iso[-1]:
                            f.write("\n")
            else:
                check_output = True
    return check_output

def import_sequence(sequencefile):
    return SeqIO.read(sequencefile, "fasta")


def start_query(myfasta, gene):
    print("Looking for homologues")
    results = blastq.qblast("blastp", "swissprot", myfasta.seq, alignments=500, word_size=6)
    blastRecord = blastparser.read(results)
    #  write hit file:
    with open(f"{gene}_hits.fasta", "w") as f:
        for alignment in blastRecord.alignments:
            f.write(f">{str(myfasta.id)}\n")
            f.write(str(myfasta.seq).strip() + "\n\n")
            for hsp in alignment.hsps:
                f.write(f">{alignment.hit_id}\n")
                f.write(str(hsp.sbjct).replace("-", "") + "\n\n")
    print("Done")


def run_bash_command(command):
    subprocess.run(command, shell=True, universal_newlines=True, check=True)


def run_cobalt(gene, cobalt_home):
    if f"{gene}_hits.fasta" != 0:
        print("Doing MSA with cobalt")
        command = f"{cobalt_home}cobalt -i {gene}_hits.fasta -outfmt mfasta " \
                  f"-end_gapopen=5 -end_gapextend=1 -gapopen=11 -gapextend=1  -blast_evalue=0.003 -norps=T " \
                  f"-treemethod=clust > aligned.fasta"
        run_bash_command(command)
        print("Done")
        return False
    else:
        return True


def build_hmm(gene, hmmer_home):
    print("Building HMM")
    command = f"{hmmer_home}hmmbuild {gene}.hmm aligned.fasta"
    run_bash_command(command)
    print("Done")


def search_templates(gene, max_pdbs):
    print("Looking for templates")
    command = "curl -L -H 'Expect:' -H 'Accept:text/xml' -F seqdb=pdb -F seq='<" + gene + ".hmm' " \
               "https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch > possible_templates.xml"
    run_bash_command(command)

    "curl -L -H 'Expect:' -H 'Accept:text/xml' -F seqdb=pdb -F seq='<" + gene + ".hmm' https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch > possible_templates.xml"

    #  extract the pdb ids
    with open("possible_templates.xml", "r") as f:
        content = f.readlines()

    check_output = "Sorry" in content
    if not check_output:
        count_pdbs = 0
        templates = open("top_templates.dat", "w")
        with open("all_pdbs.dat", "w") as all_pdbs:
            for line in content:
                if "hits name" in line and not line[16:22].isnumeric():
                    if count_pdbs < max_pdbs:
                        templates.write(str(line[16:22]) + "\n")
                    all_pdbs.write(str(line[16:22]) + "\n")
                    count_pdbs += 1
            templates.close()
        print("Done")
    return check_output

def download_pdb(pdb):
    right_pdbname = f'{pdb}.pdb'
    wrong_pdbname = f"pdb{pdb}.ent"
    pdbl = PDBList()
    filename = pdbl.retrieve_pdb_file(pdb, pdir=".", file_format="pdb", overwrite=True)
    os.rmdir("obsolete")
    if os.path.exists(filename):
        os.rename(wrong_pdbname, right_pdbname)
        return True
    else:
        return None



def get_templates():
    with open("top_templates.dat", "r") as f:
        content = f.readlines()
    pdbs = []
    chains = []
    for pdb, chain in zip(content, content):
        pdbs.append(pdb[:-3])
        chains.append(chain[-2])
    return pdbs, chains



class Select_Atoms(Select):
    def accept_atom(self, atom):
        return 1 if not atom.is_disordered() else 0


def get_structure(pdbname):
    parser_pdb = PDBParser(PERMISSIVE=1)
    structure_id = "pratchett"
    return parser_pdb.get_structure(structure_id, pdbname)


def manage_pdb(pdb, chain):
    structure = get_structure(f'{pdb}.pdb')
    io = PDBIO()
    io.set_structure(structure)
    first_model = structure[0]

    if not structure.header["resolution"]:
        return None

    mychain = Selection.unfold_entities(first_model, 'C')
    residues = Selection.unfold_entities(mychain, 'R')

    unwanted = ["UNK", "MSE", "M3L"]
    non_standard_aa = [
        residue.get_full_id()
        for residue in residues
        if residue.get_id()[0] != " "
        or residue.get_resname() in unwanted
        or residue.get_full_id()[2] != chain
        or residue.get_full_id()[3][1] < 1
    ]

    if non_standard_aa:
        for resid in non_standard_aa:
            structure[resid[1]][resid[2]].detach_child(resid[3])

    io.save(f'{pdb}_{chain}_clean.pdb', Select_Atoms())

    return "OK"

def write_fasta_file(residues, chain, pdb, f):
    acceptable_residues = [
        residue
        for residue in residues
        if is_aa(residue) and residue.get_id()[1] > 0
    ]

    f.write(f">P1;{pdb}" + "\n")
    for counter, res in enumerate(acceptable_residues, start=1):
        f.write(Polypeptide.three_to_one(res.get_resname()))
        if counter % 50 == 0:
            f.write("\n")
    f.write("\n")


def get_fasta(pdb_file, chain_id):
    # gets the PDB
    pdbname = f'{pdb_file}_{chain_id}_clean.pdb'
    structure = get_structure(pdbname)

    residues = Selection.unfold_entities(structure, 'R')
    filename = pdbname[:-4]

    with open(f'{filename}_fasta.dat', "w") as f:
        write_fasta_file(residues, chain_id, filename, f)


def merge_fastas(pdbs, chains, gene):

    for pdb, chain in zip(pdbs, chains):
        command = f"cat {pdb}_{chain}_clean_fasta.dat >>{gene}_hits.fasta"
        run_bash_command(command)


def align_templates_to_msa(gene, cobalt_home):

    print("Doing MSA with cobalt")
    command = f"{cobalt_home}cobalt -i " \
              f"{gene}_hits.fasta -outfmt mfasta -end_gapopen=5 -end_gapextend=1 -gapopen=11 " \
              f"-gapextend=1  -blast_evalue=0.003 -norps=T -treemethod=clust > final_msa.fasta"
    run_bash_command(command)
    print("Done")


def onesequence(block):
    fullchain = []
    for line in range(1, len(block)):
        fullchain.extend(block[line][aa] for aa in range(len(block[line])))
    return fullchain


def sorting_function(beginnings, ends, deltas, pdbs, chains):  # sorts the interval from the shortest(s) to the
    # longest(s) and retains the information of which interval is related to which pdb. I also removes the chains that
    # are in the same pdb and cover the same intervals
    sorted_beginnings = []
    sorted_ends = []
    sorted_deltas = []
    sorted_pdbs = []
    sorted_chains = []
    intervals = len(beginnings)

    for i in range(intervals):
        if i == 0:
            sorted_beginnings.append(beginnings[i])
            sorted_ends.append(ends[i])
            sorted_deltas.append(deltas[i])
            sorted_pdbs.append(pdbs[i])
            sorted_chains.append(chains[i])
        else:
            if deltas[i] >= sorted_deltas[i - 1]:
                sorted_beginnings.append(beginnings[i])
                sorted_ends.append(ends[i])
                sorted_deltas.append(deltas[i])
                sorted_pdbs.append(pdbs[i])
                sorted_chains.append(chains[i])
            if deltas[i] < sorted_deltas[i - 1]:
                counter = 0
                for j in sorted_deltas:
                    if deltas[i] < j:
                        sorted_deltas.insert(counter, deltas[i])
                        sorted_ends.insert(counter, ends[i])
                        sorted_beginnings.insert(counter, beginnings[i])
                        sorted_pdbs.insert(counter, pdbs[i])
                        sorted_chains.insert(counter, chains[i])
                        break
                    counter += 1
    copy_beginnings = sorted_beginnings.copy()
    copy_ends = sorted_ends.copy()
    copy_deltas = sorted_deltas.copy()
    copy_pdbs = sorted_pdbs.copy()
    for i in range(intervals):
        if i < int(len(beginnings)) - 1:
            for j in range(intervals):
                if j > i:
                    if copy_beginnings[i] == sorted_beginnings[j] and copy_ends[i] == sorted_ends[j] \
                            and copy_deltas[i] == sorted_deltas[j] and copy_pdbs[i] == sorted_pdbs[j]:
                        sorted_beginnings[j] = "-"
                        sorted_ends[j] = "-"
                        sorted_deltas[j] = "-"
                        sorted_pdbs[j] = "-"
                        sorted_chains[j] = "-"
    i = 0
    while i < intervals:
        if sorted_beginnings[i] == "-":
            sorted_beginnings.remove(sorted_beginnings[i])
            sorted_ends.remove(sorted_ends[i])
            sorted_deltas.remove(sorted_deltas[i])
            sorted_pdbs.remove(sorted_pdbs[i])
            sorted_chains.remove(sorted_chains[i])
            intervals -= 1
            continue
        i += 1
    return sorted_beginnings, sorted_ends, sorted_deltas, sorted_pdbs, sorted_chains


def filter_function(sorted_beginnings, sorted_ends, sorted_deltas, sorted_pdbs, sorted_chains, slash):

    intervals = range(len(sorted_pdbs))
    reduced_beginnings = []
    reduced_ends = []
    reduced_pdbs = []
    reduced_chains = []
    reduced_deltas = []
    counter = 0

    with open("covered_aa.txt", "w") as o:
        while counter < intervals[-2]:
            testpoints = np.arange(int(sorted_beginnings[counter]), int(sorted_ends[counter]))
            globalintervals = []
            for pdb in intervals:
                if pdb > counter:
                    pdbpoints = np.arange(int(sorted_beginnings[pdb]), int(sorted_ends[pdb]))
                    globalintervals.append(pdbpoints)
            global_intervals = np.unique(np.hstack(globalintervals))
            intersection = np.intersect1d(testpoints, global_intervals)
            testpoints_size = testpoints.size
            if testpoints_size != intersection.size:
                reduced_beginnings.append(sorted_beginnings[counter])
                reduced_ends.append(sorted_ends[counter])
                reduced_deltas.append(sorted_deltas[counter])
                reduced_chains.append(sorted_chains[counter])
                reduced_pdbs.append(sorted_pdbs[counter])
                o.write(f"{sorted_pdbs[counter]}\t{sorted_beginnings[counter]}\t{sorted_ends[counter]}\t{sorted_deltas[counter]}\t{sorted_chains[counter]}\n")
            counter += 1
        reduced_beginnings.append(sorted_beginnings[-1])
        reduced_ends.append(sorted_ends[-1])
        reduced_deltas.append(sorted_deltas[-1])
        reduced_chains.append(sorted_chains[-1])
        reduced_pdbs.append(sorted_pdbs[-1])
        o.write(f"{sorted_pdbs[-1]}\t{sorted_beginnings[-1]}\t{sorted_ends[-1]}\t{sorted_deltas[-1]}\t{sorted_chains[-1]}")
    return reduced_beginnings, reduced_ends, reduced_deltas, reduced_chains, reduced_pdbs


def refine_function(reduced_beginnings, reduced_ends, reduced_deltas, reduced_chains, reduced_pdbs, slash):
    copy_beginnings = reduced_beginnings.copy()
    copy_ends = reduced_ends.copy()
    copy_deltas = reduced_deltas.copy()
    copy_chains = reduced_chains.copy()
    copy_pdbs = reduced_pdbs.copy()
    o = open("testrefined.txt", "w")
    interval = 0
    while interval < len(reduced_pdbs):
        testpoints = np.arange(int(reduced_beginnings[interval]),
                               int(reduced_ends[interval]))  # fixed interval to be tested (i.e. 3-9)
        for i in range(len(copy_pdbs)):  # run all possible intervals
            for j in range(len(copy_pdbs)):  # run all possible intervals
                pair = []
                if i != interval and j != interval and i != j and i < j:  # select with i only the intervals smaller than the one being tested and with j only the bigger
                    lower = np.arange(reduced_beginnings[i], reduced_ends[
                        i])  # get all the number between smaller begin and smaller end (i.e. 2-3)
                    higher = np.arange(reduced_beginnings[j], reduced_ends[
                        j])  # get all the number between bigger begin and bigger end (i.e 5-12)
                    pair.append(lower)
                    pair.append(higher)
                    points = np.unique(np.hstack(
                        pair))  # get all the unique points between higher and lower (but not missing ones i.e. 2,3,5,6...)
                    intersection = np.intersect1d(points,
                                                  testpoints)  # get the points that are both in our test and the other 2

                    # if the intersection is equal to the other 2 it means it can be represented by a smaller (not previously eliminated interval) and a bigger one which mean that the smaller one covers a part of the protein not covered from our tested one, so we can eliminate the tested interval.

                    if (testpoints.size == intersection.size):
                        z = 0
                        while z < len(copy_pdbs):
                            vector_to_pop = np.arange(int(reduced_beginnings[z]), int(reduced_ends[z]))
                            if testpoints.size == vector_to_pop.size:
                                copy_beginnings.pop(z)
                                copy_ends.pop(z)
                                copy_deltas.pop(z)
                                copy_chains.pop(z)
                                copy_pdbs.pop(z)
                                break
                            z += 1
        interval += 1
    for i in range(len(copy_pdbs)):
        o.write(str(copy_pdbs[i]) + "\t" + str(copy_beginnings[i]) + "\t" +
                str(copy_ends[i]) + "\t" + str(copy_deltas[i]) + "\t" + str(
            copy_chains[i]) + "\n")


def run_templates_analysis(res_cutoff, seqid_cut, slash):
    with open("aligned_oneline.fasta", "r") as f:
        content = f.read()

    newblock = []
    for block in content.split("\n\n"):
        newblock.append(block.splitlines())
    o = open("covered_intervals.txt", "w")


    for block in newblock:
        score = 0  # the number of aligned residues (when there's a mutation it doesn't consider if it's conservative or semiconservative)
        seqid = 0  # % of the number of identical residues / the total amount of residues in the template sequence
        counter_templ = 0
        if block != newblock[-1]:  # block = template sequence. newblock[-1] = target
            # if "clean" in block[0]:
            pdbname = block[0][-12:-8] + "_chain_" + block[0][-7]

            # resolution = 100

            structure = get_structure(pdbname[0:4] + ".pdb")
            io = PDBIO()
            io.set_structure(structure)
            resolution = structure.header["resolution"]

            # if resolution > res_cutoff:
                # pdbfile.close()

            if resolution <= res_cutoff:
                o.write(pdbname + "\n")
                line = 1
                begins = []
                ends = []
                fullchain = onesequence(block)
                counter_aa = 0
                while line < len(block):
                    aa = 0
                    while aa < len(block[line]):
                        if block[line][aa] != "-":
                            score += 1  # checks if there is an aligned AA
                            counter_templ += 1
                        if block[line][aa] == newblock[-1][line][aa] and block[line][aa] != "-":
                            seqid += 1  # checks if the AA is identical to the one in the reference seq

                        if counter_aa == 0 and fullchain[counter_aa] != "-":
                            begins.append(counter_aa + 1)
                        if counter_aa > 0 and counter_aa < len(fullchain):
                            if fullchain[counter_aa - 1] == "-" and fullchain[counter_aa] != "-":
                                begins.append(counter_aa + 1)
                            if fullchain[counter_aa - 1] != "-" and fullchain[counter_aa] == "-":
                                ends.append(counter_aa)
                            if fullchain[counter_aa] != "-" and aa == len(block[line]) - 1:
                                ends.append(counter_aa)
                        counter_aa += 1
                        aa += 1
                    line += 1
                percentage = "{:.2f}".format((score / counter_aa * 100))
                seqid = "{:.2f}".format((seqid / counter_templ * 100))
                o.write("Score = " + str(score) + "\n")
                o.write("Coverage = " + str(percentage) + "\n")  # % of the reference sequence covered by the template
                o.write("Sequence identity = " + str(seqid) + "\n")
                for i in range(len(begins)):
                    o.write(str(begins[i]) + "-" + str(ends[i]) + "\n")
                o.write("\n")
    o.close()
    f.close()
    counter = 0
    pdb = ""
    ch = None
    data = {}  # creates the empty dictionary
    # del line
    del f
    del o

    with open("covered_intervals.txt", "r") as f:
        line = f.readline()
        while line:
            line = line.strip()
            if line.strip() == "":
                counter = 0
                x = data.get(pdb, [])  # creates inside the dictionary the key (pdb) and the value (ch). ch is stored as a list and will become the chain
                x.append(ch)
                data[pdb] = x
            elif counter == 0:
                x = line.split("_", 1)
                pdb = x[0]
                chainname = x[1]
                ch = Chain(chainname)
                counter = counter + 1
            elif counter == 1:
                ch.add_score(float(line.split("=")[1]))
                counter = counter + 1
            elif counter == 2:
                ch.add_coverage((float(line.split("=")[1])))
                counter += 1
            elif counter == 3:
                ch.add_seqid(float(line.split("=")[1]))
                counter += 1
            else:
                ch.add_range(range1(line))
            line = f.readline()

    if counter != 0:
        x = data.get(pdb, [])
        x.append(ch)
        data[pdb] = x
    beginnings = []
    ends = []
    deltas = []
    pdbs = []
    allchain = []
    for pdb in data:
        for chain in range(len(data[pdb])):
            for interval in range(len(data[pdb][chain].ranges)):
                if data[pdb][chain].seqid > seqid_cut:
                    beginnings.append(int(data[pdb][chain].ranges[interval].start))
                    ends.append(int(data[pdb][chain].ranges[interval].end))
                    deltas.append(
                        int(data[pdb][chain].ranges[interval].end) - int(data[pdb][chain].ranges[interval].start) + 1)
                    pdbs.append(str(pdb))
                    allchain.append(str(data[pdb][chain].chainid))

    if not pdbs:
        print("There are no suitable PDBS")

        used_beginnings, used_ends, used_chains, used_pdbs = 42, 42, 42, 42

    else:
        sorted_beginnings, sorted_ends, sorted_deltas, sorted_pdbs, sorted_chains = sorting_function(beginnings, ends,
                                                                                                     deltas, pdbs,
                                                                                                     allchain)
        output = open("test_sorting.txt", "w")
        for i in range(len(sorted_deltas)):
            output.write(str(sorted_pdbs[i]) + "\t" + str(sorted_beginnings[i]) + "\t" + str(sorted_ends[i]) + "\t"
                         + str(int(sorted_deltas[i])) + "\t" + str(sorted_chains[i]) + "\n")

        reduced_beginnings, reduced_ends, reduced_deltas, reduced_chains, reduced_pdbs = filter_function(
            sorted_beginnings, sorted_ends, sorted_deltas,
            sorted_pdbs, sorted_chains, slash)

        refine_function(reduced_beginnings, reduced_ends, reduced_deltas, reduced_chains, reduced_pdbs, slash)
        used_beginnings, used_ends, used_chains, used_pdbs = 8, 8, 8, 8

    return used_beginnings, used_ends, used_chains, used_pdbs


def chain_list(dct_name, chains):
    # This functions takes all the chain IDs and after playing a bit with the text returns
    # a dictionary containing some lists called with the correct chain ID (A, B...)
    # The text is not default for biopython, that's why it's necessary to modify it.

    dct_name = {}
    # count the position used to extract the chain_ID
    position = 0
    k = 0
    for i in chains:
        test = str(chains[position])
        index1 = test[test.find('='):]
        index2 = index1[:index1.find('>')]
        index3 = index2[index2.find('Waldo')]
        if position == 0:
            chain_ID = [index3]
        if position != 0:
            chain_ID.append(index3)
        position = position + 1
        dct_name['chain_%s' % chain_ID[k]] = []
        k = k + 1
    return dct_name, chain_ID


def write_numbers(numb, f):
    # This function writes the number above the AAs ( 1     10     20...) and the spaces every 10 AAs
    # This is the function you want to change if you want to increase the number of AA for each chain (atm it's <1000)
    # and still have a nicely formatted output

    if numb == 0:
        f.write(str(numb + 1) + "       " + str(numb + 10) + "         " + str(numb + 20) + "         " + str(
            numb + 30) + "         " + str(numb + 40) + "         " + str(numb + 50))
        f.write("\n")
    if numb > 0 and numb < 100 and numb % 50 == 0:
        f.write(
            "\n" + "        " + str(numb + 10) + "         " + str(numb + 20) + "         " + str(
                numb + 30) + "         " + str(numb + 40) + "        " + str(numb + 50))
        f.write("\n")
    if numb >= 100 and numb % 50 == 0:
        f.write("\n" + "       " + str(numb + 10) + "        " + str(numb + 20) + "        " + str(
            numb + 30) + "        " + str(numb + 40) + "        " + str(numb + 50))
        f.write("\n")
    if numb > 0 and numb % 10 == 0 and not numb % 50 == 0:
        f.write(" ")


def pdb_not_complete(first_res, numb, f):
    # If a PDB doesn't start with the AA number 1 this will fill the gap

    while numb < first_res - 1:
        write_numbers(numb, f)
        f.write("-")
        numb = numb + 1
    return numb


def fill_holes(shadow_count, res_number, numb, f):
    # If a PDB has missing residues this will fill the gaps

    if shadow_count < res_number:
        while shadow_count < res_number - 1:
            write_numbers(numb, f)
            f.write("-")
            numb = numb + 1
            shadow_count = shadow_count + 1
    return shadow_count, numb


def input_modeller(numb, pdbname, first_res, last_res, count, chain_ID, dct_chains, count_res, shadow_count, f, keys,
                   residues, AA_dict):
    # This function straight up gives you the input that can be used with modeller. It contains a slightly different version
    # of the function "fill the holes". This version doesn't use the "write_numbers" function since is not necessary.

    if numb == 0:
        f.write(">P1;" + pdbname + "\n")
        f.write(
            "structureX:" + pdbname + ":" + str(first_res) + ":" + chain_ID[count] + ":" + str(last_res) + ":" +
            chain_ID[count] + ":: : :\n")
    if numb % 50 == 0 and numb > 0:
        f.write("\n")
    res_number = int(dct_chains[keys[count]][count_res])
    resname = residues.get_resname()
    if shadow_count < res_number:
        while shadow_count < res_number - 1:
            if numb % 50 == 0:
                f.write("\n")
            f.write("-")
            numb = numb + 1
            shadow_count = shadow_count + 1
    safe = 0
    try:
        AA_dict[resname]
    except KeyError:
        f.write("-")
        safe = 1
    if safe == 0 and res_number >= 1:
        f.write(AA_dict[resname])
    numb = numb + 1
    shadow_count = shadow_count + 1
    count_res = count_res + 1
    return numb, shadow_count, count_res


# class NonHetSelect(Select):  # use to remove residues like MSE, anisotropies and others. Doesn't remove DNA
#     def accept_residue(self, residue):
#         return 1 if residue.id[0] == " " else 0


def prepare_modeller(sequence_name, gene):

    refinedfile = str(sequence_name) + "/testrefined.txt"
    filein = open(refinedfile, "r")
    lines = filein.readlines()
    used_beginnings = []
    used_ends = []
    used_chains = []
    used_pdbs = []
    for x in lines:
        used_beginnings.append(int(x.split('\t')[1]))
        used_ends.append(int(x.split('\t')[2]))
        used_chains.append(str(x.split('\t')[4]).strip())
        used_pdbs.append(str(x.split('\t')[0]))

    # prepares the modeller files with the aligned sequences
    # sequence_file = str(sequence_name) + "/aligned_oneline_withgaps.fasta"
    # directory = "modeller_" + str(sequence_name)
    # os.mkdir(directory)
    # os.mkdir(directory + "/fastas")
    # gene = str(gene.lower())
    # first_aa = min(used_beginnings)
    # last_aa = max(used_ends)
    pdbs_and_chains = {}
    count = 0
    dct = {}
    while count < len(used_pdbs):
        dct[used_pdbs[count]] = []
        count += 1
    count = 0
    while count < len(used_pdbs):
        dct[used_pdbs[count]].append(used_chains[count])
        count += 1
    for pdb in dct:
        pdbs_and_chains[pdb] = []
        appo_dct = {}
        for chain in dct[pdb]:
            appo_dct[chain] = []
        for i in appo_dct:
            pdbs_and_chains[pdb].append(i)

    return pdbs_and_chains


def prepare_alignments(pdbs_and_chains, winning_iso, gene, cutoff_slash):

    mydir = winning_iso
    os.chdir(mydir)
    with open("alignment_pdbs.fasta", "r") as f:
        content = f.read()


    with open("for_modeller.fasta", "w") as f:
        for block in content.split("\n\n"):
            for pdb in pdbs_and_chains:
                if pdb in block:
                    towrite = block.split("\n")
                    f.write(towrite[0] + "\n")
                    removed = towrite.pop(0)
                    towrite = ''.join(towrite)
                    f.write(towrite + "\n\n")
            if "1 sp" in block:
                towrite = block.split("\n")
                f.write(towrite[0] + "\n")
                removed = towrite.pop(0)
                towrite = ''.join(towrite)
                f.write(towrite)

    with open("for_modeller.fasta", "r") as f:
        content = f.read()

    sequences = []
    headers = []
    for block in content.split("\n\n"):
        for pdb in pdbs_and_chains:
            if pdb in block:
                headers.append(block.split("\n")[0])
                sequences.append(block.split("\n")[1])
        if block == content.split("\n\n")[-1]:
            target_head = block.split("\n")[0]
            target_sequence = block.split("\n")[1]
    positions = []
    counter_aa = 0

    aa_positions = []
    while counter_aa < len(target_sequence):
        for sequence in sequences:
            covered = 0
            if sequence[counter_aa] != "-":
                covered = 1
                aa_positions.append(counter_aa)
                break
        if covered == 0 and target_sequence[counter_aa] != "-":
            positions.append(counter_aa)
        counter_aa += 1
    min_aa = min(aa_positions)
    max_aa = max(aa_positions)
    to_remove = []
    for pos in positions:
        if pos < min_aa or pos > max_aa:
            to_remove.append(pos)
    for i in to_remove:
        positions.remove(i)

    appo_array = np.split(np.array(positions), np.where(np.diff(np.array(positions)) != 1)[0]+1)
    continuos_intervals = []
    for interval in appo_array:
        if len(interval) > cutoff_slash:
            continuos_intervals.append(list(interval))
    final_target_sequence = [target_head]
    counter_aa = 0
    aa_list = ""
    while counter_aa < len(target_sequence):
        if counter_aa < min_aa or counter_aa > max_aa:
            aa_list += "-"
        elif any(counter_aa in interval for interval in continuos_intervals):
            for interval in continuos_intervals:
                if counter_aa == interval[0]:
                    aa_list += "/"
            aa_list += "-"
        else:
            aa_list += target_sequence[counter_aa]
        counter_aa += 1
    final_target_sequence.append(aa_list)

    pdbs = []
    with open("input_modeller.dat", "w") as f:
        count_seq = 0
        while count_seq < len(sequences):

            pdbs.append(headers[count_seq][-12:])
            structure = get_structure(headers[count_seq][-12:] + ".pdb")
            res_list = Selection.unfold_entities(structure, "R")
            first = str(res_list[0].get_id()[1])
            last = str(res_list[-1].get_id()[1])

            sequence = list(sequences[count_seq])
            counter = 0
            for interval in continuos_intervals:
                first_element = interval[0]
                sequence.insert(first_element + counter, "/")
                counter += 1
            sequence = ''.join(sequence)

            f.write(">" + headers[count_seq][-15:] + "\n")
            f.write("structureX:" + headers[count_seq][-12:] + ":" + first + ":" + headers[count_seq][-7] + ":" + last +
                    ":" + headers[count_seq][-7] + ":: : :\n")
            f.write(str(sequence) + "*\n\n")
            count_seq +=1
        f.write(">P1;" + gene.upper() + "\n")
        f.write("sequence:" + gene.upper() + ":xxx:x:xxx:x::::\n")
        f.write(final_target_sequence[1] + "*")


    return pdbs


def write_modeller(pdbs, gene, number_of_models):
    # local_slave may need an update to local_worker
    directory = os.getcwd()
    # newdir = directory + slash + "modeller_" + str(isoform) + slash
    # myfiles = [f for f in os.listdir(newdir) if os.path.isfile(os.path.join(newdir, f))]
    # for i in myfiles:
    #     if "reduced.pdb" in i:
    #         shutil.copyfile(newdir + slash + i, newdir + slash + i[:4])
    o = open("run_modeller.py", "w")
    o.write(
        "from modeller import *\nfrom modeller.automodel import *\nfrom modeller.scripts import complete_pdb\nfrom modeller.parallel import *\nimport os\n\n")
    o.write("log.verbose()\nj=job()\ncores_number = 8\nthread = 0\n")
    o.write("while thread < cores_number:\n    j.append(local_slave())\n    thread += 1\n")
    o.write("env = environ()\nenv.io.atom_files_directory = ['" + directory + "']\nenv.io.hetatm = True\n")
    o.write("a = automodel(env,\n    alnfile =\"input_modeller.dat\",\n    knowns = (" + str(
        pdbs) + "),\n    sequence = \"" + str(gene.upper()) + "\",\n    assess_methods=(assess.DOPE, assess.GA341))\n")
    o.write("a.starting_model= 1\na.ending_model  = " + number_of_models + "\na.use_parallel_job(j)\na.make()\n")
    o.write(
        "ok_models = filter(lambda x: x['failure'] is None, a.outputs)\ntoscore = 'DOPE score'\nok_models = sorted(ok_models, key=lambda k: k[toscore])\nm = ok_models[0]\nmyout = open(\"MYOUT.dat\", \"w\")")
    o.write("\nmyout.write(\"Top model: \" + str(m['name']) + \" (DOPE SCORE: %.3f)\" % (m[toscore]))\n")
    o.write(
        "env.libs.topology.read(file='$(LIB)/top_heav.lib')\nenv.libs.parameters.read(file='$(LIB)/par.lib')\nmdl = complete_pdb(env, m['name'])\ns = selection(mdl)\n")
    o.write("s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=\"" + str(
        gene.upper()) + "_profile.dat\", normalize_profile=True, smoothing_window=15)")
    o.close()


def run_modeller(gene, slash):

    # uses a bash command to run modeller.
    cmd = "python run_modeller.py > modeller_output.txt"
    run_bash_command(cmd)

    filename = "MYOUT.dat"
    waiter = 0
    while waiter == 0:
        doesit = str(path.exists(filename))
        if doesit == "True" and os.stat(filename).st_size != 0:
            waiter = 1
            time.sleep(10)
        else:
            time.sleep(10)


def oneliner(seq_name):
    lines = []
    heading = []
    with open(seq_name, "r") as f:
        content = f.readlines()
    f.close()
    count = 0
    for line in content:
        if count < 2:
            heading.append(line.strip())
        else:
            line.strip("*")
            lines.append(line.strip())
        count += 1
    line = "".join(lines)
    return heading, line


def make_dir(dirname):

    try:
        os.mkdir(dirname)
    except FileExistsError:
        pass


def clean_templates_search(gene, line):

        replace_files("possible_templates.xml", "template_search/possible_templates.xml")
        replace_files(gene + "_hits.fasta", "template_search/" + line.strip()[:-6] + "_hits.fasta")
        replace_files(gene + ".hmm", "template_search/" + gene + ".hmm")
        replace_files("aligned.fasta", "template_search/homologous_aligned.fasta")
        replace_files("top_templates.dat", "template_search/top_templates.dat")
        replace_files("all_pdbs.dat", "template_search/all_possible_templates.dat")

def move_pdbs(my_files, pdbs):

    for f in my_files:
        if ".pdb" in f and len(f) == 8:
            replace_files(f, f"pdbs/{f}")
        if "_clean.pdb" in f and str(f[:-4]) not in pdbs:
            replace_files(f, f"pdbs/{f}")


def remove_files(filename):
    try:
        os.remove(filename)
    except FileNotFoundError:
        pass


def replace_files(filein, fileout):
    try:
        os.replace(filein, fileout)
    except FileNotFoundError:
        pass


def clean_analysis(my_files):

    file_to_remove = ["temp_align.fasta", "cleaned_align.fasta", "aligned_oneline.fasta", "test_sorting.txt",
                   "testrefined.txt", "covered_aa.txt"]
    for f in file_to_remove:
        remove_files(f)

    replace_files("covered_intervals.txt", "template_analysis/covered_intervals")

    for f in my_files:
        if "_clean_fasta.dat" in f:
            replace_files(f, "template_analysis/" + f[:-10] + ".fasta")


def cleaning(gene, pdbs):

    with open("master_isoform.txt", "r") as f:
        content = f.readlines()
    for line in content:
        os.chdir(line.strip()[:-6])
        my_files = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]

        dirs = ["template_search", "pdbs", "template_analysis", "modeller_files_wt"]
        if "run_modeller.py" in my_files:
            for mydir in dirs:
                make_dir(mydir)
            replace_files("MYOUT.dat", "modeller_files_wt/MYOUT.dat")
            for filename in my_files:
                if "slave" in filename:
                    remove_files(filename)
                elif "worker" in filename:
                    remove_files(filename)
                if gene in filename and ".pdb" in filename:
                    replace_files(filename, "modeller_files_wt/" + filename)
                if gene in filename and len(filename) == len(gene) + 4 and ".hmm" not in filename and "pdb" not in filename:
                    remove_files(filename)
                if gene in filename and len(filename) == len(gene) + 10:
                    remove_files(filename)
                if "_profile" in filename:
                    remove_files(filename)



        else:
            make_dir(dirs[0])
            make_dir(dirs[1])
            make_dir(dirs[2])

        move_pdbs(my_files, pdbs)
        clean_templates_search(gene, line)
        clean_analysis(my_files)
        os.chdir("..")

def clean_modeller():
    os.mkdir("models")
    with open("MYOUT.dat", "r") as f:
        content = f.readlines()
    mymodel = content[0][11:content[0].find("(")-1]
    for model in os.listdir():
        if mymodel not in model and "MYOUT.dat" not in model and "models" not in model:
            replace_files(model, f"models/{model}")

def clean_mutant(gene):

    for filename in os.listdir():
        if "slave" in filename:
            remove_files(filename)
        elif "worker" in filename:
            remove_files(filename)
        elif ".ini" in filename or ".rsr" in filename or "sch" in filename:
            remove_files(filename)
        elif gene in filename and (".V" or ".D") in filename:
            remove_files(filename)
        elif gene and ".V" in filename or gene and ".D" in filename:
            remove_files(filename)
    clean_modeller()


def final_cleaning(gene, master_directory):
    os.chdir(gene)
    allfiles = os.listdir()
    mydirs = [thefile for thefile in allfiles if "isoform" in thefile and "." not in thefile]
    for dir in mydirs:
        if "modeller_files_wt" in os.listdir(dir):
            os.chdir(f"{dir}/modeller_files_wt")
            clean_modeller()
            os.chdir("..")

            for mutant in os.listdir():
                if "mutation_" in mutant:
                    os.chdir(mutant)
                    print(os.getcwd())
                    clean_mutant(gene)
                    os.chdir("..")
            os.chdir("..")
    os.chdir(master_directory)


def problematic_genes(all_genes):

    command = "grep Sorry */isoform*/template_search/possible_templates.xml"
    result = subprocess.run(command, shell=True, universal_newlines=True, check=True, capture_output=True)
    partial_results = result.stdout.split("\n")
    partial_results.pop()
    # genes and isoforms that went wrong
    problematic_genes = []
    problematic_isoforms = []
    for i in partial_results:
        problematic_genes.append(i.split("/")[0])
        problematic_isoforms.append(i.split("/")[1])

    no_models = []
    command = "ls isoform*/ | grep modeller_files_wt"
    for gene in all_genes:
        os.chdir(gene)
        result = subprocess.run(command, shell=True, universal_newlines=True, capture_output=True)
        if len(result.stdout) == 0:
            no_models.append(gene)
        os.chdir("..")

    with open("Errors.txt", "w") as f:
        f.write("For the following genes and isoforms, hmmer didn't return any template, you may want to rerurn these genes:\n")
        for gene, isoform in zip(problematic_genes, problematic_isoforms):
            f.write(f"\n{gene}\t\t{isoform}")
        f.write("\n\nFor the following genes, no models were produced, you may want to take a look and rerun or do it by hand:\n")
        for unmodeled in no_models:
            f.write(f"\n{unmodeled}")
def final_report():

    command = "touch final_report.log"
    run_bash_command(command)

    myfolders = next(os.walk('.'))[1]
    for folder in myfolders:
        myfiles = os.listdir(folder)
        for f in myfiles:
            if ".log" in f:
                command = f"cat {folder}/{f} >> final_report.log"
                run_bash_command(command)