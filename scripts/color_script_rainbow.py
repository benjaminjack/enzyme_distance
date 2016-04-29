import csv
import gzip
# from spectrumany import *

def make_structure(path, pdb, pdb_name):

    cmd.load(path + pdb + '.pdb')
    cmd.hide('all')

    cor_min = 0
    cor_max = 2.2

    resi_list = []

    cmd.alter("%s and n. CA"%pdb, "b=0.0")
    cmd.color("gray50", pdb)
    cmd.show("cartoon")
    cmd.cartoon("tube")
    cmd.set('cartoon_tube_radius', '0.7')
    cmd.set("cartoon_transparency", "0.5")

    with gzip.open('master_data_table.csv.gz') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['pdb'] == pdb_name:
                resi = row['Residue']
                gam_rate = row['gam_rate']
                cmd.alter("%s and resi %s and n. CA"%(pdb,resi), "b=%s"%gam_rate)
                if int(row['ACTIVE_SITE']) == 1:
                    cmd.show("spheres", "res %s"%resi)
                    cmd.color("white", "res %s"%resi)
                resi_list.append(int(resi))

    cmd.select("to_color", "res "+str(min(resi_list))+":"+str(max(resi_list)))
    cmd.create("to_color2","to_color")
    cmd.remove(pdb + " in to_color2")
    cmd.spectrum("b", "rainbow", "to_color2 and n. CA", minimum=cor_min, maximum=cor_max)
    cmd.set("cartoon_transparency", "0.0", "to_color2")
    cmd.select("all")
    cmd.set('opaque_background', 'off')
    cmd.select("none")


    # clear out the old B Factors
    # cmd.alter("%s and n. CA"%pdb, "b=0.0")

    # update the B Factors with new properties
    # cmd.alter("%s and n. CA"%pdb, "b=stored.newB.pop(0)")

    # color the protein based on the new B Factors of the alpha carbons
    '''
    cmd.spectrum("b", "rainbow", "%s and n. CA"%pdb, minimum=cor_min, maximum=cor_max)
    cmd.show("cartoon")
    cmd.cartoon("tube")
    cmd.set('cartoon_tube_radius', '0.7')
    cmd.set('opaque_background', 'off')
    for pos in active_sites:
        cmd.show("spheres", "resi " + str(pos))
    '''
    # cmd.save('pymol_sessions_pcor/'+protein+'.pse')
    #cmd.ray("1200", "1200")
    #cmd.png("%s.png"%pdb)
    # cmd.delete(pdb)

    # print(cor_max, cor_min)

cmd.extend("makeColors", make_structure)
