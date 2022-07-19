import pandas as pd
import pyteomics
from pyteomics import mzid
from collections import defaultdict
import uuid
def mod_name(modification):
    '''Function that returns list of modifications'''
    mod_types = []
    for i in modification:
        mod_name = i['name']
        mod_types.append(mod_name)
    return mod_types
def location(modification):
    '''Function that returns positions of the modifications'''
    locations = []
    for i in modification:
        locs = i['location']
        locations.append(locs)
    return locations

def massdelta(modification):
    massdelta = []
    for i in modification:
        mass = i['monoisotopicMassDelta']
        massdelta.append(mass)
    return massdelta


def loc_mod(series):
    """ Function for returning the location:modification combined value from Series"""

    mod = series[0]  # modification type column
    loc = series[1]  # positions column
    new = {}
    for i in range(len(mod)):
        new[loc[i]] = mod[i]
    return new

def mod_loc(loc_mod):
    '''Returns dictionary of modifications:locations'''
    mod_pos = dict()
    for key, value in loc_mod.items():
        mod_pos.setdefault(value, list()).append(key)
    return mod_pos

def mod_seq(series):

    pep_seq = series[0]  # peptidesequence
    mod_loc = series[1]  # modification_location
    mod_name = series[2]  # modification_name

    loc_mod = dict(zip(mod_loc, mod_name))
    mod_loc.sort(reverse=True)

    for i in range(len(mod_loc)):
        pos = int(mod_loc[i])
        new_seq = pep_seq[:pos] + "(" + loc_mod[pos] + ")" + pep_seq[pos:]
        pep_seq = new_seq
    return pep_seq

def relational_data(series):
    a = series[0]
    rel_data = ""
    for i in range(len(a)):
        value = "<null>"
        rel_data += value
    return rel_data
def intensity(series):
    a = series[0]
    int_value = ""
    for i in range(len(a)):
        value = None

    return value

def file_type(series):
    '''Default return mzIdentML'''
    a = series[0]
    f_type = ""
    for i in range(len(a)):
        value = "mzIdentML"
        f_type += value
    return f_type
def relational(series):
    a = series[0]
    rel_data = ""
    for i in range(len(a)):
        value = "0"
        rel_data += value
    return rel_data
def protein(row):
    '''Returns Protein names'''
    acc_nos = []
    for i in row:
        i = i.split("|")[1]
        acc_nos.append(i)
    return acc_nos

def modification_number(series):
    '''Returns a list with the total number of modifications'''
    l = series[0]
    total_number = 0
    for i in range(len(l)):
        total_number += 1
    return total_number

def mass_errors(series):
    '''Returns difference in mass values'''
    ex_mass = series[0]
    cal_mass = series[1]
    count = series[2]
    mass_errors = {}
    for i in range(len(count)):
        diff = round(abs(ex_mass - cal_mass), 5)
        mass_errors["mass-error[ppm]"] = diff
    return mass_errors

def modified_sequence(row):
    '''Returns the modified sequence, the first three letters of the modification are represented at the modified positions'''
    seq = row['peptidesequence']
    mods = row['modifications']
    modlocs = row['modification_locations']
    seqd = {-1: ''}
    for i in range(len(seq)):
        seqd[i] = seq[i]

    out = ""
    for i,j in zip(mods, modlocs):
        if j == 0: continue
        seqd[j-1] = f"{seqd[j-1]}({i[:3]})"
    out = "".join([seqd[i] for i in range(len(seq))])
    #print(seqd, out)
    return out


def mzdf(mzid_file):
    read_mzid = mzid.read(mzid_file)
    mzid_df = pyteomics.mzid.DataFrame(read_mzid)
    mzid_df.columns = [col.lower() for col in mzid_df]
    mzid_df = mzid_df[
        ['peptidesequence', 'accession','protein description', 'name', 'modification', 'experimentalmasstocharge',
        'calculatedmasstocharge', 'chargestate']]
    # print(mzid_df)
    mzid_df.rename(columns={'experimentalmasstocharge':'m/z'}, inplace=True)
    mzid_df.rename(columns={'chargestate': 'charge'}, inplace=True)
    mzid_df['modifications'] = mzid_df.modification.apply(mod_name)
    # mzid_df['residues'] = mzid_df.modification.apply(residues)
    mzid_df['modification_locations'] = mzid_df.modification.apply(location)
    mzid_df['mass'] = mzid_df.modification.apply(massdelta)
    mzid_df['location:modification'] = mzid_df[['modifications', 'modification_locations']].apply(loc_mod, axis=1)
    mzid_df['modification:location'] = mzid_df['location:modification'].apply(mod_loc)
    #mzid_df['Modified sequence'] = mzid_df[['peptidesequence', 'modification']].apply(modified_sequence, axis=1)
    mzid_df['Modified sequence'] = mzid_df.apply(modified_sequence, axis=1)
    mzid_df['uniprot'] = mzid_df.accession.apply(protein)
    mzid_df['Modification number'] = mzid_df[['modification_locations']].apply(modification_number, axis=1)
    mzid_df['mass errors'] = mzid_df[
        ['m/z', 'calculatedmasstocharge', 'uniprot']].apply(mass_errors, axis=1)
    mzid_df['Relational-data'] = mzid_df[['mass errors']].apply(relational_data, axis=1)
    mzid_df['Relational'] = mzid_df[['mass errors']].apply(relational, axis=1)
    mzid_df['intensities'] = mzid_df[['mass errors']].apply(intensity, axis=1)
    mzid_df['File-type'] = mzid_df[['mass errors']].apply(file_type, axis=1)
    mzid_df.rename(columns={'peptidesequence':'Sequence'}, inplace = True )
    mzid_df.drop(['accession'], axis=1, inplace = True)
    mzid_df.drop(['calculatedmasstocharge'], axis = 1, inplace = True)
    mzid_df.drop(['modification'], axis=1, inplace=True)
    return mzid_df

def processor(name):
    mzid_file = name
    mzid_df = mzdf(mzid_file)
    mzid_df.insert(0, 'uuid', '')
    mzid_df['uuid'] = mzid_df.apply(lambda _: uuid.uuid4(), axis=1)
    #print(mzid_df.iloc[38])
    #mzid_df.to_csv("check_output.csv") # PXD_filename.csv
if __name__ == '__main__':
    main()

if __name__ == '__main__':
    main()
