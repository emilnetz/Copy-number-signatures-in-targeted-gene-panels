import numpy as np
import pandas as pd

# assign meta sheet with CN sigs and treatment information
comb = pd.read_excel('/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/MH_clin2_wMeds_filtered_wAge4_folfadjust20260212.xlsx')


# everything after the 10 most comon will be OTHER 

entcounts = comb['Entität'].value_counts()
entcounts = entcounts[0:10]

for i in range(len(comb)):
    if comb.loc[i, 'Entität'] not in entcounts.index:
        comb.loc[i, 'Entität'] = 'Other'

ents = comb['Entität'].unique() 




# threshhold for combining treatments
threshhold = 80

for ent in ents: 
    comb2 = comb[comb['Entität'] == ent]

    #define med variable 
    pos1 = comb.columns.get_loc('Entität')
    pos2 = comb.columns.get_loc('age')
    meds = comb.columns[pos1+1:pos2]

    # drop meds if less than 10 patients were treated with this med
    for med in meds: 
        medcount = sum(comb2[med])
        if medcount < 10: 
            comb2 = comb2.drop(columns=[med])

    # define new meds variable with meds that were kept 
    pos1 = comb2.columns.get_loc('Entität')
    pos2 = comb2.columns.get_loc('age')
    meds = comb2.columns[pos1+1:pos2]

    # create cotreatment df with meds for cols and index 
    df = pd.DataFrame(columns = meds, index= meds)
    
    #reset index of comb2 
    co_comb = comb2.copy()
    co_comb.reset_index(drop=True, inplace=True)


    for med in meds: 
        if sum(comb2[med]) > 0: 
            comb3 = comb2[comb2[med] == 1] 
            for med2 in meds: 
                if sum(comb3[med2]) > 0: 
                    perc = (sum(comb3[med2]) / sum(comb3[med])) * 100
                    perc = f'{perc:.2f}'
                    df.loc[med, med2] = perc
                else: 
                    df.loc[med, med2] = np.nan

    df = df.dropna(axis=0, how='all')
    df = df.dropna(axis=1, how='all')
    df = df.fillna(0)



    # empty list to store meds that will be combined 
    new_names = []
    # loop over all rows of df (meds of interest)
    for med in df.index:
        thresh_meds = []
        # loop over all columns of df 
        for med2 in df.columns: 
            # proceed if med and med2 are not the same
            if med2 != med: 
                # if med2 is given in over 80% of cases where med is given, add to thresh_meds
                if (float(df.loc[med, med2]) > threshhold): 
                    thresh_meds.append(med2)

        # proceed if there is at least one med2 that reached the threshhold
        if len(thresh_meds) > 0: 
            # create a list with all med and all med2 that reached threshhold and add to new_names
            new_name = [med]
            for i in range(len(thresh_meds)):
                name = thresh_meds[i]
                new_name.append(name)
            new_names.append(new_name)



    # remove doubles where only order of meds is different (i.e. [Carboplatin, Paclitaxel], [Paclitaxel, Carboplatin])
    new_names2 = new_names.copy()
    for names in new_names2:
        for names2 in new_names2: 
            if names != names2: 
                if (set(names) == set(names2)): 
                    new_names2.remove(names2)
            


    # loop over list of lists with meds that should be combined
    for names in new_names2:
        co_var = ''
        # loop over each med in sublist 
        for i in range(len(names)):
            var = names[i]
            # create the combined med variable name with underscores 
            co_var = f'{co_var}_{var}'
        # remove leading underscore 
        co_var = co_var.lstrip('_')


        # add combined variable to meta sheet
        co_comb[co_var] = 0
        # loop over every patient (row)
        for x in range(len(co_comb)):
            vals = []
            # loop over every single med that will get combined 
            for i in range(len(names)):
                if co_comb[names[i]][x] == 1: 
                    vals.append(1)
                else: 
                    vals.append(0)
            # change value of combined variable to 1 only if every single med is 1, else: value stays 0
            if 0 not in vals: 
                co_comb[co_var][x] = 1
                
    #co_comb.to_excel('/Users/emilnetz/Desktop/testdf.xlsx')



    # remove single meds that are now in combined variable from meta sheet, they will not be tested seperatly
    for names in new_names: 
        co_comb = co_comb.drop(columns = [names[0]])



    co_comb.to_excel(f'/Users/emilnetz/Desktop/cnv_sigs/MH/MH_therapy_outcome_pred/cotreatment/cotreatment_adjfolf_20260220/MH_{ent}_clin_wCotreatmentGroups.xlsx', index=False)    