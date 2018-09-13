import re

class countsColumnsNaming():
    
    def name_is_one_of_11_ht_reps(self, x):
        if (re.search('fbf\d', x) and (
            re.search('oo', x) or re.search('sp', x))):
            return True
        #elif (re.search('control.*lane.*', x)):
        #    return True
        else:
            return False

    def name_is_a_combined_ht_control(self, x):
        # Combined: control_sp_counts.txt
        # control_oo_counts.txt
        # control_n2_counts.txt
        if (re.search('control_oo_counts', x) or \
            re.search('control_sp_counts', x)):
            return True
        return False

    def name_is_a_combined_lt_control(self, x):
        if re.search('control_n2_counts', x):
            return True
        return False
    
    def shorten_names_and_subset_columns(self):
        def clean(n):

            n = re.sub('_counts.txt', '', n)
            n = re.sub('exp_fbf1_sp', 'SP FBF1', n)
            n = re.sub('exp_fbf2_sp', 'SP FBF2', n)
            n = re.sub('exp_fbf1_oo', 'OO FBF1', n)
            n = re.sub('exp_fbf2_oo', 'OO FBF2', n)
            n = re.sub('exp_fbf_sp', 'SP FBF', n)
            n = re.sub('exp_fbf_oo', 'OO FBF', n)
            n = re.sub('contol_', 'c_', n)
            
            #n = re.sub('fbf', 'F', re.sub('rt.*and.*', '', re.sub('exp_', '',
            #        re.sub('control_', 'c_', re.sub('_counts.txt', '',
            #        re.sub('[A-Z]{4}', 'i', n)))
            #              )))


            #n = re.sub('F_sp_', 'SP FBF', n)
            #n = re.sub('F_oo_', 'OO FBF', n)
            n = re.sub('exp_fbf1_TGGC', 'LT FBF1_1', n)
            n = re.sub('exp_fbf1_GGTT', 'LT FBF1_2', n)
            n = re.sub('exp_fbf1_CGGA', 'LT FBF1_3', n)  
            n = re.sub('exp_fbf2_CGGA', 'LT FBF2_1', n)
            n = re.sub('exp_fbf2_GGTT', 'LT FBF2_2', n)
            n = re.sub('exp_fbf2_TGGC', 'LT FBF2_3', n)

            # Replicates for FBF1 and FBF2 at 20C:
#exp_bed1: ./bed_collapsed/exp_fbf1_TGGC.bed
#exp_bed2: ./bed_collapsed/exp_fbf1_GGTT.bed
#exp_bed3: ./bed_collapsed/exp_fbf1_CGGA.bed

#exp_bed1: ./bed_collapsed/exp_fbf2_CGGA.bed
#exp_bed2: ./bed_collapsed/exp_fbf2_GGTT.bed
#exp_bed3: ./bed_collapsed/exp_fbf2_TGGC.bed

            return n  # 12 FBF + 12 controls
        
        known = []
        for i, col in enumerate(self.counts_df.columns):
                  
            if col == 'ave_neg':
                known.append(x)
                continue
                
            #rep = 1
            name = clean(col) #+ '_' + str(rep)
            
            #while (name in known):
            #    rep += 1
            #    if rep > 9:
            #        break
            #    name = name[:-1] + str(rep)
            known.append(name)
        
        print('>', known)
        self.counts_df.columns = known