import pandas as pd
import numpy as np
import scipy.stats as st

sheet_id = '1aqHnDaiRXntyrbj6qaqzR-4XnLbvtKU84FB_JLu5x5A'
sheet_name = 'EXP_8_second_main_factory'
url = f'https://docs.google.com/spreadsheets/d/{sheet_id}/gviz/tq?tqx=out:csv&sheet={sheet_name}'
exp8fac = pd.read_csv(url)

exp8fac = exp8fac.iloc[9:13,].reset_index(drop=True)

code = []
mut = []
temp = []
time = []
fert = []
swim = []

for row in exp8fac.index:
    total_count = int(exp8fac.loc[row, 'TOTAL'])
    fert_count = int(exp8fac.loc[row, 'FERT.1'])
    swim_count = int(exp8fac.loc[row, 'SWIM'])
    
    code += [exp8fac.loc[row, 'CODE']] * total_count
    mut += [exp8fac.loc[row, 'MUT']] * total_count
    temp += [exp8fac.loc[row, 'TEMP']] * total_count
    time += [exp8fac.loc[row, 'TIME']] * total_count
    
    fert += [1] * fert_count
    fert += [0] * (total_count - fert_count)
    swim += [1] * swim_count
    swim += [0] * (total_count - swim_count)

ind_data = pd.DataFrame({'CODE': code, 'MUT': mut, 'TEMP': temp, 'TIME': time, 'FERT': fert, 'SWIM': swim})

print('Start bootstrapping...')
bs_cc = st.bootstrap((ind_data[ind_data['CODE'] == 'CC']['SWIM'], ), np.count_nonzero, n_resamples=1000)
print(bs_cc)
bs_c38 = st.bootstrap((ind_data[ind_data['CODE'] == 'C38']['SWIM'], ), np.count_nonzero, n_resamples=1000)
print(bs_c38)
bs_mc = st.bootstrap((ind_data[ind_data['CODE'] == 'MC']['SWIM'], ), np.count_nonzero, n_resamples=1000)
print(bs_mc)
bs_m38 = st.bootstrap((ind_data[ind_data['CODE'] == 'M38']['SWIM'], ), np.count_nonzero, n_resamples=1000)
print(bs_m38)
    