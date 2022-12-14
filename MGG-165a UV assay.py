import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import easygui
#pip install --upgrade easygui
import csv
import os
print('deez')
#############################

# file_path_UV = easygui.fileopenbox(msg="Choose a UV file", filetypes= "*.txt", multiple=True)
# file_path_UV_lac = easygui.fileopenbox(msg="Choose a UV file for laccase", filetypes= "*.txt", multiple=True)
# file_path_UV
# file_path_UV = sorted(file_path_UV)
# file_path_UV_lac = sorted(file_path_UV_lac)
############################

def reduce_columns(data_raw_input, ll = 250):
    
    ''' Reduces columns (delete columns that are wavelenth duplicates) and defines working wavelength range. Default range is 200-700nm'''
    
    dataframe = data_raw_input.copy(deep = True)
    for i in range(2,dataframe.shape[1],2):                 
        dataframe.drop(i, axis=1, inplace=True)  # remove nm column duplicates
        dataframe = dataframe[dataframe[0] > 220]  # lower limit (from 220nm onwards as default)
        
    return dataframe

def norm_min(dataframe):
    
    ''' Normalise by minimum value '''

    #norm_data = dataframe/dataframe.min() # normalise each column with the lowest value
    
    #norm_data[0] = dataframe[0]           # keep wavelengths columns unaffected
    
    
#     return norm_data
    wavelength = dataframe[0]
    dataframe.drop(dataframe.columns[0], axis=1, inplace = True) # remove wavelegnth column
    base = dataframe.iloc[0]                                     # selected lowest values of the spectra @ 800nm, i.e. baseline
    dataframe = dataframe-base  
    dataframe.insert(0,0,wavelength)

        
    return dataframe

def get_abs_at(data, wavelength = 576):
    
    '''Gets the values of a peak defined at the given wavelenth for each scan'''
    
    peaks = data[round(data[0],0) == wavelength].T[1:]
    
    return peaks.copy(deep= True)


def abs_to_conc(absorbance, coef):

    ''' Beer-Lambert eq. Converts absorbance to concentration. coef = molar extintion coeficient (in M-1.cm-1)'''
    
    length = 1                        # path lentgh = 1 cm
   
    conc = absorbance/(coef*length)   # molar concentration (M)
    
    
    return conc
    
    

def remove_status(dataframe):
    
    '''Removes the status report that the uv-vis or fluorimter give for each scan'''
    
    dataframe.dropna(axis=1, how='all', inplace=True) # removes padding NaN columns at the end of the data
    dataframe.dropna(axis=0, how='any', inplace=True) # actually removes the status reports
    dataframe = dataframe.astype(float)                  # converts each value to numeric (when imported each element is a string of numbers)
    
    return dataframe

def samples_name(path_list):
    
    ''' Extracts the name/code of the sample from the file's name'''
    
    samples_name = []
    
    for _ in range(len(file_path_UV)):
        selected_file = os.path.basename(os.path.normpath(file_path_UV[_]))
        file_name     = selected_file[5:21]
        samples_name.append(file_name)
        
    return samples_name


########################


coef_Resorufin = 54000
Laccase_ext_coef = 16953.47034 # M.cm-1                  e
time = list(range(40))

y    = pd.DataFrame({'time' : time})
lac    = []

for i in range(len(file_path_UV)):
    
    raw_data = pd.read_csv(file_path_UV[i], header=None, sep=',', skiprows=2)
    raw_data    = remove_status(raw_data)  
    data        = reduce_columns(raw_data)   
    norm_data   = norm_min(data)
    abs_RR      = get_abs_at(norm_data, 576)
    
    conc_RR     = abs_to_conc(abs_RR,coef_Resorufin)*10**6  # in ??M
    conc_RR.reset_index(drop = True, inplace = True)
    conc_RR     = conc_RR-conc_RR.iloc[0]
    #Change data range
    conc_RR     = conc_RR[:40]
    #######


    
    raw_data_lac  = pd.read_csv(file_path_UV_lac[i], header=None, sep=',', skiprows=2)
    raw_data_lac  = remove_status(raw_data_lac)
    abs_lac       = get_abs_at(raw_data_lac, 270)
    conc_lac      = abs_to_conc(abs_lac,Laccase_ext_coef)*10**6 # in ??M
    conc_lac      = conc_lac.values.tolist()
    
    y   = pd.concat([y,conc_RR], axis = 1)
    lac.append(conc_lac[0][0])             # [0][0] because of of nested list result from.values.tolit()
  
legend = samples_name(file_path_UV)
#legend = ['MGG165a-4-2 AT 95','MGG-165a-5-2 DT 95','MGG-165a-6-2 H 95']

legend_lac = samples_name(file_path_UV_lac)

y.columns = ['time'] + legend # tidy up columns
y.drop(y.columns[0], axis=1, inplace = True)

y_final = []
for _ in range(y.shape[1]): 
    
    z = y[legend[_]]/lac[_]
    z = z.to_numpy()
    y_final.append(z)



c = ['tab:blue','tab:green','tab:red','k','tab:gray']
m = ['s','d','^','o','o']
mksize = 8 #65
ftsize = 22

fig, ax = plt.subplots(figsize=(6.4, 4.8))
for _ in range(len(file_path_UV)):
    if _ < 3:
        plt.plot(time,y_final[_], c = c[_], marker = m[_], ms = mksize, alpha = 0.7)#, mfc = 'None', mew = 2)
        #plt.scatter(time,y_final[_], c = c[_], marker = m[_], s=70, facecolors='none', alpha = 0.65)
    else:
        plt.plot(time,y_final[_], c = c[_], marker = m[_], ms = mksize, alpha = 0.65, mfc = 'None', mew = 2)

    
#plt.rcParams.update({'font.size': 10}) 
#plt.legend(legend)
plt.xlabel('time / $min$', fontsize = ftsize)
plt.ylabel('[Resorufin] / $??M$', fontsize = ftsize)
plt.ylim(-0.004,0.05)
plt.title('Laccase activity')

plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)

ax.annotate(
        annot, xy=(0.9, 0.85), xycoords='axes fraction',
        va="center", ha="center",
        bbox=dict(boxstyle='circle', pad=0.8, fc="w"),
        fontsize = 16
        )

############################### MGG-165a vesciles 95% ############### 

# file_path_UV =normal_5
# file_path_UV_lac =normal_5_lac

normal_5 = [

'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-1 AT5 +AR_t0-to-60min.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-2-1 DT5 +AR_t0-to-60min.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-3 H5 +AR_t0-to-60min.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-7 LACC 2mgmL+AR_t0-to-60min.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-8 AR1mgml ONLY_t0-to-60min.csv'
]
normal_5_lac = [
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-1 AT5 _100uLin200uL.csv',
 'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-2-1 DT5 _100uLin100uL.csv',
 'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-3 H5 _100uLin200uL.csv',
 'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-7 LACC 2mgmL 100ulin100ul.csv',
 'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-8 AR1mgml ONLY 5ulin100ul.csv'
 ]

 ############################### MGG-165a vesciles 95% ###############

# file_path_UV =normal_95
# file_path_UV_lac =normal_95_lac
# annot = '95%'

normal_95 = [
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-4-2 AT 95  _ unstirred_+AR_t0-to-60min.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-5-2 DT 95  _ unstirred_+AR_t0-to-60min.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-6-2 H 95  _ unstirred_+AR_t0-to-60min.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-7 LACC 2mgmL+AR_t0-to-60min.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-8 AR1mgml ONLY_t0-to-60min.csv'
]
normal_95_lac = [
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-4-2 AT 95  _unstirred_100ulin100uL.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-5-2 DT 95  _unstirred_100ulin100uL.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-6-2 H 95  _unstirred_100ulin100uL.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-7 LACC 2mgmL 100ulin100ul.csv',
'C:\\Users\\drb18182\\Desktop\\Docs\\Data\\uv-vis\\MGG-165 [GPSomes, 5% & 95%, TF, LAc nanoreactors]\\[MGG-165a] [UV] assay\\[UV] MGG-165a-8 AR1mgml ONLY 5ulin100ul.csv'
]

# file_path_UV =normal_95
# file_path_UV_lac =normal_95_lac