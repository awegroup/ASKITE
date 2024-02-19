import numpy as np

def fetch_surfplan_file(SURFPLAN_FILENAME,sys_path):
    file_path_surfplan = sys_path[-1]+'/data/surfplan_files/kitepower_confidential/'+str(SURFPLAN_FILENAME)+'.obj'
    surfplan_file = np.load(file_path_surfplan,allow_pickle=True)
    return surfplan_file