import os
import sys
import glob
import pandas as pd




def main():
    
    druntime_s = {}
    druntime_g = {}
    
    for log_o_fpath in sorted(glob.glob('*')):
        if '.o' not in log_o_fpath:
            continue
        log_e_fpath = log_o_fpath.replace('.o', '.e')
        
        data_type = log_o_fpath.split('.')[0].split('_')[2]
        if data_type[-1] == 'g':
            data_type = data_type[:-1]
        arch_type = None
        data_lv = None
        with open(log_e_fpath) as infh:
            for buf in infh:
                if 'model:' in buf:
                    arch_type = buf.split(';')[2].split(':')[1]
                if 'n_class:' in buf:
                    if 'n_class:86' in buf:
                        data_lv = 'g'
                    elif 'n_class:204' in buf:
                        data_lv = 's'
                
                if arch_type is not None and data_type is not None:
                    k = arch_type + '__' + data_type
                if 'user' in buf:
                    if k not in druntime_s:
                        druntime_s[k] = []
                    if k not in druntime_g:
                        druntime_g[k] = []
                    
                    runtime = buf.replace(' ', '').replace('user', '')
                    runtime_m = float(runtime.split('m')[0])
                    runtime_s = float(runtime.split('m')[1].replace('s', ''))
                    if data_lv == 's':
                        druntime_s[k].append(runtime_m / 60 + runtime_s / 60 / 60)
                    elif data_lv == 'g':
                        druntime_g[k].append(runtime_m / 60 + runtime_s / 60 / 60)
                
                if 'sys ' in buf:
                    runtime = buf.replace(' ', '').replace('sys', '')
                    runtime_m = float(runtime.split('m')[0])
                    runtime_s = float(runtime.split('m')[1].replace('s', ''))
                    if data_lv == 's':
                        druntime_s[len(druntime_s) - 1] += runtime_m / 60 + runtime_s / 60 / 60
                    elif data_lv == 'g':
                        druntime_g[len(druntime_s) - 1] += runtime_m / 60 + runtime_s / 60 / 60
 

   
    d = []
    m = []
    t = []
    n_t = {}
    for k, v in druntime_s.items():
        if k not in n_t:
            n_t[k] = 0
        for _ in v:
            if n_t[k] < 10:
                d.append(k.split('__')[1])
                m.append(k.split('__')[0])
                t.append(_)
            n_t[k] += 1
    
    dfs = pd.DataFrame({'data': d, 'model': m, 'time': t})


    d = []
    m = []
    t = []
    n_t = {}
    for k, v in druntime_g.items():
        if k not in n_t:
            n_t[k] = 0
        for _ in v:
            if n_t[k] < 10:
                d.append(k.split('__')[1])
                m.append(k.split('__')[0])
                t.append(_)
            n_t[k] += 1
    
    dfg = pd.DataFrame({'data': d, 'model': m, 'time': t})
    
    dfs.to_csv('species.log', sep='\t', header=True, index=False)
    dfg.to_csv('genus.log', sep='\t', header=True, index=False)
    
            
        
if __name__ == '__main__':
    main()
    
    


