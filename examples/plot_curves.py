SIMPLE = 1
LAST_GRAPH_SPECIAL = 0
LAST_GRAPH_STYLE = ',g'

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

REL_ERROR = 0.05

profiles = sys.argv[1:]
labels = profiles
fun = lambda x: pd.DataFrame(pd.read_table(x, header=None, sep='\s+'));
dfs = [fun(x) for x in profiles]

def chi_score(df1, df2):
	left  = max(df1[0][0], df2[0][0])
	right = min(df1[0].iloc[-1], df2[0].iloc[-1])
	
	master = df2
	slave  = df1
	
	flag   = True
	
	qnum1 = sum((df1[0] >= left) & (df1[0] <= right))
	qnum2 = sum((df2[0] >= left) & (df2[0] <= right))
	
	if qnum1 > qnum2:
		master = df1
		slave  = df2
		flag   = False	
		
	if qnum1 == qnum2:
		if df1.shape[1] == 3: 
			master = df1 
			slave = df2 
			flag = False
			
	master = master.copy()
	slave  = slave.copy()
	
	slave  = slave[(slave[0] >= left) & (slave[0] <= right)]
	master = master[(master[0] >= slave[0].iloc[0]) & (master[0] <= slave[0].iloc[-1])]

	if slave.shape[1] == 2: slave['err'] = slave[1]*REL_ERROR
	if master.shape[1] == 2: master['err'] = master[1]*REL_ERROR
	
	m_int = []
	m_err = []
	s_int = []
	q     = []
	m_id  = 0
	
	for s_id in range(slave.shape[0]-1):
		while master.iloc[m_id, 0] < slave.iloc[s_id+1, 0]:
			q.append(master.iloc[m_id, 0])
			m_int.append(master.iloc[m_id, 1])
			m_err.append(master.iloc[m_id, 2])
			s_tan = (slave.iloc[s_id+1, 1]-slave.iloc[s_id, 1])/(slave.iloc[s_id+1, 0]-slave.iloc[s_id, 0])
			s_int.append(s_tan*(master.iloc[m_id, 0]-slave.iloc[s_id, 0])+slave.iloc[s_id, 1])
			m_id+=1
			if m_id >= master.shape[0]:
				break;
			
	numer = 0.
	denom = 0.
	
	for i in range(len(q)):
		numer += s_int[i] * m_int[i] / (m_err[i] * m_err[i])
		denom += s_int[i] * s_int[i] / (m_err[i] * m_err[i])
	c = numer/denom
	
	chi = 0.
	for i in range(len(q)):
		chi += pow((m_int[i]-c*s_int[i])/m_err[i], 2)
		
	if flag:
		c = 1/c
		
	#input()
		
	return np.sqrt(chi/len(q)), c
	
	
coefs = [0.]*len(dfs)
		
def score_table(labels, dfs):
    table = pd.DataFrame(columns=labels, index=labels)
    for ind in range(len(labels)):
        for col in range(len(labels)):
            table.iloc[ind, col], coefs[col] = chi_score(dfs[ind], dfs[col])
    return table
            
print('Score table')
print(score_table(labels, dfs).to_string(index=False))
print()

def plot_curves(labels, dfs):
    for i in range(len(dfs)-1): 
        coef = coefs[i]
        plt.plot(dfs[i][0], list(map(lambda x: np.log(x * coef), dfs[i][1])), label=labels[i])
        
    if (LAST_GRAPH_SPECIAL):
    	plt.plot(dfs[-1][0], list(map(lambda x: np.log(x), dfs[-1][1])), LAST_GRAPH_STYLE, label=labels[-1])
    else:
    	plt.plot(dfs[-1][0], list(map(lambda x: np.log(x), dfs[-1][1])), label=labels[-1])
    
    plt.legend()
    plt.xlabel('q')
    plt.ylabel('Intencity')
    plt.title('SAXS')
    plt.show()
    
def plot_curves_simple(labels, dfs):
    q_start = max([df[0][0] for df in dfs])
    
    for i in range(len(dfs)-1): 
        id_start = np.argmin(list(map(lambda x: np.abs(x - q_start), dfs[i][0])))
        plt.plot(dfs[i][0], list(map(lambda x: np.log(x * dfs[-1][1][0] / dfs[i][1][id_start]), dfs[i][1])), label=labels[i])
        
    if (LAST_GRAPH_SPECIAL):
    	plt.plot(dfs[-1][0], list(map(lambda x: np.log(x), dfs[-1][1])), LAST_GRAPH_STYLE, label=labels[-1])
    else:
    	plt.plot(dfs[-1][0], list(map(lambda x: np.log(x), dfs[-1][1])), label=labels[-1])
        	
    plt.legend()
    plt.xlabel('q')
    plt.ylabel('Intencity')
    plt.title('SAXS')
    plt.show()
    
def plot_curves_super_simple(labels, dfs):
	for i in range(len(dfs)-1): 
		plt.plot(dfs[i][0], list(map(lambda x: np.log(x), dfs[i][1])), label=labels[i])

	if (LAST_GRAPH_SPECIAL):
		plt.plot(dfs[-1][0], list(map(lambda x: np.log(x), dfs[-1][1])), LAST_GRAPH_STYLE, label=labels[-1])
	else:
		plt.plot(dfs[-1][0], list(map(lambda x: np.log(x), dfs[-1][1])), label=labels[-1])
			
	plt.legend()
	plt.xlabel('q')
	plt.ylabel('Intencity')
	plt.title('SAXS')
	plt.show()

plt.figure(figsize=(12, 9))
if SIMPLE == 2:
	plot_curves_super_simple(labels, dfs)
if SIMPLE == 1:
	plot_curves_simple(labels, dfs)
if SIMPLE == 0:
	plot_curves(labels, dfs)


