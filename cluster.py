# Script that takes normalized data from each tool and makes swarmplot with variances in coexpression clusters

import numpy as np
import pandas as pd
import glob
import pdb
from collections import defaultdict
from sklearn.metrics import adjusted_rand_score
import matplotlib.pyplot as plt

def swarmplot(yList, xCenter, markersize, pW, pH, xrang, yrang):
	

	minDist = 1/72*markersize
	minYdist = (minDist*(yrang[1]-yrang[0]))/pH
	minXdist = (minDist*(xrang[1]-xrang[0]))/pW

	shift = ((minDist/10)/pW)*(xrang[1]-xrang[0])

	placedPoints = []
	stop = False
	for yPos in yList:
		if stop:
			break
		xLeft = xCenter-shift
		xRight = xCenter+shift
		if not placedPoints:
			placedPoints.append((xCenter, yPos))
		else:
			placed = False
			conflict = False
			for point in placedPoints:
				if np.abs(point[1]-yPos) < minYdist and np.abs(point[0]-xCenter) < minXdist:
					conflict = True
					break
			if not conflict:
				placedPoints.append((xCenter, yPos))
				placed = True
			while not placed:	
				left_conflict = False
				for point in placedPoints:
					if np.abs(point[1]-yPos) < minYdist and np.abs(point[0]-xLeft) < minXdist:
						left_conflict = True
						break
				right_conflict = False
				for point in placedPoints:
					if np.abs(point[1]-yPos) < minYdist and np.abs(point[0]-xRight) < minXdist:
						right_conflict = True
						break
				placed = True
				if not left_conflict and not right_conflict:
					placedPoints.append(([xLeft, xRight][np.random.randint(0,2)], yPos))
				elif not left_conflict:
					placedPoints.append((xLeft, yPos))
				elif not right_conflict:
					placedPoints.append((xRight, yPos))
				else:
					placed = False
					xLeft -= shift
					xRight += shift
					if xCenter-xLeft > 0.45:
						stop=True
						break
	return placedPoints



figureWidth = 8
figureHeight = 8

plt.figure(figsize=(figureWidth, figureHeight))


panelWidth = 7
panelHeight = 7

relativePanelWidth = panelWidth/figureWidth
relativePanelHeight = panelHeight/figureHeight

panel1=plt.axes([1/10, 1/10, relativePanelWidth, relativePanelHeight])

panel1.set_xlim(0, 4)
panel1.set_xticks(np.arange(1,4))
panel1.set_xticklabels(['HTseq', 'Kallisto', 'MISO'], size=12)
panel1.set_xlabel('Quantification Tool', size=14)

panel1.set_ylim(0,15)
panel1.set_yticks(np.arange(0, 14, 2))
panel1.set_yticklabels(np.arange(0, 14, 2), size=12)
panel1.set_ylabel(r'Log$_{2}$(Cluster Standard Deviation)', size=14)



scores = {}
truth = pd.read_csv('ensembl_truth.tab', sep='\t')
cluster_genes = {}
for cluster in range(60):
	cluster_genes[cluster] = truth[truth['cluster'] == cluster]['target_id']
for i,tool in enumerate(['htseq', 'new_kallisto', 'miso']):
	print(tool)
	data = pd.read_csv(f'{tool}/merged_tpms.tab', sep='\t')
	data = data[data['target_id'].isin(truth['target_id'])]
	#data['tpm'] = np.log2(data['tpm']+np.finfo(float).eps)
	stds = []
	for cluster in range(60):
		stds.append(np.std(data[data['target_id'].isin(cluster_genes[cluster])]['tpm']))
	stds = np.log2(np.asarray(stds)+1)

	points = swarmplot(stds, i+1, 8, panelWidth, panelHeight, (0,4), (0,15))
	plt.plot([x[0] for x in points], [x[1] for x in points], marker='o', linewidth=0, markersize=8, color='black')
	avg_std = np.median(stds)
	plt.plot([i+0.5, i+1.5], [avg_std, avg_std], color='red', linewidth=2)
	scores[tool] = avg_std
plt.show()

