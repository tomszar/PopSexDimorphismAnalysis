'''
Definition of set of functions to create admix plots
First, organize figure in terms of number of Qfiles. 
Then, for each Qfile, create its corresponding axis, already sorted 
Try to keep only one sorting options across ks
Finally, add text to figures
'''

def sort_qfile(Qfile_unsorted):
	'''
	Function to sort Qfiles. It can be used with arbitrary populations. If pop column exist, it will sort first by pop, super_pop, and then by ks. 
	This will also return the position of pop and super pop labels. If pop column does not exist, then sort only by highest K
	Usage
	Input:
		- Qfile_unsorted: the Qfile to be sorted. It can contain a super_pop and pop column to sort first by those. 
						  If those are not present, the sorting will be made only through the highest K
	Output:
		- Qfile_sorted: the sorted Qfile
		- Qfile_sort_indices: index sequence to sort the Qfiles. It can be applied to other Qfiles to maintain the order

	'''
	
	import pandas as pd
	import numpy as np
	
	Qfile_sorted       = pd.DataFrame()
	Qfile_sort_indices = pd.Series()
	
	if 'pop' in Qfile_unsorted.columns:
		sortpops = ['super_pop', 'pop']
		Qfile_sorted_initial = Qfile_unsorted.sort_values(by = sortpops)
		list_superpop = Qfile_sorted_initial['super_pop'].unique()
		list_pop      = Qfile_sorted_initial['pop'].unique()
		loc_superpop  = []
		loc_pop       = []
		for l in list_superpop:
			#Getting location of superpop with respect to data
			pop_bool = Qfile_sorted_initial['super_pop'] == l
			pop_ind  = [i for i, x in enumerate(pop_bool) if x]
			loc_lab  = ( min(pop_ind) + max(pop_ind) ) / 2
			loc_superpop.append(loc_lab)
		for l in list_pop:
			#Getting location of pop with respect to data
			pop_bool = Qfile_sorted_initial['pop'] == l
			pop_ind  = [i for i, x in enumerate(pop_bool) if x]
			loc_lab  = ( min(pop_ind) + max(pop_ind) ) / 2
			loc_pop.append(loc_lab)
				  
			#Sorting each pop by highest k, and then the remaining ones
			group      = Qfile_sorted_initial.iloc[pop_ind,:]
			sortingby  = [i[0] for i in sorted(enumerate(group.mean()), reverse=True, key=lambda x:x[1])] #Sort by most to least relevant cluster
			indices    = pd.Series(group.sort_values(by = sortingby).index)
			group      = group.sort_values(by = sortingby).reset_index(drop=True)
			Qfile_sorted = pd.concat([Qfile_sorted, group])
			Qfile_sort_indices = Qfile_sort_indices.append(indices)
			
		return(Qfile_sorted, Qfile_sort_indices, loc_superpop, loc_pop, list_superpop, list_pop)
	
	elif 'pop' not in Qfile_unsorted.columns:
		ks = Qfile_unsorted.shape[1] #Number of groups
		cols_to_group = Qfile_unsorted.idxmax(axis=1) #Create list of column with highest value per row
		for i in range(ks):
			#Sorting individuals by separating them in their k with highest proportion, and sorting by the other ks next
			#We'll paste them together to create a new sorted dataframe
			group      = Qfile_unsorted[cols_to_group == i]
			sortingby  = [i[0] for i in sorted(enumerate(group.mean()), reverse=True, key=lambda x:x[1])] #Sort by most to least relevant cluster
			indices    = pd.Series(group.sort_values(by = sortingby).index)
			group      = group.sort_values(by = sortingby).reset_index(drop=True)
			Qfile_sorted = pd.concat([Qfile_sorted, group])
			Qfile_sort_indices = Qfile_sort_indices.append(indices)
		
		return(Qfile_sorted, Qfile_sort_indices)
		
def plot_admix(Qfiles, filename="plot.png", popfile=None, sort_by_k=None, font_sizes=(10,12,10)):
	'''
	Function to plot admixture files
	Qfiles can be either a list of pandas Dataframes containing the Q proportions or a single file
	Usage
	Input:
		- Qfiles: list or dataframe of admix proportions
		- filename: name to give the output
		- popfile: pop and super_pop file to order the samples
		- sort_by_k: sort by a specific number of K. If None, the highest will be used
		- font_sizes: font sizes for K label, upper, and lower label respectively
	Output:
		- Save the file externally
	'''
	
	import pandas as pd
	import numpy as np
	import matplotlib.pyplot as plt
	
	width = 1
	fontsize_k  = font_sizes[0]
	fontsize_up = font_sizes[1]
	fontsize_lo = font_sizes[2]
	
	if isinstance(Qfiles, list): #if there is a list of pandas dataframes
		qs     = len(Qfiles)
		fig    = plt.figure(figsize=(20, 1*qs), dpi=300)
		loc_ks = 1/qs/2
		#Doing sorting in here, keep the same sort across ks
		if popfile is not None:
			#Sort samples based on known populations and ks from last list
			Qfile = pd.concat([Qfiles[qs-1], popfile], axis = 1)
			Qfile_sorted, Qfile_sort_indices, loc_superpop, loc_pop, list_superpop, list_pop = sort_qfile(Qfile)
		if popfile is None:
			if sort_by_k is None:
				sort_by_k = qs - 1                   
			Qfile_sorted, Qfile_sort_indices = sort_qfile(Qfiles[sort_by_k])
			
		for i in range(qs):
			Qfile_toplot = Qfiles[i].loc[Qfile_sort_indices,:].reset_index(drop=True)
			ks    = Qfile_toplot.shape[1]
			ind   = np.arange(Qfile_toplot.shape[0]) #Number of samples
			ax    = fig.add_subplot(qs, 1, i+1)
			for l in range(ks):
				this_stack = Qfile_toplot[l]
				if l > 0:
					bottom_stack = Qfile_toplot.iloc[:,0:l].sum(axis=1)
				else:
					bottom_stack = None
				ax.bar(ind, this_stack, bottom = bottom_stack, width = width) 
				ax.axis('off')
				ax.autoscale(tight=True)
				
			if i == 0 and popfile is not None:
				#if first plot, draw the super_pop legends
				for l in range(len(loc_superpop)):
					plt.text(loc_superpop[l], 1.2, list_superpop[l], fontsize = fontsize_up)
			if i == (qs-1) and popfile is not None:
				for l in range(len(loc_pop)):
					plt.text(loc_pop[l], -1, list_pop[l], fontsize = fontsize_lo, rotation=90)
			
			#Annotating number of ks
			trans = ax.get_yaxis_transform() # x in data untis, y in axes fraction
			k_text = "k=" + str(ks)
			ann   = ax.annotate(k_text, xy=(-0.03, 0.4 ), xycoords=trans, fontsize= fontsize_k)
					
	else:
		if popfile is not None:
			#Sort samples based on known populations and ks from last list
			Qfile = pd.concat([Qfiles, popfile], axis = 1)
			Qfile_sorted, Qfile_sort_indices, loc_superpop, loc_pop, list_superpop, list_pop = sort_qfile(Qfile)
		if popfile is None:
			Qfile_sorted, Qfile_sort_indices = sort_qfile(Qfiles)
	   
		Qfile_toplot = Qfiles[i].loc[Qfile_sort_indices,:].reset_index(drop=True)
		fig   = plt.figure(figsize=(20, 1), dpi=300)
		ks    = Qfile_toplot.shape[1]
		ind   = np.arange(Qfile_toplot.shape[0]) #Number of samples
		ax    = fig.add_subplot(111)
		for i in range(ks):
			this_stack = Qfile_toplot[i]
			if i > 0:
				bottom_stack = Qfile_toplot.iloc[:,0:i].sum(axis=1)
			else:
				bottom_stack = None
			ax.bar(ind, this_stack, bottom = bottom_stack,width = width)
			ax.axis('off')
			ax.autoscale(tight=True)
		
		#Annotating number of ks
		trans = ax.get_yaxis_transform() # x in data untis, y in axes fraction
		k_text = "k=" + str(ks)
		ann   = ax.annotate(k_text, xy=(-0.03, 0.4 ), xycoords=trans, fontsize=fontsize_k)
		if popfile is not None:
			for i in range(len(loc_superpop)):
				plt.text(loc_superpop[i], 1.2, list_superpop[i], fontsize=fontsize_up)
			for i in range(len(loc_pop)):
				plt.text(loc_pop[i], -1, list_pop[i], fontsize=fontsize_lo, rotation=90)    
		
	#plt.subplots_adjust(left=0.12)
	
	plt.savefig(filename, dpi = 300, bbox_inches='tight')
	#return(fig, Qfile_sort_indices)
	
	
def plot_admix2(Qfiles, ax, popfile=None, sort_by_k=None, horizontal=True):
	'''
	Function to plot admixture files in horizontal bars. Qfiles can be either a list of pandas Dataframes containing the Q proportions or a single file
	Usage
		- Qfiles: list or dataframe of admix proportions
		- ax: matplotlib ax
		- popfile: pop and super_pop file to order the samples
		- sort_by_k: sort by a specific number of K. If None, the highest will be used	
	'''
	
	import pandas as pd
	import numpy as np
	import matplotlib.pyplot as plt
	
	width = 1
	
	if popfile is not None:
		#Sort samples based on known populations and ks from last list
		Qfile = pd.concat([Qfiles, popfile], axis = 1)
		Qfile_sorted, Qfile_sort_indices, loc_superpop, loc_pop, list_superpop, list_pop = sort_qfile(Qfile)
	if popfile is None:
		Qfile_sorted, Qfile_sort_indices = sort_qfile(Qfiles)
	   
	Qfile_toplot = Qfiles.loc[Qfile_sort_indices,:].reset_index(drop=True)
	ks    = Qfile_toplot.shape[1]
	ind   = np.arange(Qfile_toplot.shape[0]) #Number of samples
	for i in range(ks):
		this_stack = Qfile_toplot[i]
		if i > 0:
			left_stack = Qfile_toplot.iloc[:,0:i].sum(axis=1)
		else:
			left_stack = None
		if horizontal==True:
			ax.barh(ind, this_stack, left = left_stack, height = width)
		else:
			ax.bar(ind, this_stack, bottom = left_stack, width = width)
		ax.axis('off')
		ax.autoscale(tight=True)

