import MDAnalysis as mda 
import numpy as np
import matplotlib.pyplot as plt
import copy

u = mda.Universe('EQ_chain.pdb','red.xtc')

protein_length=len(u.atoms)
protein = u.select_atoms('all')
print(len(u.trajectory))

master_array = np.zeros(protein_length)


def determine_pairs(input_value):
	list_1=[]
	list_2=[]
	for i in range(protein_length):
		for j in range(i+1,protein_length):

			if abs(i-j) == input_value:
				list_1.append(i)
				list_2.append(j)
	return list_1, list_2

pairlist_1=[]
pairlist_2=[]
for i in range(protein_length):
	pair_1,pair_2 = determine_pairs(i+1)
	pairlist_1.append(copy.deepcopy(pair_1))
	pairlist_2.append(copy.deepcopy(pair_2))




for ts in u.trajectory:
	for i in range(protein_length):
		arr_distances=[]
		pairs_1 = pairlist_1[i]
		pairs_2 = pairlist_2[i]
		for j,k in zip(pairs_1,pairs_2):
			diff = protein.positions[j]-protein.positions[k]
			distance_vec = np.sqrt(sum([coord**2 for coord in diff]))
			arr_distances.append(distance_vec)
		master_array[i]+=np.mean(arr_distances)
	print(ts.frame)


master_array = [i/len(u.trajectory) for i in master_array]
dist_values = [i+1 for i in range(protein_length)]



ax = plt.subplot(111)
ax.plot(dist_values, master_array, 'r--', lw=2, label=r"$r_{ij}$")
ax.set_xlabel("distance |i-j|")
ax.set_ylabel(r"internal scaling distance $R_{ij}$ ($\AA$)")
ax.set_xlim(0,protein_length)
ax.figure.savefig("int_dist.pdf")
plt.draw()



#for ts in u.trajectory:
