from Bio import Phylo
import sys

''' 
Dette program tager to sys argumenter, som begge er træer 
i newick format:
	sys argument 1 = tree1
	sys argument 2 = tree2
Programmet printer DF-afstanden imellem de to træer.
Det er muligt, at få programmer til at printe "en tegning" af begge 
træer og få den til at printe DF-intervallerne for træerne.
(Så skal man ændre i koden forneden)
''' 

def rename (tree1):
	''' 
	rename() tager et træ som argument (tree1), bestemmer roden af 
	træet og retunerer df_rename() af dette rodet træ.
	Funktionen df_rename(), tager en clade som argument
	(her roden af tree1) og retunerer en liste med to lister. 
	Den første liste indeholder bladenes originale navne og den
	anden liste indeholder "Depth-First numbering". 
	''' 
	return df_rename(tree1.clade)

def df_rename(clade):
	if clade.is_terminal():	
		first_name.append(clade.name)
		new_name.append(len(new_name)+1)
	else:
		for child in clade.clades:
			df_rename(child)
		return [first_name,new_name]


def RF_distance(tree1,tree2):
	''' 
	Rf_distance() tager tree1 og tree2 som argument og retunerer 
	RF-afstanden imellem disse (Antallet af splits i tree1, som ikke 
	er i tree2 + antallet af splits i tree2, som ikke er i tree1)
	intervals() tager et træ som argument og retunerer en liste med 
	længden to. På index [0] er en sorteret liste med DF-intervallerne 
	for hver node. Numrene laves ud fra listen new_name, som laves med 
	funktionen rename(). Hvis to forskellige træer køres i rename() 
	og intervals() sorteres de intervaller fra, som ikke opfylder 
	"max-min+1=size". Antallet af intervaller der er sorteret fra står
	på listens index [1].
	intervals() kalder funktionen df_intervals() og giver et rodet 
	træ som argument. df_intervals() retunerer en liste med alle subtrees.
	df_intervals() kalder funktionen df_interval_subtree() og giver den
	en node som argument. df_interval_subtree() retunerer en liste med
	det subtree, som ligger under denne node.
	''' 
	l_tree1=intervals(tree1)
	l_tree2=intervals(tree2)
	x=0
	for i in l_tree1[0]:
		if i not in l_tree2[0]:
			x=x+1
	distance = x+l_tree2[1]
	return distance

def intervals (tree):
	subtrees.clear()
	liste_substrees = df_intervals(tree.clade)
	liste_intervals = []
	x=0
	for i in liste_substrees:
		if (max(i)-min(i)+1==len(i)):
			liste_intervals.append([min(i),max(i)])
		else:
			x=x+1
	return [sorted(liste_intervals),x]

def df_intervals(clade):		
		
	if clade.is_terminal():	
		i=first_name.index(clade.name)
		return new_name[i]
	else:
		l=(df_interval_subtree(clade))
		for child in clade.clades:
			df_intervals(child)
	subtrees.append(l)	
	return subtrees

def df_interval_subtree(clade):
	
	if clade.is_terminal():	
		i=first_name.index(clade.name)
		return [new_name[i]]
	else:
		l=[]
		for child in clade.clades:
			l.extend(df_interval_subtree(child)) 
	return l

''' 
''' 

if __name__ == "__main__":
	first_name=[]
	new_name=[]
	subtrees=[]

	tree1=Phylo.read(sys.argv[1],"newick")
	tree2=Phylo.read(sys.argv[2],"newick")

	#Phylo.draw_ascii(tree1)
	#Phylo.draw_ascii(tree2)

	rename_trees=rename(tree1)
	#print(rename_trees)
	interval_tree1=intervals(tree1)
	interval_tree2=intervals(tree2)
	#print(interval_tree1)
	#print(interval_tree2)
	print(RF_distance(tree1,tree2))

			
	