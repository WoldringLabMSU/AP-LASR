#These are unused functions, in storage. If you find a repuropsed use for them, be my guest. Significat portions are written by Pattrick Finneran.
# Written by James VanAntwerp in September 2021 - vanantj@udel.edu
# Written for the Woldring Lab, Michigan State University in East Lansing, Michigan, USA.

def Select_Ancestor_Nodes_OLD(dirname,finname,cutoff=85):
    #This will evaluate the tree, which is in Newick format, and select nodes of high enough confidence to reconstruct an ancestral library.
    #The .treefile has all of the confidence values, and the node names.
    with open(f'{dirname}/{finname}') as treefin:
        treefile = treefin.readlines()[0]
    Good_Nodes=[] #these will be the nodes with poor values. This isn't quite the same value between *.contree and *.trefile, but close
    nodes=treefile.split(')')#Let's split off just the information for each node, which is stored after every close parenthiesis.
    nodes.pop(0)#The above will split off a first section without a node. This does not cause any nodes to be lost.
    for i,node in enumerate(nodes): #rearanging the newick format a little. Each node in nodes will now be stored as "name,int,float"
        try:
            node_info=node.split(',')[0]
            #nodes[i]=f"{node_info.split('/')[0]},{(node_info.split('/')[1]).split(':')[0]},{node_info.split(':')[1]}"
            nodes[i]=node_info.split('/')[0]
            if int((node_info.split('/')[1]).split(':')[0]) >= cutoff: #Select the nodes which have a value above the cutoff.
                #tup = ((node_info.split('/')[0]),(node_info.split('/')[1]).split(':')[0])
                Good_Nodes.append(nodes[i])
        except:
            if node.strip() == 'Node1;':
                pass #The root node doesn't have info - this prevents an error in handling that from causing problems
            else:
                print(f'{node} caused problems :(')
                print("There was an error handling the Newick Tree returned by IQTree - This likely indicates a bug with IQTree.")
                raise ValueError("There was an error handling the Newick Tree returned by IQTree - This likely indicates a bug with IQTree.")
    #return [n.split(',')[0] for n in nodes if n.split(',')[0] not in stinky_nodes]
    return Good_Nodes #return the list of node names whose value is above the cutoff.

def Trim_N_C_Termini_Percent (fasta_dict,terminus_cutoff=0.20): #Trim termini based on % of alignment that has a sequence. Currently not used.
    sequences_list = [i for i in fasta_dict.values()]
    gap_percentage_by_pos=[]#Represents the % of each position that is a gap.
    for position in range(len(sequences_list[0])): #For every position in the sequences,
        pos_gap_tally=0
        for i in range(len(fasta_dict)): #For every sequence
            if sequences_list[i][position] =='-': #If that position has a gap
                pos_gap_tally+=1
        gap_percentage_by_pos.append(pos_gap_tally/len(sequences_list))
    #Count N-terminus trim length
    len_to_remove_N=0
    len_to_remove_C=0
    for pos in gap_percentage_by_pos:
        if pos > terminus_cutoff:
            len_to_remove_N+=1
        else:
            break
    #Count C-terminus trim length 
    for i in range(len(gap_percentage_by_pos)-1,-1,-1):
        if gap_percentage_by_pos[i]>terminus_cutoff:
            len_to_remove_C+=1
        else:
            break
    for name,seq in fasta_dict.items():
        fasta_dict.update({name : (seq[len_to_remove_N:-(len_to_remove_C)]) })

def Statefile_to_Dict_AAs(dirname,finname): #{NodeX:[AA,AA,AA,...]}
    statefile_dict={} #This is a dicitonary made out of the statefile - its keys are the node names, 
    #and its values are a list of tuples with the amino acid and list of amino acid distributions at each position.
    node=[]
    #This will parse the .state file for the desired nodes to get their AA distribution
    with open(f'{dirname}/{finname}') as statefin: #Read in each line, skipping the header.
        statelines=[]
        for i,line in enumerate(statefin):
            if i>8:
                statelines.append(line)
    #Now let's pull the data from each line
    working_node=statelines[0].split()[0] #prime the working node
    for line in statelines: # For every line in the state file
        line_list = line.split() #Break up the stuff
        if working_node == line_list[0]: #If we're still working on the same node
            node.append(line_list[2]) #record the amino acid at that position
        else: #If we've come to the end of a node,
            statefile_dict[working_node]=node #Add a key-value pair to the statefile dictionary that is the node name and the node's list
            working_node = line_list[0] #update the working_node value
            node=[] #Clear the working node list
            node.append((line_list[2],line_list[3:])) #add to the node list a touple of the amino acid and the distribution of amino acids at that position.
    statefile_dict[working_node]=node #Be sure to add the last node into the dictionary too!
    return(statefile_dict)
