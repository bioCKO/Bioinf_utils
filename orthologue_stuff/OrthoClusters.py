import sys

infile = open(sys.argv[1])
outfile = open(sys.argv[2], 'w')

args = len(sys.argv)
if args > 3:
    mat = True
    outmat = open(sys.argv[3], 'w')

def get_species(str):
    if str.startswith('Y'):
        return 'SC'
    else:
        return str[0:2]

def search_clusters(clusters, gene):
    for cluster in sorted(clusters):
        if gene in cluster:
            return clusters.index(cluster)
            
def clusteradd(clusters, index, newgene):
    return clusters[index].add(newgene)

def trim(string):
    l = len(string)
    if string[l-2:l] == 'T0':
        return string[0:l-2]
    else:
        return string
    
clusters = [set()]

species = set()

for line in infile:
    cols = line.split()
    srcGene = trim(cols[0])
    ortho = trim(cols[2])
    species.add(get_species(ortho))
    species.add(get_species(srcGene))
    if search_clusters(clusters, srcGene) > 0:
        clusteradd(clusters, search_clusters(clusters, srcGene), ortho)
    elif search_clusters(clusters, ortho) > 0:
        clusteradd(clusters, search_clusters(clusters, ortho), srcGene)
    else:
        newcluster = set()
        newcluster.add(srcGene)
        newcluster.add(ortho)
        clusters.append(newcluster)

print 'Species:'
print species

clusterdic = {}
for cluster in clusters:
    if len(cluster) in clusterdic:
        clusterdic[len(cluster)] += 1
    else:
        clusterdic[len(cluster)] = 1

print 'cluster size distribution:'
print clusterdic
print 'total number of clusters:'
print sum(clusterdic.values())

    
def cluster_to_spdic(cluster):
    cludic = {}
    for gene in cluster:
        sp  = get_species(gene)
        if sp in cludic:
            cludic[sp].append(gene)
        else:
            cludic[sp] = [gene]
    return cludic



def print_by_sp(species, spdic):
    output = []
    for sp in sorted(species):
        if sp in spdic:
            output.append(', '.join(spdic[sp]))
        else:
            output.append('-')
    return '\t'.join(output) + '\n'

def print_matrix(species, spdic):
    output = []
    for sp in sorted(species):
        if sp in spdic:
            output.append(str(len(spdic[sp])))
        else:
            output.append('0')           
    return '\t'.join(output) + '\n'

            
clusters[0] = species


for cluster in clusters:
    outfile.write(print_by_sp(species, cluster_to_spdic(cluster)))

if mat:
    outmat.write(print_by_sp(species, cluster_to_spdic(clusters[0])))
    for cluster in clusters[1:len(clusters)]:
        outmat.write(print_matrix(species, cluster_to_spdic(cluster)))

    
outfile.close()
if mat:
    outmat.close()
