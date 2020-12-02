import csv
import pandas as pd

import sys

print('\n\nArgument List:', str(sys.argv))


## The individual we are trying to match
#ref_individual = 1355
ref_individual = int(sys.argv[1])

print("   Reference individual = ",ref_individual)



def read_genotypes(fname):
    print(" -> Reading csv file",fname)
    df = pd.read_csv(fname)
    df = df.drop(columns = df.columns[0])
    print(df.head())
    genotype = df.to_numpy()
    return df, genotype


def get_genotypes_for_marker(i):
    return genotype[i]

def get_individual_gt(i):
    return genotype[:, i]


def count_genotypes(marker):
    for i in gt:
        cnt[i] = cnt[i]+1
    
    print("   -> Genotype count of marker %d   n0: %d  n1: %d  n2: %d" % ( marker, cnt[0],cnt[1],cnt[2]))


##
## Read genotypes
## 

fname = "gt-1000genomes-chr2-2000markers.csv"
#fname = "chr1-2504-dense-2000genotypes.csv"

df, genotype = read_genotypes(fname)
#genotype = genotype[:,0:1000]

npersons = len(genotype[1,:])
nmarker = len(genotype)
print("\n\n\n ->   Number of markers =", nmarker," persons =", npersons)



## Retrieve the genotype of the individual we are trying to match
ref_gt = get_individual_gt(ref_individual)

print("    - Reference individual genotypes:", ref_gt)



##
## Init partitions
## 


def print_partitions(partitions):
    n_notempty = 0

    for i in partitions:
        if len(partitions[i])>0:
            n_notempty = n_notempty + 1
    
    print("       ** partitions: num not empty = " ,n_notempty, " total=",len(partitions))
    #for i in partitions:
    #    print("        id = " + str(i), ":", partitions[i])


partitions = {}
partitions['S'] = list(range(0,npersons))

print_partitions(partitions)


##
## Partition marker i
## 

   

def split_partitions(position):

    G = ref_gt[position]
    gt = get_genotypes_for_marker(position)

    print("\n   -> (at position %d) SPLITTING the partitions (ref genotype is %d) " % (position, G) )

    new_partitions = {}

    for i in partitions:
        #print("      looking at part = " + str(i), ":", partitions[i])
        
        x = partitions[i]
        p = {0: [], 1: [], 2: []}
        
        for j in x:
            #print("        at individual:",j, "gt=",gt[j])
            if gt[j] == 0:
                p[0].append(j)
            if gt[j] == 1:
                p[1].append(j)
            if gt[j] == 2:
                p[2].append(j)
        
        #print("        finished p=",p)
        
        if G==0:
            new_partitions[ i + '0' ] = p[0]
            new_partitions[ i + '1' ] = p[1]
        if G==1:
            new_partitions[ i + '0' ] = p[0]
            new_partitions[ i + '1' ] = p[1]
            new_partitions[ i + '2' ] = p[2]
        if G==2:
            new_partitions[ i + '1' ] = p[1]
            new_partitions[ i + '2' ] = p[2]

    empty_keys = []

    for i in new_partitions:
        if len(new_partitions[i]) == 0:
            empty_keys.append(i)
    
    for i in empty_keys:
        del new_partitions[i]
    
    return new_partitions




## Iterate over 500 markers
for i in range(0,500):
    partitions = split_partitions(i)
    print_partitions(partitions)

print("\n\n\n")
print(partitions)


print("\n\n\n")
nmatching = 0
matches = []
for i in partitions:
    for j in partitions[i]:
        matches.append(j)

    nmatching = nmatching + len(partitions[i])

print("Matches = ", matches)

cols = df.columns.tolist()

matches.sort()
for i in matches:
    print ( cols[ i ], "  index=",i)

print("Number of individuals that have a match with " , cols[ref_individual], " = ",  nmatching)
