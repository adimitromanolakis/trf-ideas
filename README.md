# Ideas for speeding up truffle

## C++ Set Intersections

How to run:

Make sure genotypes.csv exists in the current folder,. run make and run ./fast-ibd1-set-intersection-scanner.

Overview:

class Set : a bit vector ( 1 bit per individual). In the algorithms the ones are compatible people (matches) with the current one.

CollectionOfPacked12s : A list of m sets S_i (m: number of markers) sets with 1 bit per individual. For each marker position and individual, a bit is 1 if they individual has a genotype of 1 or 2. That means that set S_i contains the individuals that are compatible with a 2 at position i.

Functions:
find_matches_single_window: find matches for a single window and a single individual.
scanner_v1: first implementation of the whole chromosome scanner.

To switch between the implementations, change function scan_matches.



## Python partitions test:

Run python3 partitions.py 3

where 3 is the index of the invididual which is the query.
