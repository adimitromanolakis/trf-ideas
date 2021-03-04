
// Compile with:
// g++ -O2 fast-ibd1-set-intersection.cpp  printf-binary-format.o -Wformat=0

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "timers.h"
#include "Colors.h"

#define debug 0
#define RR __restrict


int nmark = 30000;

int nind = 3000;



// a: Genotype matrix a[individual][marker]
int ** a;

// SampleNames: IDs of individuals 
char **SampleNames;


typedef uint64_t TYPE;


#define PACKED_BITS 64

TYPE PACKED_MASK = -1;
//#define PACKED_MASK ( ((1ull)<<PACKED_BITS)-1  )

//#define PACKED_BITS 15
//#define PACKED_MASK ( (1<<PACKED_BITS)-1  )


int * ibs0;


double allele_freq = 0.15;

int debug_nops = 0;

TYPE **G;




void * malloc_align64(int size) {
    void *t;
    int e = posix_memalign(&t,64,size);

    return t;

}

template<typename TYPE>
TYPE * alloc_1d(int n1) {
    printf("alloc vector %d : %6d \n",(int)sizeof(TYPE), n1);

    TYPE *a = new TYPE[n1];
    return a;
}


template <typename TYPE>
TYPE** alloc_2d(int n1, int n2)
{
    uint64_t size = sizeof(TYPE) * (uint64_t)n1 * n2;
    size = size / 1024 / 1024;

    printf("alloc_2d: alloc array size=%.2fMB (%d): %6d * %6d\n", (double)size, sizeof(TYPE), n1, n2);

#if 0
    TYPE ** a = new TYPE*[n1];
    for(int i=0;i<n1;i++) a[i] = new TYPE[n2];
#endif

#if 1
    TYPE** a = (TYPE**)malloc_align64(sizeof(TYPE*) * n1 + 32);
    for (int i = 0; i < n1; i++)
        a[i] = (TYPE*)malloc(sizeof(TYPE) * n2 + 1024);
#endif

#if 0
    TYPE ** a = (TYPE**)malloc_align64(sizeof(TYPE*)*n1 + 32);
    TYPE *b = (TYPE*)malloc(sizeof(TYPE)*n2*n1 + 1024);
    int linesize = sizeof(TYPE)*n2;
    for(int i=0;i<n1;i++) a[i] = b + n2*i;
#endif

    return a;
}



// Helper functions


void print_matrix(int **a, int nx = 100, int ny = 30)
{
    int i,j,k;

    printf("%s -- [[ Genotype matrix ]] -- %s\n\n",Colors.bold_red,Colors.white);


    for(int i=0 ; i < nx ; i++) {
        printf("%s %5d : %s",Colors.yellow,i+1,Colors.white);
        
        for(int j=0 ; j < ny ; j++) {
            if(j % 10 == 0) printf("%s| %s",Colors.blue,Colors.white);

            printf("%d ", a[i][j] );

        }
        printf("\n");
    }

    printf("\n\n");
}



void generate_random_genotypes()
{
    srand48(23987);

    printf(" allocate genotype array: ");
    a = alloc_2d<int> (nind,nmark);


    if(1) for(int i=0;i<nind;i++) {

        for(int j=0;j<nmark;j++) {
            int x = 0;

            if( drand48() < allele_freq ) x = x + 1;
            if( drand48() < allele_freq ) x = x + 1;

            a[i][j] = x;
        }
    }

}





void pack_genotypes_v1()
{
    Timer t1;
    t1.start();


    // line size
    int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

    G = alloc_2d<TYPE> (nmark,n2+1);

    for(int i=0 ; i < nmark ; i++) {

        int pos = 0;
        TYPE v = 0;

        for (int k = 0; k < n2-1; k++) {

            v = 0;

            for (int z = 0; z < PACKED_BITS; z++) {
                int g = a[pos][i];

                g = (g == 0) ? 0 : 1;

                v = (v | g) << 1;
                pos++;
            }

            //printf("STORE (%4d , %4d): %x %x\n",  i, k ,(uint32_t)(v>>32), (uint32_t)(v) );
            G[i][k] = v;

        }
    }

    printf("pack_genotypes V1: time %s%.2f%s\n", Colors.blue, t1.duration() , Colors.white);
}





class TableReader
{
    public:

    int num_individuals;
    char buf[1024 * 1024 + 1];
    FILE *f;

    TableReader() {}

    void get_num_individuals()
    {
        char *p;
        int pos = 0;
        p = fgets(buf, sizeof(buf), f);

        char *strtok_data;

        p = strtok_r(buf, " ", &strtok_data);

        while (p != NULL)
        {
            pos++;
            p = strtok_r(NULL, " ", &strtok_data);
        }

        printf("First line : read %d tokens\n", pos);

        num_individuals = pos;
        fseek(f, 0L, SEEK_SET);

        SampleNames = alloc_2d<char>(num_individuals, 500);

    }

    void read_file(char *fname, int **dest)
    {

        f = fopen(fname, "r");

        printf("Reading file with genotypes : %s\n", fname);

        if (f == NULL)
        {
            printf("cannot open file %s\n", fname);
            exit(1);
        }

        get_num_individuals();

        char *p;
        int nline = 0;
        int pos = 0;

        // Read sample names
        p = fgets(buf, sizeof(buf), f);
   

        char *strtok_data;

        p = strtok_r(buf, " ", &strtok_data);
        pos = 0;

        while (p != NULL)
        {
            //printf(" token: \"%s\" N=%d\n",p ,(int)atoi(p));

            strncpy(SampleNames[pos], p, 500 - 1);
            pos++;

            p = strtok_r(NULL, " ", &strtok_data);
        }

        // Read genotype data
        while (1)
        {
            p = fgets(buf, sizeof(buf), f);

            if (p == NULL)
                break;

            char *strtok_data;

            char *p = strtok_r(buf, " ", &strtok_data);
            pos = 0;

            while (p != NULL)
            {
                //printf(" token: \"%s\" N=%d\n",p ,(int)atoi(p));

                a[pos][nline] = atoi(p);
                pos++;

                p = strtok_r(NULL, " ", &strtok_data);
            }

            //printf(" line %d read %d tokens\n", nline, pos);
            nline++;
        }

        printf("Read %d markers with %d individuals\n", nline, pos );
        nind = pos;
        nmark = nline;

        fclose(f);
    }
};


void print_error_matrix(int dual = 0)
{
    int e[102][102];
    int N = nind;

    if(N>15) N = 15;

    for (int i1 = 0; i1 < N; i1++) {

        for (int i2 = 0; i2 < N; i2++) {
            int err = 0;

            for (int j = 0; j < nmark; j++) {
                int g1 = a[i1][j];
                int g2 = a[i2][j];

                if(dual) if(g1==0 && g2==2) err++;
                if(g1==2 && g2==0) err++;
            }

            e[i1][i2] = err;
                
        }
    }

    printf("%s -- [[ Error matrix ]] -- %s\n\n",Colors.bold_red,Colors.white);

    for(int i=0 ; i < N; i++) {
        printf("%s%4d - %s",Colors.blue, i+1 , Colors.white );
        int num_0s = 0;

        for(int j=0 ; j < N ; j++) {
            //if(e[i][j]==0) printf("%s",Colors.bold_red);
            printf("%2d ", e[i][j] );
            if(e[i][j] == 0) num_0s ++;
            //if(e[i][j]==0) printf("%s",Colors.white);
        }

        printf("   num_0s=%d", num_0s);
        printf("\n");
    }

    printf("\n\n");

}



// -----------------------------------------------
//            Main algorithm below
// -----------------------------------------------



int set_sizes = 0;


class CollectionOfPacked12s {

public:

    TYPE **packed12s;

    CollectionOfPacked12s() {
    }

    void compute_compatible_sets_12()
    {
        //printf("compute_compatible_sets_12: BITS=%d MASK=%x\n", PACKED_BITS, PACKED_MASK);
        packed12s = alloc_2d<TYPE>(nmark, nind/PACKED_BITS + 2);


        set_sizes = (nind+PACKED_BITS-1)/PACKED_BITS;



        int p0 = 0;

        for (int j = 0; j < nmark; j++)
        {
            int numpacked = 0, p0 = 0;
            TYPE v = 0;

            for (int i = 0; i < nind; i++)
            {

                int x1 = 0;

                if (a[i][j] != 0)
                    v = v | 1;

                numpacked++;
                if (numpacked == PACKED_BITS)
                {
                    //printf("long int: %llx\n", v); // long long decimal

                    packed12s[j][p0] = v;
                    //printf("fill [%d,%d] = %0x %0x  -- " , j,p0, (uint32_t)(v>>32),(uint32_t)(v & 0xffffffff) );
                    //printBits(4, (void*)&packed12s[j][p0]);

                    v = 0;
                    numpacked = 0;
                    p0++;
                }

                v = v << 1;
            }

            v = v << (PACKED_BITS-numpacked-1);
            packed12s[j][p0] = v;

            //printf("fill [%d,%d] = %0x %0x  -- " , j,p0, (uint32_t)(packed12s[j][p0]>>32),(uint32_t)(packed12s[j][p0] & 0xffffffff) );
            //printBits(4, (void*)&packed12s[j][p0]);
            //printf("P0=%d round=%f %d set_sizes=%d overflow numpacked=%d\n",p0, (double)nind/PACKED_BITS, (int)(nind/PACKED_BITS) , set_sizes, numpacked);
            //printf("  overflow numpacked=%d\n",numpacked);
        }
    }

    void counts()
    {
        unsigned int t = 0x808080;

        for (int j = 0; j < nmark; j++) {

            printf("Marker %d : ", j);

            int n2 = nind / PACKED_BITS;
            //n2 = n2+1;

            int cnt = 0;
            for (int i = 0; i < n2; i++)
            {
                cnt = cnt + __builtin_popcountll(packed12s[j][i]);
            }

            printf(" num 1s in set = %d\n", cnt);
        }
    }

};



CollectionOfPacked12s packed12s;


int set_id = 0;

int total_intersections = 0;


class Set {
    public:
    TYPE *S;

    int id;

    Set() {
        void *t;

        int e = posix_memalign(&t,64,8024*sizeof(S[0]) );
        S = (TYPE*) t;

        // printf("create set id %d\n",set_id);
        id = set_id;
        set_id++;
    }

    ~Set() {
        // printf("free set id %d  : %x\n",id, S);
        free(S);
    }

    void copy_from(CollectionOfPacked12s * __restrict K, int marker) 
    {
        TYPE * __restrict p = K->packed12s[marker];
        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

        //("COPY_FROM: %d\n", marker);

        for (int i = 0; i < n2; i++)
        {
            S[i] = p[i];
            //S2[i] = 0;
        }
        
    }


    void intersect(CollectionOfPacked12s * __restrict K, int marker) 
    {
        TYPE * __restrict p = K->packed12s[marker];
        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

        for (int i = 0; i < n2; i++)
        {
   
            S[i] &= p[i];
        }
    }



    void copy_from(Set *s2) 
    {
        TYPE * __restrict p = s2->S;
        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

        //("COPY_FROM: %d\n", marker);

        for (int i = 0; i < n2; i++)
        {
            S[i] = p[i];
        }
    }


    void intersect(Set *s2) 
    {
        TYPE * __restrict p = s2->S;
        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

        for (int i = 0; i < n2; i++)
        {
           
            S[i] &= p[i];
        }
    }  


    void intersect(Set *s1, Set *s2) 
    {
        TYPE * __restrict p1 = s1->S;
        TYPE * __restrict p2 = s2->S;

        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;
        //printf("intersect1 %d %d\n",s1->id,s2->id);

        for (int i = 0; i < n2; i++)
        {
            S[i] = p1[i] & p2[i];
        }
    }      

  void intersect2(CollectionOfPacked12s * __restrict K, int marker) 
    {
        TYPE * __restrict p = K->packed12s[marker];
        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

        for (int i = 0; i < n2; i++)
        {
            if(S[i] != 0)
                S[i] &= p[i];
        }
    }


    void update_start_positions(int *match_start_pos, int genome_pos)
    {
        //printf("(update start positions)\n");

        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

        int pos = 0;
        for (int i = 0; i < n2; i++)
        {
            TYPE k = S[i];
            TYPE mask = 1ULL << (PACKED_BITS-1);

            //if(k!=0)
            for(int j=0;j<PACKED_BITS;j++) {

                if( (k & mask) != 0) {

                    // Set at position "pos" is 1
                    if(match_start_pos[pos] == -1) {
                        match_start_pos[pos] = genome_pos;
                    }

                } else {
                    // Set is 0
                    if(match_start_pos[pos] == -1) {
                        //match_start_pos[pos] = genome_pos;
                    } else {
                        //printf("match POS=%d IND=%d\n", genome_pos, pos);
                        match_start_pos[pos] = -1;
                    }

                }

                mask >>= 1;

                pos ++;
            }
        }
    }

    void print_matches(int *match_start_pos, int individual_id, int genome_pos)
    {
        //printf("(print_matches)\n");

        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

        int pos = 0;
        for (int i = 0; i < n2; i++)
        {
            TYPE k = S[i];
            TYPE mask = 1ULL << (PACKED_BITS-1);

            //if(k!=0)
            for(int j=0;j<PACKED_BITS;j++) {

                if( (k & mask) != 0) {

                    // Set at position "pos" is 1
                    
                } else {
                    // Set is 0
                    if(match_start_pos[pos] == -1) {
                        //match_start_pos[pos] = genome_pos;
                    } else {
                        int left_location = match_start_pos[pos];
                        
                        if(genome_pos - left_location > 800)
                        if(individual_id != pos)
                            printf("MATCH %s %s %d %d\n", SampleNames[individual_id], SampleNames[pos] , match_start_pos[pos], genome_pos );

                        //match_start_pos[pos] = -1;
                    }

                }

                mask >>= 1;

                pos ++;
            }
        }
    }


    void print_1s()
    {
        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

        int pos = 0;
        for (int i = 0; i < n2; i++)
        {
            TYPE k = S[i];
            TYPE mask = 1ULL << (PACKED_BITS-1);
           // printf("i=%d mask=%x S=%x %x\n",i, mask>>32,k>>32,k);

            for(int j=0;j<PACKED_BITS;j++) {
                //printf("pos=%d mask=%x i=%d has1=%d\n", pos,mask,i, (int)((k&mask)!=0) );
                if( (k & mask) != 0) {
                    printf("%d ", pos+1);
                }
                mask >>= 1;

                pos ++;
            }
        }

        printf("\n");

    }

    void print_counts()
    {
        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

        int pos = 0;

        printf("    %sC= %s",Colors.blue,Colors.yellow);
        for (int i = 0; i < n2; i++)
        {
            TYPE k = S[i];
            int cnt = __builtin_popcountll(k);
            printf("%2d ", cnt);
        }
        
        printf("%s\n",Colors.white);

    }

    
    void print_brief()
    {
       int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

        int pos = 0;
        for (int i = 0; i < n2; i++)
        {
            TYPE k = S[i];
            TYPE mask = 1ULL << (PACKED_BITS-1);

            for(int j=0;j<PACKED_BITS;j++) {
                //printf("pos=%d mask=%x i=%d has1=%d\n", pos,mask,i, (int)((k&mask)!=0) );
                
                char c = '0';
                if( (k & mask) != 0) {
                    c = '1';
                }

                printf("%c",c);

                mask >>= 1;
                pos ++;

                if(pos>=100) break;
            }
        }

        printf("\n");

    }

    int count()
    {
        int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;
        //n2 = n2+1;

        int cnt = 0;
        for (int i = 0; i < n2; i++)
        {
            TYPE x = S[i];

            if(x)
                cnt = cnt + __builtin_popcountll(x);
        }

        return cnt;
    }

};



Set s2;

int find_matches_single_window(int individual, int start_marker, int end_marker) {

	// Get genotypes of individual and find the position of the 2s
    int * P = a[individual];

    //printf("%s Process individual %s %d %s\n", Colors.green,Colors.blue, individual,Colors.white);

    int j;
    int pos = 0;

    for (j = start_marker; j < end_marker; j++) {
        if (P[j] == 2) {
            
            s2.copy_from(&packed12s, j);

            total_intersections++;

            pos++;
            j++;
            break;
        }
    }

    for (; j < end_marker; j++) {
        if (P[j] == 2) {
            
            s2.intersect(&packed12s, j); 

            total_intersections++;
            pos++;
            //if(pos == 5 && s2.count() < 5 ) break;
        }
    }

    int total_compatible = s2.count();

    //printf("Matches for %d: ", individual+1);
    //s2.print_1s();

    //printf("num compatibles = %d\n", s2.count() );

    return total_compatible;
}

int match_start_positions[50000];

int scanner_v1(int individual_id) 
{
    // Get the genotypes of the individual
    int *P = a[individual_id];

    //printf("%s Process individual %s %d %s\n", Colors.green,Colors.blue, individual,Colors.white);

    int j;

    // Create 2 sets:
    // S1: markers 0 to 199
    // S2: markers 200 to 400

    int first = 1;

    Set matches;


    for(int i=0;i<nind+128;i++) match_start_positions[i] = -1;
    
    int pos = 0;
    int end = pos + 200;

    //printf("1 Scan from %d-%d\n", pos, end);

    
    // find the compatible individuals for positions pos..pos+200
    for (j = pos; j < end; j++) {
        if (P[j] == 2) {
            if(first==1)
                {
                    //printf("%d copy from %d\n", individual_id, j);
                    matches.copy_from(&packed12s, j); total_intersections++;
                    first = 0;
                }
            else 
                {
                    //printf("%d intersect from %d\n", individual_id, j);
                    matches.intersect(&packed12s, j); total_intersections++;
                }
            }
    }


    matches.update_start_positions(match_start_positions, 0);

    //current_set->print_brief();

    //printf("Match start positions:\n");
    //for(int i=0;i<100;i++) { printf("%c", (match_start_positions[i]==-1)?'.':'M'  ); }
    //printf("\n");

    for(int i=0;i<nind;i++) {
     //printf("MS(%d) = %d\n", i, match_start_positions[i]);
    }

    pos += 200;

    while(1) {
    
        end = pos + 200;
        //printf("2 Scan from %d-%d\n", pos, end);
        
        first = 1;

        // find the compatible individuals for positions pos..pos+200
        for (j = pos; j < end; j++) {
            if (P[j] == 2) {
                if(first==1)
                    {
                        //printf("%d copy from %d\n", individual_id, j);
                        matches.copy_from(&packed12s, j);total_intersections++;
                        first = 0;
                    }
                else 
                    {
                        //printf("%d intersect from %d\n", individual_id, j);
                        matches.intersect(&packed12s, j);total_intersections++;
                    }
                }
        }

        //matches.print_matches(match_start_positions, individual_id, pos);

        matches.update_start_positions(match_start_positions, pos);

        pos += 200;
        if(pos > nmark-200 ) break;
    }

    for(int i=0;i<nind;i++) {
        if(match_start_positions[i] != -1) {
            if(end - match_start_positions[i] >= 800)
            if(individual_id != i)
            printf("Match at end: %s %s %d %d\n", SampleNames[individual_id], SampleNames[i], match_start_positions[i], end);
        }        
    }

    int total_compatible = s2.count();

    return total_compatible;
}

void print_matches(int individual_id, Set *s)
{
    int n2 = (nind+PACKED_BITS-1) / PACKED_BITS;

    int pos = 0;
    for (int i = 0; i < n2; i++)
    {
        TYPE k = s->S[i];

        TYPE mask = 1ULL << (PACKED_BITS-1);
        // printf("i=%d mask=%x S=%x %x\n",i, mask>>32,k>>32,k);

        if(k==0) { pos += PACKED_BITS; continue; }

        for(int j=0;j<PACKED_BITS;j++) {
            // printf("pos=%d mask=%x i=%d has1=%d\n", pos,mask>>32,i, (int)((k&mask)!=0) );
            if( (k & mask) != 0) {
                
                //if(pos != individual_id)  printf("SEG %s %s\n", SampleNames[individual_id], SampleNames[pos]);
                if(pos != individual_id)  printf("SEG %s %s %d %d\n", SampleNames[individual_id], SampleNames[pos] , individual_id, pos ) ;

            }
            mask >>= 1;

            pos ++;
        }
    }

    //printf("\n");
}



int 
scan_matches()
{
    int total = 0;

    printf("nmark=%d\n", nmark);

    for (int i = 0; i < nind; i++) {
            total += scanner_v1(i);
            //total += find_matches_single_window(i,0,nmark);
        }

    return total;
}


void test_method_intersections()
{
    printf("\n");
    printf("%s   Intersections method   %s\n", Colors.reverse_red, Colors.white);
    printf("\n");
 
    // For debugging we can print the genotype matrix and the computed error matrix
    //print_matrix();
    //print_error_matrix();

    Timer t2;
    t2.start();

    packed12s.compute_compatible_sets_12();
    printf("sets were filled time=%s %.2f %s\n",Colors.blue, t2.duration() , Colors.white);

    printf("\n\n");

    int print_sets = 0;

    if(print_sets) {
        printf("set packes12s:\n");

        for(int i=0;i<10;i++) {
            for(int j=0;j<10;j++) {
                printf("%9x ", (uint32_t)(packed12s.packed12s[i][j]>>32) );
            }
            printf("\n");

        }
        printf("\n\n");

        printf("set G:\n");

        for(int i=0;i<10;i++) {
            for(int j=0;j<10;j++) {
                printf("%9x ", (uint32_t)(G[i][j]>>32) );
            }
            printf("\n");

        }

    }

    printf("\n\n");

    for (int rep = 0; rep < 5; rep++) {
        Timer t1;

        int total = 0;

        total_intersections = 0;
       
    
        total += scan_matches();

        printf("total=%d\n", total);
        printf("total intersections=%d\n",total_intersections);
        printf("time = %s%f%s\n", Colors.blue, t1.duration(), Colors.white);
    }
}



int
main()
{
    printf("%s   Algorithms for speeding up TRUFFLE   %s\n", Colors.reverse_gray, Colors.white);
    printf("\n");
 	generate_random_genotypes();
    
    // A file written from R as : write.table(gt[1:400,1:1000], row=F,col=F,quote=F,"/tmp/genotypes.csv")
    // Columns are individuals.

    TableReader R;
	R.read_file("genotypes.csv", a);

    test_method_intersections();
}

