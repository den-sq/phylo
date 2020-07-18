/*
int is_prime(const int n);
struct Clade;
typedef struct Clade Clade;
Clade* re_export_rt(char *file_name, char *file_type);
*/
struct PAMLRecord;
struct Alignment;
struct Clade;
struct Per_Site_Annotation;
typedef struct PAMLRecord PAMLRecord;
typedef struct Alignment Alignment;
typedef struct Clade Clade;
typedef struct Per_Site_Annotation Per_Site_Annotation;

typedef struct {
    double mean;
    double median;
    double max;
} Ranges;

typedef struct {
    int x;
    int y;
    int z;
} QuickPair;

typedef struct {
    int length;
    int capacity;
} ArrayMeta;

typedef struct {
	char* name;
	int id;
	int ancestor;
	double dNdS;
} PAMLSelection;

// Basic Error Handling
int last_error_length();
int last_error_message(char *buffer, int length);

// Read and Process PAML Data
void process_paml_file(char* file_name);
PAMLSelection* get_paml_selection(char* file_name, ArrayMeta* meta);
void free_paml_selection(PAMLSelection* ptr, ArrayMeta* meta);
Alignment* get_paml_alignment(char* filename, ArrayMeta* meta);
void free_paml_alignment(Alignment* ptr, ArrayMeta* meta);
void free_meta(ArrayMeta* meta);

// Detailed Site Data Calculation & Data Fetching
Per_Site_Annotation* calculate_site_data(Alignment* sequences, char* pdb);
long* get_secondary_structure_distribution(Per_Site_Annotation* annotation);
long* get_regional_distribution(Per_Site_Annotation* annotation);
long* get_unique_sites(Per_Site_Annotation* annotation);

// Read Tree & Get Species Data
Clade* read_tree(char* file_name, char* file_type);
char* get_arctic(Clade* tree, char* name_list);
char* get_non_arctic(Clade* tree, char* name_list);
char** get_nonarctic_species(Clade* tree, int size);

// Tests
int quick_test(int test_input);
int* int_array_test(int test_input);
long* long_array_test(long test_input);
char** char_array_test(int test_input);
void preallocate_array_test(int value, int* preallocate, ArrayMeta* meta);

QuickPair* pair_test(int input_one, int input_two, int input_three);
QuickPair* pair_array_test(int input_one, int input_two, int input_three, int input_four, ArrayMeta* meta);
void free_pair_array(QuickPair* ptr, ArrayMeta* meta);

void print_pair(QuickPair* ptr);
void print_pair_array(QuickPair* ptr, ArrayMeta* meta);

ArrayMeta* array_meta();
void bootstrap(char* path, char* counts, char* related_counts, char* data, char* matrix, char* lookup);
void basic_analysis(char* path, char* counts, char* data, char* matrix);
void related_basic_analysis(char* path, char* counts, char* data, char* matrix, char* lookup);

// ??
char** get_genes(Clade* trees, int len, int size, int* gene_len);

Ranges* get_protein_identity(char* alignment_path, char* protein_path);
Ranges* get_identity(char* alignment_path, char* alignment_format);
void process(char* db_path, char* db_data);
