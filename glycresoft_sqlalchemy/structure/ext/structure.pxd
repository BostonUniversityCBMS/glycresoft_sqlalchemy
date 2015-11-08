cdef struct AminoAcidStruct:
    char* name
    char* symbol
    double mass


cdef struct AminoAcidStructArray:
    AminoAcidStruct* residues
    size_t size


cdef struct ModificationStruct:
    char* name
    double mass


cdef struct ModificationStructArray:
    size_t size
    ModificationStruct* modifications


cdef struct SequencePositionStruct:
    AminoAcidStruct* residue
    ModificationStructArray* modifications
    double mass


cdef struct SequenceStruct:
    size_t size
    SequencePositionStruct* sequence
    double mass
    int n_glycan_sequon_sites
    SimpleGlycanStruct* glycan


cdef struct SimpleGlycanStruct:
    double mass
    char* name


cdef struct FragmentIonStruct:
    char* series
    double mass
    int index
    AminoAcidStruct* n_term
    AminoAcidStruct* c_term


cdef struct FragmentIonStructArray:
    FragmentIonStruct* fragments
    size_t size


cdef AminoAcidStructArray* AminoAcids
cdef ModificationStructArray* Modifications


cdef AminoAcidStruct* get_residue_by_symbol(char* symbol, AminoAcidStructArray* residues) nogil
cdef AminoAcidStruct* get_residue_by_name(char* name, AminoAcidStructArray* residues) nogil
cdef inline AminoAcidStruct* _get_residue_by_symbol(char* symbol) nogil
cdef inline AminoAcidStruct* _get_residue_by_name(char* name) nogil

cdef ModificationStruct* get_modification(char* name, ModificationStructArray* modifications) nogil
cdef inline ModificationStruct* _get_modification(char* name) nogil

cdef double sequence_position_mass(SequencePositionStruct* position) nogil
cdef char* sequence_to_string(SequenceStruct* seq) nogil
cdef FragmentIonStructArray* fragment_series(SequenceStruct* sequence, char* series, int direction, double offset) nogil

cdef AminoAcidStructArray* create_amino_acid_table(dict residue_map)
cdef ModificationStructArray* create_modification_table(dict modification_map)

cdef SequenceStruct* sequence_from_object(object obj)
cdef SimpleGlycanStruct* glycan_from_object(object obj)


cdef void free_sequence(SequenceStruct* seq) nogil
cdef void free_sequence_position(SequencePositionStruct* pos) nogil
cdef void free_modification_array(ModificationStructArray* mod_array) nogil
cdef void free_amino_acid_array(AminoAcidStructArray* aa_array) nogil
cdef void free_fragment_ion_array(FragmentIonStructArray* frag_array) nogil
cdef void free_simple_glycan(SimpleGlycanStruct* g) nogil
