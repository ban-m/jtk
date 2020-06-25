import json

class DataSet:
    """
    A dataset of HLA decomposer.
    
    Attributes 
    --------------------
    raw_reads: list of `RawRead`
    selected_chunks: list of `Unit`
    encoded_reads: list of `EncodedRead`
    assginments: List of `Assignments`
    hic_pairs: list of `HiCPair`
    hic_edges: list of `HiCEdge`
    """
    def __init__(self, raw_reads, selected_chunks, encoded_reads, assignments,
                 hic_pairs, hic_edges):
        """
        The initialization method.
        Usually, users should not call this method directlly.
        Rather, it is recommended to parse this type from a JSON file.
        See `as_dataset` for more details.
        """
        self.raw_reads = raw_reads
        self.selected_chunks = selected_chunks
        self.encoded_reads = encoded_reads
        self.assignments = assignments
        self.hic_pairs = hic_pairs
        self.hic_edges = hic_edges

class RawRead:
    """
    A raw reads. In other words, it is a fasta record. Period.
    
    Attributes
    ------------------
    name: str
    desc: str
        A description of the original fasta.
    id: int
        Entry program(`jtk entry`) gives a unique identification number to each of reads. It is recorded as self.id.
    seq: str
    """
    def __init__(self, name, desc, identification, seq):
        """
        Initialization method. Users should not call this method directlly.
        """
        self.name = name
        self.desc = desc
        self.id = identification
        self.seq = seq

class Unit:
    """
    A unit sequence. It is sometimes called a `segment`, `section`, or `chunk`.
    
    Attributes
    ----------------
    id: int
        Unit encoder should give a unique identification number of each of units. It is recorded as self.id
    seq: str
    """
    def __init__(self, identification, seq):
        """
        Initialization method. Users should not call this method directlly.
        """
        self.id = identification
        self.seq = seq

class EncodedRead:
    """
    A encoded read. Note that each `EncodedRead` should have one and only one 
    `RawRead` having the same identification number.
    Attributes
    ------------------
    id: int
        The id number of the corresponding raw reads.
    original_length: int
        The length of the corresponding raw reads.
    leading_gap:int
        The length of un-encoded segment at the beggining in the read.
    trailing_gap:int 
        The length of un-encoded segment at the end of the read.
    nodes: list of `Node`
    edges: list of `Edge`
    """
    def __init__(self, idn, original_length, leading_gap, trailing_gap,
                 nodes, edges):
        self.id = idn;
        self.original_length = original_length
        self.leading_gap = leading_gap
        self.trailing_gap = trailing_gap
        self.nodes = nodes
        self.edges = edges

class Node:
    """
    A node in encoded reads.
    
    Attributes
    -----------------
    position_from_start: int
        The start position where the encoding started.
    unit: int
        The ID number of the unit which encods the part of the read.
    cluster: int
        The (local) cluster assigned to this node. Note that the default value is 0.
    seq: str
        The raw sequence of the **read** (not unit sequence, which can be 
        obtained by checking unit ID number). 
        Note that the seq is rev-cmped if needed.
    is_forward: bool
    cigar: list of :obj:
        Each :obj: is one of the following:
        {"Match":int}, {"Del": int}, or {"Ins":int}
    """
    def __init__(self, position_from_start, unit, cluster, seq, is_forward, cigar):
        self.position_from_start = position_from_start
        self.cluster = cluster
        self.unit = unit
        self.seq = seq
        self.is_forward = is_forward
        self.cigar = cigar

class Edge:
    """
    An edge between two nodes.
    In attributes, let (`u`, `v`) be the edge we represent.
    
    Attributes
    ----------------------
    frm: int,
        The ID number of the unit `u`
    to: int,
        The ID number of the unit `v`
    offset: int
        The distance between `u` and `v` in **the raw read**. 
        Thus, if the two units overlap, the offset would be negative.
    label: str
        The sequence between the end position of `u` and the start position of `v`.
        If the two units overlap in the raw read, label would be emtpy string.
    """
    def __init__(self, frm, to, offset, label):
        self.frm = frm
        self.to = to
        self.offset = offset
        self.label = label


class HiCPair:
    """
    A raw HiC pair. It can be changed in the future.
    
    Attributes 
    --------------
    pair1: int
        The ID number of pair1.
    pair2: int
        The ID number of pair2.
    pair_id: int
    seq1: str
    seq2: str
    """
    def __init__(self, pair1, pair2, pair_id, seq1, seq2):
        self.pair1 = pair1
        self.pair2 = pair2
        self.pair_id = pair_id
        self.seq1 = seq1
        self.seq2 = seq2

class HiCEdge:
    """
    A pair of HiC reads. It can be changed in the future.
    
    Attributes
    --------------
    pair_id: int
    pair1: int
        The ID number of the unit which one of the hi_c pair hits.
    pair2: int
    """
    def __init__(self, pair_id, pair1, pair2):
        self.pair_id = pair_id
        self.pair1 = pair1
        self.pair2 = pair2

class Assignments:
    """
    An assignment of raw read into cluster.
    
    Attributes
    --------------
    id: int
        The ID number of the raw read.
    cluster: int
    """
    def __init__(self, idn, cluster):
        self.id = idn
        self.cluster = cluster

def as_dataset(dct):
    """
    A decorder function of the JSON file encoding a DataSet.
    
    Parameters
    ---------------
    dct: Dictonary
    
    Examples 
    ------------------
    >>> import json
    >>> f = open("dataset.json", 'r')
    >>> dataset = json.load(f, object_hook = as_dataset)
    """
    if 'raw_reads' in dct and \
       'selected_chunks'in dct and\
       'encoded_reads'in dct and\
       'assignments' in dct and\
       'hic_pairs' in dct and \
       'hic_edges' in dct:
        return DataSet(raw_reads = dct['raw_reads'],
                       selected_chunks = dct['selected_chunks'],
                       encoded_reads = dct['encoded_reads'],
                       assignments = dct['assignments'],
                       hic_pairs = dct['hic_pairs'],
                       hic_edges = dct['hic_edges'])
    elif 'name' in dct and \
         'desc' in dct and \
         'id' in dct and \
         'seq'in dct :
        return RawRead(name = dct['name'],
                       desc = dct['desc'],
                       identification = dct['id'],
                       seq = dct['seq'])
    elif 'id' in dct and\
         'seq' in dct:
        return Unit(identification = dct['id'], seq = dct['seq'])
    elif 'original_length' in dct :
        return EncodedRead(idn = dct['id'],
                           original_length = dct['original_length'],
                           leading_gap = dct['leading_gap'],
                           trailing_gap = dct['trailing_gap'],
                           nodes = dct['nodes'],
                           edges = dct['edges'])
    elif 'position_from_start' in dct:
        return Node(position_from_start = dct['position_from_start'],
                    unit = dct['unit'],
                    cluster = dct['cluster'],
                    seq = dct['seq'],
                    is_forward = dct['is_forward'],
                    cigar = dct['cigar'])
    elif 'offset' in dct:
        return Edge(frm = dct['from'], to = dct['to'], offset = dct['offset'],
                    label = dct['label'])
    elif 'pair_id' in dct and 'seq1' in dct and 'seq2' in dct:
        return HiCPair(pair1 = dct['pair1'],
                       pair2 = dct['pari2'],
                       pair_id = dct['pair_id'],
                       seq1 = dct['seq1'],
                       seq2 = dct['seq2'])
    elif 'pair_id' in dct and 'pair1' in dct and 'pair2' in dct:
        return HiCEdge(pair_id = dct['pair_id'],
                       pair1 = dct['pair1'],
                       pair2 = dct['pair2'])
    elif 'id' in dct and 'cluster'in dct:
        return Assignments(idn = dct['id'], cluster = dct['cluster'])
    else:
        return dct


class DataSetEncoder(json.JSONEncoder):
    """
    A encoder for DataSet object.
    
    Examples
    ---------------------
    >>> import json
    >>> dataset = DataSet()
    >>> json.dumps(dataset, cls = DataSetEncoder)
    """
    def default(self, obj):
        if isinstance(obj, DataSet):
            return {
                'raw_reads':list(map(self.default, obj.raw_reads)),
                'selected_chunks':list(map(self.default, obj.selected_chunks)),
                'encoded_reads':list(map(self.default, obj.encoded_reads)),
                'assignments':list(map(self.default, obj.assignments)),
                'hic_pairs':list(map(self.default, obj.hic_pairs)),
                'hic_edges':list(map(self.default, obj.hic_edges)),
            }
        elif isinstance(obj, RawRead):
            return {
                'name':obj.name,
                'desc':obj.desc,
                'id':obj.id,
                'seq':obj.seq
                }
        elif isinstance(obj, Unit):
            return {
                'id': obj.id,
                'seq': obj.seq
            }
        elif isinstance(obj, EncodedRead):
            return {
                'id': obj.id,
                'original_length': obj.original_length,
                'leading_gap': obj.leading_gap,
                'trailing_gap': obj.trailing_gap,
                'nodes': list(map(self.default, obj.nodes)),
                'edges':list(map(self.default, obj.edges)),
            }
        elif isinstance(obj, Node):
            return {
                'position_from_start': obj.position_from_start,
                'unit': obj.unit,
                'cluster': obj.cluster,
                'seq': obj.seq,
                'is_forward':obj.is_forward,
                'cigar': obj.cigar,
            }
        elif isinstance(obj, Edge):
            return {
                'from': obj.frm,
                'to': obj.to,
                'offset': obj.offset,
                'label': obj.label,
            }
        elif isinstance(obj, HiCPair):
            return {
                'pair1': obj.pair1,
                'pair2': obj.pair2,
                'pair_id': obj.pair_id,
                'seq1': obj.seq1,
                'seq2': obj.seq2,
            }
        elif isinstance(obj, HiCEdge):
            return {
                'pair_id': obj.self_id,
                'pair1': obj.pair1,
                'pair2': obj.pair2,
            }
        elif isinstance(obj, Assignments):
            return {
                'id': obj.id,
                'cluster': obj.cluster,
            }
        else:
            print(obj)
            return json.JSONEncoder.default(self, obj)
