import jtk
import sys
import json
ARGS = sys.argv
with open(ARGS[1], 'r') as f:
    dataset = json.load(f, object_hook=jtk.as_dataset)
    total_nodes = sum(map(lambda read: len(read.nodes), dataset.encoded_reads))
    print(total_nodes)
