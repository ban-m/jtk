import jtk
import sys
import json
ARGS = sys.argv
with open(ARGS[1], 'r') as f:
    dataset = json.load(f, object_hook=jtk.as_dataset)
