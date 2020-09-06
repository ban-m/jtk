import definitions
import json
import sys
if __name__ == "__main__":
    ARGV = sys.argv;
    FILE = ARGV[1]
    with open(FILE, 'r') as f:
        DATASET = json.load(f, object_hook=definitions.as_dataset)
        print("DESERIALIZE OK")
        serialize = json.dumps(DATASET, cls=definitions.DataSetEncoder)
        print("SERIALIZE OK")
        deserialize = json.loads(serialize, object_hook=definitions.as_dataset)
        print("SER->DE OK")
        print("NUMBER OF READS:", len(deserialize.raw_reads))
        print("NUMBER OF UNITS:", len(deserialize.selected_chunks))



