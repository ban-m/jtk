import definitions
import json
FILE = "../result/CCS_reads.15000.1M.entry.units.encode.json"
if __name__ == "__main__":
    print("test")
    # with open(FILE, 'r') as f:
    #     DATASET = json.load(f)
    #     print(len(DATASET['raw_reads']))
    #     print(len(DATASET['selected_chunks']))
    with open(FILE, 'r') as f:
        DATASET = json.load(f, object_hook=definitions.as_dataset)
        serialize = json.dumps(DATASET, cls=definitions.DataSetEncoder)
        print("OK")
        deserialize = json.loads(serialize, object_hook=definitions.as_dataset)
        print(len(deserialize.raw_reads))
        print(len(deserialize.selected_chunks))



