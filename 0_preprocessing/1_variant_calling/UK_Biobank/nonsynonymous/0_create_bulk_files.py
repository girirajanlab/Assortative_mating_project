# creates multiple bulk files and store them in appropriate directories
# directory name should be the final two digits of an eid
import pandas as pd
import numpy as np
import os

### GLOBALS ###
META_BULK_FILE="./ukb48799_23152.bulk"
ROOTDIR="./"


def get_last_two(row):
    value = str(row)
    return value[-2:]

def main():
    df = pd.read_csv(META_BULK_FILE, sep=" ", header=None, engine="python")
    # get the field name from bulk file 
    field_name = os.path.splitext(os.path.basename(META_BULK_FILE))[0].split("_")[1]
    # verify all rows have the same field name
    assert df[1].all() == f"{field_name}_0_0"
    # get the last two digits from eid
    df[2] = df[0].apply(get_last_two)
    # group the df by the last two digits
    grp_obj = df.groupby(2)
    for name, grp in grp_obj:
        # create directory
        dirname = os.path.join(ROOTDIR, name)
        print(f"Creating directory: {dirname}")
        os.makedirs(dirname, exist_ok=True)
        # create bulk file
        basename = os.path.splitext(os.path.basename(META_BULK_FILE))[0]
        bulkfilename = os.path.join(ROOTDIR, name,  basename+"_"+name+".bulk")
        print(f"Creating bulk file: {bulkfilename}")
        grp.iloc[:, [0,1]].to_csv(bulkfilename, sep = " ", header=False, index=False)
    return

if __name__ == "__main__":
    main()