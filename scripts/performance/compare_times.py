
regions = ["ConnectivityMap_computeConnectivity", "CheapRK2"]

def compare_times(manager):
    for test in manager.testlist:
        run_dir = t.directory
        cfile = os.path.join(run_dir, t.options["cali_file"])
        r = cr.CaliperReader()
        r.read(cali_file)
        records = r.records
        for i in records:
            if ("region" in i):
                fname = i["region"]
                if (type(fname) is list):
                    fname = fname[-1]
                if (fname in regions):
                    print(f"{
