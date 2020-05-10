import mysql.connector
import json

USER = "tcrd"
HOST = "tcrd.kmc.io"


def connect(db):
    return mysql.connector.connect(user=USER, host=HOST, database=db, buffered=True)


def disease_and_compound(db):
    cursor = db.cursor(dictionary=True, buffered=True)
    diseases = []
    # first pass collect all the disease indications
    row = {"disease": ""}
    query = (
        "select * from disease where dtype in "
        "('DisGeNET', 'DrugCentral Indication',"
        "'Monarch', 'UniProt Disease')  order by binary name,dtype"
    )
    cursor.execute(query)
    for result in cursor:
        name = result["name"].encode("utf8")
        if row["disease"] != name:
            if row["disease"] != "":
                # print "{}:".format(len(diseases)),
                # print row
                diseases.append(row)
            row = {"disease": name, "source": result["dtype"].encode("utf8")}

        drug = result["drug_name"]
        if drug != None:
            if "drug" not in row:
                row["drugs"] = set()
            row["drugs"].add(drug.encode("utf8"))
        else:
            if "targets" not in row:
                row["targets"] = set()
            row["targets"].add(result["target_id"])
    diseases.append(row)
    cursor.close()

    # print "{} diseases!".format(len(diseases))
    query = (
        "select * from protein a, target b, t2tc c, drug_activity d "
        "where c.target_id = b.id "
        "and c.protein_id = a.id "
        "and d.target_id = c.target_id "
        "and b.id = %(id)s"
    )
    cursor = db.cursor(dictionary=True)
    n = 0
    print("[", end=" ")
    for d in diseases:
        drugs = []
        if "drugs" in d:
            drugs = list(d["drugs"])
        targets = []
        if "targets" in d:
            for id in d["targets"]:
                cursor.execute(query, {"id": id})
                t = cursor.fetchone()
                if t != None and t["tdl"] == "Tclin":
                    # print "{}: {}".format(t['uniprot'],t['name'])
                    x = {
                        "uniprot": t["uniprot"].encode("utf8"),
                        "name": t["name"].encode("utf8"),
                        "drug": t["drug"].encode("utf8"),
                    }
                    if "action_type" in t and t["action_type"] != None:
                        x["moa"] = t["action_type"].encode("utf8")
                    #                    if 'nlm_drug_info' in t and t['nlm_drug_info'] != None:
                    #                        x['description'] = t['nlm_drug_info'].encode('utf8')
                    targets.append(x)
            if len(targets) > 0:
                d["targets"] = targets

        if len(drugs) > 0 or len(targets) > 0:
            if len(drugs) > 0:
                d["drugs"] = drugs
            elif "drugs" in d:
                del d["drugs"]

            if len(targets) > 0:
                d["targets"] = targets
            elif "targets" in d:
                del d["targets"]

            if n > 0:
                print(",")
            print(json.dumps(d, indent=4, separators=(",", ": ")), end=" ")
            n = n + 1
    print("]")
    sys.stderr.write("{} diseases!\n".format(n))


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage:", sys.argv[0], "DATABASE")
        sys.exit(1)
    db = connect(sys.argv[1])

    # print "Successfully connect to DATABASE",sys.argv[1]
    disease_and_compound(db)
    db.close()
