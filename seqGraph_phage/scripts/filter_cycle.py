import re
import os
import sys
res = set()
ignore_len = int(sys.argv[2])
with open(sys.argv[1]) as r:
    for line in r.readlines():
        if "loop" in line or "iter" in line:
            continue
        line_len = 0
        splited=re.split(r'[+-]',line.strip())
        for v in splited:
            if v == "" or v ==" ":
                continue
            if ignore_len == 0:
                line_len = line_len + int(v.split('_')[3])
        if line_len >= 10000 and ignore_len == 0:
        #print(line,"xxx\n")
            liner = line.replace("cycle","").replace("score","").replace("self","").replace("gene","").replace("ref","")
            res.add(liner.strip("\n"))
        else:
            pass
            #liner = line.replace("cycle","").replace("score","").replace("self","").replace("gene","").replace("ref","")
            #res.add(liner.strip("\n"))
for item in res:
    print(item.replace("+","+\t").replace("-","-\t"))
