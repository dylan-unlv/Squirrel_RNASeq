import json as js

fjson = js.load(open('data/fastqc/BLAST_gc_fail.json','r'))
wjson = js.load(open('data/fastqc/BLAST_gc_warn.json','r'))

out = []
for i in range(0,len(fjson['BlastOutput2'])):
    if len(fjson['BlastOutput2'][i]['report']['results']['search']['hits'])>0:
        sname=fjson['BlastOutput2'][i]['report']['results']['search']['query_title']
        description=fjson['BlastOutput2'][i]['report']['results']['search']['hits'][0]['description'][0]['title']
        out.append([sname,'fail',description])
    else:
        sname=fjson['BlastOutput2'][i]['report']['results']['search']['query_title']
        out.append([sname,'fail', 'No hits'])


for i in range(0,len(wjson['BlastOutput2'])):
    if len(wjson['BlastOutput2'][i]['report']['results']['search']['hits'])>0:
        sname=wjson['BlastOutput2'][i]['report']['results']['search']['query_title']
        description=wjson['BlastOutput2'][i]['report']['results']['search']['hits'][0]['description'][0]['title']
        out.append([sname,'warn', description])
    else:
        sname=fjson['BlastOutput2'][i]['report']['results']['search']['query_title']
        out.append([sname,'warn', 'No hits'])
    
    
with open('data/fastqc/blast_summ.txt','w') as f:
    for j in out:
        f.write(','.join(j)+'\n')


