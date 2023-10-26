import json as js


with open("data/genomic.gff") as f:
  json_obj = "\n".join(f.readlines()[5:])
  jsonparse = json.loads(json_obj)
  print(jsonparse)


gff = js.load(open('data/genomic.gff','r'))
print(gff.keys()[0:10])
