http://doc.yonyoucloud.com/doc/Biopython-cn/cn/chr06.html#id8 <br/>

```
from Bio.Align.Applications import ClustalwCommandline
cline = ClustalwCommandline("clustalw2", infile="opuntia.fasta")
print cline
stdout, stderr = clustalw_cline()
```
