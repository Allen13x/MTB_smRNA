#!/usr/bin/env nextflow

//-------------------------channels----------------------------//
params.train=false
params.outdir = "results"
params.list=''
params.sc = 8

scan = Channel.value("${params.sc}")

outdir = params.outdir
list= file("${params.list}")





//------------------------Process_1---------------------------//
process Export_list{
conda 'conda-forge::numpy conda-forge::pandas'
when:
params.train

input:
file list

output:
file("*.csv") into p1
file("*.csv") into p2a
"""
#!/usr/bin/env python
import numpy as np
import pandas as pd
df=pd.read_csv('$list',delimiter=';',names=['srna','mrna','spec','top','s','sf','mf'])
df['ss']=df[['s','srna']].agg('_'.join,axis=1)
df=df[['ss','mrna','sf','mf','spec']]
df.to_csv('tp_couple.csv',sep=';',header=False,index=False)
"""

}


//------------------------Process_2---------------------------//
process Blast{
conda 'bioconda::blast'
when:
params.train

input:
file(p1)

output:
file("blast/blasted_tp.csv") into p2
"""
#!/bin/bash
filename='$p1'
rm -r blast
mkdir blast
while read p
do
echo "\$p" > out
mrna=`awk -F ';' '{print \$2}' out `
sp=`awk -F ';' '{print \$5}' out `
mkdir blast/\${sp}

awk -F ';' '{print ">"\$2}' out> blast/\${sp}/\${mrna}.fasta
awk -F ';' '{print \$4}' out>> blast/\${sp}/\${mrna}.fasta
done < \$filename

cd blast
for i in */

do
cd \$i
s=`basename \$i`
mkdir db
cd db
echo \$s
esearch -db nuccore -query \$s | efetch -format gene_fasta >all.fasta
makeblastdb -in all.fasta -dbtype nucl -out ppe_all_db
cd ../..
done
cd ..

cd blast
rm blasted_tp.csv
for i in */

do
cd \$i
t=`basename \$i`
echo \$t
cat db/all.fasta | awk '/^>/ {printf("\n%s\n",\$0);next; } { printf("%s",\$0);}  END {printf("\n");}' | tail -n +2 > db/all.fa
for k in *.fasta
do
sc=`blastn -query \$k -db db/ppe_all_db -task megablast -outfmt 6 | head -n 1 | awk '{print \$3}'`
if [[ \$sc == "100.000" ]]
then
f=`blastn -query \$k -db db/ppe_all_db -task megablast -outfmt 6 | head -n 1 | awk -F '\t' '{print \$2}'`
m=`blastn -query \$k -db db/ppe_all_db -task megablast -outfmt 6 | head -n 1 | awk -F '\t' '{print \$1}'`
s=`grep -A 1 -w \${f} db/all.fa | tail -n +2`
c=1
else
m=`head -n 1 \$k | sed "s/>//g"`
s=`tail -n +2 \$k`
c=0
fi
echo "\$m;\$s;\$t;\$c" >> ../blasted_tp.csv
done
cd ..
done
cd ..


"""

}

//------------------------Process_3---------------------------//
process Blasted_list{
conda 'conda-forge::numpy conda-forge::pandas'
when:
params.train


input:
file(p2) from p2
file(p2a) from p2a

output:
file("*.csv") into p3
"""
#!/usr/bin/env python
import pandas as pd
import numpy as np
df=pd.read_csv("$p2a",delimiter=';',names=['ss','mrna','sf','mf','spec'])
df2=pd.read_csv("$p2",delimiter=';',names=['mrna','fasta','spec','boh'])
df3=pd.merge(df,df2,on=['mrna','spec'],how='outer')
df3['fasta']=df3['fasta'].str.lower()
df3['y.l']= df3['fasta'].str.len()
df3['x.l']=df3['mf'].str.len()
df3['fasta']=np.where(df3['y.l']>6000,df3['mf'],np.where(df3['y.l']>df3['x.l'],df3['fasta'],df3['mf']))
df3=df3[['ss','mrna','sf','fasta','spec']]
df3=df3.sort_values(by=['ss'])
df3.to_csv('btp_couple.csv',sep=',',header=True,index=False)
"""

}

//------------------------Process_4---------------------------//
process Inta{
conda 'bioconda::intarna'
when:
params.train


input:
set srna, mrna, sf, fasta, spec from p3.splitCsv(header:true).map{ row-> tuple(row.ss, row.mrna, row.sf, row.fasta,row.spec)}


output:
tuple mrna,val('t'),  srna,file("input/$srna/mrna/${mrna}.fasta"), file("input/$srna/srna/${srna}.fasta") into p4i
tuple srna, mrna, file("false/$srna/srna/${srna}.fasta"), file("false/$srna/mrna/${mrna}.fasta"), file("false/res.csv") into p4f
"""
#!/bin/bash

mkdir input


mkdir input/${srna}
mkdir input/${srna}/mrna
mkdir input/${srna}/srna

echo ">"${srna} > input/${srna}/srna/${srna}.fasta
echo $sf >> input/${srna}/srna/${srna}.fasta
echo ">"$mrna >> input/${srna}/mrna/${mrna}.fasta
echo $fasta >> input/${srna}/mrna/${mrna}.fasta


cd input
echo "id1;id2;start1;end1;start2;end2;E" > res.csv
for d in */; do
cat \${d}mrna/* > \${d}mrna.1fasta
cat \${d}srna/* > \${d}srna.1fasta
echo "at \$d"
IntaRNAsTar -q \${d}srna.1fasta -t \${d}mrna.1fasta --threads 0 | tail -n +2 >> res.csv
rm \${d}mrna.1fasta
rm \${d}srna.1fasta
done
cd ..
mkdir false
cp -r input/* false/
"""

}

//------------------------Process_5---------------------------//
process Cut_false{
conda 'conda-forge::biopython conda-forge::numpy conda-forge::pandas'
when:
params.train


input:
tuple srna, mrna, file(srnaf), file(mrnaf), file(res) from p4f

output:

tuple mrna, val('f'), srna,  file("${mrna}.fasta"), file("${srna}.fasta") into p5f
"""
#!/usr/bin/env python
import numpy as np
from Bio import SeqIO
import pandas as pd
import os
fasta = SeqIO.read(open("$mrnaf"),'fasta')
a = list(fasta.seq)
os.remove("$mrnaf")
ofile = open("$mrnaf",'w')
b = a.copy()
res2=pd.read_csv("$res",sep=';')
del b[int(res2.start1)-1:int(res2.end1)]
b = "".join(b)
id = fasta.id
ofile.write(">"+ id +"\\n"+b+"\\n")
ofile.close
"""

}




(dlist) = (params.train
    ? [Channel.empty()] 
    : [Channel
    .fromPath(params.list)
    .splitCsv(header:true)
    .map{ row-> tuple(row.mrna,row.i, row.srna, file(row.pm),file(row.ps))}])
dset=p4i.concat(p5f).mix(dlist)
dset.into{ mrna_dataset;srna_dataset }


//------------------------Process_6---------------------------//
process Create_modified_seqs{


conda 'conda-forge::biopython conda-forge::numpy conda-forge::pandas'


input:
tuple mrnaID,i,srna,file(pm),file(sp) from mrna_dataset
val sc from scan

output:
set mrnaID,i, file("${mrnaID}.fa") into processed1
script:
"""
#!/usr/bin/env python
import numpy as np
import os
from Bio import SeqIO
fasta = SeqIO.read(open("$pm"),'fasta')
dic = {'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c'}
a = list(fasta.seq)
ofile=open("${mrnaID}.fa",'w')
ofile.write(">" + fasta.id +'_'+'1'+'_'+'0'+ '\\n' + str(fasta.seq) + '\\n')
for j in range($sc):
    for i in range(len(a)):
        b=a.copy()
        b[i:i+j+1]=(list(map(dic.get,b[i:i+j+1])))
        b ="".join(b)
        id = fasta.id + '_' +str(i+1)+'_'+str(j+1)
        ofile.write(">"+ id +'\\n'+b+'\\n')

ofile.close
"""
}

//------------------------Process_7---------------------------//
process Target_pred{

conda 'bioconda::intarna'
joined = srna_dataset.join(processed1,by: [0,1])

input:
tuple mrnaID,i, srnaID, file(pt), file(ps),file(pm) from joined

output:
set mrnaID,srnaID, file("res_${mrnaID}_${srnaID}.csv"),i into processed3

script:
"""
#!/bin/bash
IntaRNAsTar -q $ps -t $pm --threads 0 --outmode=C --out res_${mrnaID}_${srnaID}.csv
"""
}

//------------------------Process_8--------------------------//
process Energy_var_matrix{
conda 'conda-forge::numpy conda-forge::pandas'
input:
tuple mrnaID, srnaID, file(datasetFile), i from processed3

output:
tuple file("*.csv"),srnaID,mrnaID, i into fres1

script:
"""
#!/usr/bin/env python 
import numpy as np
import pandas as pd


df=pd.read_csv('$datasetFile',delimiter=';')
df[['id','mut','q']]=df['id1'].str.extract("(.*)_([^_]+)_([^_]+)\$")
df['mut']=pd.to_numeric(df['mut'])
df['q']=pd.to_numeric(df['q'])
df['id2']='$i'+'_'+df['id2']
ref=df[df['q'] == 0]
ref['mut']=0
if ref.shape[0] == 0:
    ref=df.sort_values(by=['E'],ascending=False).iloc[[0],:].assign(mut=0).assign(q=0).assign(id2='fake'+'_'+df['id2'])



seqa=pd.Series(range(1,df['mut'].max()+1))
#seq1=pd.Series(range(ref.start))
seq2=pd.Series(range(ref['start1'].max(),ref['end1'].max()+1))
seq1=pd.Series(range(ref['start1'].max()-10,ref['end1'].max()+11))

conf_mut = df[df.start1.isin(seq1) & df.end1.isin(seq1)].pop('id1')

f_db=pd.concat([pd.Series(np.tile(range(0,9),len(seqa)),name='q'),pd.Series(np.repeat(range(1,df['mut'].max()+1),9),name='mut')],axis=1)

c_db=pd.merge(f_db,df,on=['q','mut'],how='left')




db_ref=c_db.loc[c_db['q'] == 0]
db_ref['E']=np.where(db_ref.mut.isin(seq2),100,ref['E'])
db_ref



db_filt=c_db.loc[c_db['q']!=0]
db_filt['E']=np.where(db_filt.id1.isin(conf_mut),db_filt['E'],100)

db_filtered=pd.concat([db_ref,db_filt])[['E','mut','q']].drop_duplicates().pivot(index='mut',columns='q',values='E')

mat=db_filtered.loc[~db_filtered.index.isin([0]),:].transpose()

name=ref['id2'][ref.index[0]]+'_'+ref['id'][ref.index[0]]+'.csv'
name


mat.to_csv(name,index=False,header=True)
"""
}

//------------------------Process_9--------------------------//
process Final_window{
conda 'conda-forge::numpy conda-forge::pandas'
publishDir "$outdir/$i", mode: 'copyNoFollow'

input:
tuple file(mat),srna,mrna, i from fres1

output:
file("*.csv") into out1
 """
#!/usr/bin/env python 
import numpy as np
import pandas as pd
import os


df=pd.read_csv("$mat",delimiter=',')

n=int(df.columns[(np.where((df.iloc[0]==100)))][0])-1


cmax=int(np.where(n+100>int(df.shape[1]),int(df.shape[1]),max(n+100,min(df.shape[1],200))))
cmin=int(np.where(n-99<1,1,min(n-99,max(df.shape[1]-199,1))))

mat=df.iloc[:,cmin:cmax+1]
name="$i"+'_'+"$srna"+'_'+"$mrna"+'_w.csv'
mat.to_csv(str(name),index=False,header=True)

"""

}






//-------------------------summary---------------------------//

workflow.onComplete {
println(
"""
Pipeline execution summary
---------------------------
Run as : ${workflow.commandLine}
Completed at: ${workflow.complete}
Duration : ${workflow.duration}
Success : ${workflow.success}
workDir : ${workflow.workDir}
exit status : ${workflow.exitStatus}
""")
}


