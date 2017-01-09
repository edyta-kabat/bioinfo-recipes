# Przygotowanie danych
### Przygotowanie odczytów
Przekopiuj katalog *metagenome* do swojego katalogu roboczego
```sh
cp -r metagenome path/to/dir/
```

Przejdź do katalogu *mothur/MiSeq_samples*. Pliki pliki wynikowe przeprowadzonych analiz będą zapisywane w tym katalogu.

# Analiza narzędziem mothur
Aby uruchomić narzędzie mothur należy wpisać w konsoli
```sh
mothur
```

### Podstawowe operacje
Wywołanie poniższego polecania utworzy dwa pliki: *fasta* zawierający sekwencje odczytów oraz *qual* zawierający ciągi znaków informujące o jakości odczytów
```sh
mothur > fastq.info(fastq=F3D0_S188_L001_R1_001.fastq)
```
Wskazówka: naciskając dwukrotnie klawisz "tab" otrzymamy podpowiedź jakie pliki znajdują się w katalogu roboczym. 

Teraz można użyć polecenia *summary.seqs*. Zwraca ono podsumowane informacje o pliku fasta. Funkcja ta pozwala kontolować kolejne etapy analizy.
```sh
mothur > summary.seqs(fasta=F3D0_S188_L001_R1_001.fasta)

Using 1 processor.

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	249	249 	0	    3   	1
2.5%-tile:	1	251	251 	0	    4   	195
25%-tile:	1	251	251	    0   	4   	1949
Median: 	1	251	251 	0   	4       3897
75%-tile:	1	251	251	    0   	5   	5845
97.5%-tile:	1	251	251	    0   	6   	7599
Maximum:	1	251	251	    25  	13	    7793
Mean:	1	250.983	250.983	0.00975234	4.63249
# of Seqs:	7793
```

### Utworzenie kontigów
Następnie należy połączyć odczyty z końca R1 oraz R2 za pomocą komendy *make.contigs*. Polecenie to wymaga specjalnego pliku składającego się z 3 następujących kolumn: 1) nazwa próbki, 2) nazwa pliku zawierającego odczyty z końca R1 3) nazwa pliku zawierającego odczyty z R2 (zobacz plik *stability.files*)
```sh
head stability.files

F3D0	F3D0_S188_L001_R1_001.fastq 	F3D0_S188_L001_R2_001.fastq
F3D141	F3D141_S207_L001_R1_001.fastq	F3D141_S207_L001_R2_001.fastq
F3D142	F3D142_S208_L001_R1_001.fastq	F3D142_S208_L001_R2_001.fastq
F3D143	F3D143_S209_L001_R1_001.fastq	F3D143_S209_L001_R2_001.fastq
F3D144	F3D144_S210_L001_R1_001.fastq	F3D144_S210_L001_R2_001.fastq
F3D145	F3D145_S211_L001_R1_001.fastq	F3D145_S211_L001_R2_001.fastq
F3D146	F3D146_S212_L001_R1_001.fastq	F3D146_S212_L001_R2_001.fastq
F3D147	F3D147_S213_L001_R1_001.fastq	F3D147_S213_L001_R2_001.fastq
F3D148	F3D148_S214_L001_R1_001.fastq	F3D148_S214_L001_R2_001.fastq
F3D149	F3D149_S215_L001_R1_001.fastq	F3D149_S215_L001_R2_001.fastq
```
Utworzenie kontigów
```sh
mothur > make.contigs(file=stability.files)
```
Plik wynikowy:
* stability.trim.contigs.fasta

```sh
mothur > summary.seqs(fasta=stability.trim.contigs.fasta)

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	248	248	    0   	3   	1
2.5%-tile:	1	252	252 	0   	3   	3810
25%-tile:	1	252	252 	0   	4   	38091
Median: 	1	252	252 	0   	4   	76181
75%-tile:	1	253	253     0   	5   	114271
97.5%-tile:	1	253	253	    6   	6   	148552
Maximum:	1	502	502	    249	    243	    152360
Mean:	1	252.811	252.811	0.70063	4.44854
# of Seqs:	152360
```
### Filtrowanie danych 
Funkcja *screen.seqs* umożliwia użytkownikowi przefiltrowanie sekwencji.
```sh
mothur > screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=280)
```
Pliki wynikowe:
* stability.trim.contigs.good.fasta - zawiera przefiltrowane sekwencje
* stability.trim.contigs.bad.accnos - zawiera informacje o odrzuconych sekwencjach (nazwa sekwencji; na podstawie jakiego kryterium została odrzucona)
* stability.contigs.good.groups - zawiera informacje o zachowanych sekwencjach (nazwa sekwencji; nazwa próbki z której pochodzi)

### Zduplikowane sekwencje
Poniższe polecenie zwraca plik *fasta* zawierający tylko unikatowe sekwencje. Informacja o pozostałych sekwencjach nie jest utracana przez utworzenie tabeli, w której pogrupowane są wszystkie odczyty o takim samym ciągu znaków.
```sh
mothur > unique.seqs(fasta=stability.trim.contigs.good.fasta)
```
Pliki wynikowe:
* stability.trim.contigs.good.names - tabela: nazwa; nazwy wszystkich takich samych odczytów
* stability.trim.contigs.good.unique.fasta - unikatowe odczyty

```sh
mothur > summary.seqs(fasta=stability.trim.contigs.good.unique.fasta, name=stability.trim.contigs.good.names)

# of unique seqs:	16439
total # of seqs:	128885
```
### Utworzenie "count table"
Kolejnym krokiem jest utworzenie tabeli zawierającej informacje czy dana próbka zawiera daną (unikatową) sekwencję czy nie.
```sh
mothur > count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups)
```
Plik wynikowy:
* stability.trim.contigs.good.count_table 

### Uliniowienie do referencyjnej bazy
Następne etapy analizy przeprowadzimy na innym zbiorze danych (zabieg ten służy jedynie celom szkoleniowym). 
Przejdź do katalogu *PRJDB2729_samples*
```sh
cd ../PRJDB2729_samples
```
Poniżej znajdziemy informacje, jakie pliki z katalogu *PRJDB2729_samples* odpowiadają plikom z katalogu *MiSeq_samples*
```sh
**MiSeq_samples**                               **PRJDB2729_samples**
stability.trim.contigs.good.unique.fasta        PRJDB2729.shhh.trim.unique.fasta
stability.trim.contigs.good.count_table         PRJDB2729.shhh.trim.count_table
```
#### Uliniowienie
```sh
mothur > align.seqs(fasta=PRJDB2729.shhh.trim.unique.fasta, reference=silva.bacteria/silva.bacteria.fasta)
```
Pliki wynikowe:
* PRJDB2729.shhh.trim.unique.align - uliniwiony fasta
* PRJDB2729.shhh.trim.unique.align.report - szczegółowe info do jakich fragmentów każdy odczyt
* PRJDB2729.shhh.trim.unique.flip.accnos - nazwy odczytów

##### Raport z uliniowionych odczytów
```sh
summary.seqs(fasta=PRJDB2729.shhh.trim.unique.align)

		    Start	End	    NBases	Ambigs	Polymer	NumSeqs
Minimum:	-1	    -1	    0	    0	    1   	1
2.5%-tile:	0	    0	    0	    0	    1   	2022
25%-tile:	9964	23961	2	    0	    1   	20215
Median: 	10352	25318	287	    0	    4   	40430
75%-tile:	13139	25495	306	    0	    5   	60645
97.5%-tile:	43115	43116	325 	0   	7   	78838
Maximum:	43116	43116	356	    0   	8       80859
Mean:	16259.6	23682.3	157.332	0	3.20662
# of Seqs:	80859

```

### Filtrowanie
Przefiltrowanie uliniowionych odczytów wg kryteriów narzuconych przed użytkownika
```sh
mothur > screen.seqs(fasta=PRJDB2729.shhh.trim.unique.align, count=PRJDB2729.shhh.trim.count_table, end=0, maxhomop=8) 
```
Pliki wynikowe:
* PRJDB2729.shhh.trim.unique.good.align
* PRJDB2729.shhh.trim.unique.bad.accnos
* PRJDB2729.shhh.trim.good.count_table


### Usunięcie redundancji
```sh
unique.seqs(fasta=PRJDB2729.shhh.trim.unique.good.align, count=PRJDB2729.shhh.trim.good.count_table)
```
Pliki wynikowe:
* PRJDB2729.shhh.trim.unique.good.count_table
* PRJDB2729.shhh.trim.unique.good.unique.align


### Usunięcie błędów sekwencjonowania
Poniższa komenda używa jednego z algorytmów klastrowania w celu usunięcia sekwencji, które prawdopodobnie powstały jako wynik błędu sekwencjonowania
```sh
mothur > pre.cluster(fasta=PRJDB2729.shhh.trim.unique.good.unique.align, count=PRJDB2729.shhh.trim.unique.good.count_table, diffs=2)
```
Output File Names: 
PRJDB2729.shhh.trim.unique.good.unique.precluster.align
PRJDB2729.shhh.trim.unique.good.unique.precluster.count_table
PRJDB2729.shhh.trim.unique.good.unique.precluster.*.map

### Wyszukanie i usunięcie chimer
Wyszukanie sekwencji będących chimerami
```sh
chimera.uchime(fasta=PRJDB2729.shhh.trim.unique.good.unique.precluster.align, count=PRJDB2729.shhh.trim.unique.good.unique.precluster.count_table, dereplicate=t)
```
Pliki wynikowe:
* PRJDB2729.shhh.trim.unique.good.unique.precluster.uchime.pick.count_table
* PRJDB2729.shhh.trim.unique.good.unique.precluster.uchime.chimeras
* PRJDB2729.shhh.trim.unique.good.unique.precluster.uchime.accnos

Usunięcie ich
```sh
remove.seqs(fasta=PRJDB2729.shhh.trim.unique.good.unique.precluster.align, accnos=PRJDB2729.shhh.trim.unique.good.unique.precluster.uchime.accnos)
```
Plik wynikowy: 
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.align

### Przypisanie sekwencjom opisu taksonomicznego
```sh
classify.seqs(fasta=PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.align, count=PRJDB2729.shhh.trim.unique.good.unique.precluster.uchime.pick.count_table, reference=silva.bacteria/silva.bacteria.fasta, taxonomy=silva.bacteria/silva.bacteria.silva.tax)
```
Pliki wynikowe: 
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.silva.wang.taxonomy
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.silva.wang.tax.summary
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.silva.wang.flip.accnos

### Usunięcie wybranych taksonów
Te taksony, które nas nie interesują (Archaea, Chloroplast, mitochondria, Eukaryota, unknown) zostają usunięte z plików do dalszych analiz.
```sh
remove.lineage(fasta=PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.align, count=PRJDB2729.shhh.trim.unique.good.unique.precluster.uchime.pick.count_table, taxonomy=PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.silva.wang.taxonomy, taxon=Archaea-Chloroplast-mitochondria-Eukaryota-unknown)
```
Pliki wynikowe: 
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.silva.wang.pick.taxonomy
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.align
* PRJDB2729.shhh.trim.unique.good.unique.precluster.uchime.pick.pick.count_table

### Sklastrowanie sekwencji do OTU
Utworzenie macierzy odległości 
```sh
dist.seqs(fasta=PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.align, cutoff=0.20)
```
To zajmuje najwięcej czasu (na tym zbiorze danych około 45 minut na 20 rdzeniach)
Plik wynikowy:
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.dist

Klastrowanie sekwencji
```sh
cluster(column=PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.dist, count=PRJDB2729.shhh.trim.unique.good.unique.precluster.uchime.pick.pick.count_table)
```
Zajęło około 100 minut
Plik wynikowy:
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.an.unique_list.list

### Tabela z ilością sekwencji na OTU na próbkę
To polecenie utworzy tabelę zawierającą liczbę sekwencji przypadającą na każde OTU dla każdej próbki.
```sh
make.shared(list=PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.an.unique_list.list, count=PRJDB2729.shhh.trim.unique.good.unique.precluster.uchime.pick.pick.count_table, label=0.01)
```
Pliki wynikowe:
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.an.unique_list.shared
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.an.unique_list.*.rabund


### Utworzenie opisu taksonomicznego dla OTU
```sh
classify.otu(list=PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.an.unique_list.list, count=PRJDB2729.shhh.trim.unique.good.unique.precluster.uchime.pick.pick.count_table, taxonomy=PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.silva.wang.pick.taxonomy, label=0.01
```
Pliki wynikowe: 
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.an.unique_list.0.01.cons.taxonomy
* PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.an.unique_list.0.01.cons.tax.summary


###
Dla ułatwienia dalszej pracy, zmienimy nazwy otrzymanych plików
```sh
cp PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.an.unique_list.shared PRJDB2729_final.shared
cp PRJDB2729.shhh.trim.unique.good.unique.precluster.pick.pick.an.unique_list.0.01.cons.taxonomy PRJDB2729_final.cons.taxonomy
```

### Normalizacja 
W celu normalizacji otrzymanych zliczeń przeskalujemy je na wielkość biblioteki. Poniższe polecenie wypisuje liczbę sekwencji dla każdej próbki
```sh
mothur > count.groups(shared=PRJDB2729_final.shared)

DRR001466 contains 6185.
DRR001467 contains 6089.
DRR001468 contains 3639.
DRR001469 contains 6963.
DRR001470 contains 3797.
DRR001471 contains 5760.
DRR001472 contains 3289.
DRR001473 contains 6364.
DRR001474 contains 3538.
DRR001475 contains 6660.
DRR001476 contains 4373.

Total seqs: 56657.

Output File Names: 
PRJDB2729_final.count.summary
```

```sh
mothur > sub.sample(shared=PRJDB2729_final.shared, size=3289)
```
Plik wynikowy:
* PRJDB2729_final.0.01.subsample.shared

## Analiza w programie R
### 
Wczytanie potrzebnych bibliotek
```r
library(magrittr)
library(ggplot2)
```
Przygotowanie tabeli ze zliczeniami
```r
countTable = read.table("../PRJDB2729_samples/PRJDB2729_final.0.01.subsample.shared")
countTable = t(countTable)
rownames(countTable) = countTable[,1]
colnames(countTable) = countTable[2,]
countTable = countTable[,-1]
countTable = countTable[-(1:3),]
class(countTable) = "numeric"
```
Przygotowanie tabeli z taksonomicznym opisem OTU
```r
taxFile = read.table('../PRJDB2729_samples/PRJDB2729_final.cons.taxonomy', header=T, sep='\t')

# utworzenie tabeli pomocniczej, w ktorej znajda sie koluny z: domena, typ, rzad, rodzina, rodzaj

taxFile.tmp = vector()
for(i in 1:nrow(taxFile)){
  tmp = strsplit(as.character(taxFile$Taxonomy[i]), split = ";", fixed = T) %>% unlist() 
  taxFile.tmp = rbind(taxFile.tmp, tmp)
}

# usuniecie nawiasow i cyfr

taxFile.tmp = gsub("\\(\\d*\\)", "", taxFile.tmp)

# utworzenie ostatecznej tabeli

colnames(taxFile.tmp) = c("Domain",	"Phylum",	"Class",	"Order",	"Family",	"Genus", "NA", "NA", "NA")
rownames(taxFile.tmp) = taxFile$OTU
taxFile = taxFile.tmp[, 1:6]
taxFile = as.matrix(taxFile)

```
Zsumowanie zliczeń dla każdego OTU
```r
otuSum = apply(countTable, 1, sum)
```
Wybranie top 17 OTU
```r
whRow = order(outSum, decreasing = T)
whRow = whRow[!grepl(pattern = "(uncultured)|(unclassified)", taxFile[whRow, 6])]
whRow = whRow[1:17]
```
Przygotownie tabeli do wykresu
```r
values = vector()
sample = vector()
genus = vector()
for(i in 1:11){
  values = c(values, 100 * countTable[tmp.wh, i] / sum(countTable[tmp.wh, i]))
  sample = c(sample, rep(colnames(countTable)[i], length(tmp.wh)))
  genus = c(genus, taxFile[tmp.wh, 6])
}

plot.data = data.frame(values, sample, genus)
```
Zobrazowanie za pomocą funkcji z pakietu *ggplot2*
```r
ggplot(plot.data, aes(x = sample, y = values, fill = genus)) + geom_bar(stat = "identity")
```








