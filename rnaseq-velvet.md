# Przygotowanie danych
### Przygotowanie odczytów
Należy utworzyć katalog *fastq* i przekopiować do niego wszystkie pliki zawierające surowe odczyty, które będziemy chcieli poddać analizie
```sh
mkdir fastq
cp path/to/fastq/* fastq/
```

### Zlepienie odczytów
Ze względu na ograniczony czas i RAM nie wykonamy tego kroku
```sh
cat fastq/sample1 fastq/sample2 > fastq.all
```

# Przygotowanie alignmentu z wykorzystaniem różnych długości k-mer
W celu przygotowania alignmentu użyjemy narzędzia Velvet:
* Velvet https://github.com/dzerbino/velvet/wiki/Manual


### Przygotowanie k-merów narzędziem velveth
Utwórz katalog *assembly*, do niego zapisywane zostaną wszystkie pliki wygenerowane przez velvetg
```sh
mkdir assembly
```
Następnie uruchom program (ciąg znaków *sample1* powinien zostać zamieniony na właściwą nazwę pliku fastq, który chcemy uliniowić)
```sh
velveth assembly/ 23,35,2 -fastq.gz -short fastq/sample1.fq.gz
```
Kolejne argumenty programu to :

* 23,35,2 - dlugosci k-mer do alignmentu od,do,krok
* assembly/ - katalog, w którym mają zostać zapisane wyniki
* short - krótkie odczyty niesparowane
* fastq/sample1.fq.gz - plik fastq, który ma zostać uliniowiony

Narzędzie velveth utworzy katalogi *assembly/23-35*, w którym znajdą się wszystkie możliwe k-mery

# Wykonanie alignmentu z wykorzystaniem różnych długości k-mer

### Wykonanie grafów de Bruijna narzędziem velvetg
Uruchom program na wszystkich długościach
```sh
for i in assembly*; do velvetg $i -read_trkg yes; done
```
Kolejne argumenty programu to :

* $i - katalog, z którego czytane są pliki
* read_trkg - assembly bogatszy w informacje

Narzęrzie velvetg utworzy w katalogach *assembly/23-35* contigi w pliku contigs.fa

### Mergowanie wyników
Przygotujemy k-mery ze zmergowanych contigów 
```sh
velveth Merged 27 -long dir*/contigs.fa
```

### Mergowanie wyników
Drugie wykonanie contigów 
```sh
velvetg Merged/ -read_trkg yes -conserverLong yes
```

# Ostateczny assembly
W tym kroku złożymy transkrypty dla otrzymanych locusów
```sh
oases Merged/ -merge -min_trans_lgth 200
```

# Zbudowanie genomu referencyjnego
```sh
bowtie2-build Merged/transcripts.fa transcripts
```

# Uliniowienie odczytów do genomu referencyjnego
```sh
bowtie2 -x Merged/transcripts.fa -U fastq/sample1.fq.gz -S sample1
```

# Obliczenie poziomu ekspresji transkryptów.
W tej części użyjemy narzędzia *samtools*
```sh 
mkdir exprs
```
Uruchamiamy program
```sh
samtools idxstats sample1.sam > exprs/sample1
```

# Przeprowadzenie analizy na pozostałych próbkach
Wszystkie wcześniej wymienione kroki (poza tworzeniem katalogów) należy powtórzyć dla wszystkich próbek.

# Analiza statystyczna EdgeR
Gdy wszystkie próbki zostały już uliniowione i zostały policzone dla nich poziomy ekspresji traksryptów, możemy wykonać analizę statystyczną przy pomocy programu edgeR w R. EdgeR można wykorzystywać do zmiennych dyskretnych o rozkłądzie bini

```sh
### jedna z kolumn wczytywanego pliku powinna stanowić annotacje transkryptów
anno <- <- read.table("sample.file1", colClasses = "character")[,1]
### wybieramy kolumnę z danymi dla każdej próbki
sample1 <- as.numeric(read.table("sample.file1", colClasses = "character")[,2])
sample2 <- as.numeric(read.table("sample.file2", colClasses = "character")[,2])
### mergujemy wszystkie kolumny
counts <- cbind(sample1, sample2, sample3)
### pr<yporządkowujemy kolumny do grup
group = as.factor(c(1,1,1,2,2,2))
### budujemy model
y <- DGEList(counts=counts,group=group)
### liczymy normalizację
y <- calcNormFactors(y)
### obliczenie macierzy prawdopodobienstwa
y <- estimateDisp(y,design)
### obliczenie statystyki
et <- exactTest(y)
### korekcja na wielokrotne testowanie
fdr <- topTags(et, n = nrow(counts), sort.by = "none")
```
W ten sposób otrzymujemy tabelę *annotation* oraz tebelę *FPKM.raw*. W obu tabelach transkrypty są uporządkowane w takiej samej kolejności.
