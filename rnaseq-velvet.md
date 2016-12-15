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

### Wykonanie grafów de Bruijna narzędziem velvetg
Uruchom program na wszystkich długościach
```sh
for i in assembly/*; do velvetg $i -read_trkg yes; done
```
Kolejne argumenty programu to :

* $i - katalog, z którego czytane są pliki
* read_trkg - assembly bogatszy w informacje

Narzęrzie velvetg utworzy w katalogach *assembly/23-35* contigi w pliku contigs.fa

### Przygotowanie zmergeowanego zestawu k-merów
Przygotujemy k-mery ze zmergowanych contigów 
```sh
velveth Merged 27 -long assembly/*/contigs.fa
```

### Przygotowanie zmergeowanego zestawu contigów
Drugie wykonanie contigów 
```sh
velvetg Merged/ -read_trkg yes -conserveLong yes
```

# Ostateczny assembly
W tym kroku złożymy transkrypty dla otrzymanych locusów. Transkrypty mają wiele izoform. Velvet tworzy osobny contig dla każdej izoformy. Potrzebne jest narzędzie, które zgrupowałyby transkrypty pochodzące z jednego locusa (genu). W tym celu użyjemy narzędzia oases. https://github.com/dzerbino/oases/tree/master
```sh
oases Merged/ -merge yes -min_trans_lgth 200
```

# Annotacja transkryptów
Aby dowiedzieć się jakie baiłko kodują nasze transkrypty należy je porównać ze znanymi białkami lub transkryptami. Do tego posłuży nam program blast (blastx do białek i blastn do sekwencji nukleotydowych)
```sh
blastx -query Trinity.fasta -db nr -out blastx-Trinity.xml -evalue 1e-6 -num_threads 24 -max_target_seqs 1 -outfmt 
```

Niredundantną bazę białek można pobrać np. z Uniprota. Program tren zwraca xml. Gdy transkryptów nie ma za wiele, można skorzystać z blasta w wersji online.

# Przygotowanie referencji
Aby zliczyć odczyty przy pomocy oprogramowania bowtie. Trzeba najpierw przygotowac referencję. Plik fasta zostanie zindeksowany, co umożliwi szybkie jego wyszukiwanie
```sh
bowtie2-build Merged/transcripts.fa Merged/transcripts
```

# Uliniowienie odczytów do genomu referencyjnego
Uliniowienie do transkryptomu, rożni się od uliniowienia do genomu. Gdy uliniawiamy odczyty pochodzące z sekwencjonowania z RNA do genomu musimy się liczyć z obecnością intronów.

```sh
for i in `ls fastq/*.gz`; do bowtie2 -x Merged/transcripts -U $i | samtools view -bS - | samtools sort - > $i.bam; samtools index $i.bam; done; mkdir bowtie; mv fastq/*.bam* bowtie/
```
Uliniawiamy odczyty do genomu refrencyjnego. Bowtie2 jako argumenty przyjmuje -x czyli indeks do genomu oraz -U ze ścieżką do pliku fastq. Bowtie2 zwraca wynik w formacie SAM. Następnie program samtools view z opcją -b zamienia SAM na BAM. Kreseczka w opcjach programu samtools oznacza wejście standardowe zamiast ścieżki do pliku. Następnie plik jest sortowany. Sortowanie jest konieczne do zrobienia indeksu. Następnie indeksujemy plik.


# Obliczenie poziomu ekspresji transkryptów.
W tej części użyjemy narzędzia *samtools*. Samtool jest pakietem, w ktorym znajduje się wiele narzędzi. Jednym z nich jest idxstats, które wypisuje informację o ilości odczytów uliniowionych do każdego contigu. 
```sh
# utworzenie katalogu
mkdir counts
# zapoczątkowanie naglowka pliku z wynikami
echo "anno" > counts/header
# wrzucenie nazw contigów do pliku z ze zliczeniami
samtools idxstats bowtie/A1_CTRL_trimmed.fq.gz.bam | cut -f1 > counts/counts
# utworzenie pliku ze zliczeniami
for i in `ls bowtie/*.bam`; do echo $i >> counts/header; samtools idxstats $i | cut -f3 | paste counts/counts - > counts.tmp; mv counts.tmp counts/counts; done
cat counts/counts | grep -v "*" > counts/counts.final
```

# Przeprowadzenie analizy na pozostałych próbkach
Wszystkie wcześniej wymienione kroki (poza tworzeniem katalogów) należy powtórzyć dla wszystkich próbek.

# Analiza statystyczna EdgeR
Gdy wszystkie próbki zostały już uliniowione i zostały policzone dla nich poziomy ekspresji traksryptów, możemy wykonać analizę statystyczną przy pomocy programu edgeR w R. EdgeR można wykorzystywać do zmiennych dyskretnych o rozkłądzie bini

```sh
### jedna z kolumn wczytywanego pliku powinna stanowić annotacje transkryptów
anno <- read.table("~/html/rnaseq/counts/counts", colClasses = "character")[,1]
### wybieramy kolumnę z danymi dla każdej próbki
counts <- (read.table("~/html/rnaseq/counts/counts", colClasses = c("character",rep("numeric",10)))[,2:11])
### pr<yporządkowujemy kolumny do grup
group <- as.factor(rep(c("ctrl", "dex"), 5))
colnames(counts) <- paste(rep(c("ctrl", "dex"), 5), rep(1:5, each = 2), sep = "")
### budujemy model
require(edgeR)
y <- DGEList(counts=counts,group=group)
### liczymy normalizację
y <- calcNormFactors(y)
### obliczenie macierzy prawdopodobienstwa
design <- model.matrix(~group)
y <- estimateDisp(y,design)
### obliczenie znormalizowanych wartości ekspresji
y <- estimateCommonDisp(y)
pseudo.counts <- y$pseudo.counts
### obliczenie statystyki
et <- exactTest(y)
### korekcja na wielokrotne testowanie
fdr <- topTags(et, n = nrow(counts), sort.by = "none")
```
W ten sposób otrzymujemy tabelę *fdr*, którą następnie możemy wykorzystać do filtrowania wyników


Pakiety można doinstalować komendami
```sh
setRepositories()
### wybrać 1 2
install.packages("edgeR")
### wybrac server cloud czyli nr 1
```
