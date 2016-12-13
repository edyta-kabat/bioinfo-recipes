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
```
bowtie2-build Merged/transcripts.fa transcripts
```

# Obliczenie poziomu ekspresji transkryptów.
W tej części użyjemy narzędzia *Cufflinks* http://cole-trapnell-lab.github.io/cufflinks/releases/v2.2.1/.
Dla porządku tworzymy katalog *cufflinks*, w którym znajdą się wynikowe pliki
```sh 
mkdir cufflinks
```
Uruchamiamy program
```sh
cufflinks -G mm10/mm10.tss.gtf -o cufflinks/sample1 bam/sample1.bam
```
Poszczególne parametry to:
* *mm10/mm10.tss.gtf* - annotacje; program obliczy poziomy ekspresji transkryptów zawartych w tym pliku
* cufflinks/sample1 - katalog, w którym zostaną zapisane wynikowe pliki
* bam/sample1.bam - uliniowiona próbka, dla której chcemy obliczyć poziomy ekspresji transkryptów

Program Cufflinks zwraca wiele plików, na potrzeby przykładu skupy się na wybranym:
* *isoforms.fpkm_tracking* - tabela ta zawiera wartości FPKM obliczone dla transkryptów.

# Przeprowadzenie analizy na pozostałych próbkach
Wszystkie wcześniej wymienione kroki (poza tworzeniem katalogów) należy powtórzyć dla wszystkich próbek.

# Identyfikacja transkryptów o zmienionej ekspresji
W tej części użyjemy narzędzia *cuffdiff* http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/.
Dla porządku tworzymy katalog *cuffdiff*, w którym znajdą się wynikowe pliki
```sh 
mkdir cuffdiff
```
Uruchamiamy program
```sh
cuffdiff -o cuffdiff mm10/mm10.tss.gtf bam/sample1.bam,bam/sample2.bam bam/sample3.bam,bam/sample4.bam
```
Poszczególne parametry to:
* *mm10/mm10.tss.gtf* - annotacje; program obliczy poziomy ekspresji transkryptów zawartych w tym pliku
* cuffdiff - katalog, w którym zostaną zapisane wynikowe pliki
* bam/sample1.bam,bam/sample2.bam - zestaw próbek z pierwszej grupy
* bam/sample3.bam,bam/sample4.bam - zestaw próbek z drugiej grupy

Program Cufflinks zwraca wiele plików, na potrzeby przykładu skupy się na wybranym:
* *isoforms.fpkm_tracking* - tabela ta zawiera wartości FPKM obliczone dla transkryptów oraz roznice w ekspresji trankryptów

# Tworzenie tabeli z wartościami FPKM oraz tabeli z annotacjami
Gdy wszystkie próbki zostały już uliniowione i zostały policzone dla nich poziomy ekspresji traksryptów, pomocne w dalszej analizie jest utworzenie zbiorczej tabeli zawierającej wartości FPKM dla każdej z próbek, oraz tabeli zawierającej annotacje. Można to zrobić za pomocą poniżeszgo skryptu.

```sh
cat cufflinks/sample1/isoforms.fpkm_tracking | head -1 | cut -f1-6 > annotation
cat cufflinks/sample1/isoforms.fpkm_tracking | tail -n+2 | sort -k1 | cut -f1-6 >> annotation
touch FPKM.tmp

for i in `ls cufflinks`
    do
        echo $i > tmp.data
        cat cufflinks/$i/isoforms.fpkm_tracking | tail -n+2 | sort -k1 | cut -f10 >> tmp.data
        paste -d"\t" FPKM.tmp tmp.data > FPKM.raw
        cat FPKM.raw > FPKM.tmp
    done
    
rm FPKM.tmp && rm tmp.data
```
W ten sposób otrzymujemy tabelę *annotation* oraz tebelę *FPKM.raw*. W obu tabelach transkrypty są uporządkowane w takiej samej kolejności.

# Analiza w programie R
### Przygotowanie danych
Dalszą analizę wygodnie jest przeprowadzić w programie **R**. W pierwszej kolejności należy wczytać wcześniej utworzone tabele
```r
annotation = read.table("path/to/file/annotation", sep = "\t", header = T)
FPKM.raw = read.table("path/to/file/FPKM.raw", sep = "\t", header = T)
```
Dobrze jest nadać wierszom tabeli *FPKM.raw* nazwy transkryptów
```r
rownames(FPKM.raw) = annotaion&tracking_id
```

### Analiza statystyczna
W celach przeprowadzenia analizy statystycznej, napierw utówrzmy tabelę zawierającą zlogarytmowane wartości FPKM. 

```r
FPKM.log = log(FPKM.raw, 2)
```

Kroki podstawowej analizy statystycznej
1. Przeprowadzenie testu T-Studenta, wybranie z wyniku analizy wartości *p-value*
```r
pValue = apply(FPKM.log, 1, function(x){
    t.test(x[indeksy dla 1. grupy], x[indeksy dla 2. grupy], var.equal = T)$p.value})
```
gdzie *indeksy dla 1. grupy* oznaczają numery kolumn, w których znajdują się próbki z pierwszej grupy ekperymentalnej, *indeksy dla 2. grupy* oznaczają numery kolumn, w których znajdują się próbki z drugiej grupy ekperymentalnej.

2. Obliczenie skorygowanej wartości *p.value*
```r
fdr = p.adjust(pValue, method = "fdr")
```
### Wizualizacja wyników
Uzyskane wyniki (traskrypty, których poziomy ekspresji różnią się między grupami) można przedstawić na tak zwanej mapie cieplnej. W tym celu należy wybrać transkrypty, których wartość *fdr* jest miejsza od jakieś zadanej, na przykład 0.01
```r
wh = which(fdr < 0.01)
```
Tak wybrane traksrypty można przedstawić na mapie cieplnej
```r
heatmap.2(
  FPKM.log[wh, ],
  hclustfun = function(x) hclust(x, method = "average"),
  distfun = function(x) as.dist(1-cor(t(x))),  
  col = bluered(50), 
  scale = 'row', 
  trace = "none",
  Colv = NA
)
```







