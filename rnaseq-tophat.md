# Przygotowanie danych
### Przygotowanie odczytów
Należy utworzyć katalog *fastq* i przekopiować do niego wszystkie pliki zawierające surowe odczyty, które będziemy chcieli poddać analizie
```sh
mkdir fastq
cp path/to/fastq/* fastq/
```
### Przygotowanie genomu referencyjnego i annotacji
Wszystkie referencyjne pliki potrzebne do przeprowadzenia analizy umieścić w katalogu *mm10*
```sh
mkdir mm10
cp -r path/to/reference/* mm10/
```

# Uliniowienie odczytów do genomu referencyjnego
W cel uliniwienia odczytów do genomu referencyjnego użyjemy narzędzia Tophat2: 
* TopHat v.2 https://ccb.jhu.edu/software/tophat/index.shtml


### Uliniowienie narzędziem TopHat
Utwórz katalog *tophat*, do niego zapisywane zostaną wszystkie pliki wygenerowane przez TopHat
```sh
mkdir tophat
```
Następnie uruchom program (ciąg znaków *sample1* powinien zostać zamieniony na właściwą nazwę pliku fastq, który chcemy uniliowić)
```sh
tophat2 -p 6 --GTF mm10/mm10.tss.gtf -o tophat/sample1 mm10/fasta/genome fastq/sample1.fq.gz
```
Kolejne argumenty programu to :

* mm10/mm10.tss.gtf - plik GTF zawierający annotajce
* tophat/sample1 - katalog, w którym mają zostać zapisane wyniki
* mm10/fasta/genome - genom referencyjny wraz z indeksem
* fastq/sample1.fq.gz - plik fastq, który ma zostać uliniowiony

Narzęrzie TopHat utworzy katalog *tophat/sample1*, w którym znajdą się wynikowe pliki. Na potrzeby przykładu skupmy się tylko na dwóch wybranych:
* plik *accepted_hits.bam* - zawiera uliniowione odczyty
* plik *unmapped.fq* - zawiera nieuliniowione odczyty

# Indeksowanie pliku bam
Utworzenie indeksu pliku bam w celu wizualizacji przez IGV
```sh 
samtools index tophat/sample/accepted_hits.bam
```
Argumenty:
index - rodzaj uruchamianego programu
tophat/sample/accepted_hits.bam - sciezka do indeksowanego pliku

# Tworzenie nowej sesji
Nową sesję możemy utworzyć programeme screen
```sh
screen -R nazwa_sesji
```
CTRL+A i zaraz D - wyjście ze screen
CTRL+A i zaraz [ - copy mode, umozliwia przewijanie ekranu

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

Przecinkami oddzielamy próbki z danej grupy. Spacją oddzielamy grupy.

Program Cufflinks zwraca wiele plików, na potrzeby przykładu skupy się na wybranym:
* *isoforms.fpkm_tracking* - tabela ta zawiera wartości FPKM obliczone dla transkryptów oraz roznice w ekspresji trankryptów

# Tworzenie tabeli z wartościami FPKM oraz tabeli z annotacjami
Gdy wszystkie próbki zostały już uliniowione i zostały policzone dla nich poziomy ekspresji traksryptów, pomocne w dalszej analizie jest utworzenie zbiorczej tabeli zawierającej wartości FPKM dla każdej z próbek, oraz tabeli zawierającej annotacje. Można to zrobić za pomocą poniżeszgo skryptu.

```sh
cat cufflinks/sample1/isoforms.fpkm_tracking | head -1 | cut -f1-6 > annotation
cat cufflinks/sample1/isoforms.fpkm_tracking | tail -n+2 | sort -k1 | cut -f1-6 >> annotation
cat annotation | cut -f1 > FPKM.tmp

for i in `ls cufflinks`
    do
        echo $i > tmp.data
        cat cufflinks/$i/isoforms.fpkm_tracking | tail -n+2 | sort -k1 | cut -f10 >> tmp.data
        paste FPKM.tmp tmp.data > FPKM.raw
        cat FPKM.raw > FPKM.tmp
    done
    
rm FPKM.tmp && rm tmp.data
```
W ten sposób otrzymujemy tabelę *annotation* oraz tebelę *FPKM.raw*. W obu tabelach transkrypty są uporządkowane w takiej samej kolejności.

# Analiza w programie R
### Przygotowanie danych
Dalszą analizę wygodnie jest przeprowadzić w programie **R**. W pierwszej kolejności należy wczytać wcześniej utworzone tabele
```r
annotation = read.table("../annotation", sep = "\t", header = T, colClasses = "character")
FPKM.raw = read.table("../FPKM.raw", sep = "\t", header = T, row.names = 1, colClasses = c("character", rep("numeric", 10)))
```

colClasses mówi nam o typie danych w kolumnie
row.names mówi nam o tym, w której kolumnie są nazwy wierszy
sep to rodzaj użytego separatora
header mówi nam czy w pliku jest nagłówek


### Analiza statystyczna
W celach przeprowadzenia analizy statystycznej, napierw utówrzmy tabelę zawierającą zlogarytmowane wartości FPKM. 

```r
FPKM.log = log(1 + FPKM.raw, 2)
```

Kroki podstawowej analizy statystycznej
1. Przeprowadzenie testu T-Studenta, wybranie z wyniku analizy wartości *p-value*
```r
pValue = apply(FPKM.log, 1, function(x){
    t.test(x[indeksy dla 1. grupy], x[indeksy dla 2. grupy], var.equal = T)$p.value})
### alternatywnie mozna utworzyć wektor z czynnikiem (przyporządkowaniem do grupy) oraz użyć go do funkcji t.test jako specyfikacji modelu statystycznej    
group <- as.factor(substr(colnames(FPKM.log), 1, 3))
pValue = apply(FPKM.log, 1, function(x){
  t.test(x~group, var.equal = T)$p.value})
```
gdzie *indeksy dla 1. grupy* oznaczają numery kolumn, w których znajdują się próbki z pierwszej grupy ekperymentalnej, *indeksy dla 2. grupy* oznaczają numery kolumn, w których znajdują się próbki z drugiej grupy ekperymentalnej.

2. Obliczenie skorygowanej wartości *p.value*
```r
fdr = p.adjust(pValue, method = "fdr")
```
### Wizualizacja wyników
Uzyskane wyniki (traskrypty, których poziomy ekspresji różnią się między grupami) można przedstawić na tak zwanej mapie cieplnej. W tym celu należy wybrać transkrypty, których wartość *fdr* jest miejsza od jakieś zadanej, na przykład 0.01
```r
wh = which(fdr < 0.1)
```
Tak wybrane traksrypty można przedstawić na mapie cieplnej
```r
require(gplots)
dev.off()
heatmap.2(
  as.matrix(FPKM.log[wh, ]),
  hclustfun = function(x) hclust(x, method = "average"),
  distfun = function(x) as.dist(1-cor(t(x))),  
  col = bluered(50), 
  scale = 'row', 
  trace = "none",
  margins = c(8,16)
)

###podejrzenie annotacji dla istotnych transkryptów
annotations[wh,]

```







