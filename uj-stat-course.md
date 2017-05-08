# Tydzień 1 - wstęp do R i statystyki
### Zadanie 1.1
Z jednego rozkładu normalnego losujemy dwie grupy, każda o liczebności 100. Porównujemy te dwie grupy przy pomocy testu t Studenta. Czynność tą wykonujemy 1000 razy. Należy naszkicować rozkład wartości P otrzymanych z testu t.
```
### Mikrowykład
Błąd pierwszego rodzaju
- P
- istotność
- próg istotności alpha
- podejście pra i post fisherowskie)
Błąd drugiego rodzaju
- moc
```
### Zadanie 1.2
Należy przeprowadzić symulację zadania 1.1 w R i narysować histogram wartości P
```
### Mikrowykład
stackoverflow.com
```
### Zadanie 1.3
Histogram z zadania 1.2 należy narysować przy pomocy biblioteki ggplot2
### Zadanie 1.4
Wylosować dwie grupy o liczebności 100 z rozkładu normalnego o parametrach odpowiednio mean=0 sd=1 dla pierwszej grupy i mean=0.2 i sd=1 dla drugiej grupy. Obliczyć P dla różnicy między grupami. Obliczyć różnicę (nazwaną potem DIFF) pomiędzy średnimi z obydwu grup. Następnie wykonać permutację wektora powstałego z połączenia dwóch grup i obliczyć róznicę między średnimi odpowiednio dla pierwszej i drugiej połowy zpermutowanego wektora. Różnicę zapisać do wektora V. Permutację przeprowadzić 20000 razy. Obliczyć prawdopodobieństwo, że V >= DIFF. Dlaczego wartość P z testu t jest różna od otrzymanego prawdopodobieństwa?
```
### Mikrowykład
parametry, kwartet Ascombiego
```
# Tydzień 2 - wielokrotne testowanie
```
### Mikrowykład
potoki w R
```
### Zadanie 2.1
Wylosować dwie grupy o liczebności 100 z rozkładu normalnego o parametrach odpowiednio mean=0 sd=1 dla pierwszej grupy i mean=0.5 i sd=1 dla drugiej grupy. Obliczyć P dla różnicy między grupami. Czynność powtórzyć 100 razy a P zapisać do wektora V. Wylosować dwie grupy o liczebności 100 z rozkładu normalnego o parametrach odpowiednio mean=0 sd=1 dla pierwszej grupy i mean=0 i sd=1 dla drugiej grupy. Obliczyć P dla różnicy między grupami. Czynność powtórzyć 9900 razy a P zapisać do wektora V. 

Obliczyć ile mamy istotnych wyników w wektorze V dla progu P < 0.05.
Jaki procent istotnych porównań jest przez przypadek?
Jaki procent nie jest przez przypadek?
Odpowiedzieć na to samo pytanie dla progu P < 0.01 i P < 0.001.
Następnie przeprowadzić korekcję na wielokrotne testowanie FDR i powiedzieć ile jest istotnych przy progu q < 0.05.
Ile wśród istotnych porównań jest w pierwszej 100 w wektorze V, ile w pozostałych 9900?
Wnioski?
Przeprowadzić korekcję Bonferroniego. Odpowiedzieć na dwa powyższe pytania.

# Tydzień 3 - data wrangling

```
### Mikrowykład
wprowadzenie do zbioru danych
```
### Zadanie 3.1
W katalogu biostat-data/sixdrugs jest plik cpp_6_lekow.xls. Plik uporządkować. Sprawdzić występowanie wyników odstających i usunąć.
Narysowac boxplot dla wszystkich lekow. Nazwa leku na wykresie to pierwsze trzy litery pisane uppercasem (np. ethanol = ETH)
Narysować boxplot przy pomocy pakietu ggplot2. Kolory slupkow pobrac z pakiet RColorBrewer (Set 1)
Do wykresu z punktu 3 dolozyc wąsy błędów
Policzyc statystyke dla cpp (t testy w stosunku do kontroli)

### Zadanie 3.2
Wczytac dane dla aktywnosci lokomotorycznej
https://github.com/marpiech/bioinfo-recipes/blob/master/biostat-data/sixdrugs/lokomotor_6_lekow.xls
Narysowac wykresy liniowe przy pomocy ggplot z wykorzystaniem colorbrewera (trzeba dane zagregowac)
Dolozyc odchylenie standardowe do wykresu
Policzyc ANOVA z uwzględnieniem punktów czasowych (repeated measures)

# Tydzień 4,5 - analiza zmian ekspresji genów w skali całego trankryptomu

### Zadanie 4.1
W katalogu znajdują się 3 pliki: 1) anno.xls to annotacje sond czyli mapowanie sondy na gen, 2) samples.xls to opis probek, 3) data.xls to pomiary ekspresji genów dla wssystkich sond i wszystkich probek. Wczytać do przestrzeni roboczej.

### Zadanie 4.2
Policzyć dwuczynnikową anovę dla czynników traktowanie (lek) i czas (od podania) z interakcją. Skorygować na wielokrotne testowanie (FDR) wartosci p dla wszystkich czynników. Ile jest istotnych wyników dla poszczególnych czynników dla FDR < 1%

### Zadanie 4.3
Sprawdzić istotność dla czynnika batch. Wystandaryzować przy pomocy mediany

### Zadanie 4.4
Policzyć dwuczynnikową anovę dla czynników traktowanie (lek) i czas (od podania) z interakcją. Skorygować na wielokrotne testowanie (FDR) wartosci p dla wszystkich czynników. Ile jest istotnych wyników dla poszczególnych czynników dla FDR < 1%

### Zadanie 4.5
Przy pomocy mapy cieplnej zwizualizować geny z istotnościa dla czynnika traktowanie.
Skala kolorów: od blue dla najniższego poziomu przez biały do czerwonego dla najwyzszego poziomu.
Wyskalować po wierszach.
Odciąć wyniki odstające na poziomie -2.5sd i 2.5sd.
Sklastrować wiersze przy pomocy metody complete.
Miara odległości - korelacja.
Kolumn nie klastrować. Kolumny poukładać wg. kolejności. NAI (1h, 2h, 4h, 8h), SAL, NIC, ETO, MOR, HER, COC, MET
Wiersze opisać symbolami genów.
Kolumny opisać wg. schematu COC_1h_rep1

# Tydzien 6,7 - analiza zmian ekspresji genow w skali calego genomu
W katalogu bioinfo-recipes/biostat-data/aneurysm znajdują się dane pochodzące z profilowania ekspresji genów w krwi osób po pęknięciu tetniaka.

### Zadanie 6.1
Wczytac dane. Uporządkować przestrzeń roboczą.

### Zadanie 6.2
Sprawdzić, czy plec podana rozni się od plci ustalonej na podstawie profilu ekspresji genow. W trakcie rozwiazywania zadanie wykonac wizualizację

### Zadanie 6.3
Wybrac geny rożniące osoby z tętniakami pękniętymi od kontroli. Wyniki zwizualizować

### Zadanie 6.4
Wybrać top20 genów o ekspresji obniżonej ekspresji w krwi osób z pęknietymi tętniakami i sprawdzić z jakiego typu komórek mogąpochodzić (http://amp.pharm.mssm.edu/Enrichr/)

# Tydzien 8 - analiza czułości i specyficzności testu

W katalogu bioinfo-recipes/biostat-data/aneurysm znajdują sie pomiar stosunku limfocytów do monocytów we krwi pacjentów po pęknięciu tętniaka i bez (lmratio.csv)

### Zadanie 8.1
Narysować wykres czułości i specyficzności wskaźnika (krzywa ROC)

# Tydzien 9 - prezentacja wyników - serwer Shiny

# Zadanie 9.1
Wyjonac aplikację liczącą statystykę chi-kwadrat i Fishera dla tabeli liczności 2x2. Zaprezentowac dla liczebności w poniższym tekście: "W badanej próbie 56 mężczyzn na 810 oraz 14 kobiet na 750 wykazywało objawy dyschromatomii"

### Zadanie do wyboru 9.1 (15 pktów)
Aplikacja licząca nadreprezentację listy genów w terminach ontologii Gene Ontology. Wejście: lista genów. Wyjście: nadreprezentowane terminy ontologiczne. Aby uzyskać maksimum prezentacja powinna być przejrzysta. Można pomyśleć o ustalenia liczby wyników do wyświetlenia

### Zadanie do wyboru 9.2 (15 pktów)
Wykreślić krzywe przeżywalności dla danych 'lung' z pakietu survival. Należy dać użytkownikowi możliwość podziału na centrum, na płeć. Należy dać możliwość ustalenia zakresu wieku. Aby uzsykać max punktóœ należy także pomyśleć nad pozostałymi zmiennymi

### Zadanie do wyboru 9.3 (10 pktów)
Narysować wykres ilości nadawanych imion w zależności od czasu. Wejście: imię, początek imienia. Można uwzględnić płeć. Aby uzyskać maksimum należy użytkownikowi pozwolić wpisać kilka roznych imion

### Zadanie do wyboru 9.4 (10 pktów)
Narysowac wykres zaleznosci oczekiwanej dlugosci zycia od PKB/mieszkańca. Input: rok wykresu. Wielkość koła powinna odwierciedlac wielkosc populacji (logarytmicznie). A kolor koła kontynent. Dane: http://www.gapminder.org/data/

### Zadanie do wyboru 9.5 (15 pktów)
Narysować mapę świata z zaznaczoną kolorem oczekiwaną długością życia w zależności od kraju. Input: rok wykresu. Dane: http://www.gapminder.org/data/

Zadanie należy udostępnić poprzez dowolny publiczny system kontroli wersji. Polecam github lub bitbucket. Init należy wykonać
