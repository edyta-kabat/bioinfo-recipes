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

```
### Mikrowykład
wprowadzenie do zbioru danych
```
### Zadanie 2.2
W katalogu biostat-data jest plik cpp_6_lekow.xls. Plik uporządkować. Sprawdzić występowanie wyników odstających.
Narysowac boxplot dla wszystkich lekow. Nazwa leku na wykresie to pierwsze trzy litery pisane uppercasem (np. ethanol = ETH)
Narysować boxplot przy pomocy pakietu ggplot2. Kolory slupkow pobrac z pakiet RColorBrewer (Set 1)
Do wykresu z punktu 3 dolozyc wąsy błędów
Policzyc statystyke dla cpp (t testy)

### Zadanie 2.3??
Wczytac dane dla aktywnosci lokomotorycznej
https://github.com/marpiech/bioinfo-recipes/blob/master/biostat-data/lokomotor_6_lekow.xls
Narysowac wykresy liniowe przy pomocy ggplot z wykorzystaniem colorbrewera (trzeba dane zagregowac)
Dolozyc odchylenie standardowe do wykresu
Policzyc ANOVA z uwzględnieniem punktów czasowych
