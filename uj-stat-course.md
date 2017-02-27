# Tydzień 1 - wstęp do R i statystyki
## Zadanie 1.1
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
## Zadanie 1.2
Należy przeprowadzić symulację zadania 1.1 w R i narysować histogram wartości P
```
### Mikrowykład
stackoverflow.com
```
## Zadanie 1.3
Histogram z zadania 1.2 należy narysować przy pomocy biblioteki ggplot2
## Zadanie 1.4
Wylosować dwie grupy o liczebności 100 z rozkładu normalnego o parametrach odpowiednio mean=0 sd=1 dla pierwszej grupy i mean=0.2 i sd=1 dla drugiej grupy. Obliczyć P dla różnicy między grupami. Obliczyć różnicę (nazwaną potem DIFF) pomiędzy średnimi z obydwu grup. Następnie wykonać permutację wektora powstałego z połączenia dwóch grup i obliczyć róznicę między średnimi odpowiednio dla pierwszej i drugiej połowy zpermutowanego wektora. Różnicę zapisać do wektora V. Permutację przeprowadzić 20000 razy. Obliczyć prawdopodobieństwo, że V >= DIFF. Dlaczego wartość P z testu t jest różna od otrzymanego prawdopodobieństwa?
```
### Mikrowykład
parametry, kwartet Ascombiego
```
