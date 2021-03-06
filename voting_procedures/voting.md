---
title: "Abstimmungsverfahren ('voting methods')"
author: "Hans W Borchers<br/>Duale Hochschule BW Mannheim"
date: "Januar 2021"
output:
  html_document:
    number_sections: false
    toc: true
    toc_depth: 3
    toc_float:
      collapsed: false
    css: style_rmd.css
    theme: cosmo
    keep_md: true
    self_contained: false
---

------------------------------------------------------------------------

## Einführung

Nehmen wir an, eine Gruppe von Teilnehmern einer R Benutzergruppe muss sich für das Schwerpunktthema ihres nächsten Treffens entscheiden. Verfügbare Themen könnten zum Beispiel sein:

    A : "R Markdown: Stunning new features"
    B : "Voting procedures: 'vote' package"
    C : "Machine Learning with 'tidymodels'"
    D : "Functional Time Series Analysis"

Es sind 24 Mitglieder anwesend, die nach ihren Präferenzen, d.h. einer Rangfolge der Themen nach ihren persönlichen Interessen, gefragt werden. Es ergibt sich die folgende Liste (die hier zufällig erzeugt wird).


```
##       A B C D
##  [1,] 3 4 2 1
##  [2,] 2 4 1 3
##  [3,] 1 3 4 2
##  [4,] 1 2 3 4
##  [5,] 3 2 1 4
##  [6,] 1 3 4 2
##  [7,] 2 1 3 4
##  [8,] 4 1 3 2
##  [9,] 1 2 3 4
## [10,] 4 2 1 3
## [11,] 2 4 3 1
## [12,] 1 4 3 2
## [13,] 2 4 1 3
## [14,] 4 3 1 2
## [15,] 3 4 1 2
## [16,] 2 4 3 1
## [17,] 2 3 4 1
## [18,] 4 3 1 2
## [19,] 1 4 3 2
## [20,] 4 3 2 1
## [21,] 3 1 2 4
## [22,] 3 1 4 2
## [23,] 4 1 2 3
## [24,] 2 1 3 4
```

Wir haben hier eine sogenannte Präferenzmatrix. `P1[1, ] = c(3, 4, 2, 1)` bedeutet, Teilnehmer 1 bevorzugt 'D' vor 'C' vor 'A' vor 'B'. Solche Listen könnten etwa aufgrund eingereichter Stimmzettel erstellt werden.

Wie kann man auf Grundlage dieser Präferenzen eine "faire" Entscheidung herbeiführen? In der Literatur zur **Entscheidungstheorie** sind viele verschiedene Verfahren mit ihren Vor- und Nachteilen beschrieben worden.

Als Verfahren zum Auszählen der Stimmen und zur Entscheidung über Alternativen werden wir einige Berechnungen von Hand durchführen und im Vergleich auch das Paket *vote* aufrufen.


```r
library(vote)
```

------------------------------------------------------------------------

## Mehrheitsverfahren (engl. 'majority voting')

### Einfache Mehrheit

Wir können einfach zählen, wie oft jede Alternative als erste Präferenz des Benutzers gewählt wird. Wie wir sehen, wird dann Thema 'C' der Gewinner sein.

Um das Paket *vote* zu verwenden, müssen wir alle Stimmen außer `1` löschen. Wir sehen, dass die Alternative 'C' die meisten Stimmen erhält.


```r
P_plurality = P1
P_plurality[P_plurality != 1] = 0
head(P_plurality)
```

```
##      A B C D
## [1,] 0 0 0 1
## [2,] 0 0 1 0
## [3,] 1 0 0 0
## [4,] 1 0 0 0
## [5,] 0 0 1 0
## [6,] 1 0 0 0
```


```r
# Count 'manually':
apply(P_plurality, 2, sum)
```

```
## A B C D 
## 6 6 7 5
```


```r
# Count with the 'vote' package
count.votes(P_plurality, method = "plurality")
```

```
## 
## Results of Plurality voting
## ===========================                           
## Number of valid votes:   24
## Number of invalid votes:  0
## Number of candidates:     4
## Number of seats:          1
## 
## 
## |    |Candidate | Total| Elected |
## |:---|:---------|-----:|:-------:|
## |1   |C         |     7|    x    |
## |2   |A         |     6|         |
## |3   |B         |     6|         |
## |4   |D         |     5|         |
## |Sum |          |    24|         |
## 
## Elected: C
```

### Absolute Mehrheit

Wieder hat jedes Mitglied eine Stimme. Wenn eine Alternative mehr als die Hälfte der Stimmen erreicht, wird sie gewählt. Andernfalls gehen die beiden (oder drei) besten Alternativen in eine Stichwahl. In *vote* wird dies als "tworound.runoff"-Methode bezeichnet.


```r
count.votes(P_plurality, method = "tworound.runoff")
```

```
## 
## Results of two-round-runoff voting
## ==================================                           
## Number of valid votes:   24
## Number of invalid votes:  0
## Number of candidates:     4
## Number of seats:          1
## 
## 
## |    |Candidate | Total| Percent| ROTotal| ROPercent| Elected |
## |:---|:---------|-----:|-------:|-------:|---------:|:-------:|
## |1   |A         |     6|    25.0|       6|      46.2|         |
## |2   |B         |     6|    25.0|       0|       0.0|         |
## |3   |C         |     7|    29.2|       7|      53.8|    x    |
## |4   |D         |     5|    20.8|       0|       0.0|         |
## |Sum |          |    24|   100.0|      13|     100.0|         |
## 
## Elected: C 
## 
## Runoff candidates chosen by a coin toss.
```

Im obigen Fall sollten 'A', 'B' und 'C' in der zweiten Runde laufen. Leider lässt *vote* diesen Fall nicht zu. Deshalb heißt es, die Kandidaten für die Stichwahl werden durch Münzwurf ausgewählt.

Stattdessen können wir wieder 'per Hand' die schwächste Alternative löschen und auf die übrig gebliebenen Alternativen die gleiche Methode anwenden.


```r
# Remove column 'D' from P1
P2 = P1[, c(1,2,3)]

# Generate the 'approval' matrix
for (n in 1:24) {
    d = which.min(P2[n, ])
    P2[n, ] = 0; P2[n, d] = 1
}
```

Wir werden die Abstimmungsmethode auf die drei Alternativen anwenden.


```r
count.votes(P_plurality, method = "tworound.runoff")
```

```
## 
## Results of two-round-runoff voting
## ==================================                           
## Number of valid votes:   24
## Number of invalid votes:  0
## Number of candidates:     4
## Number of seats:          1
## 
## 
## |    |Candidate | Total| Percent| ROTotal| ROPercent| Elected |
## |:---|:---------|-----:|-------:|-------:|---------:|:-------:|
## |1   |A         |     6|    25.0|       6|      46.2|         |
## |2   |B         |     6|    25.0|       0|       0.0|         |
## |3   |C         |     7|    29.2|       7|      53.8|    x    |
## |4   |D         |     5|    20.8|       0|       0.0|         |
## |Sum |          |    24|   100.0|      13|     100.0|         |
## 
## Elected: C 
## 
## Runoff candidates chosen by a coin toss.
```

Leider stellt sich in diesem Fall wieder ein 'Patt' und abermals muss ein Gewinner durch "Münzwurf" entschieden werden.

------------------------------------------------------------------------

## Positionsverfahren (engl. 'positional rules')

### Borda

Die Mehrheitsverfahrenerfahren berücksichtigen nicht, wie oft die Themen auf den Rängen 2, 3 oder 4 gewählt werden? Stattdessen können Punkte vergeben werden, bei `k` Alternativen `k-1` Punkte für den ersten Platz, usw., einen Punkt für den vorletzten Platz und `0` Punkte für den letzten (oder `k` bis 1 Punkte).

Dieses Verfahren wird oft als Borda-Verfahren bezeichnet und zum Beispiel beim "European Song Contest" (ESC) oder bei Formel-1 Rennen angewandt.


```r
P3 = 4 - P1
apply(P3, 2, sum)
```

```
##  A  B  C  D 
## 37 32 38 37
```

Der Gewinner wäre 'C' und 'B' wäre die schlechteste aller Alternativen. Der Grund ist offensichtlich, dass 'C' öfters als 'A' oder 'D' in der Präferenzliste der Teilnehmer auftaucht.

Die Borda-Regel kommt im *vote* Paket leider nicht vor.


### Pairwise comparisons

Wir können uns ansehen, wie oft Alternativen in paarweisen Vergleichen bevorzugt werden. Wir schauen uns also an, wie viele Benutzer 'A' gegenüber 'B' bevorzugen und wie viele 'B' gegenüber 'A' bevorzugen, unabhängig davon, wo in ihrer Präferenzliste dies geschieht.

Die Ausrechnung ist ein wenig komplizierter.


```r
library(hash)
```

```
## hash-2.2.6.1 provided by Decision Patterns
```

```r
H = hash()
for (n in 1:24) {
    for (i in 1:3) {
        for (j in (i+1):4) {
            k = paste(alternatives[P1[n,i]], alternatives[P1[n,j]])
            if (has.key(k, H)) {
                H[k] = values(H, k) + 1
            } else {
                H[k] = 1
            }
        }
    }
}

for (cp in c('A B', 'A C', 'A D', 'B C', 'B D', 'C D')) {
    x = values(H, cp)
    cat(cp, ':\t', x, ':', (24-x), '\n')
}
```

```
## A B :	 13 : 11 
## A C :	 13 : 11 
## A D :	 11 : 13 
## B C :	 11 : 13 
## B D :	 11 : 13 
## C D :	 12 : 12
```

Wenn wir mit einem Vergleich 'A' vs. 'B' und 'C' beginnen, dann wird 'D' gewinnen. Wenn wir zum Beispiel mit 'C' gegen 'D' beginnen, dann kann 'A' gewinnen, da es 'B' und 'C' dominiert. Aber wenn wir mit 'D' gegen 'A' und 'B' beginnen, dann wird 'C' gewinnen.

Das Gleiche können wir mit `count.votes()` und der *Condorcet*-Abstimmungsmethode sehen.


```r
count.votes(P2, method = "condorcet")
```

```
## 
## Results of Condorcet voting
## ===========================                           
## Number of valid votes:   24
## Number of invalid votes:  0
## Number of candidates:     3
## Number of seats:          0
## 
## 
## |   |  A|  B|  C| Total| Loser |
## |:--|--:|--:|--:|-----:|:-----:|
## |A  |  0|  1|  0|     1|       |
## |B  |  0|  0|  0|     0|   x   |
## |C  |  0|  1|  0|     1|       |
## 
## There is no condorcet winner (no candidate won over all other candidates).
## Condorcet loser: B
```

------------------------------------------------------------------------

## Zustimmungsregel

### 'Approval voting'

Beim 'approval voting' hat jeder Teilnehmer (höchstens) so viele Stimmen, wie es Alternativen gibt. Er muss (sollte?) nicht alle Stimmen vergeben, aber er kann nicht Stimmen auf bestimmten Alternativen anhäufen.

Wir simulieren eine solche Abstimmung.


```r
B = matrix(0, 24, 4, dimnames = list(1:24, c("A", "B", "C", "D")))
for (i in 1:24) {
    k = sample(c(1,2,3), 1, prob = c(0.1, 0.6, 0.3))
    a = P1[i, 1:k]
    B[i, a] = 1
}
head(B)
```

```
##   A B C D
## 1 0 0 1 1
## 2 0 1 0 1
## 3 1 0 1 1
## 4 1 1 1 0
## 5 0 1 1 0
## 6 1 0 1 0
```


```r
# apply(B, 2, sum)

count.votes(B, method = "approval")
```

```
## 
## Results of Approval voting
## ==========================                           
## Number of valid votes:   24
## Number of invalid votes:  0
## Number of candidates:     4
## Number of seats:          1
## 
## 
## |    |Candidate | Total| Elected |
## |:---|:---------|-----:|:-------:|
## |1   |D         |    14|    x    |
## |2   |A         |    12|         |
## |3   |B         |    12|         |
## |4   |C         |    12|         |
## |Sum |          |    50|         |
## 
## Elected: D
```

In vielen Abstimmungen wird das Häufen von Stimmen auf bestimmten Alternativen erlaubt. Davon ist dringend abzuraten, weil es die Möglichkeit der Manipulation von Wahlen durch Vorabsprachen eröffnet.


### 'Scoring'

Beim 'scoring' werden die Teilnehmer aufgefordert, alle Alternativen auf einer Skala von 0 bis `k` zu bewerten. Gewinner ist dann die Alternative mit dem höchsten aggregierten Bewertung.


```r
C = matrix(0, 24, 4, dimnames = list(1:24, c("A", "B", "C", "D")))
for (n in 1:24) {
    k = sample(c(2, 3), 1)              # how many
    j = sample(1:4, k)
    C[n, j] = sample(0:8, k)
}
head(C)
```

```
##   A B C D
## 1 5 0 0 1
## 2 0 4 0 7
## 3 0 0 7 0
## 4 0 3 0 1
## 5 8 7 0 0
## 6 0 0 0 5
```

Es werden einfach die Summen in den Spalten gebildet.


```r
count.votes(C, method = "score",
            max.score = 8, larger.wins = TRUE)
```

```
## 
## Results of Score voting
## =======================                           
## Number of valid votes:   24
## Number of invalid votes:  0
## Number of candidates:     4
## Number of seats:          1
## 
## 
## |    |Candidate | Total| Elected |
## |:---|:---------|-----:|:-------:|
## |1   |C         |    72|    x    |
## |2   |A         |    68|         |
## |3   |B         |    52|         |
## |4   |D         |    52|         |
## |Sum |          |   244|         |
## 
## Elected: C
```

Die Benotungen sollten möglichst objektiv sein und nicht mit den Interessen der Teilnehmer auf ein gutes Ranking ihrer Favoriten kollidieren.

------------------------------------------------------------------------

## Weitere Verfahren

### STV: 'Single Transferable Vote'

Beim STV Verfahren (deutsch "übertragbaren Einzelstimmgebung") soll das Problem der unwirksamen Stimmen bei der reinen Mehrheitswahl beheben. Es werden mehrere Sieger ermittelt, daher dient es vor allem bei der Wahl von mehreren Peronen, etwa bei Wahlen zu Ausschüssen.

Von jedem Wähler wird eine Rangfolge aller Kandidaten erstellt. Die Wahlzettel werden abgearbeitet; ist ein Kandidat bereits gewählt, kommt diese Stimme dem nächsten Kandidaten auf der persönlichen Rangliste des Wählers zugute.

Das Verfahren ist komplizierter als die bisherigen soll hier nur mithilfe des *vote* Paketes berechnet werden. Als Beispiel nehmen wir die Präferenzmatrix vom Anfang -- aber beachte, dass für STV auch weniger Stimmen als 4 pro Wähler möglich sind.

Es sollen zwei Themen ausgewählt werden, daher `mcan=2`. Die zu erfüllende Quote ist $Q = \lfloor\frac{v}{s+1}\rfloor + 1$, $v$ die abgegebenen Stimmen, $s$ die Anzahl der Sitze.


```r
count.votes(P1, method = "stv", mcan = 2, verbose =TRUE)
```

```
## 
## Single transferable vote count
## ===================================
## Number of votes cast is 24 
## 
## List of 1st preferences in STV counts: 
## 
## Count: 1 
##   QUOTA A B C D
## 1 8.001 6 6 7 5
## Candidate D eliminated 
## 
## Count: 2 
##   QUOTA A B C
## 2 8.001 9 6 9
## Candidate C elected using forwards tie-breaking method 
## 
## Count: 3 
##   QUOTA     A     B
## 3 8.001 9.444 6.555
## Candidate A elected 
## 
## Results of Single transferable vote
## ===================================                           
## Number of valid votes:   24
## Number of invalid votes:  0
## Number of candidates:     4
## Number of seats:          2
## 
## 
## |           |     1| 2-trans|     2| 3-trans|     3|
## |:----------|-----:|-------:|-----:|-------:|-----:|
## |Quota      | 8.001|        | 8.001|        | 8.001|
## |A          | 6.000|       3| 9.000|   0.444| 9.444|
## |B          | 6.000|       0| 6.000|   0.555| 6.555|
## |C          | 7.000|       2| 9.000|  -0.999|      |
## |D          | 5.000|      -5|      |        |      |
## |Tie-breaks |      |        |     f|        |      |
## |Elected    |      |        |     C|        |     A|
## |Eliminated |     D|        |      |        |      |
## 
## Elected: C, A
```

Thema `D` wird eliminiert wegen zu weniger Erststimmen -- das hatten wir früher schon gesehen.


### Coombs Wahl

Es wird ein einzelner Sieger bestimmt. Ähnlich wie beim STV werden Stimmen übertragen. Allerdings werden so lange Kandidaten eliminiert und ihre Stimmen umverteilt, bis ein Kandidat die absolute Mehrheit erreicht.

Diese Methode ist leider (bisher) nicht im *vote* Paket implementiert.

------------------------------------------------------------------------

## Ein Unmöglichkeitssatz

Welches dieser Verfahren ist nun gerechter bzw. weniger manipulierbar?

Vernünftigerweise sollte ein gerechtes Abstimmungsverfahren die folgenden Bedingungen erfüllen:

* (U) Uneingeschränkter Definitionsbereich
* (P) Pareto-Bdingung
* (I) Unabhängig von irrelevanten Alternativen
* (D) Kein Diktator möglich

Nun gilt:

**Arrows Unmöglichkeitssatz**: *Ist die Menge der Alternativen grösser als 2, dann gibt es kein allgemeines Abstimmungsverfahren, dass die Beingungen (U)-(D) alle erfüllt.*

Dieser Satz wird auch *Arrow's paradox* genannt. Kenneth Arrow erhielt für seine Arbeiten zur "Theorie kollektiver Entscheidungen" (engl. 'social choice theory') im Jahr 1972 den Wirtschafts-Nobelpreis.

------------------------------------------------------------------------

## Nachweise

Deutsche Wikipedia:

* [Präferenzwahlsysteme](https://de.wikipedia.org/wiki/Vorzugswahl)

* [Arrow Theorem](https://de.wikipedia.org/wiki/Arrow-Theorem)
