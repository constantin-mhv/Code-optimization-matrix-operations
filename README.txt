Pentru rezultate auxiliare am alocat doar doi vectori dinamic: temp
si left_term. Temp s-a folosit de cateva ori. Rezultatul este returnat
folosind vector temp.

C = A * B * B' + A' * A

===== Blas =====
Ordinea de operatii
1) A * B -> temp
2) res * B' -> left_term
3) A' * A -> temp
4) left_term + temp -> temp
Am folosit functie cblas_dtrmm pentru a inmulti matrici triunghiulare,
cblas_dgemm pentru inmultire matrici simple si cblas_daxpy pentru adunare
matrici. Fosoind aceste functii de inmultire s-a realizat si transpunerea
matricelor specificand parametru CblasTrans. Functia cblas_daxpy face
a*x + y, am specificat a = 1, astfel a s-a realizat adunare simpla de matrici.

===== Neopt si opt =====
Ordinea de operatii
1) B * B' -> temp
2) A * temp -> left_term
3) A' * A -> temp
4) left_term + temp -> temp

Transpunerea matricii nu a fost realizata in bulca suplementara, doar am
accesat linii si coloane invers facand inmultirea. X[i][j] = X'[j][i].
Am observat ca rezultatul la X * X' si  X' * X va fi o matrice simetrica.
Acest lucru mi-a permis sa calculez doar jumatatea de matrice inclusiv
diagonala la operatii 1) si 3). A doua jumatatea se completa folosind prima
jumatatea schimband cu locuri i si j (i != j) : X[i][j] = X[j][i], i != j.
Din aceasta cauza prima operatii a fost B * B' si nu A * B.
Pentru a lua in seama proprietatea matricei triunghiulare a fost "prescurtat"
for loop in cazul operatiilor cu A.
Din cauza ca matricea este pusa in memorie ca un vector la adunarea nu a fost
nevoie de doua loop-uri, s-au adunat elemente la i = 0 : N * N.

===== Inbunatatiri la opt =====
In primul rand am detectat constanta de bucla. Acum resultatul inmultirelor
se face in register sum.
Acum accesul la zona din memorie pentru a scrie rezultat se va face dupa
sfarsitul ultimului for loop si nu in timpul lui.
Dupa am modificat accesul la vectori. Pentru a nu face multe inmultiri si
adunari  pentru a accesa zona din memorie, am folosit pointeri care se
incrementeaza / se inmultesc cu N. Astfel am redus semnificativ timpul pentru
adunari si inmultiri la accesarea zonelor de memorie. Aici foarte mult m-a
ajutat laboaratorul 5. Aceasta ambunatarie a fost folosita si la adunare.
La completarea matricei simetrice (cand se dubleaza rezultatul caclulat la
prima jumatate) se faceau mai putine  adunari/inmultiri daca se facea
accesul la memorie ca la neopt (am si comparat timpii de executie).
Astfel am reusit sa scad timpul de executie pentru N = 1200 de la 22-23 secunde
pana la 8-9 secunde.

===== Analiza comparativa =====
Se poate vedea ca la toate metode timpul de executie creste exponential cu
cresterea N. Totusi se pot vedea timpii de
executie mult mai mici la metoda blas. Eu la metoda optimizata nu am reusit
sa ating performante obtinute la blas. La N = 1600 blas calculeaza in
2 secunde, dar opt_m in 23.