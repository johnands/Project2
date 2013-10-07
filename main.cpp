#include <iostream>
#include <armadillo>
#include <math.h>
#include <time.h>
#include <fstream>

using namespace std;
using namespace arma;

// erklærer funksjoner
double maxoffdiag(mat &A, int* k, int* l, int n);
void rotate(mat &A, int k, int l, int n);

// main fyller matrise A og finner egenverdier vha. både armadillo-funksjonen eig_sym og Jacobis algoritme
// eig_sym finner også egenvektorer
// skriver deretter ut de tre første egenvektorene, dvs. de som korresponderer til lambda 1, 2 og 3
int main()
{
    int n = 150;                // antall steg
    int k, l;                   // indeks til største ikke-diag. element
    double omega_r = 5.0;       // styrken til potensialet
    double epsilon = 1.0e-8;    // toleranse
    int iterations = 0;         // teller antall iterasjoner
    mat A(n,n);                 // inputmatrise
    vec rho(n);                 // rho = r/alfa

    // fyller vektor rho
    double rho_max = 4.5;
    double h = rho_max/(n+1);   // steglengde
    double rho_min = h;         // fordi rho_min = 0 er utelatt her, tas med senere
    rho = linspace<vec>(rho_min, rho_max-h, n);     // rho_max utelates

    // lager en potensialvektor
    //vec V = rho % rho;
    // % er komponentvis ganging
    vec V = (omega_r*omega_r) * (rho%rho) + (1.0/rho);

    // fyller matrise A
    vec d = 2.0/(h*h) + V;       // diagonalen
    vec e(n-1);                  // super- og subdiagonal
    e.fill(-1/(h*h));
    A.diag() = d;
    // A.diag(-1) og A.diag(1) aksesserer hhv. subdiagonalen og superdiagonalen til A
    A.diag(-1) = e;
    A.diag(1) = e;

    // finner egenverdier og -vektorer vha. armadillo-funksjon og ser på tidsbruk
    mat B = A;
    vec eigval;                 // vektor som inneholder egenverdiene til A
    mat eigvec;                 // den j-te kollonen til eigvec er den j-te egenvoktoren til A
    clock_t start, finish;
    start = clock();
    eig_sym(eigval, eigvec, B);
    finish = clock();
    // skriver ut de tre minste egenverdiene og tid brukt
    //cout << "n = " << n << endl;
    cout << "omega_r = " << omega_r << endl;
    //cout << "Time elapsed Arma: " << ((finish-start)/CLOCKS_PER_SEC) << endl;
    //cout << "Eigenvalues Arma:" << endl << eigval(0) << endl << eigval(1) << endl << eigval(2) << endl;

    // finner egenverdier vha. Jacobis algoritme
    double maxoff = maxoffdiag(A, &k, &l, n); // finner maks ikke-diag. el. før løkken begynner
    // finner det største ikke-diag. elementet (maxoff) og roterer A ved theta funnet ut i fra indeksene til maxoff helt til
    // Frobenius-normen til det største ikke-diag. elementet er mindre enn epsilon
    start = clock();
    while (fabs(maxoff) > epsilon) {
        maxoff = maxoffdiag(A, &k, &l, n);
        rotate(A, k, l, n);
        iterations++;
    }
    finish = clock();
    cout << "Time elapsed Jacobi: " << ((finish-start)/CLOCKS_PER_SEC) << endl;

    // printer de tre minste egenverdiene og antall iterasjoner brukt
    vec lambda = A.diag();
    lambda = sort(lambda);
    cout << "Eigenvalues Jacobi:" << endl << lambda(0) << endl << lambda(1) << endl << lambda(2) << endl;
    cout << "Number of iterations: " << iterations << endl;

    // skriver ut de tre første egenvektorene til fil data.dat
    fstream outFile;
    outFile.open("data.dat", ios::out);

    outFile << n << endl;
    outFile << rho_max << endl;
    outFile << omega_r << endl;
    for (int i=0; i<n; i++) {
        outFile << eigvec(i,0) << " " << eigvec(i,1) << " " << eigvec(i,2) << endl;
    }
    outFile.close();


    return 0;
}


// funksjon som finner ikke-diag. element med størst absoluttverdi, lagrer indeksene i k, l
double maxoffdiag(mat &A, int* k, int* l, int n)
{
    double maxoff = 0.0;
    for ( int i = 0; i < n; i++ ) {
        for ( int j = i + 1; j < n; j++ ) {
            if ( fabs(A(i,j)) > maxoff ) {
                maxoff = fabs(A(i,j));
                *l = i;
                *k = j;
            }
        }
    }
    return maxoff;
}


// funksjon som finner korrekt sin og cos ut i fra k og l, for deretter å rotere A
void rotate(mat &A, int k, int l, int n)
{
    double s, c;            // s = sin(theta), c = cos(theta)
    if (A(k,l) != 0.0) {
        double t, tau;
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        // sjekke om tau er negativ eller ikke for å få minst mulig t, slik at c blir størst mulig
        // uttrykket for t er omskrevet for å forhindre tap av numerisk presisjon
        if (tau > 0) {
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        c = 1/sqrt(1 + t*t);
        s = c*t;
    }
    else {
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    // oppdaterer matriseel. med indekser k og l
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    // de største elementene settes så til null
    A(k,l) = 0.0;
    A(l,k) = 0.0;
    // oppdaterer så de andre elementene
    for ( int i = 0; i < n; i++ ) {
        if ( i != k && i != l ) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
    }
    return;
}




