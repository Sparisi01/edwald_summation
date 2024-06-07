#ifndef CONSTANTS_H
#define CONSTANTS_H

#define PI 3.14159265358
#define SQR_2 1.41421356237  // Square root of 2
#define SQR_3 1.73205080757  // Square root of 3
#define SQR_PI 1.77245385091 // Square root of PI

// Random seed for reproducibility
#define RAND_SEED 1

const unsigned int _N_ELEMENTS = 118;

// from https://www.qmul.ac.uk/sbcs/iupac/AtWt/index.html#02 (2019)
const double atomic_masses[118] =
    {1.008, 4.002, 6.94, 9.012, 10.81, 12.011, 14.007, 15.999,
     18.998, 20.1797, 22.989, 24.305, 26.981, 28.085, 30.973, 32.06,
     35.45, 39.948, 39.0983, 40.078, 44.955, 47.867, 50.9415, 51.9961,
     54.938, 55.845, 58.933, 58.6934, 63.546, 65.38, 69.723, 72.630,
     74.921, 78.971, 79.904, 83.798, 85.4678, 87.62, 88.905, 91.224,
     92.906, 95.95, 97., 101.07, 102.905, 106.42, 107.8682, 112.414,
     114.818, 118.710, 121.760, 127.60, 126.904, 131.293, 132.905, 137.327,
     138.905, 140.116, 140.907, 144.242, 145., 150.36, 151.964, 157.25,
     158.925, 162.500, 164.930, 167.259, 168.934, 173.045, 174.9668, 178.486,
     180.947, 183.84, 186.207, 190.23, 192.217, 195.084, 196.966, 200.592,
     204.38, 207.2, 208.980, 209., 210., 222., 223., 226.,
     227., 232.0377, 231.035, 238.028, 237., 244., 243., 247.,
     247., 251., 252., 257., 258., 259., 262., 267.,
     270., 269., 270., 270., 278., 281., 281., 285.,
     286., 289., 289., 293., 293., 294.};

const char elements_symbol[118][2] =
    {" H", "He", "Li", "Be", " B", " C", " N", " O", " F", "Ne", "Na", "Mg", "Al", "Si", " P", " S",
     "Cl", "Ar", " K", "Ca", "Sc", "Ti", " V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
     "As", "Se", "Br", "Kr", "Rb", "Sr", " Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
     "In", "Sn", "Sb", "Te", " I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
     "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", " W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
     "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", " U", "Np", "Pu", "Am", "Cm",
     "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
     "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

double getElementMass(int n)
{
    if (n > _N_ELEMENTS)
    {
        printf("ERROR: %d must be in [1,%d]", n, _N_ELEMENTS);
        exit(EXIT_FAILURE);
    }

    return atomic_masses[n - 1];
}

void printElementSymbol(int n)
{
    if (n > _N_ELEMENTS)
    {
        printf("ERROR: n must be in [1,%d]", _N_ELEMENTS);
        exit(EXIT_FAILURE);
    }

    printf("%c%c", elements_symbol[n - 1][0], elements_symbol[n - 1][1]);
}

#endif