/** 
    Software to simulate a wavepacket in a potential V
    based on 
    [1] -   http://www.physics.buffalo.edu/phy410-505
    [2] -   physics/0701150 
            Accurate numerical solutions of the time-dependent Schrödinger equation
            W. van Dijk, F. M. Toyama
    [3] -   Shiff - Quantum Mechanics

**/

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

#ifdef __APPLE__
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

const double pi = 4*atan(1.0);  // pi
int N = 600*10;                    // number of interior points
double L = 100.0*10;               // system extends from x = 0 to x = L
double dx = L / double(N+1);    // grid spacing
double h_bar = 1;               // natural units
double mass = 1;                // natural units
double dt = 0.1;                // time step
vector<double> x;               // coordinates of grid points
bool periodic = true;           // false = oo potential, true = periodic

// Potential V(x)
double V_0 = 0.5;               // height of potential well
double V_width = 10.0;          // width of potential well
double V_center = 0.75 * L;     // center of potential well
bool gaussian = false;          // true = Gaussian, false = step potential

bool harmonic = true;            //true= harmonic potential, false= one of the previous potentials
double Kappa =  0.0002;// 0.002; // constant in the harmonic potential V(x)=1/2 Kappa (x-x0_harmonic)^2
double x0_harmonic = 0.5 * L;


double V(double x) {
    double half_width = abs(0.5 * V_width);
  if (harmonic){
	double d_harmonic = (x-x0_harmonic);
	return  0.5 * Kappa * d_harmonic * d_harmonic;

  }else{
    if (gaussian) {
        double d = x - V_center;
        return V_0 * exp( - d * d / (2 * half_width * half_width));
    } else {
        if (abs(x - V_center) <= half_width)
            return V_0;
        else
            return 0.0;
    }
  }
}

// define names for complex numbers and vectors
typedef complex<double> cdouble;
typedef vector<cdouble> cvector;
const cdouble i(0.0, 1.0);


cvector psi, chi;               // complex wavefunction
cvector a, b, c;                // to represent tridiagonal elements of Q
cdouble alpha, beta;            // corner elements of Q
cvector psi_exact_tmp;           // exact psi

// Initial wave packet
//double x_0 = L / 4.0;           // location of center
double x_0 = L / 2.5;           
double E = 1.0;                 // average energy
double velocity = 1.0;          // average velocity
double k_0;                     // average wavenumber
//double sigma_0 = L / 10.0;      // width of wave packet
double alp  =  sqrt(sqrt(mass* Kappa/( h_bar * h_bar)));  // is the parameter in (5.2) of [2]
double sigma_0 = 1 / alp;
double psi_norm;                // normalization of psi

double t;                       // time


ofstream error_file;        // file where to write error

// initialize waveform at time t=0
void init_psi0 ( vector<double> x, cvector & psi0 ) {
    // inititalize the packet
    //k_0 = mass * velocity / h_bar;
    k_0 = 0 ;
    double psi_norm = 1 / sqrt(sigma_0 * sqrt(pi));
    
    for (int j = 0; j < N; j++) {
        double d = x[j] - x_0;
        double exp_fac = exp(- d * d / (2 * sigma_0 * sigma_0));
        psi0[j] = psi_norm * exp(i * k_0 * x[j]) * exp_fac;
    }
 
    
}


// Define exact solution as a function of t. from (13.22) pag.75 of [3] 
cvector psi_exact( vector<double> x, double t ){
    
    cvector psi_exa = cvector(N);
    double psi_norm = 1 / sqrt(sigma_0 * sqrt(pi)); // normalization factor
    double omega_c = sqrt( Kappa / mass);
    double xi_0 = alp * (x_0 -  x0_harmonic);
    
    for (int j = 0; j < N; j++) {
        double xi   = alp * (x[j] - x0_harmonic);

            
        psi_exa[j] = psi_norm * 
          exp( -0.5 * pow(xi - xi_0 * cos(omega_c * t ) ,2.0 )
            - i *(0.5* omega_c * t + xi * xi_0 * sin(omega_c * t) 
            - 0.25 * (xi_0*xi_0 * sin(2 * omega_c *t) ) )
            );
    }
    return psi_exa;
}

// evaluate the difference between two complex vectors function (5.4) of [2]
double eval_err(cvector & phi_1, cvector & phi_2) {
    double error = 0.0;
    for (int j = 1; j < N; j++) {
        error += pow(abs(phi_1[j] - phi_2[j]),2);         
        //double quadr = pow(deltaphi.real(),2) + pow(deltaphi.imag(),2); 
    }

    return sqrt(error);
}








void initialize () {

    // create and initialize vectors to zero
    chi = psi = psi_exact_tmp = cvector(N);
    a = b = c = cvector(N);

    // reset the lattice and time
    x = vector<double>(N);
    dx = L / double(N+1);
    for (int j = 0; j < N; j++)
        x[j] = (j + 1) * dx;
    t = 0.0;

    init_psi0 ( x, psi ); // see function above 
    psi_exact_tmp =  psi_exact(x,0.0);

    // elements of tridiagonal matrix Q = (1/2)(1 + i dt H / (2 hbar))
    for (int j = 0; j < N; j++) {
        a[j] = - i * dt * h_bar / (8 * mass * dx * dx);
        b[j] = 0.5 + 1j * dt / (4 * h_bar) *
              (V(x[j]) + h_bar * h_bar / (mass * dx * dx));
        c[j] = a[j];
    }
    alpha = c[N-1];
    beta = a[0];
    
    // initialize the file
    error_file.open ("error.csv");
    error_file<<"t    err"<<endl;
    error_file<<t<<" "<<0.0<<endl;
}

cvector solve_tridiagonal(cvector& a, cvector& b, cvector& c, cvector& r) {
    // solve Ax = r where A is tridiagonal with diagonals (a, b, c) and return x
    int n = r.size();
    cvector x(n), gama(n);
    cdouble beta = b[0];
    x[0] = r[0] / beta;
    for (int j = 1; j < n; j++) {
        gama[j] = c[j-1] / beta;
        beta = b[j] - a[j] * gama[j];
        x[j] = (r[j] - a[j] * x[j-1]) / beta;
    }
    for (int j = n-2; j >= 0; j--)
        x[j] -= gama[j+1] * x[j+1];
    return x;
}

cvector solve_tridiagonal_cyclic(
    cvector& a, cvector& b, cvector& c, cvector& r,
    complex<double> alpha, complex<double> beta)
{
    // solve Ax = r where A is tridagonal with corner elements alpha, beta
    int n = r.size();
    cvector x(n), b_prime(n), u(n), z(n);
    cdouble gama = -b[0];
    b_prime[0] = b[0] - gama;
    b_prime[n-1] = b[n-1] - alpha * beta / gama;
    for (int j = 1; j < n; j++)
        b_prime[j] = b[j];
    x = solve_tridiagonal(a, b_prime, c, r);
    u[0] = gama;
    u[n-1] = alpha;
    for (int j = 1; j < n-1; j++)
        u[j] = 0;
    z = solve_tridiagonal(a, b_prime, c, u);
    cdouble fact = x[0] + beta * x[n-1] / gama;
    fact /= 1.0 + z[0] + beta * z[n-1] / gama;
    for (int j = 0; j < n; j++)
        x[j] -= fact * z[j];
    return x;
}



void takeStep () {              // time step using sparse matrix routines
    if (periodic)
        chi = solve_tridiagonal_cyclic(a, b, c, psi, alpha, beta);
    else
        chi = solve_tridiagonal(a, b, c, psi);
    for (int j = 0; j < N; j++)
        psi[j] = chi[j] - psi[j];
    t += dt;
    
    psi_exact_tmp = psi_exact(x,t);
    double errore = eval_err(psi,psi_exact_tmp);
    
    error_file<<t<<" "<<errore<<endl;
    glutPostRedisplay();
}



void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    // draw probability exact
    glColor3ub(200, 200, 0);
    glBegin(GL_LINE_STRIP);
        double pOld_exact = psi_exact_tmp[1].real() * psi_exact_tmp[1].real() +
                      psi_exact_tmp[1].imag() * psi_exact_tmp[1].imag();
        for (int j = 1; j < N; j++) {
            double p_exact = psi_exact_tmp[j].real() * psi_exact_tmp[j].real() +
                       psi_exact_tmp[j].imag() * psi_exact_tmp[j].imag();
            glVertex2d(x[j-1], -4 * pOld_exact);
            glVertex2d(x[j], -4 * p_exact);
            pOld_exact = p_exact;
        }
    glEnd();


    // draw real part of psi
    glColor3ub(0, 0, 255);
    glBegin(GL_LINE_STRIP);
        for (int j = 1; j < N; j++) {
            glVertex2d(x[j-1], psi[j-1].real());
            glVertex2d(x[j], psi[j].real());
        }
    glEnd();

    // draw imaginary part of psi
    glColor3ub(0, 255, 0);
    glBegin(GL_LINE_STRIP);
        for (int j = 1; j < N; j++) {
            glVertex2d(x[j-1], psi[j-1].imag());
            glVertex2d(x[j], psi[j].imag());
        }
    glEnd();

    // draw probability
    glColor3ub(255, 0, 0);
    glBegin(GL_LINE_STRIP);
        double pOld = psi[1].real() * psi[1].real() +
                      psi[1].imag() * psi[1].imag();
        for (int j = 1; j < N; j++) {
            double p = psi[j].real() * psi[j].real() +
                       psi[j].imag() * psi[j].imag();
            glVertex2d(x[j-1], 4 * pOld);
            glVertex2d(x[j], 4 * p);
            pOld = p;
        }
    glEnd();
    


    // draw potential well
    glColor3ub(255, 0, 255);
    glBegin(GL_LINE_STRIP);
        double Vold = V(x[1]);
        for (int j = 1; j < N; j++) {
            double Vnew = V(x[j]);
            glVertex2d(x[j-1], 0.2 * Vold);
            glVertex2d(x[j], 0.2 * Vnew);
            Vold = Vnew;
        }
    glEnd();
   
   /* 
    // draw error in time 
    glColor3ub(255, 0, 255);
    glBegin(GL_LINE_STRIP);
        double Vold = V(x[1]);
        for (int j = 1; j < N; j++) {
            double Vnew = V(x[j]);
            glVertex2d(x[j-1], 0.2 * Vold);
            glVertex2d(x[j], 0.2 * Vnew);
            Vold = Vnew;
        }
    glEnd();*/

    glutSwapBuffers();
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, L, -0.3, 0.3);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

bool running;

void mouse(int button, int state, int x, int y) {
    switch (button) {
    case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN) {
            if (running) {
                glutIdleFunc(NULL);
                running = false;
                error_file.close();
            } else {
                glutIdleFunc(takeStep);
                running = true;
            }
        }
        break;
    default:
        break;
    }
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    initialize();
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(600*2, 500);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Wave packet motion in Schroedinger's equation");
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glShadeModel(GL_FLAT);
    glutDisplayFunc(display);
    glutMouseFunc(mouse);
    glutReshapeFunc(reshape);
    
    /*glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(600*2, 500);
    glutInitWindowPosition(100, 400);
    glutCreateWindow("Seconda"); */
    glutMainLoop();
}
