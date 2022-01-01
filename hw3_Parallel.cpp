// Programmer: Taylor Grubbs
// Start Date: Tuesday October 1, 2019
// HW 3
//

#include "GL/freeglut.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <random>
#include <algorithm>
// #include <omp.h>
// #include <iostream> //breaks when compiling in atom for some reason

using namespace std;

//global variables for all functions to easily access
//x and v represent the current values of the particle positions and velocities
//f stores the current value of force acting on all particles
//posList keeps tracks of all particle positions over time
double ** x;
double ** v;
double ** f;
double *** posList;
double ** realPosList;//stores actual particle positions
double currentPotential = 0.;

int numParticles = 16;

//timestep
double dt = .005;
double currtime = 0.;

//actually the minimum number of ms to wait between frames
unsigned int framerate = 1;

//used for dynamically changing timestep. Copied method from Ryan Wilmington
double sigma = .0034;
int step = 1;
int numSteps = 100000;

//half of Side length Falseof "box". L0/2. Makes it easier for distributions
double boxSize = 10.0;

//Heating factor
double g = .99;

//files for recording stuff
FILE *energyFile;
FILE *tempFile;
FILE *pos2File; //stores r^2
FILE *timeFile; //records time. Important if using variable timestep
FILE *comFile;
FILE *speedDistFile;

bool runGraphics = true;

//graphics functions
void drawPoints();
void update(int);
void updateNoGraphics();

//Initialization functions
double ** create2dArray(int,int);
double *** create3dArray(int,int,int);
double ** createRandomPositions(int);
double ** createRandomVelocities(int);
double ** createRandomVFromDist(int);

//copies x to certain time in posList
void recordPositions(int);

//physics calculations
double findMinDistance(int,int);
inline double forceCalc(double);
void calculateAllForces();

int main(int argc, char **argv) {

  auto begin = chrono::high_resolution_clock::now();

  //initializing positions and velocities
  energyFile = fopen("energyOverTime.csv", "w+");
  tempFile = fopen("temp.csv", "w+");
  pos2File = fopen("rSquared.csv", "w+");
  timeFile = fopen("time.csv", "w+");
  comFile = fopen("centerOfMassVelocity.csv", "w+");
  speedDistFile = fopen("speeddistributionfile.csv", "w+");
  x = createRandomPositions(numParticles);
  v = createRandomVelocities(numParticles);
  f = create2dArray(numParticles, 2);
  posList = create3dArray(numSteps+1, numParticles, 2);
  realPosList = create2dArray(numParticles, 2);
  recordPositions(0);

  //using trick to get second position
  for(int i=0; i<numParticles; i++) {
    x[i][0] = x[i][0] + v[i][0]*dt;
    x[i][1] = x[i][1] + v[i][1]*dt;
    realPosList[i][0] = x[i][0];
    realPosList[i][1] = x[i][1];
  }
  recordPositions(1);

  //Graphics stuff
  if(runGraphics) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE);

    //initial window size and position
    glutInitWindowSize(1000, 700);
    glutInitWindowPosition(500, 100);

    //Window title and declaration of draw function
    glutCreateWindow("Lennard Jones");
    glutDisplayFunc(drawPoints);

    //calls update function every "framerate" milliseconds
    //still limited by the speed of the algorithm though
    glutTimerFunc(framerate, update, 0);

    //returns you back to main() after simlation is over
    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);
    glutMainLoop();
  }//else should run without graphics
  else for(int i=0; i<numSteps; i++) update(0);

  //closing things, freeing memory
  fclose(energyFile); fclose(tempFile); fclose(pos2File); fclose(timeFile); fclose(comFile); fclose(speedDistFile);
  free(x); free(v); free(posList); free(f);

  auto end = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(end - begin);
  printf("Done in %lf seconds\n", duration.count()/1000.);
  getchar();
  return 0;
}

double ** create2dArray(int xdim, int ydim) {
  double ** v;
  v = (double **) malloc(xdim * sizeof(double *));
  for(int i=0; i<xdim; i++) {
    v[i] = (double *) malloc(ydim * sizeof(double));
  }
  return v;
}

double *** create3dArray(int xdim, int ydim, int zdim) {
  double *** v;
  v = (double ***) malloc(xdim * sizeof(double **));
  for(int i=0; i<xdim; i++) {
    v[i] = (double **) malloc(ydim * sizeof(double *));
    for(int j=0; j<ydim; j++) {
      v[i][j] = (double*) malloc(zdim*sizeof(double));
    }
  }
  return v;
}

double ** createRandomPositions(int numParticles) {
  double ** xNew = create2dArray(numParticles, 2);

  //creates evenly spaced grid of particles.
  //bin divides separation
  //xi and yi are starting points for lattice
  double bin = 1.;
  double xi = -10.;
  double yi = -10.;
  int j = 0;
  int k = 0;
  for(int i=0; i<numParticles; i++) {
    xNew[i][0] = xi + k*bin;
    xNew[i][1] = yi + j*bin;
    k++;
    //reached end of row
    if((xi+k*bin) >= 2.) {
      k = 0;
      j++;
    }
  }
  return xNew;
}

double ** createRandomVelocities(int numParticles) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-1.5, 1.5);
  double ** vNew = create2dArray(numParticles, 2);
  double vxTot = 0.;
  double vyTot = 0.;
  for(int i=0; i<numParticles-1; i++) {
    vNew[i][0] = dis(gen);
    vNew[i][1] = dis(gen);
    vxTot += vNew[i][0];
    vyTot += vNew[i][1];
  }
  //fixes center of mass velocity
  vNew[numParticles-1][0] = -vxTot;
  vNew[numParticles-1][1] = -vyTot;
  return vNew;
}

//update function that runs for every frame
void update(int value) {
  //calculating forces on all particles
  currentPotential = 0.;
  calculateAllForces();
  // printf("%d\n", step);
  for(int i=0; i<numParticles; i++) {

    //allows for heating or cooling of system
    if(step % 100 == 0) {
      // printf("%d\n", step);
      double rpxPrime = posList[step][i][0] - g*(posList[step][i][0] - posList[step-1][i][0]);
      double rpyPrime = posList[step][i][1] - g*(posList[step][i][1] - posList[step-1][i][1]);
      posList[step-1][i][0] = rpxPrime;
      posList[step-1][i][1] = rpyPrime;
    }
    //updates particle positions
    x[i][0] = 2.*posList[step][i][0] - posList[step-1][i][0] + f[i][0]*pow(dt,2);
    x[i][1] = 2.*posList[step][i][1] - posList[step-1][i][1] + f[i][1]*pow(dt,2);

    //if particle is outside of box, wrap it to other side
    //also need to translate previous steps in order for velocity to be stable
    if(x[i][0] > boxSize) {
      x[i][0] = x[i][0] - 2.*boxSize;
      posList[step][i][0] = posList[step][i][0] - 2.*boxSize;
      posList[step-1][i][0] = posList[step-1][i][0] - 2.*boxSize;
    } else if(x[i][0] < -boxSize) {
      x[i][0] = x[i][0] + 2.*boxSize;
      posList[step][i][0] = posList[step][i][0] + 2.*boxSize;
      posList[step-1][i][0] = posList[step-1][i][0] + 2.*boxSize;
    }
    if(x[i][1] > boxSize) {
      x[i][1] = x[i][1] - 2.*boxSize;
      posList[step][i][1] = posList[step][i][1] - 2.*boxSize;
      posList[step-1][i][1] = posList[step-1][i][1] - 2.*boxSize;
    } else if(x[i][1] < -boxSize) {
      x[i][1] = x[i][1] + 2.*boxSize;
      posList[step][i][1] = posList[step][i][1] + 2.*boxSize;
      posList[step-1][i][1] = posList[step-1][i][1] + 2.*boxSize;
    }
    //incrementing real position
    realPosList[i][0] += x[i][0] - posList[step][i][0];
    realPosList[i][1] += x[i][1] - posList[step][i][1];
  }

  // printf("%lf %lf\n", x[0][0], x[0][1]);
  //lets me look at the initial positions
  // if(step==2) getchar();
  //records positions
  //does other calculations
  if(step < numSteps) {
    step+=1;
    recordPositions(step);
    double totalE = 0.;
    double totalVx = 0.;
    double totalVy = 0.;
    double totalr2 = 0.;
    double totalSpeed = 0.;

    for(int i=0; i<numParticles; i++) {
      v[i][0] = (posList[step][i][0] - posList[step-2][i][0]) / (2.*dt);
      v[i][1] = (posList[step][i][1] - posList[step-2][i][1]) / (2.*dt);
      double E = (pow(v[i][0], 2.) + pow(v[i][1], 2.))*.5;

      //technically starting timer at 2nd step
      double rx = realPosList[i][0] - posList[1][i][0];
      double ry = realPosList[i][1] - posList[1][i][1];
      double r2 = pow(rx, 2.) + pow(ry, 2.);

      totalVx += v[i][0];
      totalVy += v[i][1];
      totalSpeed += sqrt(pow(v[i][0], 2.) + pow(v[i][1], 2.));
      // printf("%lf\n", pow(v[i][0], 2.) + pow(v[i][1], 2.));
      fprintf(speedDistFile, "%lf\n", sqrt(pow(v[i][0], 2.) + pow(v[i][1], 2.)));
      totalE += E;
      totalr2 += r2;
    }
    currtime+=dt;
    // printf("%lf\n", currtime);
    fprintf(energyFile, "%lf\n", (totalE+currentPotential)/numParticles);
    fprintf(pos2File, "%lf\n", totalr2/numParticles);
    fprintf(comFile, "%lf\n", sqrt(pow(totalVx/numParticles, 2.) + pow(totalVy/numParticles, 2.)));
    fprintf(timeFile, "%lf\n", currtime);
    fprintf(tempFile, "%lf\n", pow(totalSpeed/numParticles, 2.)/(2.));

      // for(int i=0; i<numParticles; i++) {
      //   double spd = sqrt(pow(v[i][0], 2.) + pow(v[i][1], 2.));
      //   // fprintf(tempFile, "%lf\n", spd);
      // }

    if(runGraphics) {
      glutPostRedisplay();
      glutTimerFunc(framerate, update, 0);
    }
  }
  else {
    if(runGraphics) glutLeaveMainLoop();
  }
}

void drawPoints() {
  //drawing functions. Not sure what they all do
  glClearColor(0.4, 0.4, 0.4, 0.4);
  glClear(GL_COLOR_BUFFER_BIT);
  glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);

  //creates circles for particles somehow
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

  //particle size
  glPointSize(20);
  glBegin(GL_POINTS);

  // double com[2] = {0., 0.};

  //draws particles
  glColor3f(1., 1., 1.);
  for(int i=0; i<numParticles; i++) {
    //need to scale particle positions by boxsize since GL window is only 1x1
    glVertex3f(x[i][0]/boxSize, x[i][1]/boxSize, 0);
    // com[0] += x[i][0];
    // com[1] += x[i][1];
  }
  glColor3f(1., 0., 0.);
  // glVertex3f(com[0]/(boxSize*numParticles), com[1]/(boxSize*numParticles), 0);

  glEnd();
  glFlush();
}

inline double forceCalc(double r) {
  // printf("%lf\n", 24. * (2./pow(r, 13.) - 1./pow(r,7.)));
  return 24. * (2./pow(r, 13.) - 1./pow(r,7.));
}

inline double potentialCalc(double r) {
  return 4. * (1./pow(r, 12.) - 1./pow(r,6.));
}

void calculateAllForces() {
  //reset force matrix to zero
  for(int i=0; i<numParticles; i++) {
    f[i][0] = 0.;
    f[i][1] = 0.;
  }

  //defines global minimum for r to control timestep
  double globalMin = 10000.;

  // #pragma omp parallel for
  for(int i=0; i<numParticles; i++) {
    for(int j=i+1; j<numParticles; j++) {

      //use fancy trick here to calculate shortest image distance
      //found at https://dasher.wustl.edu/ from Jay Ponder's lecture notes
      double rx = x[i][0] - x[j][0];
      double ry = x[i][1] - x[j][1];
      rx = rx - 2.*boxSize*floor(rx/(2.*boxSize) + .5);
      ry = ry - 2.*boxSize*floor(ry/(2.*boxSize) + .5);
      double rMag = sqrt(pow(rx, 2.) + pow(ry, 2.));

      if(rMag < globalMin) globalMin = rMag;

      //don't count force if particle is too far away
      if(rMag > 3.) continue;

      //calculate force normally
      double fmag = forceCalc(rMag);
      currentPotential += potentialCalc(rMag);
      double force[2];
      force[0] = fmag * rx/rMag;
      force[1] = fmag * ry/rMag;
      f[i][0] += force[0];
      f[i][1] += force[1];
      f[j][0] += -1.*force[0];
      f[j][1] += -1.*force[1];
    }
  }
  // dt = log(1+globalMin*sigma);
  // printf("%lf\n", dt);
}

void recordPositions(int t) {
  for(int i=0; i<numParticles; i++) {
    posList[t][i][0] = x[i][0];
    posList[t][i][1] = x[i][1];
  }
}
