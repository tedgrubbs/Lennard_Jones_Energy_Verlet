// Programmer: Taylor Grubbs
// Start Date: Tuesday October 1, 2019
// HW 3
//

#include "GL/freeglut.h"
#include "GL/gl.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <chrono>
#include <random>
#include <algorithm>
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

int numParticles = 16;

//timestep
double dt = .01;

//actually the minimum number of ms to wait between frames
double framerate = .1;

//used for dynamically changing timestep. Copied method from Ryan Wilmington
double sigma = .0034;
int step = 1;
int numSteps = 1000000;

//half of Side length of "box". L0/2. Makes it easier for distributions
double boxSize = 5.0;

//file for energy
FILE *energyFile;

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
  x = createRandomPositions(numParticles);
  v = createRandomVelocities(numParticles);
  f = create2dArray(numParticles, 2);
  posList = create3dArray(numSteps+1, numParticles, 2);
  recordPositions(0);

  //using trick to get second position
  for(int i=0; i<numParticles; i++) {
    x[i][0] = x[i][0] + v[i][0]*dt;
    x[i][1] = x[i][1] + v[i][1]*dt;
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
  fclose(energyFile);
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
  double xi = -5.;
  double yi = 0.;
  int j = 0;
  int k = 0;
  for(int i=0; i<numParticles; i++) {
    xNew[i][0] = xi + k*bin;
    xNew[i][1] = yi + j*bin;
    k++;
    if((xi+k*bin) >= 5.) {
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
  for(int i=0; i<numParticles; i++) {
    vNew[i][0] = dis(gen);
    vNew[i][1] = dis(gen);
  }
  return vNew;
}

//update function that runs for every frame
void update(int value) {
  //calculating forces on all particles
  calculateAllForces();

  for(int i=0; i<numParticles; i++) {
    //updates particle positions
    x[i][0] = 2*posList[step][i][0] - posList[step-1][i][0] + f[i][0]*pow(dt,2);
    x[i][1] = 2*posList[step][i][1] - posList[step-1][i][1] + f[i][1]*pow(dt,2);

    //if particle is outside of box, wrap it to other side
    //also need to translate previous steps in order for velocity to be stable
    if(x[i][0] > boxSize) {
      x[i][0] = -boxSize;
      posList[step][i][0] = posList[step][i][0]-2*boxSize;
      posList[step-1][i][0] = posList[step-1][i][0]-2*boxSize;
    }
    if(x[i][0] < -boxSize) {
      x[i][0] = boxSize;
      posList[step][i][0] = posList[step][i][0]+2*boxSize;
      posList[step-1][i][0] = posList[step-1][i][0]+2*boxSize;
    }
    if(x[i][1] > boxSize) {
      x[i][1] = -boxSize;
      posList[step][i][1] = posList[step][i][1]-2*boxSize;
      posList[step-1][i][1] = posList[step-1][i][1]-2*boxSize;
    }
    if(x[i][1] < -boxSize) {
      x[i][1] = boxSize;
      posList[step][i][1] = posList[step][i][1]+2*boxSize;
      posList[step-1][i][1] = posList[step-1][i][1]+2*boxSize;
    }
  }
  // printf("%lf %lf\n", x[0][0], x[0][1]);

  //records positions
  //does other calculations
  if(step < numSteps) {
    step+=1;
    recordPositions(step);
    double totalE = 0.;

    for(int i=0; i<numParticles; i++) {
      v[i][0] = (posList[step][i][0] - posList[step-2][i][0]) / (2*dt);
      v[i][1] = (posList[step][i][1] - posList[step-2][i][1]) / (2*dt);
      double E = pow(v[i][0], 2) + pow(v[i][1], 2);
      totalE += E;
    }
    fprintf(energyFile, "%lf\n", totalE);

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

  double com[2] = {0., 0.};

  //draws particles
  glColor3f(1., 1., 1.);
  for(int i=0; i<numParticles; i++) {
    //need to scale particle positions by boxsize since GL window is only 1x1
    glVertex3f(x[i][0]/boxSize, x[i][1]/boxSize, 0);
    com[0] += x[i][0];
    com[1] += x[i][1];
  }
  glColor3f(1., 0., 0.);
  // glVertex3f(com[0]/(boxSize*numParticles), com[1]/(boxSize*numParticles), 0);

  glEnd();
  glFlush();
}

inline double forceCalc(double r) {
  return 24. * (2./pow(r, 13) - 1./pow(r,7));
}

void calculateAllForces() {
  //reset force matrix to zero
  for(int i=0; i<numParticles; i++) {
    f[i][0] = 0.;
    f[i][1] = 0.;
  }

  //defines global minimum for r to control timestep
  double globalMin = 10000.;

  for(int i=0; i<numParticles; i++) {
    for(int j=i+1; j<numParticles; j++) {

      //don't calculate a self-force
      // if(i == j) continue;

      //stores all possible images of r2
      double r2List[9][2];
      double distances[9];
      int dindex = 0;

      for(int m=-1; m<2; m++) {
        for(int n=-1; n<2; n++) {

          double result[2];
          result[0] = x[i][0] - x[j][0] + m*boxSize*2;
          result[1] = x[i][1] - x[j][1] + n*boxSize*2;

          r2List[dindex][0] = result[0];
          r2List[dindex][1] = result[1];

          distances[dindex] = sqrt((pow(result[0], 2) + pow(result[1], 2)));
          dindex++;
        }
      }
      //gets index of smallest distance
      long int minIndex = std::min_element(distances, distances+9)-distances;
      double rMag = distances[minIndex];
      if(rMag < globalMin) globalMin = rMag;

      //don't count force if particle is too far away
      if(rMag > 3.) continue;

      //calculate force normally
      double fmag = forceCalc(rMag);

      double disp[2] = {r2List[minIndex][0], r2List[minIndex][1]};
      double force[2];
      force[0] = fmag * disp[0]/rMag;
      force[1] = fmag * disp[1]/rMag;
      f[i][0] += force[0];
      f[i][1] += force[1];
      f[j][0] += -1.*force[0];
      f[j][1] += -1.*force[1];
    }
  }
  dt = std::log(1+globalMin*sigma);
}

void recordPositions(int t) {
  for(int i=0; i<numParticles; i++) {
    posList[t][i][0] = x[i][0];
    posList[t][i][1] = x[i][1];
  }
}
