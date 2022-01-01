#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>

using namespace std;

int ** create2dArray(int,int);
void parTest(int **,int,int);
void serTest(int **,int,int);

int main(int argc, char const *argv[]) {

  int w = 10000;
  int h = 10000;
  int ** arr = create2dArray(w, h);
  for (int i = 0; i < w; i++) {
    for (int j = 0; j < h; j++) {
      arr[i][j] = 0;
    }
  }
  printf("Start\n");
  //parallel test
  auto begin = chrono::high_resolution_clock::now();
  parTest(arr, w, h);
  auto end = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::milliseconds>(end - begin);
  printf("Parallel done in %lf seconds\n", duration.count()/1000.);
  getchar();

  // check results
  // for (int i = 0; i < w; i++) {
  //   for (int j = 0; j < h; j++) {
  //     printf("%d, %d has %d\n", i, j, arr[i][j]);
  //   }
  // }
  //
  //serial time
  begin = chrono::high_resolution_clock::now();
  serTest(arr, w, h);
  end = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::milliseconds>(end - begin);
  printf("Serial done in %lf seconds\n", duration.count()/1000.);

  // for (int i = 0; i < w; i++) {
  //   for (int j = 0; j < h; j++) {
  //     printf("%d, %d has %d\n", i, j, arr[i][j]);
  //   }
  // }

  return 0;
}

int ** create2dArray(int xdim, int ydim) {
  int ** v;
  v = (int **) malloc(xdim * sizeof(int *));
  for(int i=0; i<xdim; i++) {
    v[i] = (int *) malloc(ydim * sizeof(int));
  }
  return v;
}

void parTest(int ** arr, int w, int h) {
  #pragma omp parallel for
  for (int i = 0; i < w*10; i++) {
    for (int j = 0; j < h*10; j++) {
      arr[i%w][j%h] = 1;
    }
  }
}

void serTest(int ** arr, int w, int h) {
  for (int i = 0; i < w*10; i++) {
    for (int j = 0; j < h*10; j++) {
      arr[i%w][j%h] = 1;
    }
  }
}
