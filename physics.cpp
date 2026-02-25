#include "jello.h"
#include "physics.h"

// External force field
// Get external force
void getForceField(struct world * jello, struct point pos, struct point &force)
{
    // Initialize as 0
    pMAKE(0, 0, 0, force);
    if (jello->resolution == 0 || jello->forceField == NULL)
        return;
    
    // Number of points on each edge
    int n = jello->resolution;
    double L = 4.0 / (n - 1);
    
    // Convert world coordinates to grid indices
    int i = (int)floor((pos.x + 2.0) / L);
    int j = (int)floor((pos.y + 2.0) / L);
    int k = (int)floor((pos.z + 2.0) / L);
    
    // Boundary checks
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (k < 0) k = 0;
    if (i >= n - 1) i = n - 2;
    if (j >= n - 1) j = n - 2;
    if (k >= n - 1) k = n - 2;
    
    double x0 = -2.0 + i * L;
    double y0 = -2.0 + j * L;
    double z0 = -2.0 + k * L;
    
    double alpha = (pos.x - x0) / L;
    double beta  = (pos.y - y0) / L;
    double gamma = (pos.z - z0) / L;
    
    // Convert grid coordinates to array index
#define FORCE_INDEX(ii,jj,kk) ((ii)*n*n + (jj)*n + (kk))
    struct point f000 = jello->forceField[FORCE_INDEX(i,   j,   k)];
    struct point f001 = jello->forceField[FORCE_INDEX(i,   j,   k+1)];
    struct point f010 = jello->forceField[FORCE_INDEX(i,   j+1, k)];
    struct point f011 = jello->forceField[FORCE_INDEX(i,   j+1, k+1)];
    struct point f100 = jello->forceField[FORCE_INDEX(i+1, j,   k)];
    struct point f101 = jello->forceField[FORCE_INDEX(i+1, j,   k+1)];
    struct point f110 = jello->forceField[FORCE_INDEX(i+1, j+1, k)];
    struct point f111 = jello->forceField[FORCE_INDEX(i+1, j+1, k+1)];
    
    // Calculate weights
    double A000 = (1-alpha) * (1-beta) * (1-gamma);
    double A001 = (1-alpha) * (1-beta) * gamma;
    double A010 = (1-alpha) * beta * (1-gamma);
    double A011 = (1-alpha) * beta * gamma;
    double A100 = alpha * (1-beta) * (1-gamma);
    double A101 = alpha * (1-beta) * gamma;
    double A110 = alpha * beta * (1-gamma);
    double A111 = alpha * beta * gamma;
    
    // Calculate force
    force.x = A000*f000.x + A001*f001.x + A010*f010.x + A011*f011.x + A100*f100.x + A101*f101.x + A110*f110.x + A111*f111.x;
    
    force.y = A000*f000.y + A001*f001.y + A010*f010.y + A011*f011.y + A100*f100.y + A101*f101.y + A110*f110.y + A111*f111.y;
    
    force.z = A000*f000.z + A001*f001.z + A010*f010.z + A011*f011.z + A100*f100.z + A101*f101.z + A110*f110.z + A111*f111.z;
#undef FORCE_INDEX
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
    /* for you to implement ... */
    // Initialize acceleration as 0
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
            {
                pMAKE(0, 0, 0, a[i][j][k]);
            }
    
    // Calculate the force
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
            {
                struct point currentPos = jello->p[i][j][k];
                struct point currentVel = jello->v[i][j][k];
                
                // Structural Springs
                // Define relative positions of neighbors
                int structural[6][3] = {
                    {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}
                };
                
                // Actual positions
                for (int n = 0; n < 6; n++)
                {
                    int ni = i + structural[n][0];
                    int nj = j + structural[n][1];
                    int nk = k + structural[n][2];
                    
                    // Assure all points are within the boundaries
                    if (ni >= 0 && ni < 8 && nj >= 0 && nj < 8 && nk >= 0 && nk < 8)
                    {
                        // Define point info
                        struct point L, direction, elasticForce, dampingForce, totalForce;
                        struct point neighborPos = jello->p[ni][nj][nk];
                        struct point neighborVel = jello->v[ni][nj][nk];
                        
                        pDIFFERENCE(neighborPos, currentPos, L);
                        double length = pLENGTH(L);
                        
                        // Calculate elastic force
                        pCPY(L, direction);
                        double len;
                        len = length;
                        pNORMALIZE(direction);
                        double elasticScalar = jello->kElastic * (len - 1.0/7.0);
                        pMULTIPLY(direction, elasticScalar, elasticForce);
                        
                        // Calculate damping force
                        struct point velDiff;
                        pDIFFERENCE(neighborVel, currentVel, velDiff);
                        double dampingScalar = jello->dElastic * pDOT(velDiff, direction);
                        pMULTIPLY(direction, dampingScalar, dampingForce);
                        
                        // Total force
                        pSUM(elasticForce, dampingForce, totalForce);
                        pSUM(a[i][j][k], totalForce, a[i][j][k]);
                    }
                }
                
                // Shear Springs
                // Diagonal of face
                // Define relative positions
                int shear[12][3] = {
                    {1,1,0}, {1,-1,0}, {-1,1,0}, {-1,-1,0}, {1,0,1}, {1,0,-1}, {-1,0,1}, {-1,0,-1}, {0,1,1}, {0,1,-1}, {0,-1,1}, {0,-1,-1}
                };
                
                // Actual positions
                for (int n = 0; n < 12; n++)
                {
                    int ni = i + shear[n][0];
                    int nj = j + shear[n][1];
                    int nk = k + shear[n][2];
                    
                    // Assure all points are within the boundaries
                    if (ni >= 0 && ni < 8 && nj >= 0 && nj < 8 && nk >= 0 && nk < 8)
                    {
                        // Define point info
                        struct point L, direction, elasticForce, dampingForce, totalForce;
                        struct point neighborPos = jello->p[ni][nj][nk];
                        struct point neighborVel = jello->v[ni][nj][nk];
                        
                        pDIFFERENCE(neighborPos, currentPos, L);
                        double length = pLENGTH(L);
                        
                        // Calculate elastic force
                        pCPY(L, direction);
                        double len;
                        len = length;
                        pNORMALIZE(direction);
                        double elasticScalar = jello->kElastic * (len - sqrt(2.0)/7.0);
                        pMULTIPLY(direction, elasticScalar, elasticForce);
                        
                        // Calculate damping force
                        struct point velDiff;
                        pDIFFERENCE(neighborVel, currentVel, velDiff);
                        double dampingScalar = jello->dElastic * pDOT(velDiff, direction);
                        pMULTIPLY(direction, dampingScalar, dampingForce);
                        
                        // Total force
                        pSUM(elasticForce, dampingForce, totalForce);
                        pSUM(a[i][j][k], totalForce, a[i][j][k]);
                    }
                }
                
                // Body diagonal
                int bodyDiag[4][3] = {
                    {1,1,1}, {1,-1,1}, {-1,1,1}, {-1,-1,1}
                };
                
                for (int n = 0; n < 4; n++)
                {
                    int ni = i + bodyDiag[n][0];
                    int nj = j + bodyDiag[n][1];
                    int nk = k + bodyDiag[n][2];
                    
                    if (ni >= 0 && ni < 8 && nj >= 0 && nj < 8 && nk >= 0 && nk < 8)
                    {
                        struct point L, direction, elasticForce, dampingForce, totalForce;
                        struct point neighborPos = jello->p[ni][nj][nk];
                        struct point neighborVel = jello->v[ni][nj][nk];
                        
                        pDIFFERENCE(neighborPos, currentPos, L);
                        double len = pLENGTH(L);
                        
                        pCPY(L, direction);
                        double length;
                        pNORMALIZE(direction);
                        double elasticScalar = jello->kElastic * (len - sqrt(3.0)/7.0);
                        pMULTIPLY(direction, elasticScalar, elasticForce);
                        
                        struct point velDiff;
                        pDIFFERENCE(neighborVel, currentVel, velDiff);
                        double dampingScalar = jello->dElastic * pDOT(velDiff, direction);
                        pMULTIPLY(direction, dampingScalar, dampingForce);
                        
                        pSUM(elasticForce, dampingForce, totalForce);
                        pSUM(a[i][j][k], totalForce, a[i][j][k]);
                    }
                }
                
                // Bend Springs
                int bend[6][3] = {
                    {2,0,0}, {-2,0,0}, {0,2,0}, {0,-2,0}, {0,0,2}, {0,0,-2}
                };
                
                for (int n = 0; n < 6; n++)
                {
                    int ni = i + bend[n][0];
                    int nj = j + bend[n][1];
                    int nk = k + bend[n][2];
                    
                    if (ni >= 0 && ni < 8 && nj >= 0 && nj < 8 && nk >= 0 && nk < 8)
                    {
                        struct point L, direction, elasticForce, dampingForce, totalForce;
                        struct point neighborPos = jello->p[ni][nj][nk];
                        struct point neighborVel = jello->v[ni][nj][nk];
                        
                        pDIFFERENCE(neighborPos, currentPos, L);
                        double len = pLENGTH(L);
                        
                        pCPY(L, direction);
                        double length;
                        pNORMALIZE(direction);
                        
                        double elasticScalar = jello->kElastic * (len - 2.0/7.0);
                        pMULTIPLY(direction, elasticScalar, elasticForce);
                        
                        struct point velDiff;
                        pDIFFERENCE(neighborVel, currentVel, velDiff);
                        double dampingScalar = jello->dElastic * pDOT(velDiff, direction);
                        pMULTIPLY(direction, dampingScalar, dampingForce);
                        
                        pSUM(elasticForce, dampingForce, totalForce);
                        pSUM(a[i][j][k], totalForce, a[i][j][k]);
                    }
                }
            }
    
    // Collision detection
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
            {
                struct point pos = jello->p[i][j][k];
                struct point vel = jello->v[i][j][k];
                
                // Calculate force
                // Right wall
                if (pos.x > 2.0)
                {
                    double length = pos.x - 2.0;
                    double elasticForce = -jello->kCollision * length;
                    double dampingForce = -jello->dCollision * vel.x;
                    a[i][j][k].x += elasticForce + dampingForce;
                }
                
                // Left wall
                if (pos.x < -2.0)
                {
                    double length = -2.0 - pos.x;
                    double elasticForce = jello->kCollision * length;
                    double dampingForce = -jello->dCollision * vel.x;
                    a[i][j][k].x += elasticForce + dampingForce;
                }
                
                // Front wall
                if (pos.y > 2.0)
                {
                    double length = pos.y - 2.0;
                    double elasticForce = -jello->kCollision * length;
                    double dampingForce = -jello->dCollision * vel.y;
                    a[i][j][k].y += elasticForce + dampingForce;
                }
                
                // Back wall
                if (pos.y < -2.0)
                {
                    double length = -2.0 - pos.y;
                    double elasticForce = jello->kCollision * length;
                    double dampingForce = -jello->dCollision * vel.y;
                    a[i][j][k].y += elasticForce + dampingForce;
                }
                
                // Up wall
                if (pos.z > 2.0)
                {
                    double length = pos.z - 2.0;
                    double elasticForce = -jello->kCollision * length;
                    double dampingForce = -jello->dCollision * vel.z;
                    a[i][j][k].z += elasticForce + dampingForce;
                }
                
                // Down wall
                if (pos.z < -2.0)
                {
                    double length = -2.0 - pos.z;
                    double elasticForce = jello->kCollision * length;
                    double dampingForce = -jello->dCollision * vel.z;
                    a[i][j][k].z += elasticForce + dampingForce;
                }
            }
    
    // Add external force field
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
            {
                struct point externalForce;
                getForceField(jello, jello->p[i][j][k], externalForce);
                a[i][j][k].x += externalForce.x;
                a[i][j][k].y += externalForce.y;
                a[i][j][k].z += externalForce.z;
            }
            
    // Calculate acceleration
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            for (int k = 0; k < 8; k++)
            {
                pMULTIPLY(a[i][j][k], 1.0 / jello->mass, a[i][j][k]);
            }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;
}

